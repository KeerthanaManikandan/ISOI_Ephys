function [tempProfile,tempProfileSuper,tempProfileDeep,tempProfileMid,lags] =...
    getCrossCorrROI(dataDir,monkeyName,expDate,runName, processedDat, probe, rawCh, badTimes,badCh,ch,...
    seedLocIn,seedRad,clipMaskROI,chSplit,refType)
% This function computes the follwing:
% 1. Cross correlations between physiology (different frequencies and layers)
% 2. Plots the temporal profile
% 3. Plots the ROI at peak negative and peak positive correlation


clear inDatTemp inDatUp gammaBand alphaBand fs bA aA bG aG szLFP szIm szMin lags ccROI ccROIAlpha...
    ccROIBeta ccROITheta ccROIRaw ccROISuper ccROIAlphaSuper ccROIBetaSuper ccROIThetaSuper...
    ccROIRawSuper ccROIDeep ccROIAlphaDeep ccROIBetaDeep ccROIThetaDeep ccROIRawDeep

% Get filter parameters...
fs = 1e3;
gammaBand  = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); % Gamma band filtering parameters
alphaBand  = [8 12];  [bA,aA] = butter(3,alphaBand./(fs/2),'bandpass'); % Alpha band filtering parameters
betaBand   = [13 30]; [bB,aB] = butter(3,betaBand./(fs/2),'bandpass');  % Beta band filtering parameters
thetaBand  = [6 8];   [bT,aT] = butter(3,thetaBand./(fs/2),'bandpass'); % Theta band filtering parameters

% Remove bad channels from LFP
probe(:,badCh) = [];
rawCh(:,badCh) = [];
szLFP          = size(probe,1);
imSize         = size(processedDat);
processedDat   = reshape(processedDat,[imSize(1)*imSize(2) imSize(3)]);

% Upsampling imaging data to 10 Hz
disp('Upsampling imaging data to 10 Hz... ');
parfor iP = 1:size(processedDat,1)
    processedDat10(iP,:) = interp(processedDat(iP,:),5);
end

% Get the ROI
processedDat10R = reshape(processedDat10,[imSize(1) imSize(2) size(processedDat10,2)]);

inDat = processedDat10R(seedLocIn(2)-round(seedRad/2):seedLocIn(2)+round(seedRad/2),...
    seedLocIn(1)-round(seedRad/2):seedLocIn(1)+round(seedRad/2),:);

inDatSize = size(inDat);
inDat = reshape(inDat,[inDatSize(1)*inDatSize(2) inDatSize(3)]);

% Make both matrices equal...
szIm = size(processedDat10,2)*100;
if ~(szLFP == szIm)
    szMin     = min([szLFP, szIm]);
    probe     = probe(1:szMin,:);
    rawCh     = rawCh(1:szMin,:);
    badTimes(badTimes>szMin) = [];
else
    szMin = szLFP;
end
probe(badTimes,:) = []; % Remove bad times from LFP
rawCh(badTimes,:) = []; % Remove bad times from raw data

% IMAGING: Upsample ROI (x100 times) and remove concurrent bad frames
tic;
parfor iP = 1:size(inDat,1)
    inDatTemp = interp(inDat(iP,:),100);
    if ~(szLFP == szIm); inDatTemp = inDatTemp(:,1:szMin); end
    inDatTemp(:,badTimes) = []; % Remove bad frames
    inDatUp(iP,:) = downsample(inDatTemp,100); % Downsample imaging data
end
toc;

% Remove bad times determined visually from spectrogram
[probe,rawCh,inDatUp] = removeBadTimesFromSpec(monkeyName,expDate,runName,probe,rawCh,inDatUp);

% Re-Reference physiology
if (strcmp(refType,'AvgRefTop5_Bottom5') || strcmp(refType,'AvgRef')) && (ch(2)-ch(1)~=0)
    probe = probe - mean(probe(:,ch(1):ch(2)),2);  % Average referencing
    rawCh = rawCh - mean(rawCh(:,ch(1):ch(2)),2); 

elseif (strcmp(refType,'BipolarRef_Top5_Bottom5')|| strcmp(refType,'BipolarRef'))&& (ch(2)-ch(1)~=0)
    % Bipolar referencing
    chCount = ch(1);  probeTemp = []; rawChTemp = []; % Bipolar referencing
    while chCount<ch(2)
        probeTemp(:,chCount) = probe(:,chCount)- probe(:,chCount+1); %#ok<AGROW>
        rawChTemp(:,chCount) = rawCh(:,chCount)- rawCh(:,chCount+1); %#ok<AGROW> 
        chCount = chCount+1;
    end
    probe(:,ch(1):ch(2)-1) = probeTemp(:,ch(1):ch(2)-1);
    rawCh(:,ch(1):ch(2)-1) = rawChTemp(:,ch(1):ch(2)-1);
    ch(2) = ch(2)-1; % Bipolar reference reduces the channel count by 1

elseif (strcmp(refType,'CSDRef_Top5_Bottom5')|| strcmp(refType,'CSDRef'))&& (ch(2)-ch(1)~=0)
    for iter = 1:2
         chCount = ch(1);  probeTemp = []; rawChTemp = [];  % CSD referencing
        while chCount<ch(2)
            probeTemp(:,chCount) = probe(:,chCount)- probe(:,chCount+1); %#ok<AGROW>
            rawChTemp(:,chCount) = rawCh(:,chCount)- rawCh(:,chCount+1); %#ok<AGROW>
            chCount = chCount+1;
        end
        probe(:,ch(1):ch(2)-1) = probeTemp(:,ch(1):ch(2)-1);
        rawCh(:,ch(1):ch(2)-1) = rawChTemp(:,ch(1):ch(2)-1);
        ch(2) = ch(2)-1; % CSD reference reduces the channel count by 2
    end
end

% Determine infraslow powers of different frequency bands (ephys)
infraEphys      = getInfraSlowPowerLFP(probe,bG,aG,ch); % Gamma band
infraEphysAlpha = getInfraSlowPowerLFP(probe,bA,aA,ch); % Alpha band
infraEphysBeta  = getInfraSlowPowerLFP(probe,bB,aB,ch); % Beta band
infraEphysTheta = getInfraSlowPowerLFP(probe,bT,aT,ch); % Theta band
infraEphysRaw   = getInfraSlowPowerLFP(rawCh,[],[],ch); % MUA

% Determine infraslow powers for superficial, middle and deep compartments
if ch(2)-ch(1)== 0 
    infraSuper      = mean(infraEphys,2,'omitnan'); % Gamma band 
    infraMid       = mean(infraEphys,2,'omitnan'); 
    infraDeep       = mean(infraEphys,2,'omitnan');

    infraAlphaSuper = mean(infraEphysAlpha,2,'omitnan'); % Alpha band
    infraAlphaMid = mean(infraEphysAlpha,2,'omitnan');
    infraAlphaDeep  = mean(infraEphysAlpha,2,'omitnan');

    infraBetaSuper  = mean(infraEphysBeta,2,'omitnan'); % Beta band
    infraBetaMid   = mean(infraEphysBeta,2,'omitnan');
    infraBetaDeep   = mean(infraEphysBeta,2,'omitnan');

    infraThetaSuper = mean(infraEphysTheta,2,'omitnan'); % Theta band
    infraThetaMid  = mean(infraEphysTheta,2,'omitnan'); 
    infraThetaDeep  = mean(infraEphysTheta,2,'omitnan');

    infraRawSuper   = mean(infraEphysRaw,2,'omitnan'); % Spiking
    infraRawMid    = mean(infraEphysRaw,2,'omitnan');
    infraRawDeep    = mean(infraEphysRaw,2,'omitnan');


else % Single channel, superficial and deep compartments are the same here.
    infraSuper      = mean(infraEphys(:,1:chSplit),2,'omitnan');
    infraDeep       = mean(infraEphys(:,(chSplit*2)+1:end),2,'omitnan');

    infraAlphaSuper = mean(infraEphysAlpha(:,1:chSplit),2,'omitnan');
    infraAlphaDeep  = mean(infraEphysAlpha(:,(chSplit*2)+1:end),2,'omitnan');

    infraBetaSuper  = mean(infraEphysBeta(:,1:chSplit),2,'omitnan');
    infraBetaDeep   = mean(infraEphysBeta(:,(chSplit*2)+1:end),2,'omitnan');

    infraThetaSuper = mean(infraEphysTheta(:,1:chSplit),2,'omitnan');
    infraThetaDeep  = mean(infraEphysTheta(:,(chSplit*2)+1:end),2,'omitnan');

    infraRawSuper   = mean(infraEphysRaw(:,1:chSplit),2,'omitnan');
    infraRawDeep    = mean(infraEphysRaw(:,(chSplit*2)+1:end),2,'omitnan');

    if chSplit~=10  && (ch(2)-ch(1)~= 0) % Get the middle channels if the split is separable
        infraMid      = mean(infraEphys(:,chSplit+1:(chSplit*2)),2,'omitnan');
        infraAlphaMid = mean(infraEphysAlpha(:,chSplit+1:(chSplit*2)),2,'omitnan');
        infraBetaMid  = mean(infraEphysBeta(:,chSplit+1:(chSplit*2)),2,'omitnan');
        infraThetaMid = mean(infraEphysTheta(:,chSplit+1:(chSplit*2)),2,'omitnan');
        infraRawMid   = mean(infraEphysRaw(:,chSplit+1:(chSplit*2)),2,'omitnan');
    else
        infraMid      = NaN(size(infraSuper));
        infraAlphaMid =  NaN(size(infraSuper));
        infraBetaMid  =  NaN(size(infraSuper));
        infraThetaMid =  NaN(size(infraSuper));
        infraRawMid   =  NaN(size(infraSuper));
    end

end

% Get average time series
infraEphys      = mean(infraEphys,2,'omitnan'); 
infraEphysAlpha = mean(infraEphysAlpha,2,'omitnan');
infraEphysBeta  = mean(infraEphysBeta,2,'omitnan');
infraEphysTheta = mean(infraEphysTheta,2,'omitnan');
infraEphysRaw   = mean(infraEphysRaw,2,'omitnan') ;

% Compute cross correlations between physiology and imaging
tic;
disp('Performing cross correlations...');
chLen = (ch(2)-ch(1));
parfor iP = 1:size(inDatUp,1)
    if iP == 1
        [ccROI(:,iP),lags(:,iP)]  = xcorr(infraEphys',inDatUp(iP,:),200,'normalized');
    else
        [ccROI(:,iP),~]  = xcorr(infraEphys',inDatUp(iP,:),200,'normalized');
    end
   
    [ccROIAlpha(:,iP),~] = xcorr(infraEphysAlpha',inDatUp(iP,:),200,'normalized'); % All channels
    [ccROIBeta(:,iP),~]  = xcorr(infraEphysBeta',inDatUp(iP,:),200,'normalized');
    [ccROITheta(:,iP),~] = xcorr(infraEphysTheta',inDatUp(iP,:),200,'normalized');
    [ccROIRaw(:,iP),~]   = xcorr(infraEphysRaw',inDatUp(iP,:),200,'normalized');

    [ccROISuper(:,iP),~]      = xcorr(infraSuper',inDatUp(iP,:),200,'normalized'); % Superficial compartment
    [ccROIAlphaSuper(:,iP),~] = xcorr(infraAlphaSuper',inDatUp(iP,:),200,'normalized');
    [ccROIBetaSuper(:,iP),~]  = xcorr(infraBetaSuper',inDatUp(iP,:),200,'normalized');
    [ccROIThetaSuper(:,iP),~] = xcorr(infraThetaSuper',inDatUp(iP,:),200,'normalized');
    [ccROIRawSuper(:,iP),~]   = xcorr(infraRawSuper',inDatUp(iP,:),200,'normalized');

    [ccROIDeep(:,iP),~]      = xcorr(infraDeep',inDatUp(iP,:),200,'normalized'); % Deep compartment
    [ccROIAlphaDeep(:,iP),~] = xcorr(infraAlphaDeep',inDatUp(iP,:),200,'normalized');
    [ccROIBetaDeep(:,iP),~]  = xcorr(infraBetaDeep',inDatUp(iP,:),200,'normalized');
    [ccROIThetaDeep(:,iP),~] = xcorr(infraThetaDeep',inDatUp(iP,:),200,'normalized');
    [ccROIRawDeep(:,iP),~]   = xcorr(infraRawDeep',inDatUp(iP,:),200,'normalized');

    if chSplit~=10 && (chLen~= 0)
        [ccROIMid(:,iP),~]      = xcorr(infraMid',inDatUp(iP,:),200,'normalized'); % Middle compartment
        [ccROIAlphaMid(:,iP),~] = xcorr(infraAlphaMid',inDatUp(iP,:),200,'normalized');
        [ccROIBetaMid(:,iP),~]  = xcorr(infraBetaMid',inDatUp(iP,:),200,'normalized');
        [ccROIThetaMid(:,iP),~] = xcorr(infraThetaMid',inDatUp(iP,:),200,'normalized');
        [ccROIRawMid(:,iP),~]   = xcorr(infraRawMid',inDatUp(iP,:),200,'normalized');
    end

end

% Reshaping cross correlations
tempProfile.ccROI      = reshape(ccROI,[401 inDatSize(1) inDatSize(2)]);
tempProfile.ccROIAlpha = reshape(ccROIAlpha,[401 inDatSize(1) inDatSize(2)]);
tempProfile.ccROIBeta  = reshape(ccROIBeta,[401 inDatSize(1) inDatSize(2)]);
tempProfile.ccROITheta = reshape(ccROITheta,[401 inDatSize(1) inDatSize(2)]);
tempProfile.ccROIRaw   = reshape(ccROIRaw,[401 inDatSize(1) inDatSize(2)]);

tempProfileSuper.ccROI      = reshape(ccROISuper,[401 inDatSize(1) inDatSize(2)]);
tempProfileSuper.ccROIAlpha = reshape(ccROIAlphaSuper,[401 inDatSize(1) inDatSize(2)]);
tempProfileSuper.ccROIBeta  = reshape(ccROIBetaSuper,[401 inDatSize(1) inDatSize(2)]);
tempProfileSuper.ccROITheta = reshape(ccROIThetaSuper,[401 inDatSize(1) inDatSize(2)]);
tempProfileSuper.ccROIRaw   = reshape(ccROIRawSuper,[401 inDatSize(1) inDatSize(2)]);

tempProfileDeep.ccROI      = reshape(ccROIDeep,[401 inDatSize(1) inDatSize(2)]);
tempProfileDeep.ccROIAlpha = reshape(ccROIAlphaDeep,[401 inDatSize(1) inDatSize(2)]);
tempProfileDeep.ccROIBeta  = reshape(ccROIBetaDeep,[401 inDatSize(1) inDatSize(2)]);
tempProfileDeep.ccROITheta = reshape(ccROIThetaDeep,[401 inDatSize(1) inDatSize(2)]);
tempProfileDeep.ccROIRaw   = reshape(ccROIRawDeep,[401 inDatSize(1) inDatSize(2)]);

if chSplit~=10 && (ch(2)-ch(1)~= 0)
    tempProfileMid.ccROI      = reshape(ccROIMid,[401 inDatSize(1) inDatSize(2)]);
    tempProfileMid.ccROIAlpha = reshape(ccROIAlphaMid,[401 inDatSize(1) inDatSize(2)]);
    tempProfileMid.ccROIBeta  = reshape(ccROIBetaMid,[401 inDatSize(1) inDatSize(2)]);
    tempProfileMid.ccROITheta = reshape(ccROIThetaMid,[401 inDatSize(1) inDatSize(2)]);
    tempProfileMid.ccROIRaw   = reshape(ccROIRawMid,[401 inDatSize(1) inDatSize(2)]);
else
    tempProfileMid.ccROI      = NaN([401 inDatSize(1) inDatSize(2)]);
    tempProfileMid.ccROIAlpha = NaN([401 inDatSize(1) inDatSize(2)]);
    tempProfileMid.ccROIBeta  = NaN([401 inDatSize(1) inDatSize(2)]);
    tempProfileMid.ccROITheta = NaN([401 inDatSize(1) inDatSize(2)]);
    tempProfileMid.ccROIRaw   = NaN([401 inDatSize(1) inDatSize(2)]);
end

toc;

% Initializing variables
tempProfileMid.profile = NaN(401,1);
tempProfileMid.magLow  = NaN;
tempProfileMid.lagLow  = NaN;

tempProfileMid.profileAlpha = NaN(401,1);
tempProfileMid.magLowAlpha  = NaN;
tempProfileMid.lagLowAlpha  = NaN;

tempProfileMid.profileBeta = NaN(401,1);
tempProfileMid.magLowBeta  = NaN;
tempProfileMid.lagLowBeta  = NaN;

tempProfileMid.profileTheta = NaN(401,1);
tempProfileMid.magLowTheta  = NaN;
tempProfileMid.lagLowTheta  = NaN;

tempProfileMid.profileRaw = NaN(401,1);
tempProfileMid.magLowRaw  = NaN;
tempProfileMid.lagLowRaw  = NaN;

% Compute the temporal profiles, peak negative amplitude for different
% compartments and frequencies
for iLayer = 1:4
    clear ccTemp ccTempAlpha ccTempBeta ccTempTheta ccTempRaw tempGamma tempMagGamma tempLagGamma...
        tempAlpha tempMagAlpha tempLagAlpha  tempMagBeta tempLagBeta tempBeta ...
        tempTheta tempMagTheta tempLagTheta tempRaw tempMagRaw tempLagRaw
    switch iLayer
        case 1 % All channels
            ccTemp      = tempProfile.ccROI;
            ccTempAlpha = tempProfile.ccROIAlpha;
            ccTempBeta  = tempProfile.ccROIBeta;
            ccTempTheta = tempProfile.ccROITheta;
            ccTempRaw   = tempProfile.ccROIRaw;
            layerLabel  = 'All';

        case 2 % Superficial compartment
            ccTemp      = tempProfileSuper.ccROI;
            ccTempAlpha = tempProfileSuper.ccROIAlpha;
            ccTempBeta  = tempProfileSuper.ccROIBeta;
            ccTempTheta = tempProfileSuper.ccROITheta;
            ccTempRaw   = tempProfileSuper.ccROIRaw;
            layerLabel  = 'Superficial';

        case 3 % Deep compartment 
            ccTemp      = tempProfileDeep.ccROI;
            ccTempAlpha = tempProfileDeep.ccROIAlpha;
            ccTempBeta  = tempProfileDeep.ccROIBeta;
            ccTempTheta = tempProfileDeep.ccROITheta;
            ccTempRaw   = tempProfileDeep.ccROIRaw;
            layerLabel  = 'Deep';

        case 4 % Middle compartment 
            if chSplit~=10 && (ch(2)-ch(1)~= 0)
                ccTemp      = tempProfileMid.ccROI;
                ccTempAlpha = tempProfileMid.ccROIAlpha;
                ccTempBeta  = tempProfileMid.ccROIBeta;
                ccTempTheta = tempProfileMid.ccROITheta;
                ccTempRaw   = tempProfileMid.ccROIRaw;
                layerLabel  = 'Mid';
            else
                continue;
            end
    end

    lowIdx = lags<0 & lags>=-80; % Set the lag limit to determine peak negative correlations
    xLow   = lags(lowIdx);

    % Get the temporal profiles for gamma band
    tempGamma                    = median(ccTemp,[2,3],'omitnan');
    [tempMagGamma,minMedcorrInd] = min(median(ccTemp(lowIdx,:,:),[2,3],'omitnan'));
    tempLagGamma                 = xLow(minMedcorrInd);
    plotLagProfiles([dataDir '\' refType],['Gamma band ' layerLabel],ccTemp,lags,monkeyName,expDate,...
        runName,clipMaskROI);

    % Get the temporal profiles for alpha band
    tempAlpha                         = median(ccTempAlpha,[2,3],'omitnan');
    [tempMagAlpha,minMedcorrIndAlpha] = min(median(ccTempAlpha(lowIdx,:,:),[2,3],'omitnan'));
    tempLagAlpha                      = xLow(minMedcorrIndAlpha);
    plotLagProfiles([dataDir '\' refType],['Alpha band ' layerLabel],ccTempAlpha,lags,monkeyName,expDate,...
        runName,clipMaskROI);

    % Get the temporal profile for beta band
    tempBeta                        = median(ccTempBeta,[2 3],'omitnan');
    [tempMagBeta,minMedcorrIndBeta] = min(median(ccTempBeta(lowIdx,:,:),[2,3],'omitnan'));
    tempLagBeta                     = xLow(minMedcorrIndBeta);
    plotLagProfiles([dataDir '\' refType],['Beta band ' layerLabel],ccTempBeta,lags,monkeyName,expDate,...
        runName,clipMaskROI);

    % Get the temporal profile for theta band
    tempTheta                         = median(ccTempTheta,[2 3],'omitnan');
    [tempMagTheta,minMedcorrIndTheta] = min(median(ccTempTheta(lowIdx,:,:),[2,3],'omitnan'));
    tempLagTheta                      = xLow(minMedcorrIndTheta);
    plotLagProfiles([dataDir '\' refType],['Theta band ' layerLabel],ccTempTheta,lags,monkeyName,expDate,runName,clipMaskROI);

    % Get the ROI and temporal profile for frequencies>250 Hz cross correlations
    tempRaw                       = median(ccTempRaw,[2 3],'omitnan');
    [tempMagRaw,minMedcorrIndRaw] = min(median(ccTempRaw(lowIdx,:,:),[2,3],'omitnan'));
    tempLagRaw                    = xLow(minMedcorrIndRaw);
    plotLagProfiles([dataDir '\' refType],['Spike band ' layerLabel],ccTempRaw,lags,monkeyName,expDate,runName,clipMaskROI);

    switch iLayer % Organize the profile, magnitude and lag for each frequency in a single structure
        case 1 % All channels
            tempProfile.profile = tempGamma; % Gamma
            tempProfile.magLow  = tempMagGamma;
            tempProfile.lagLow  = tempLagGamma;

            tempProfile.profileAlpha = tempAlpha; % Alpha
            tempProfile.magLowAlpha  = tempMagAlpha;
            tempProfile.lagLowAlpha  = tempLagAlpha;

            tempProfile.profileBeta = tempBeta; % Beta
            tempProfile.magLowBeta  = tempMagBeta;
            tempProfile.lagLowBeta  = tempLagBeta;

            tempProfile.profileTheta = tempTheta; % Theta
            tempProfile.magLowTheta  = tempMagTheta;
            tempProfile.lagLowTheta  = tempLagTheta;

            tempProfile.profileRaw = tempRaw; % Spiking
            tempProfile.magLowRaw  = tempMagRaw;
            tempProfile.lagLowRaw  = tempLagRaw;

        case 2 % Superficial compartment
            tempProfileSuper.profile = tempGamma;% Gamma
            tempProfileSuper.magLow  = tempMagGamma;
            tempProfileSuper.lagLow  = tempLagGamma;

            tempProfileSuper.profileAlpha = tempAlpha;% Alpha
            tempProfileSuper.magLowAlpha = tempMagAlpha;
            tempProfileSuper.lagLowAlpha  = tempLagAlpha;

            tempProfileSuper.profileBeta = tempBeta;% Beta
            tempProfileSuper.magLowBeta  = tempMagBeta;
            tempProfileSuper.lagLowBeta  = tempLagBeta;

            tempProfileSuper.profileTheta = tempTheta;% Theta
            tempProfileSuper.magLowTheta  = tempMagTheta;
            tempProfileSuper.lagLowTheta  = tempLagTheta;

            tempProfileSuper.profileRaw = tempRaw;% Spiking
            tempProfileSuper.magLowRaw  = tempMagRaw;
            tempProfileSuper.lagLowRaw  = tempLagRaw;

        case 3 % Deep compartment
            tempProfileDeep.profile = tempGamma;% Gamma
            tempProfileDeep.magLow  = tempMagGamma;
            tempProfileDeep.lagLow  = tempLagGamma;

            tempProfileDeep.profileAlpha = tempAlpha;% Alpha
            tempProfileDeep.magLowAlpha  = tempMagAlpha;
            tempProfileDeep.lagLowAlpha  = tempLagAlpha;

            tempProfileDeep.profileBeta = tempBeta;% Beta
            tempProfileDeep.magLowBeta  = tempMagBeta;
            tempProfileDeep.lagLowBeta  = tempLagBeta;

            tempProfileDeep.profileTheta = tempTheta;% Theta
            tempProfileDeep.magLowTheta  = tempMagTheta;
            tempProfileDeep.lagLowTheta  = tempLagTheta;

            tempProfileDeep.profileRaw = tempRaw; % Spiking
            tempProfileDeep.magLowRaw  = tempMagRaw;
            tempProfileDeep.lagLowRaw  = tempLagRaw;
    end

    if chSplit~=10 && (ch(2)-ch(1)~= 0) && iLayer == 4 % Middle compartment
        tempProfileMid.profile = tempGamma;% Gamma
        tempProfileMid.magLow  = tempMagGamma;
        tempProfileMid.lagLow  = tempLagGamma;

        tempProfileMid.profileAlpha = tempAlpha;% Alpha
        tempProfileMid.magLowAlpha  = tempMagAlpha;
        tempProfileMid.lagLowAlpha  = tempLagAlpha;

        tempProfileMid.profileBeta = tempBeta;% Beta
        tempProfileMid.magLowBeta  = tempMagBeta;
        tempProfileMid.lagLowBeta  = tempLagBeta;

        tempProfileMid.profileTheta = tempTheta;% Theta
        tempProfileMid.magLowTheta  = tempMagTheta;
        tempProfileMid.lagLowTheta  = tempLagTheta;

        tempProfileMid.profileRaw = tempRaw; % Spiking
        tempProfileMid.magLowRaw  = tempMagRaw;
        tempProfileMid.lagLowRaw  = tempLagRaw;

    end

end
end