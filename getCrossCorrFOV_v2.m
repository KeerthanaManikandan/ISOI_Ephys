function getCrossCorrFOV_v2(monkeyName,expDate,runName,dataDir,processedDat,probe,rawCh,badTimes,badTimeThresh,...
    badCh,ch,timeStamp,chSplit,refType)
% This function calculates the cross correlation between physiology and
% imaging data for the entire FOV. Refer ephysImaging_compilation

clear inDatTemp processedDat10 ccFull ccFullAlpha gammaBand alphaBand fs bA aA bG aG szLFP szIm szMin 

% Get filter parameters...
fs = 1e3; 
gammaBand = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); % Gamma band filtering parameters
alphaBand = [8 12];  [bA,aA] = butter(3,alphaBand./(fs/2),'bandpass'); % Alpha band filtering parameters
betaBand  = [13 30]; [bB,aB] = butter(3,betaBand./(fs/2),'bandpass');  % Beta band filtering parameters
thetaBand = [6 8];   [bT,aT] = butter(3,thetaBand./(fs/2),'bandpass'); % Theta band filtering parameters

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

% Make both matrices equal...
szIm = size(processedDat10,2)*100;

szMin          = min([szLFP, szIm,size(rawCh,1)]);
probe          = probe(1:szMin,:);
rawCh          = rawCh(1:szMin,:);
processedDat10 = processedDat10(:,1:floor(szMin/100));

badTimes(badTimes>szMin)           = [];
badTimeThresh(badTimeThresh>szMin) = [];

% Upsampling full FOV and removing bad frames...
% timeStamp = timeStamp-timeStamp(1);
% timeStamp = int32(floor(timeStamp.*1e3));
% badTimes  = int32(badTimeThresh);
% 
% for iT = 1: length(badTimes)
%     [~,badTimeIm2(iT)] = (min(abs(timeStamp - badTimes(iT))));
% end
% badTimeIm2 = unique(badTimeIm2);

timeStampSorted = timeStamp- timeStamp(1);
badTimes10Hz    = unique(badTimeThresh./1000);
badTimeIm       = [];

% Identifying frames to be removed from RS-ISOI
for iT = 1: length(badTimes10Hz)
    badTimeIm(iT) = find((floor(abs(timeStampSorted - badTimes10Hz(iT))*100)./100)<=0.05,1,'first'); %#ok<*AGROW> 
end

badTimeIm = unique(badTimeIm);
badTimeIm(badTimeIm>size(processedDat10,2)) = [];
processedDat10(:,badTimeIm) = [];
probe(badTimes,:) = []; % Remove bad times from LFP
rawCh(badTimes,:) = [];


% Remove bad times determined visually from spectrogram
[probe,rawCh,processedDat10] = removeBadTimesFromSpec(monkeyName,expDate,runName,probe,rawCh,processedDat10);
[bH,aH] = butter(3,250./(1e3/2),'high'); % High pass filtering >250 Hz

% Get infraslow power fluctuations for alpha and gamma LFP
% Re-Reference physiology
if (strcmp(refType,'BipolarRef_BL'))&& (ch(2)-ch(1)~=0) 

    % Bipolar referencing
    chCount = ch(1);  probeTemp = []; rawTemp = [];
    while chCount<ch(2)
        probeTemp(:,chCount) = probe(:,chCount)- probe(:,chCount+1);
        rawTemp(:,chCount) = rawCh(:,chCount)- rawCh(:,chCount+1);
        chCount = chCount+1;
    end
    probe(:,ch(1):ch(2)-1) = probeTemp(:,ch(1):ch(2)-1);
    rawCh(:,ch(1):ch(2)-1) = rawTemp(:,ch(1):ch(2)-1);
    ch(2) = ch(2)-1; % Bipolar reference reduces the channel count by 1
    
    clear rawTemp
    % Filter and downsample
    rawCh = single(filtfilt(bH,aH,double(rawCh))); % Highpass frequencies beyond 250 Hz
    

    % Determine infraslow powers of different frequency bands (ephys)
    infraEphys      = getInfraSlowPowerLFP(probe,bG,aG,ch(1):ch(2)); % Gamma band
    infraEphysAlpha = getInfraSlowPowerLFP(probe,bA,aA,ch(1):ch(2)); % Alpha band
    infraEphysBeta  = getInfraSlowPowerLFP(probe,bB,aB,ch(1):ch(2)); % Beta band
    infraEphysTheta = getInfraSlowPowerLFP(probe,bT,aT,ch(1):ch(2)); % Theta band
    infraEphysRaw   = getInfraSlowPowerLFP(rawCh,[],[],ch(1):ch(2)); % MUA

elseif (strcmp(refType,'BipolarBL_Ref'))&& (ch(2)-ch(1)~=0)
    % Bandlimit the signals 
    theta = single(filtfilt(bT,aT,double(probe)));
    alpha = single(filtfilt(bA,aA,double(probe)));
    beta = single(filtfilt(bB,aB,double(probe)));
    gamma = single(filtfilt(bG,aG,double(probe)));
    rawCh = single(filtfilt(bH,aH,double(rawCh))); % Highpass frequencies beyond 250 Hz 

    % Re-reference all time courses
    [theta,~] = bipolarRef(theta,ch);
    [alpha,~] = bipolarRef(alpha,ch);
    [beta,~] = bipolarRef(beta,ch);
    [gamma,~] = bipolarRef(gamma,ch);
    [rawCh,ch] = bipolarRef(rawCh,ch);

    % Get infraslow powers
    infraEphys = infraEphys_v2(gamma,ch);
    infraEphysAlpha = infraEphys_v2(alpha,ch); % Alpha band
    infraEphysBeta  = infraEphys_v2(beta,ch); % Beta band
    infraEphysTheta = infraEphys_v2(theta,ch); % Theta band
    infraEphysRaw   = infraEphys_v2(rawCh,ch); % MUA  

elseif (strcmp(refType,'CSDRef_BL'))&& (ch(2)-ch(1)~=0) 
    % CSD referencing
    for iter = 1:2
        chCount = ch(1);  probeTemp = []; rawTemp = [];
        while chCount<ch(2)
            probeTemp(:,chCount) = probe(:,chCount)- probe(:,chCount+1);
            rawTemp(:,chCount) = rawCh(:,chCount)- rawCh(:,chCount+1);
            chCount = chCount+1;
        end
        probe(:,ch(1):ch(2)-1) = probeTemp(:,ch(1):ch(2)-1);
        rawCh(:,ch(1):ch(2)-1) = rawTemp(:,ch(1):ch(2)-1);
        ch(2) = ch(2)-1;
    end
    % CSD reference reduces the channel count by 2

    % Filter and downsample
   rawCh = single(filtfilt(bH,aH,double(rawCh))); % Highpass frequencies beyond 250 Hz

    % Determine infraslow powers of different frequency bands (ephys)
    infraEphys      = getInfraSlowPowerLFP(probe,bG,aG,ch(1):ch(2)); % Gamma band
    infraEphysAlpha = getInfraSlowPowerLFP(probe,bA,aA,ch(1):ch(2)); % Alpha band
    infraEphysBeta  = getInfraSlowPowerLFP(probe,bB,aB,ch(1):ch(2)); % Beta band
    infraEphysTheta = getInfraSlowPowerLFP(probe,bT,aT,ch(1):ch(2)); % Theta band
    infraEphysRaw   = getInfraSlowPowerLFP(rawCh,[],[],ch(1):ch(2)); % MUA


elseif (strcmp(refType,'CSDBL_Ref'))&& (ch(2)-ch(1)~=0)
    % Bandlimit the signals 
    theta = single(filtfilt(bT,aT,double(probe)));
    alpha = single(filtfilt(bA,aA,double(probe)));
    beta = single(filtfilt(bB,aB,double(probe)));
    gamma = single(filtfilt(bG,aG,double(probe)));
    rawCh = single(filtfilt(bH,aH,double(rawCh))); % Highpass frequencies beyond 250 Hz

    % Re-reference all time courses
    [theta,~]  = csdRef(theta,ch);
    [alpha,~]  = csdRef(alpha,ch);
    [beta,~]   = csdRef(beta,ch);
    [gamma,~]  = csdRef(gamma,ch);
    [rawCh,ch] = csdRef(rawCh,ch);

    % Get infraslow powers
    infraEphys      = infraEphys_v2(gamma,ch);
    infraEphysAlpha = infraEphys_v2(alpha,ch); % Alpha band
    infraEphysBeta  = infraEphys_v2(beta,ch); % Beta band
    infraEphysTheta = infraEphys_v2(theta,ch); % Theta band
    infraEphysRaw   = infraEphys_v2(rawCh,ch); % MUA  

else
    % Filter 
    rawCh = single(filtfilt(bH,aH,double(rawCh))); % Highpass frequencies beyond 250 Hz

    infraEphys      = getInfraSlowPowerLFP(probe,bG,aG,ch(1):ch(2)); % Gamma band
    infraEphysAlpha = getInfraSlowPowerLFP(probe,bA,aA,ch(1):ch(2)); % Alpha band
    infraEphysBeta  = getInfraSlowPowerLFP(probe,bB,aB,ch(1):ch(2)); % Beta band
    infraEphysTheta = getInfraSlowPowerLFP(probe,bT,aT,ch(1):ch(2)); % Theta band
    infraEphysRaw   = getInfraSlowPowerLFP(rawCh,[],[],ch(1):ch(2)); % MUA
end

infraMid      = []; infraMidAlpha = []; infraMidBeta  = []; infraMidTheta = []; infraMidRaw   = []; 

if ch(2)-ch(1)== 0
    infraSuper      = mean(infraEphys,2,'omitnan');
    infraDeep       = mean(infraEphys,2,'omitnan');
    infraAlphaSuper = mean(infraEphysAlpha,2,'omitnan');
    infraAlphaDeep  = mean(infraEphysAlpha,2,'omitnan');
    infraBetaSuper  = mean(infraEphysBeta,2,'omitnan');
    infraBetaDeep   = mean(infraEphysBeta,2,'omitnan');
    infraThetaSuper = mean(infraEphysTheta,2,'omitnan');
    infraThetaDeep  = mean(infraEphysTheta,2,'omitnan');
    infraRawSuper   = mean(infraEphysRaw,2,'omitnan');
    infraRawDeep    = mean(infraEphysRaw,2,'omitnan');

else

    
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
        infraMidAlpha = mean(infraEphysAlpha(:,chSplit+1:(chSplit*2)),2,'omitnan');
        infraMidBeta  = mean(infraEphysBeta(:,chSplit+1:(chSplit*2)),2,'omitnan');
        infraMidTheta = mean(infraEphysTheta(:,chSplit+1:(chSplit*2)),2,'omitnan');
        infraMidRaw   = mean(infraEphysRaw(:,chSplit+1:(chSplit*2)),2,'omitnan');
    end 
end

% Infra ephys powers for all channels 
infraEphys      = mean(infraEphys,2,'omitnan');
infraEphysAlpha = mean(infraEphysAlpha,2,'omitnan');
infraEphysBeta  = mean(infraEphysBeta,2,'omitnan');
infraEphysTheta = mean(infraEphysTheta,2,'omitnan');
infraEphysRaw   = mean(infraEphysRaw,2,'omitnan') ;

% Check size of timeseries of both modalities
szIm = size(processedDat10,2);
szLFP = size(infraEphys,1);
szMin           = min([szLFP, szIm, size(infraEphysRaw,1)]);

infraEphys      = infraEphys(1:szMin,:);
infraEphysAlpha = infraEphysAlpha(1:szMin,:);
infraEphysBeta  = infraEphysBeta(1:szMin,:);
infraEphysTheta = infraEphysTheta(1:szMin,:);
infraEphysRaw   = infraEphysRaw(1:szMin,:);

infraSuper = infraSuper(1:szMin,:);
infraDeep  = infraDeep(1:szMin,:);

infraAlphaSuper = infraAlphaSuper(1:szMin,:);
infraAlphaDeep  = infraAlphaDeep(1:szMin,:);

infraBetaSuper = infraBetaSuper(1:szMin,:);
infraBetaDeep  = infraBetaDeep(1:szMin,:);

infraThetaSuper = infraThetaSuper(1:szMin,:);
infraThetaDeep = infraThetaDeep(1:szMin,:);

infraRawSuper = infraRawSuper(1:szMin,:);
infraRawDeep = infraRawDeep(1:szMin,:);

processedDat10  = processedDat10(:,1:szMin);

if chSplit~=10 && (ch(2)-ch(1)~= 0)
    infraMid      = infraMid(1:szMin,:);
    infraMidAlpha = infraMidAlpha(1:szMin,:);
    infraMidBeta  = infraMidBeta(1:szMin,:);
    infraMidTheta = infraMidTheta(1:szMin,:);
    infraMidRaw   = infraMidRaw(1:szMin,:);
end


disp('Performing cross correlations...')

tic;
[~,lagFull]  = xcorr(infraEphys',processedDat10(1,:),200,'normalized');

chLen = ch(2) - ch(1);
parfor iP = 1:size(processedDat10,1)

    [ccFull(:,iP),~]      = xcorr(infraEphys',processedDat10(iP,:),200,'normalized');
    [ccFullAlpha(:,iP),~] = xcorr(infraEphysAlpha',processedDat10(iP,:),200,'normalized');
    [ccFullBeta(:,iP),~]  = xcorr(infraEphysBeta',processedDat10(iP,:),200,'normalized');
    [ccFullTheta(:,iP),~] = xcorr(infraEphysTheta',processedDat10(iP,:),200,'normalized');
    [ccFullRaw(:,iP),~]   = xcorr(infraEphysRaw',processedDat10(iP,:),200,'normalized');

    [ccFullSuper(:,iP),~]      = xcorr(infraSuper',processedDat10(iP,:),200,'normalized');
    [ccFullAlphaSuper(:,iP),~] = xcorr(infraAlphaSuper',processedDat10(iP,:),200,'normalized');
    [ccFullBetaSuper(:,iP),~]  = xcorr(infraBetaSuper',processedDat10(iP,:),200,'normalized');
    [ccFullThetaSuper(:,iP),~] = xcorr(infraThetaSuper',processedDat10(iP,:),200,'normalized');
    [ccFullRawSuper(:,iP),~]   = xcorr(infraRawSuper',processedDat10(iP,:),200,'normalized');

    [ccFullDeep(:,iP),~]      = xcorr(infraDeep',processedDat10(iP,:),200,'normalized');
    [ccFullAlphaDeep(:,iP),~] = xcorr(infraAlphaDeep',processedDat10(iP,:),200,'normalized');
    [ccFullBetaDeep(:,iP),~]  = xcorr(infraBetaDeep',processedDat10(iP,:),200,'normalized');
    [ccFullThetaDeep(:,iP),~] = xcorr(infraThetaDeep',processedDat10(iP,:),200,'normalized');
    [ccFullRawDeep(:,iP),~]   = xcorr(infraRawDeep',processedDat10(iP,:),200,'normalized');

    if chSplit~=10 && (chLen~= 0)
        [ccFullMid(:,iP),~]      = xcorr(infraMid',processedDat10(iP,:),200,'normalized');
        [ccFullAlphaMid(:,iP),~] = xcorr(infraMidAlpha',processedDat10(iP,:),200,'normalized');
        [ccFullBetaMid(:,iP),~]  = xcorr(infraMidBeta',processedDat10(iP,:),200,'normalized');
        [ccFullThetaMid(:,iP),~] = xcorr(infraMidTheta',processedDat10(iP,:),200,'normalized');
        [ccFullRawMid(:,iP),~]   = xcorr(infraMidRaw',processedDat10(iP,:),200,'normalized');
    end

end



% Reshaping cross correlations 
spatialProfile.ccFull      = reshape(ccFull,[401 imSize(1) imSize(2)]);  
spatialProfile.ccFullAlpha = reshape(ccFullAlpha,[401 imSize(1) imSize(2)]);
spatialProfile.ccFullBeta  = reshape(ccFullBeta,[401 imSize(1) imSize(2)]);
spatialProfile.ccFullTheta = reshape(ccFullTheta,[401 imSize(1) imSize(2)]);
spatialProfile.ccFullRaw   = reshape(ccFullRaw,[401 imSize(1) imSize(2)]);

spatialProfileSuper.ccFull      = reshape(ccFullSuper,[401 imSize(1) imSize(2)]);  
spatialProfileSuper.ccFullAlpha = reshape(ccFullAlphaSuper,[401 imSize(1) imSize(2)]);
spatialProfileSuper.ccFullBeta  = reshape(ccFullBetaSuper,[401 imSize(1) imSize(2)]);
spatialProfileSuper.ccFullTheta = reshape(ccFullThetaSuper,[401 imSize(1) imSize(2)]);
spatialProfileSuper.ccFullRaw   = reshape(ccFullRawSuper,[401 imSize(1) imSize(2)]);

spatialProfileDeep.ccFull      = reshape(ccFullDeep,[401 imSize(1) imSize(2)]);  
spatialProfileDeep.ccFullAlpha = reshape(ccFullAlphaDeep,[401 imSize(1) imSize(2)]);
spatialProfileDeep.ccFullBeta  = reshape(ccFullBetaDeep,[401 imSize(1) imSize(2)]);
spatialProfileDeep.ccFullTheta  = reshape(ccFullThetaDeep,[401 imSize(1) imSize(2)]);
spatialProfileDeep.ccFullRaw   = reshape(ccFullRawDeep,[401 imSize(1) imSize(2)]);

if chSplit~=10 && (ch(2)-ch(1)~= 0)
    spatialProfileMid.ccFull      = reshape(ccFullMid,[401 imSize(1) imSize(2)]);
    spatialProfileMid.ccFullAlpha = reshape(ccFullAlphaMid,[401 imSize(1) imSize(2)]);
    spatialProfileMid.ccFullBeta  = reshape(ccFullBetaMid,[401 imSize(1) imSize(2)]);
    spatialProfileMid.ccFullTheta = reshape(ccFullThetaMid,[401 imSize(1) imSize(2)]);
    spatialProfileMid.ccFullRaw   = reshape(ccFullRawMid,[401 imSize(1) imSize(2)]);
end 
toc;

disp('Saving cross correlations for the entire FOV...');
if chSplit~=10 && (ch(2)-ch(1)~= 0)
    save([dataDir '\crossCorrFOV_' num2str(chSplit) '_' refType '.mat'],'spatialProfile','spatialProfileSuper',...
        'spatialProfileDeep','spatialProfileMid','lagFull','-v7.3');
else
    save([dataDir '\crossCorrFOV_' num2str(chSplit) '_' refType '.mat'],'spatialProfile','spatialProfileSuper',...
        'spatialProfileDeep','lagFull','-v7.3');
end
toc;

end

% Functions used
function [timeCourse,ch] = bipolarRef(timeCourse,ch) % Get bipolar referencing
chCount = ch(1);  timeCourseTemp = []; 
while chCount<ch(2)
    timeCourseTemp(:,chCount) = timeCourse(:,chCount)- timeCourse(:,chCount+1);
    chCount = chCount+1;
end
timeCourse(:,ch(1):ch(2)-1) = timeCourseTemp(:,ch(1):ch(2)-1);
ch(2) = ch(2)-1; % Bipolar reference reduces the channel count by 1
end

function  [timeCourse,ch] = csdRef(timeCourse,ch)  % Get CSD referencing
% CSD referencing
for iter = 1:2
    chCount = ch(1);  timeCourseTemp = [];
    while chCount<ch(2)
        timeCourseTemp(:,chCount) = timeCourse(:,chCount)- timeCourse(:,chCount+1);
        chCount = chCount+1;
    end
    timeCourse(:,ch(1):ch(2)-1) = timeCourseTemp(:,ch(1):ch(2)-1);
    ch(2) = ch(2)-1;
    % CSD reference reduces the channel count by 2
end        
end

function infraEphys = infraEphys_v2(probeDat,chInCortex)  % Get infraslow ephys

% Get the envelope of the signal
envelopeDat = envelope(abs(probeDat(:,chInCortex(1):chInCortex(2))),5); %envelope(probeDat,5,'peak');

clear z p k sos g enSize
fs      = 1e3;

% Bandpass - 0.01 Hz - 0.1 Hz
[z,p,k] = butter(3,[0.01 0.1]./(fs/2),'bandpass');
[sos,g] = zp2sos(z,p,k);
enSize  = size(envelopeDat);
envelopeFiltered = filtfilt(sos,g,double([envelopeDat; envelopeDat ;envelopeDat ]));
envelopeFiltered = envelopeFiltered(enSize(1)+1:(end-enSize(1)),:);
infraEphys       = single(downsample(envelopeFiltered,100));

end