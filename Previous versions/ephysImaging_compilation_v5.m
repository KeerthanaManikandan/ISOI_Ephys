% ephysImaging_compilation_v5
% This script compiles the data procesed from ephysImaging_v3 and saves the
% processed variables so that combining across animals become simpler
% See ephysImaging_v3 for reference
% see ephysImaging_compilation.m and ephysImaging_compilation_v4 for
% previous versions of this script
% October 7, 2024 - KM
% Set paths
clc; clear;
commonDir = 'C:\Users\kem294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\Imaging']));
addpath(genpath([commonDir '\Codes\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\fly_track\videoIO']));

%% Initialize variables and get monkey data
hemisphere = 'Left'; spatialBin = 3;
iM = 2;
switch iM
    case 1
        monkeyName = 'CharlieSheen';
        goodRuns  = ([1 1 NaN NaN NaN;... % Date x run
            1 1 NaN NaN NaN; ...
            1 1 1 0 1; ...
            1 1 1 1 NaN]);

        singleChFlag  = ([0 0 NaN NaN NaN;... % Date x run
            0 0 NaN NaN NaN; ...
            0 0 0 0 0; ...
            0 0 0 0 NaN]);

        goodRunsSpatial  = ([1 1 NaN NaN NaN;... % Date x run
            1 1 NaN NaN NaN; ...
            1 1 1 0 1; ...
            1 1 1 1 NaN]);

        smFlag  = (['S' 'M' '#' '#' '#';... % Date x run
            'S' 'M' '#' '#' '#'; ...
            'S' 'S' 'M' 'M' 'M'; ...
            'S' 'M' 'M' 'M' '#']); 

        isoLevel = ([0.75 0.75 NaN NaN NaN;... % Date x run
            0.75 0.75 NaN NaN NaN; ...
            0.9 0.9 0.9 0.9 0.9; ...
            0.75 0.8 0.8 0.9 NaN]); 

    case 2
        monkeyName = 'Whiskey';
        goodRuns   = ([1 1 1 NaN NaN NaN NaN; ... % Date x run
            1 1 1 1 1 1 1 ; ...
            1 0 1 1 NaN NaN NaN; ...
            1 1 0 1 1 1 0; ...
            1 0 0 1 1 1 1]);

        singleChFlag = ([1 1 1 NaN NaN NaN NaN; ... % Date x run
            0 0 0 0 0 0 0 ; ...
            0 0 0 0 NaN NaN NaN; ...
            0 0 0 0 0 0 0; ...
            0 0 0 0 0 0 0]);

        goodRunsSpatial = ([1 1 1 NaN NaN NaN NaN; ... % Date x run
            1 1 1 1 1 1 1 ; ...
            1 0 1 1 NaN NaN NaN; ...
            1 1 0 1 1 1 1; ...
            1 0 0 1 1 1 0]);

         smFlag = (['M' 'M' 'M' '#' '#' '#' '#'; ... % Date x run; # - in place of NaNs
            'M' 'M' 'M' 'M' 'S' 'S' 'S' ; ...
            'S' 'M' 'M' 'S' '#' '#' '#'; ...
            'S' 'M' 'M' 'S' 'M' 'M' 'M'; ...
            'S' 'S' 'M' 'M' 'M' 'S' 'M']);
         
         isoLevel = ([1 1.1 1.25 NaN NaN NaN NaN; ... % Date x run
            1.2 1 1 1.1 1.1 1.1 1.1 ; ...
            0.8 1.05 1.75 1.75 NaN NaN NaN; ...
            0.9 1 1 1 1 1 1; ...
            0.7 0.8 1 1.1 1.3 1.3 1.3]);
end

goodRuns = reshape(goodRuns,[size(goodRuns,1)*size(goodRuns,2) 1]);
goodRuns(isnan(goodRuns)) = []; goodRuns = logical(goodRuns);

singleChFlag = reshape(singleChFlag,[size(singleChFlag,1)*size(singleChFlag,2) 1]);
singleChFlag(isnan(singleChFlag)) = []; singleChFlag = logical(singleChFlag); 

goodRunsSpatial = reshape(goodRunsSpatial,[size(goodRunsSpatial,1)*size(goodRunsSpatial,2) 1]);
goodRunsSpatial(isnan(goodRunsSpatial)) = []; goodRunsSpatial = logical(goodRunsSpatial);

smFlag = reshape(smFlag,[size(smFlag,1)*size(smFlag,2) 1]);
smFlag((smFlag == '#')) = [];
sensoryGoodRuns = (smFlag =='S') & goodRuns; 
sensoryGoodSpatialRuns = (smFlag == 'S') & goodRunsSpatial;
motorGoodRuns = (smFlag == 'M') & goodRuns; 
motorGoodSpatialRuns = (smFlag == 'M') & goodRuns; 

isoLevel = reshape(isoLevel,[size(isoLevel,1)*size(isoLevel,2) 1]);
isoLevel(isnan(isoLevel)) = []; 
isoLevelGoodRuns = isoLevel; isoLevelGoodRuns(~goodRuns) = [];
isoLevelSpatial = isoLevel; isoLevelSpatial(~goodRunsSpatial) = [];

% Get monkey data and parameters
[allDates,allRuns, refDate, refDir,lensCombo, roiSize, ephysFileNameAll, serverPath,probeLabel,...
    chInCortexNotes, greenMapRef] = getMonkeyParams_Imaging_Ephys(monkeyName,commonDir, hemisphere);

% Get monkey data....
[processedDat,greenIm,probe,badCh,badTimesLFP,badTimeThresh,estChInCortex] = ...
    getAllData_Imaging_Ephys(monkeyName,hemisphere,commonDir,serverPath,allDates,allRuns,...
    ephysFileNameAll,greenMapRef,chInCortexNotes,probeLabel,spatialBin);

clc; disp(['All physiology and imaging data for ' monkeyName ' loaded']);

%% Correlate LFP between compartments
tic;
fs = 1e3;
gammaBand = [30 90]; % Gamma band filtering parameters
alphaBand = [8 12];
betaBand  = [13 30];
thetaBand = [6 8];

params.Fs     = fs;
params.fpass  = [1 120];
params.pad    = -1;
params.tapers = [3 5];

% paramsRaw       = params;
% paramsRaw.fpass = [250 500];
if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\refLFPCorr.mat'],'file')
    superMid = NaN(size(probe,2),size(probe,1),5,3); % LFP
    midDeep = NaN(size(probe,2),size(probe,1),5,3);
    superDeep = NaN(size(probe,2),size(probe,1),5,3);

    superMidPow = NaN(size(probe,2),size(probe,1),5,3); % Powers
    midDeepPow = NaN(size(probe,2),size(probe,1),5,3);
    superDeepPow = NaN(size(probe,2),size(probe,1),5,3);

    infraSM = NaN(size(probe,2),size(probe,1),5,3); % Infraslow powers
    infraMD = NaN(size(probe,2),size(probe,1),5,3);
    infraSD = NaN(size(probe,2),size(probe,1),5,3);

    for iDate = 1: size(allDates,1)
        clear expDate
        expDate = allDates(iDate,:);

        for iRun = 1:size(allRuns{iDate,1})
            clear runName dataDir
            runName = allRuns{iDate,1}(iRun,:);
            dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

            clc; disp([monkeyName ' '  expDate ' run: ' runName]);
            probeCh    = probe{iRun,iDate}.probeCh;
            rawChTemp  = rawCh{iRun,iDate}.rawChTemp;
            timeStamp  = probe{iRun,iDate}.timeStamp;
            badTimes   = badTimesLFP{iDate,iRun};
            chInCortex = estChInCortex{1,iDate}(iRun,:);

            % Removing bad channels
            probeCh(:,badCh{iDate,iRun}) = [];
            rawChTemp(:,badCh{iDate,iRun})   = [];

            % Remove bad times
            probeCh(badTimes,:) = [];
            rawChTemp(badTimes,:)   = [];
            szLFP   = size(probeCh,1);

            if chInCortex(1)-chInCortex(2) ~= 0
                clear probeChCortex spec aLim bLim gLim spec timeValsSpec freqValsSpec
                for iRef = 1:3
                    clear chInCortex probeChCortex rawChCortex
                    chInCortex = estChInCortex{1,iDate}(iRun,:);

                    switch iRef
                        case 1 % No ref
                            probeChCortex = probeCh;
                            rawChCortex   = rawChTemp;

                        case 2 % Bipolar ref
                            probeChCortex = probeCh; rawChCortex = rawChTemp;
                            chCount = chInCortex(1);  probeTemp = []; rawTemp = []; % Average referencing

                            while chCount<chInCortex(2)
                                probeTemp(:,chCount) = probeChCortex(:,chCount)- probeChCortex(:,chCount+1);
                                rawTemp(:,chCount)   = rawChCortex(:,chCount)- rawChCortex(:,chCount+1);
                                chCount = chCount+1;
                            end

                            probeChCortex(:,chInCortex(1):chInCortex(2)-1) = probeTemp(:,chInCortex(1):chInCortex(2)-1);
                            rawChCortex(:,chInCortex(1):chInCortex(2)-1) = rawTemp(:,chInCortex(1):chInCortex(2)-1);
                            chInCortex(2) = chInCortex(2)-1; % Bipolar reference reduces the channel count by 1

                        case 3 % CSD ref
                            probeChCortex = probeCh; rawChCortex = rawChTemp;
                            for iter = 1:2
                                chCount = chInCortex(1);  probeTemp = []; rawTemp = []; % Average referencing

                                while chCount<chInCortex(2)
                                    probeTemp(:,chCount) = probeChCortex(:,chCount)- probeChCortex(:,chCount+1);
                                    rawTemp(:,chCount) = rawChCortex(:,chCount)- rawChCortex(:,chCount+1);
                                    chCount = chCount+1;
                                end
                                probeChCortex(:,chInCortex(1):chInCortex(2)-1) = probeTemp(:,chInCortex(1):chInCortex(2)-1);
                                rawChCortex(:,chInCortex(1):chInCortex(2)-1) = rawTemp(:,chInCortex(1):chInCortex(2)-1);
                                chInCortex(2) = chInCortex(2)-1; % CSD reference reduces the channel count by 2
                            end
                    end

                    for iBand = 1:5
                        switch iBand
                            case 1
                                [bCoeff,aCoeff] = butter(3,thetaBand./(fs/2),'bandpass');
                                bandName = 'Theta';
                            case 2
                                [bCoeff,aCoeff] = butter(3,alphaBand./(fs/2),'bandpass');
                                bandName = 'Alpha';
                            case 3
                                [bCoeff,aCoeff] = butter(3,betaBand./(fs/2),'bandpass');
                                bandName = 'Beta';
                            case 4
                                [bCoeff,aCoeff] = butter(3,gammaBand./(fs/2),'bandpass');
                                bandName = 'Gamma';
                            case 5
                                [bCoeff,aCoeff] = butter(3,250./(fs/2),'high');
                                bandName = 'Spiking';
                        end

                        clear timeCourse infraPow
                        if iBand~=5
                            timeCourse = single(filtfilt(bCoeff,aCoeff,double(probeChCortex(:,chInCortex(1):chInCortex(2)))));
                            infraPow = getInfraSlowPowerLFP(probeChCortex,bCoeff,aCoeff,chInCortex);

                        else
                            timeCourse = single(filtfilt(bCoeff,aCoeff,double(rawChCortex(:,chInCortex(1):chInCortex(2)))));           
                            infraPow = getInfraSlowPowerLFP(single(filtfilt(bCoeff,aCoeff,double(rawChCortex))),[],[],chInCortex);
                        end

                        % Time series
                        clear super mid deep superInfra midInfra deepInfra
                        envelopeDat = envelope(abs(timeCourse),5);

                        super = mean(timeCourse(:,1:6),2,'omitnan');
                        mid   = mean(timeCourse(:,7:12),2,'omitnan');
                        deep  = mean(timeCourse(:,13:end),2,'omitnan');

                        superMid(iDate,iRun,iBand,iRef)  = corr(super,mid,'rows','complete');
                        midDeep(iDate,iRun,iBand,iRef)   = corr(mid,deep,'rows','complete');
                        superDeep(iDate,iRun,iBand,iRef) = corr(super,deep,'rows','complete');
                        
                        % Powers
                        superPow = mean(envelopeDat(:,1:6),2,'omitnan');
                        midPow   = mean(envelopeDat(:,7:12),2,'omitnan');
                        deepPow  = mean(envelopeDat(:,13:end),2,'omitnan');

                        superMidPow(iDate,iRun,iBand,iRef)  = corr(superPow,midPow,'rows','complete');
                        midDeepPow(iDate,iRun,iBand,iRef)   = corr(midPow,deepPow,'rows','complete');
                        superDeepPow(iDate,iRun,iBand,iRef) = corr(superPow,deepPow,'rows','complete');

                        % Infraslow powers
                        superInfra = mean(infraPow(:,1:6),2,'omitnan');
                        midInfra   = mean(infraPow(:,7:12),2,'omitnan');
                        deepInfra  = mean(infraPow(:,13:end),2,'omitnan');

                        infraSM(iDate,iRun,iBand,iRef) = corr(superInfra,midInfra,'rows','complete');
                        infraMD(iDate,iRun,iBand,iRef) = corr(midInfra,deepInfra,'rows','complete');
                        infraSD(iDate,iRun,iBand,iRef) = corr(superInfra,deepInfra,'rows','complete');


                    end
                end
            else
                superMid(iDate,iRun,1:5,1:3) = 0;
                midDeep(iDate,iRun,1:5,1:3)  = 0;
                superDeep(iDate,iRun,1:5,1:3) = 0;

                superMidPow(iDate,iRun,1:5,1:3) = 0;
                midDeepPow(iDate,iRun,1:5,1:3)  = 0;
                superDeepPow(iDate,iRun,1:5,1:3) = 0;

                infraSM(iDate,iRun,1:5,1:3) = 0;
                infraMD(iDate,iRun,1:5,1:3)  = 0;
                infraSD(iDate,iRun,1:5,1:3) = 0;

            end
        end
    end

    toc;

    szVar = size(superMid);
    superMid  = reshape(superMid,[szVar(1)*szVar(2) szVar(3) szVar(4)]);
    midDeep   = reshape(midDeep,[szVar(1)*szVar(2) szVar(3) szVar(4)]);
    superDeep = reshape(superDeep,[szVar(1)*szVar(2) szVar(3) szVar(4)]);

    superMidPow  = reshape(superMidPow,[szVar(1)*szVar(2) szVar(3) szVar(4)]);
    midDeepPow   = reshape(midDeepPow,[szVar(1)*szVar(2) szVar(3) szVar(4)]);
    superDeepPow = reshape(superDeepPow,[szVar(1)*szVar(2) szVar(3) szVar(4)]);

    infraSM = reshape(infraSM,[szVar(1)*szVar(2) szVar(3) szVar(4)]);
    infraMD = reshape(infraMD,[szVar(1)*szVar(2) szVar(3) szVar(4)]);
    infraSD = reshape(infraSD,[szVar(1)*szVar(2) szVar(3) szVar(4)]);


    nanRows = find(isnan(superMid(:,1)));
    superMid(nanRows,:,:) = []; midDeep(nanRows,:,:) = [];  superDeep(nanRows,:,:) = [];
    superMidPow(nanRows,:,:) = []; midDeepPow(nanRows,:,:) = [];  superDeepPow(nanRows,:,:) = [];
    infraSM(nanRows,:,:) = []; infraMD(nanRows,:,:) = [];  infraSD(nanRows,:,:) = [];

    % Time series
    sensorySM = superMid; sensorySM(~(~singleChFlag & sensoryGoodSpatialRuns),:,:) = [];
    sensoryMD = midDeep; sensoryMD(~(~singleChFlag & sensoryGoodSpatialRuns),:,:) = [];
    sensorySD = superDeep; sensorySD(~(~singleChFlag & sensoryGoodSpatialRuns),:,:) = [];

    motorSM = superMid; motorSM(~(~singleChFlag & motorGoodSpatialRuns),:,:) = [];
    motorMD = midDeep;  motorMD(~(~singleChFlag & motorGoodSpatialRuns),:,:) = [];
    motorSD = superDeep; motorSD(~(~singleChFlag & motorGoodSpatialRuns),:,:) = [];

    superMid(~(~singleChFlag & goodRunsSpatial),:,:) = [];  midDeep(~(~singleChFlag & goodRunsSpatial),:,:) = [];
    superDeep(~(~singleChFlag & goodRunsSpatial),:,:) = [];

    % Powers
    sensoryPowSM = superMidPow; sensoryPowSM(~(~singleChFlag & sensoryGoodSpatialRuns),:,:) = [];
    sensoryPowMD = midDeepPow; sensoryPowMD(~(~singleChFlag & sensoryGoodSpatialRuns),:,:) = [];
    sensoryPowSD = superDeepPow; sensoryPowSD(~(~singleChFlag & sensoryGoodSpatialRuns),:,:) = [];

    motorPowSM = superMidPow; motorPowSM(~(~singleChFlag & motorGoodSpatialRuns),:,:) = [];
    motorPowMD = midDeepPow;  motorPowMD(~(~singleChFlag & motorGoodSpatialRuns),:,:) = [];
    motorPowSD = superDeepPow; motorPowSD(~(~singleChFlag & motorGoodSpatialRuns),:,:) = [];

    superMidPow(~(~singleChFlag & goodRunsSpatial),:,:) = [];  midDeepPow(~(~singleChFlag & goodRunsSpatial),:,:) = [];
    superDeepPow(~(~singleChFlag & goodRunsSpatial),:,:) = [];

    % Infraslow powers
    sensorySMInfra = infraSM; sensorySMInfra(~(~singleChFlag & sensoryGoodSpatialRuns),:,:) = [];
    sensoryMDInfra = infraMD; sensoryMDInfra(~(~singleChFlag & sensoryGoodSpatialRuns),:,:) = [];
    sensorySDInfra = infraSD; sensorySDInfra(~(~singleChFlag & sensoryGoodSpatialRuns),:,:) = [];

    motorSMInfra = infraSM; motorSMInfra(~(~singleChFlag & motorGoodSpatialRuns),:,:) = [];
    motorMDInfra = infraMD;  motorMDInfra(~(~singleChFlag & motorGoodSpatialRuns),:,:) = [];
    motorSDInfra = infraSD; motorSDInfra(~(~singleChFlag & motorGoodSpatialRuns),:,:) = [];

    infraSM(~(~singleChFlag & goodRunsSpatial),:,:) = [];  infraMD(~(~singleChFlag & goodRunsSpatial),:,:) = [];
    infraSD(~(~singleChFlag & goodRunsSpatial),:,:) = [];


    bandLabels = {'Theta';'Alpha';'Beta';'Gamma';'Spiking'};
    save(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\refLFPCorr.mat'],'superMid','midDeep','superDeep',...
        'sensorySM','sensorySD','sensoryMD','motorSM','motorSD','motorMD','superMidPow','midDeepPow','superDeepPow',...
        'sensoryPowSM','sensoryPowSD','sensoryPowMD','motorPowSM','motorPowSD','motorPowMD','infraSM','infraSD','infraMD',...
        'sensorySMInfra','sensoryMDInfra','sensorySDInfra','motorSMInfra','motorSDInfra','motorMDInfra');
else
    load(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\refLFPCorr.mat'],'superMid','midDeep','superDeep',...
        'sensorySM','sensorySD','sensoryMD','motorSM','motorSD','motorMD','superMidPow','midDeepPow','superDeepPow',...
        'sensoryPowSM','sensoryPowSD','sensoryPowMD','motorPowSM','motorPowSD','motorPowMD','infraSM','infraSD','infraMD',...
        'sensorySMInfra','sensoryMDInfwinLenra','sensorySDInfra','motorSMInfra','motorSDInfra','motorMDInfra');
end

%% Obtain the correlations 
allMonkeyVars(2) =  load(['X:\Data\Whiskey_SqM\' hemisphere ' Hemisphere\refLFPCorr.mat']);
allMonkeyVars(1) =  load(['X:\Data\CharlieSheen_SqM\' hemisphere ' Hemisphere\refLFPCorr.mat']);
monkeys = {'Charlie Sheen'; 'Whiskey';'Combined'};

iBand = 4; iRef = 1;

% Focus on gamma band and no re-referencing
timeSeriesGammaSM = [allMonkeyVars(1).superMid(:,iBand,iRef); allMonkeyVars(2).superMid(:,iBand,iRef)];
timeSeriesGammaMD = [allMonkeyVars(1).midDeep(:,iBand,iRef); allMonkeyVars(2).midDeep(:,iBand,iRef)];
timeSeriesGammaSD = [allMonkeyVars(1).superDeep(:,iBand,iRef); allMonkeyVars(2).superDeep(:,iBand,iRef)];

powGammaSM = [allMonkeyVars(1).superMidPow(:,iBand,iRef); allMonkeyVars(2).superMidPow(:,iBand,iRef)];
powGammaMD = [allMonkeyVars(1).midDeepPow(:,iBand,iRef); allMonkeyVars(2).midDeepPow(:,iBand,iRef)];
powGammaSD = [allMonkeyVars(1).superDeepPow(:,iBand,iRef); allMonkeyVars(2).superDeepPow(:,iBand,iRef)];

infraGammaSM = [allMonkeyVars(1).infraSM(:,iBand,iRef); allMonkeyVars(2).infraSM(:,iBand,iRef)];
infraGammaMD = [allMonkeyVars(1).infraMD(:,iBand,iRef); allMonkeyVars(2).infraMD(:,iBand,iRef)];
infraGammaSD = [allMonkeyVars(1).infraSD(:,iBand,iRef); allMonkeyVars(2).infraSD(:,iBand,iRef)];

% Comparing betwween timescales
timeSeries = [timeSeriesGammaSM; timeSeriesGammaMD; timeSeriesGammaSD];
powAll     = [powGammaSM; powGammaMD; powGammaSD];
infraAll   = [infraGammaSM; infraGammaMD; infraGammaSD];

figure; boxplot([timeSeries powAll infraAll],{'Time series';'Power';'Infraslow power'});
box off; ylim([-1 1]); title('Time scale comparison');

% Comparing between compartments
superMid = [timeSeriesGammaSM; powGammaSM;infraGammaSM];
midDeep  = [timeSeriesGammaMD; powGammaMD; infraGammaMD];
superDeep = [timeSeriesGammaSD; powGammaSD; infraGammaSD];

figure; boxplot([superMid midDeep superDeep],{'Superficial-Middle';'Middle-Deep';'Superficial-Deep'});
box off; ylim([-1 1]); title('Compartment comparison');


figure; 
subplot(131); boxplot([timeSeriesGammaSM timeSeriesGammaMD timeSeriesGammaSD],{'Superficial-Middle';'Middle-Deep';'Superficial-Deep'});
ylim([-0.5 1]); box off; title('Time series'); yticks(-0.5:0.1:1);
subplot(132); boxplot([powGammaSM powGammaMD powGammaSD],{'Superficial-Middle';'Middle-Deep';'Superficial-Deep'});
ylim([-0.5 1]); box off; title('Power');yticks(-0.5:0.1:1);
subplot(133); boxplot([infraGammaSM infraGammaMD infraGammaSD],{'Superficial-Middle';'Middle-Deep';'Superficial-Deep'});
ylim([-0.5 1]); box off; title('Infra slow power');yticks(-0.5:0.1:1);
sgtitle('Gamma band- no ref-compartment correlations');

figure;
xVar = [1 2 3];
yVarTS = [median(timeSeriesGammaSM) median(timeSeriesGammaMD) median(timeSeriesGammaSD)];
yVarP = [median(powGammaSM) median(powGammaMD) median(powGammaSD)];
yVarIS = [median(infraGammaSM) median(infraGammaMD) median(infraGammaSD)];
plot(xVar,yVarTS,'LineWidth',1); hold on; plot(xVar,yVarP,'LineWidth',1); plot(xVar,yVarIS,'LineWidth',1);
legend('Time series', 'Power','Infraslow power','Location','southwest'); ylim([0 1]); box off;
xticks(1:3);xticklabels({'Superficial-middle', 'Middle-deep','Superficial-deep'}); 

%% Plot all variables
for iVar = 1:3
    clear SM MD SD
    switch iVar
        case 1
            SM = [allMonkeyVars(1).superMid; allMonkeyVars(2).superMid];
            MD = [allMonkeyVars(1).midDeep; allMonkeyVars(2).midDeep]; 
            SD = [allMonkeyVars(1).superDeep; allMonkeyVars(2).superDeep];

            SMSensory = [allMonkeyVars(1).sensorySM; allMonkeyVars(2).sensorySM];
            MDSensory = [allMonkeyVars(1).sensoryMD; allMonkeyVars(2).sensoryMD];
            SDSensory = [allMonkeyVars(1).sensorySD; allMonkeyVars(2).sensorySD];

            SMMotor = [allMonkeyVars(1).motorSM; allMonkeyVars(2).motorSM];
            MDMotor = [allMonkeyVars(1).motorMD; allMonkeyVars(2).motorMD];
            SDMotor = [allMonkeyVars(1).motorSD; allMonkeyVars(2).motorSD]; 
           
            varName = 'Time series';
            
        case 2

            SM = [allMonkeyVars(1).superMidPow; allMonkeyVars(2).superMidPow];
            MD = [allMonkeyVars(1).midDeepPow; allMonkeyVars(2).midDeepPow];
            SD = [allMonkeyVars(1).superDeepPow; allMonkeyVars(2).superDeepPow];

            SMSensory = [allMonkeyVars(1).sensoryPowSM; allMonkeyVars(2).sensoryPowSM];
            MDSensory = [allMonkeyVars(1).sensoryPowMD; allMonkeyVars(2).sensoryPowMD];
            SDSensory = [allMonkeyVars(1).sensoryPowSD; allMonkeyVars(2).sensoryPowSD];

            SMMotor = [allMonkeyVars(1).motorPowSM; allMonkeyVars(2).motorPowSM];
            MDMotor = [allMonkeyVars(1).motorPowMD; allMonkeyVars(2).motorPowMD];
            SDMotor = [allMonkeyVars(1).motorPowSD; allMonkeyVars(2).motorPowSD];

            varName = 'Power';

        case 3

            SM = [allMonkeyVars(1).infraSM; allMonkeyVars(2).infraSM];
            MD = [allMonkeyVars(1).infraMD; allMonkeyVars(2).infraMD];
            SD = [allMonkeyVars(1).infraSD; allMonkeyVars(2).infraSD];

            SMSensory = [allMonkeyVars(1).sensorySMInfra; allMonkeyVars(2).sensorySMInfra];
            MDSensory = [allMonkeyVars(1).sensoryMDInfra; allMonkeyVars(2).sensoryMDInfra];
            SDSensory = [allMonkeyVars(1).sensorySDInfra; allMonkeyVars(2).sensorySDInfra];

            SMMotor = [allMonkeyVars(1).motorSMInfra; allMonkeyVars(2).motorSMInfra];
            MDMotor = [allMonkeyVars(1).motorMDInfra; allMonkeyVars(2).motorMDInfra];
            SDMotor = [allMonkeyVars(1).motorSDInfra; allMonkeyVars(2).motorSDInfra];

            varName = 'Infraslow powers';
    end 
    for iRef = 1:3
        switch iRef
            case 1
                refName = 'No ref';
            case 2
                refName = 'Bipolar';
            case 3
                refName = 'CSD';
        end
        figure;
        for iBand = 1:5
            subplot(2,3,iBand);
            boxplot([SM(:,iBand,iRef)'; MD(:,iBand,iRef)'; SD(:,iBand,iRef)']',{'S_M';'M_D';'S_D'});
            title(bandLabels{iBand}); ylim([-1 1]); box off;
        end
        sgtitle([refName ' - ' varName]);

        figure;
        for iBand = 1:5
            subplot(2,3,iBand);
            boxplot([SMSensory(:,iBand,iRef)'; MDSensory(:,iBand,iRef)'; SDSensory(:,iBand,iRef)']',{'S_M';'M_D';'S_D'});
            title(bandLabels{iBand}); ylim([-1 1]); box off;
        end
        sgtitle([refName ' Sensory - ' varName]);

        figure;
        for iBand = 1:5
            subplot(2,3,iBand);
            boxplot([SMMotor(:,iBand,iRef)'; MDMotor(:,iBand,iRef)'; SDMotor(:,iBand,iRef)']',{'S_M';'M_D';'S_D'});
            title( bandLabels{iBand}); ylim([-1 1]); box off;
        end
        sgtitle([ refName ' Motor - ' varName]);
    end
end

%% Get powers for LFP - clean up this and the next section....
% fs = 1e3;
% gammaBand = [30 90]; 
% alphaBand = [8 12];  
% betaBand  = [13 30]; 
% thetaBand = [6 8];
% params.Fs     = fs;
% params.fpass  = [1 120];
% params.pad    = -1;
% params.tapers = [3 5];
% 
% clear meanSpecAlpha meanSpecBeta meanSpecGamma
% meanEnvelopeAlpha = NaN(5,7,21);meanEnvelopeBeta = NaN(5,7,21);meanEnvelopeGamma = NaN(5,7,21);
% for iDate = 1: size(allDates,1)
%     clear expDate;
%     expDate = allDates(iDate,:);
% 
%     for iRun = 1:size(allRuns{iDate,1})
%         clear runName dataDir
%         runName = allRuns{iDate,1}(iRun,:);
%         dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
% 
%         clc; disp([monkeyName ' '  expDate ' run: ' runName]);
% 
%         % Obtaining ephys data
%         clear probeCh timeStamp badTimes badTimesThIm chInCortex szLFP...
%             szIm szMin badTimeIm badTimes10Hz
%         probeCh      = probe{iRun,iDate}.probeCh;
%         timeStamp    = probe{iRun,iDate}.timeStamp;
%         badTimes     = badTimesLFP{iDate,iRun};
%         chInCortex   = estChInCortex{1,iDate}(iRun,:);
% 
%         % Removing bad channels
%         probeCh(:,badCh{iDate,iRun}) = [];
%         probeCh(badTimes,:) = []; 
%         szLFP   = size(probeCh,1);
% 
%         if chInCortex(1)-chInCortex(2) ~= 0
%             clear probeDat envelopeDat
%             for iBand = 1:5
%                 switch iBand
%                     case 1
%                         [bCoeff,aCoeff] = butter(3,thetaBand./(fs/2),'bandpass'); 
%                     case 2
%                         [bCoeff,aCoeff] = butter(3,alphaBand./(fs/2),'bandpass');
%                     case 3
%                         [bCoeff,aCoeff] = butter(3,betaBand./(fs/2),'bandpass');
%                     case 4
%                         [bCoeff,aCoeff] = butter(3,gammaBand./(fs/2),'bandpass');
%                     case 5
%                         bCoeff = []; aCoeff = []; 
%                 end
%                 
%                 % Filter the LFP recorded from channels in cortex and rectify it
%                 if ~isempty(bCoeff) && ~isempty(aCoeff) % Rectify the signal
%                     probeDat =  abs(single(filtfilt(bCoeff,aCoeff,double(probeCh(:,chInCortex(1):chInCortex(2))))));
%                 else
%                     probeDat =  abs(probeDat(:,chInCortex));
%                 end
% 
%                 % Get the envelope of the signal
%                 envelopeDat = envelope(probeDat,5); %envelope(probeDat,5,'peak');
% 
%                 switch iBand
%                     case 1
%                          meanEnvelopeAlpha(iDate,iRun,:) = median(envelopeDat,1,'omitnan');
%                     case 2
%                          meanEnvelopeBeta(iDate,iRun,:) = median(envelopeDat,1,'omitnan');
%                     case 3
%                         meanEnvelopeGamma(iDate,iRun,:) = median(envelopeDat,1,'omitnan');
%                 end
%             end
%         else
%             meanEnvelopeAlpha(iDate,iRun,1:21) = 0;
%             meanEnvelopeBeta(iDate,iRun,1:21) = 0;
%             meanEnvelopeGamma(iDate,iRun,1:21) = 0;
%         end
%     end
% end
% 
% meanEnvelopeAlphaN = reshape(meanEnvelopeAlpha,[size(meanEnvelopeAlpha,1)*size(meanEnvelopeAlpha,2) size(meanEnvelopeAlpha,3)]);
% meanEnvelopeBetaN = reshape(meanEnvelopeBeta,[size(meanEnvelopeBeta,1)*size(meanEnvelopeBeta,2) size(meanEnvelopeBeta,3)]);
% meanEnvelopeGammaN = reshape(meanEnvelopeGamma,[size(meanEnvelopeGamma,1)*size(meanEnvelopeGamma,2) size(meanEnvelopeGamma,3)]);
% nanRow = (isnan(meanEnvelopeAlphaN(:,1)));
% 
% meanEnvelopeAlphaN(nanRow,:) = []; meanEnvelopeBetaN(nanRow,:) = []; meanEnvelopeGammaN(nanRow,:)  = []; 
% 
% motorEnvelopeGamma = meanEnvelopeGammaN; motorEnvelopeAlpha = meanEnvelopeAlphaN; motorEnvelopeBeta = meanEnvelopeBetaN; 
% motorEnvelopeGamma(~motorGoodSpatialRuns,:)= []; motorEnvelopeAlpha(~motorGoodSpatialRuns,:)= []; motorEnvelopeBeta(~motorGoodSpatialRuns,:)= [];
% 
% sensoryEnvelopeGamma = meanEnvelopeGammaN; sensoryEnvelopeAlpha = meanEnvelopeAlphaN; sensoryEnvelopeBeta = meanEnvelopeBetaN; 
% sensoryEnvelopeGamma(~sensoryGoodSpatialRuns,:) = []; sensoryEnvelopeAlpha(~sensoryGoodSpatialRuns,:)= []; sensoryEnvelopeBeta(~sensoryGoodSpatialRuns,:)= [];
% 
% meanEnvelopeAlphaN(~(~singleChFlag & goodRunsSpatial),:) = []; meanEnvelopeBetaN(~(~singleChFlag & goodRunsSpatial),:) = []; meanEnvelopeGammaN(~(~singleChFlag & goodRunsSpatial),:)  = []; 
% 
% % Plotting
% for iType = 1:3
%     switch iType
%         case 1
%             alphaSpec = meanEnvelopeAlphaN;
%             betaSpec  = meanEnvelopeBetaN;
%             gammaSpec = meanEnvelopeGammaN; 
%             specTitle = 'All';
%         case 2
%             alphaSpec = sensoryEnvelopeAlpha;
%             betaSpec   = sensoryEnvelopeBeta;
%             gammaSpec  = sensoryEnvelopeGamma;
%              specTitle = 'Sensory';
%         case 3
%             alphaSpec = motorEnvelopeAlpha;
%             betaSpec   = motorEnvelopeBeta;
%             gammaSpec  = motorEnvelopeGamma;
%              specTitle = 'Motor';
%     end
% 
%     figure;subplot(131); plot(median(alphaSpec,1),1:21); set(gca,'YDir','reverse');
%     semVal = (mad(alphaSpec,1,1))./sqrt(size(alphaSpec,1)); yVar = [(1:21) fliplr((1:21))];
%     med = median(alphaSpec,1); ylim([0 22]); title('Alpha band'); box off;
%     patch([med-2*semVal fliplr(med+2*semVal)],yVar,'b','FaceAlpha',0.3,'EdgeColor','none');
% 
%     subplot(132); plot(median(betaSpec,1),1:21);  set(gca,'YDir','reverse');
%     semVal = (mad(betaSpec,1,1))./sqrt(size(betaSpec,1)); yVar = [(1:21) fliplr((1:21))];
%     med = median(betaSpec,1);ylim([0 22]); title('Beta band');box off;
%     patch([med-2*semVal fliplr(med+2*semVal)],yVar,'b','FaceAlpha',0.3,'EdgeColor','none');
% 
%     subplot(133); plot(median(gammaSpec,1),1:21);  set(gca,'YDir','reverse');
%     semVal = (mad(gammaSpec,1,1))./sqrt(size(gammaSpec,1)); yVar = [(1:21) fliplr((1:21))];
%     med = median(gammaSpec,1);ylim([0 22]); title('Gamma band');box off;
%     patch([med-2*semVal fliplr(med+2*semVal)],yVar,'b','FaceAlpha',0.3,'EdgeColor','none');
%     
%     sgtitle(specTitle);
% 
%     figure; subplot(131); boxplot(alphaSpec); title('Alpha band'); box off;
%     subplot(132); boxplot(betaSpec);title('Beta band');box off;
%     subplot(133); boxplot(gammaSpec);title('Gamma band');box off;
%     sgtitle(specTitle);
% 
%     % Grouping it by compartment
%     superAlpha = median(alphaSpec(:,1:6),2);
%     superBeta  = median(betaSpec(:,1:6),2);
%     superGamma = median(gammaSpec(:,1:6),2);
% 
%     midAlpha = median(alphaSpec(:,7:12),2);
%     midBeta  = median(betaSpec(:,7:12),2);
%     midGamma = median(gammaSpec(:,7:12),2);
% 
%     deepAlpha = median(alphaSpec(:,13:end),2);
%     deepBeta  = median(betaSpec(:,13:end),2);
%     deepGamma = median(gammaSpec(:,13:end),2);
% 
%     figure; subplot(131); boxplot([superAlpha superBeta superGamma]);
%     subplot(132); boxplot([midAlpha midBeta midGamma]);
%     subplot(133); boxplot([deepAlpha deepBeta deepGamma]);sgtitle(specTitle);
% 
%     figure; subplot(131); boxplot([superAlpha midAlpha deepAlpha],{'Super';'Mid';'Deep'}); title('Alpha band')
%     subplot(132); boxplot([superBeta midBeta deepBeta],{'Super';'Mid';'Deep'}); title('Beta band')
%     subplot(133); boxplot([superGamma midGamma deepGamma],{'Super';'Mid';'Deep'}); title('Gamma band');sgtitle(specTitle);
% 
% end
% 
% [pSpCorr,tblSpCorr,statsSpCorr] = anova1([superAlpha midAlpha deepAlpha],{'Superficial' ; 'Middle'; 'Deep'},'off');
% [rSpCorr,mSpCorr,~,gnamesSpCorr] = multcompare(statsSpCorr,"Alpha",0.01,"CriticalValueType","bonferroni");
% 

%% Get powers for LFP
fs = 1e3;
gammaBand = [30 90]; % Gamma band filtering parameters
alphaBand = [8 12];
betaBand  = [13 30];
thetaBand = [6 8];

params.Fs     = fs;
params.fpass  = [1 120];
params.pad    = -1;
params.tapers = [3 5];

paramsRaw       = params;
paramsRaw.fpass = [250 500];

bandLabels = {'Theta'; 'Alpha'; 'Beta'; 'Gamma'; 'Spiking'};

if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\specValsAll.mat'],'file')
    meanSpecVals = NaN(5,7,5,21);
    for iDate = 1: size(allDates,1)
        clear expDate;
        expDate = allDates(iDate,:);

        for iRun = 1:size(allRuns{iDate,1})
            clear runName dataDir
            runName = allRuns{iDate,1}(iRun,:);
            dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

            clc; disp([monkeyName ' '  expDate ' run: ' runName]);


            % Obtaining ephys data
            clear probeCh timeStamp badTimes badTimesThIm chInCortex szLFP...
                szIm szMin badTimeIm badTimes10Hz
            probeCh    = probe{iRun,iDate}.probeCh;
            rawCh      = probe{iRun,iDate}.rawCh;
            timeStamp  = probe{iRun,iDate}.timeStamp;
            badTimes   = badTimesLFP{iDate,iRun};
            chInCortex = estChInCortex{1,iDate}(iRun,:);

            % Removing bad channels
            probeCh(:,badCh{iDate,iRun}) = [];
            rawCh(:,badCh{iDate,iRun})   = [];

            % Remove bad times
            probeCh(badTimes,:) = [];
            rawCh(badTimes,:)   = [];
            szLFP   = size(probeCh,1);

            if chInCortex(1)-chInCortex(2) ~= 0
                clear probeChCortex spec aLim bLim gLim spec timeValsSpec freqValsSpec
                probeChCortex = probeCh(:,chInCortex(1):chInCortex(2));
                rawChCortex   = rawCh(:,chInCortex(1):chInCortex(2));
                [spec,timeValsSpec,freqValsSpec] = mtspecgramc(probeChCortex,[5 2],params);
                [specRaw,timeValsSpecR,freqValsSpecR] = mtspecgramc(rawChCortex,[5 2],paramsRaw);

                for iBand = 1:5
                    switch iBand
                        case 1
                            freqLim = freqValsSpec>= thetaBand(1) & freqValsSpec<=thetaBand(2);
                            specVal = spec;
                        case 2
                            freqLim = freqValsSpec>=alphaBand(1) & freqValsSpec<=alphaBand(2);
                            specVal = spec;
                        case 3
                            freqLim = freqValsSpec>=betaBand(1) & freqValsSpec<=betaBand(2);
                            specVal = spec;
                        case 4
                            freqLim = freqValsSpec>=gammaBand(1) & freqValsSpec<= gammaBand(2);
                            specVal = spec;
                        case 5
                            freqLim = freqValsSpecR>=250 & freqValsSpecR<= 500;
                            specVal = specRaw;
                    end

                    meanSpecVals(iDate,iRun,iBand,1:size(specVal,3)) = (squeeze(mean(specVal(:,freqLim,:),[1,2],'omitnan')));
                    meanSpecVals(iDate,iRun,iBand,1:size(specVal,3)) = (meanSpecVals(iDate,iRun,iBand,1:size(specVal,3)))./max(meanSpecVals(iDate,iRun,iBand,1:size(specVal,3)));
                end
                %             if ~exist([dataDir '\LFP_Laminar_Powers.png'],'file')
                %                 figure; plot(squeeze(meanSpecAlpha(iDate,iRun,:)),1:size(probeChCortex,2)); hold on; box off;
                %                 plot(squeeze(meanSpecBeta(iDate,iRun,:)),1:size(probeChCortex,2));
                %                 plot(squeeze(meanSpecGamma(iDate,iRun,:)),1:size(probeChCortex,2));
                %                 xlim([0 1]); ylim([0 size(probeChCortex,2)+1]); yticks((1:size(probeChCortex)));
                %                 set(gca,'YDir','reverse'); legend('Alpha','Beta','Gamma','Location','northeast');
                %                 sgtitle(strrep(['Average power (dB) vs channels for ' monkeyName ' Date: ' expDate ...
                %                     ' - Run: ' runName(end)],'_','\_'));
                %                 f = gcf; exportgraphics(f,[dataDir '\LFP_Laminar_Powers.png'],'Resolution',300); close gcf;
                %             end
            else
                meanSpecVals(iDate,iRun,1:5,1:21) = 0;
            end
        end
    end

    meanSpecValsN = reshape(meanSpecVals,[size(meanSpecVals,1)*size(meanSpecVals,2) size(meanSpecVals,3) size(meanSpecVals,4)]);
    nanRow = isnan(meanSpecValsN(:,1,1));
    meanSpecValsN(nanRow,:,:) = [];

    motorSpecVals   = meanSpecValsN; motorSpecVals(~(~singleChFlag & motorGoodSpatialRuns),:,:)     = [];
    sensorySpecVals = meanSpecValsN; sensorySpecVals(~(~singleChFlag & sensoryGoodSpatialRuns),:,:) = [];

    meanSpecValsN(~(~singleChFlag & goodRunsSpatial),:,:) = [];
    save(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\specValsAll.mat'],'meanSpecValsN','motorSpecVals','sensorySpecVals');

else
    monkeyVars = load(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\specValsAll.mat']);
    meanSpecValsN = monkeyVars.meanSpecValsN;
    motorSpecVals = monkeyVars.motorSpecVals;
    sensorySpecVals = monkeyVars.sensorySpecVals;
end


%% Plotting
for iType = 1:3
    switch iType
        case 1
            specType  = meanSpecValsN;
            specTitle = 'All';
        case 2
            specType  = sensorySpecVals;
            specTitle = 'Sensory';
        case 3
            specType  = motorSpecVals;
            specTitle = 'Motor';
    end

    figure;
    for iBand = 1:5
        clear bandSpec semVal med
        bandSpec = squeeze(specType(:,iBand,:));
        subplot(1,5,iBand); plot(median(bandSpec,1),1:21); xlim([0.1 1]); set(gca,'YDir','reverse');
        semVal = (mad(bandSpec,1,1))./sqrt(size(bandSpec,1)); yVar = [(1:21) fliplr((1:21))];
        med = median(bandSpec,1); ylim([0 22]); title(bandLabels{iBand});  box off;
        patch([med-2*semVal fliplr(med+2*semVal)],yVar,'b','FaceAlpha',0.3,'EdgeColor','none');
    end
    sgtitle(specTitle);

    % Grouping it by compartment
    figure;
    superLayer = squeeze(median(specType(:,:,1:6),3,'omitnan'));
    midLayer   = squeeze(median(specType(:,:,7:12),3,'omitnan'));
    deepLayer  = squeeze(median(specType(:,:,13:end),3,'omitnan'));

    for iBand = 1:5
        subplot(2,3,iBand);
        boxplot([superLayer(:,iBand) midLayer(:,iBand) deepLayer(:,iBand)],{'Super';'Mid';'Deep'});
        ylim([0 1]);title(bandLabels{iBand}); box off;
    end
    sgtitle(specTitle);
end

%% Combined plotting for both animals
clear monkeyVars
monkeyVars = load(['X:\Data\CharlieSheen_SqM\' hemisphere ' Hemisphere\specValsAll.mat']);
meanSpecCh = monkeyVars.meanSpecValsN; 
motorSpecCh = monkeyVars.motorSpecVals;
sensorySpecCh = monkeyVars.sensorySpecVals; 
for iType = 1:3
    switch iType
        case 1
            specType  = [meanSpecCh;meanSpecValsN];
            specTitle = 'All';
        case 2
            specType  = [sensorySpecCh;sensorySpecVals];
            specTitle = 'Sensory';
        case 3
            specType  = [motorSpecCh;motorSpecVals];
            specTitle = 'Motor';
    end

    figure;
    for iBand = 1:5
        clear bandSpec semVal med
        bandSpec = squeeze(specType(:,iBand,:));
        subplot(1,5,iBand); plot(median(bandSpec,1),1:21); xlim([0.1 1]); set(gca,'YDir','reverse');
        semVal = (mad(bandSpec,1,1))./sqrt(size(bandSpec,1)); yVar = [(1:21) fliplr((1:21))];
        med = median(bandSpec,1); ylim([0 22]); title(bandLabels{iBand});  box off;
        patch([med-2*semVal fliplr(med+2*semVal)],yVar,'b','FaceAlpha',0.3,'EdgeColor','none');
    end
    sgtitle(specTitle);

    % Grouping it by compartment
    figure;
    superLayer = squeeze(median(specType(:,:,1:6),3,'omitnan'));
    midLayer   = squeeze(median(specType(:,:,7:12),3,'omitnan'));
    deepLayer  = squeeze(median(specType(:,:,13:end),3,'omitnan'));

    for iBand = 1:5
        subplot(2,3,iBand);
        boxplot([superLayer(:,iBand) midLayer(:,iBand) deepLayer(:,iBand)],{'Super';'Mid';'Deep'});
        ylim([0 1]);title(bandLabels{iBand}); box off;
    end
    sgtitle(specTitle);
end
%%
bandVals = [2 4 5]; colorVals = ['r';'b';'y'];
figure;
for iBand = 1:3
    clear bandSpec semVal med
    bandSpec = squeeze(specType(:,bandVals(iBand),:));
    plot(median(bandSpec,1),1:21,colorVals(iBand)); xlim([0.1 1]); set(gca,'YDir','reverse');
    semVal = (mad(bandSpec,1,1))./sqrt(size(bandSpec,1)); yVar = [(1:21) fliplr((1:21))];
    med = median(bandSpec,1); ylim([0 22]);   box off;
    patch([med-2*semVal fliplr(med+2*semVal)],yVar,colorVals(iBand),'FaceAlpha',0.3,'EdgeColor','none'); hold on;
end


%% Cross correlating physiology and imaging at ROI 
clear tempProfileNoRef tempProfileSuperNoRef tempProfileDeepNoRef tempProfileMidNoRef...
    tempProfile_6Ch_NoRef tempProfileSuper_6Ch_NoRef tempProfileDeep_6Ch_NoRef tempProfileMid_6Ch_NoRef
% tempProfileNoRef           = NaN(iM,size(probe,2), size(probe,1));
% tempProfileSuperNoRef      = NaN(iM,size(probe,2), size(probe,1));
% tempProfileDeepNoRef       = NaN(iM,size(probe,2), size(probe,1));
% tempProfile_6Ch_NoRef      = NaN(iM,size(probe,2), size(probe,1));
% tempProfileSuper_6Ch_NoRef = NaN(iM,size(probe,2), size(probe,1));
% tempProfileDeep_6Ch_NoRef  = NaN(iM,size(probe,2), size(probe,1));

for iDate = 1:size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);
    for iRun = 1:size(allRuns{iDate,1})
        clear runName dataDir clipMask elecMask clipMaskCortex corrMask x negIdx lowIdx
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

        clc; disp(['Analyzing data for ' monkeyName ' '  expDate ' run: ' runName]);

        if ~exist([dataDir '\ROIAllVars.mat'],'file')
            % IMAGING: Load the appropriate masks for the imaging data
            % Load clipmask
            if exist([dataDir '\clipMask0' runName(end) '.BMP'],'file') == 0
                clipMask = imread([dataDir '\clipMask0' runName(end) '.png']);
            else
                clipMask = imread([dataDir '\clipMask0' runName(end) '.bmp']);
            end

            % This mask includes vessels and electrode(s)
            if exist([dataDir '\maskSkull0' runName(end) '.bmp'],'file') == 0
                elecMask = imread([dataDir '\maskSkull0' runName(end) '.png']);
            else
                elecMask = imread([dataDir '\maskSkull0' runName(end) '.bmp']);
            end

            clipMask       = imresize(clipMask,1/3); % Resize mask
            clipMask       = clipMask(:,:,1)>0; % Converting to 0s and 1s
            elecMask       = imresize(elecMask-100,1/3); % Binarize image
            elecMask       = elecMask(:,:,1)>0;
            clipMaskCortex = clipMask & ~elecMask; % Mask of cortex sans vessels

            % Mask for comparing the imaging and hybrid map (after removing
            % areas beyond lateral sulcus)
            if exist([dataDir '\corrMask0' runName(end) '.BMP'],'file')
                corrMask = imread([dataDir '\corrMask0' runName(end) '.bmp']);
                corrMaskFlag = 1;
            elseif exist([dataDir '\corrMask0' runName(end) '.PNG'],'file')
                corrMask = imread([dataDir '\corrMask0' runName(end) '.png']);
                corrMaskFlag = 1;
            else
                corrMaskFlag = 0;
            end

            if corrMaskFlag
                corrMask = imresize(corrMask,1/3); % Resize the mask
                corrMask = corrMask(:,:,1)>0;
                corrMask = corrMask & ~elecMask;
            else
                corrMask = clipMaskCortex;
            end

            % Get the cross correlations for the ROI and FOV
            % 1. Get the temporal profile for all recordings...
            % Pick the ROI (1mm x 1mm) - 1mm diameter/ 500um radius;
            disp('Getting temporal profile...');
            seedRad = round(roiSize{iDate}(iRun)*2./spatialBin);
            imSize = size(imresize(greenIm{iDate,iRun},1/spatialBin));
            greenFig = imresize(greenIm{iDate,iRun},1/spatialBin,'OutputSize',[imSize(1) imSize(2)]);

            % Check if ROI center and location of probe are saved
            clear varInfo;
            if exist(fullfile(dataDir,'roiCenterLoc.mat'),'file')
                varInfo = who('-file', fullfile(dataDir,'roiCenterLoc.mat'));
            else
                varInfo =[];
            end
            if find(ismember(varInfo,'seedLocProbe')); probeLocFlag = 1; else; probeLocFlag = 0; end

            % Save ROI center and probe location if not saved already
            if ~exist(fullfile(dataDir,'roiCenterLoc.mat'),'file') || ~probeLocFlag
                if ~exist(fullfile(dataDir,'roiCenterLoc.mat'),'file')
                    disp('Picking ROI near the electrode...');
                    figure; imagesc(greenFig); colormap gray; axis image off;
                    title('Pick a seed to get the ROI');
                    seedLocIn = ginput(1); seedLocIn = (round(seedLocIn)); close gcf;
                else
                    seedLocIn = load([dataDir '\roiCenterLoc.mat'],'seedLocIn');
                    seedLocIn = seedLocIn.seedLocIn;
                end

                if ~probeLocFlag % Check if probe location has been saved
                    figure; imagesc(greenFig); colormap gray; axis image off; hold on;
                    title ('Pick the location of the probe...')
                    seedLocProbe = ginput(1); seedLocProbe = (round(seedLocProbe)); close gcf;
                end

                save([dataDir '\roiCenterLoc.mat'],'seedLocProbe','seedLocIn');

            else
                disp('Seed for ROI already picked...');

                seedLocIn = load([dataDir '\roiCenterLoc.mat'],'seedLocIn');
                seedLocIn = seedLocIn.seedLocIn;

                seedLocProbe = load([dataDir '\roiCenterLoc.mat'],'seedLocProbe');
                seedLocProbe = seedLocProbe.seedLocProbe;
            end

            % Show the ROI on the green map
            if ~exist(fullfile(dataDir,'ROI.png'),'file')
                figure; imagesc(greenFig); axis image off; colormap gray;
                rectangle('Position',[seedLocIn(1)-round(seedRad/2),seedLocIn(2)-round(seedRad/2),...
                    seedRad,seedRad],'EdgeColor','r','LineWidth',2);
                title('ROI near the electrode');
                f = gcf; exportgraphics(f,[dataDir '\ROI.png'],'Resolution',300); close gcf;
            end

            if ~exist(fullfile(dataDir,'FCMap_ROI.png'),'file')
                clear pDatTemp;
                circleRad = round(roiSize{1}(1)/(spatialBin));
                pDatTemp  = processedDat{iDate,iRun}.tempBandPass;
                seedSigT  = calculateSeedSignal(greenFig,corrMask,...
                    seedLocIn,circleRad,pDatTemp); % Get Gaussian weighted seed signal
                corrMapT   = plotCorrMap(seedSigT,pDatTemp,0);
                greenImRGB = ind2rgb(greenFig,gray(256));
                figure; imagesc(greenImRGB);axis image; colormap jet; axis image off;
                hold on; imagesc(corrMapT,'AlphaData',corrMapT.*0.8);clim([0 1]);colorbar;
                f = gcf; exportgraphics(f,[dataDir '\FCMap_ROI.png'],'Resolution',300); close gcf;
            end

            % Get the ROI for gamma cross correlations
            clipMaskROI = clipMaskCortex(seedLocIn(2)-round(seedRad/2):seedLocIn(2)+round(seedRad/2),...
                seedLocIn(1)-round(seedRad/2):seedLocIn(1)+round(seedRad/2));

            clear tempProfileNoRefRun tempProfileSuperNoRefRun tempProfileDeepNoRefRun tempProfileMidNoRefRun...
                tempProfile_6Ch_NoRefRun tempProfileSuper_6Ch_NoRefRun tempProfileDeep_6Ch_NoRefRun tempProfileMid_6Ch_NoRefRun

            [tempProfileNoRefRun,tempProfileSuperNoRefRun,... % 10-channel split No-Ref
                tempProfileDeepNoRefRun,tempProfileMidNoRefRun,~] = ...
                getCrossCorrROI(dataDir,monkeyName,expDate,runName,...
                processedDat{iDate,iRun}.tempBandPass,probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,...
                badTimesLFP{iDate,iRun},badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),seedLocIn,...
                seedRad,clipMaskROI,10,'NoRef');

            [tempProfile_6Ch_NoRefRun,tempProfileSuper_6Ch_NoRefRun,... % 6-channel split No-Ref
                tempProfileDeep_6Ch_NoRefRun,tempProfileMid_6Ch_NoRefRun,~] = ...
                getCrossCorrROI(dataDir,monkeyName,expDate,runName,...
                processedDat{iDate,iRun}.tempBandPass,probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,...
                badTimesLFP{iDate,iRun},badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),seedLocIn,...
                seedRad,clipMaskROI,6,'NoRef');

            save([dataDir '\ROIAllVars.mat'],'tempProfileNoRefRun','tempProfileSuperNoRefRun',...
                'tempProfileDeepNoRefRun','tempProfile_6Ch_NoRefRun','tempProfileDeep_6Ch_NoRefRun',...
                'tempProfileSuper_6Ch_NoRefRun','tempProfileMidNoRefRun','tempProfileMid_6Ch_NoRefRun');

            tempProfileNoRef(iDate,iRun)           = tempProfileNoRefRun; %#ok<*SAGROW> 
            tempProfileSuperNoRef(iDate,iRun)      = tempProfileSuperNoRefRun;
            tempProfileMidNoRef(iDate,iRun)        = tempProfileMidNoRefRun;
            tempProfileDeepNoRef(iDate,iRun)       = tempProfileDeepNoRefRun;
            tempProfile_6Ch_NoRef(iDate,iRun)      = tempProfile_6Ch_NoRefRun;
            tempProfileSuper_6Ch_NoRef(iDate,iRun) = tempProfileSuper_6Ch_NoRefRun;
            tempProfileMid_6Ch_NoRef(iDate,iRun)   = tempProfileMid_6Ch_NoRefRun;
            tempProfileDeep_6Ch_NoRef(iDate,iRun)  = tempProfileDeep_6Ch_NoRefRun;

        else
            clear allVars
            allVars                                 = load([dataDir '\ROIAllVars.mat']);
            tempProfileNoRef(iDate,iRun)            = allVars.tempProfileNoRefRun;
            tempProfileSuperNoRef(iDate,iRun)       = allVars.tempProfileSuperNoRefRun;
            tempProfileMidNoRef(iDate,iRun)         = allVars.tempProfileMidNoRefRun;
            tempProfileDeepNoRef(iDate,iRun)        = allVars.tempProfileDeepNoRefRun;
            tempProfile_6Ch_NoRef(iDate,iRun)       = allVars.tempProfile_6Ch_NoRefRun;
            tempProfileSuper_6Ch_NoRef(iDate,iRun)  = allVars.tempProfileSuper_6Ch_NoRefRun;
            tempProfileMid_6Ch_NoRef(iDate,iRun)    = allVars.tempProfileMid_6Ch_NoRefRun;
            tempProfileDeep_6Ch_NoRef(iDate,iRun)   = allVars.tempProfileDeep_6Ch_NoRefRun;          
        end
    end
end


%% Compile the ROI level variables 
clear bandLabels tempProfilesAll medTempProfileAll allChCorr allChLagVal...
    allChCorrSuper allChCorrMid allChCorrDeep
bandLabels = {'Theta'; 'Alpha'; 'Beta'; 'Gamma'; 'Spiking'};

% Temporal profiles for all frequency bands
tempProfilesAll(:,:,1)  = [tempProfileNoRef.profileTheta]; 
tempProfilesAll(:,:,2)  = [tempProfileNoRef.profileAlpha]; 
tempProfilesAll(:,:,3)  = [tempProfileNoRef.profileBeta]; 
tempProfilesAll(:,:,4)  = [tempProfileNoRef.profile];      
tempProfilesAll(:,:,5)  = [tempProfileNoRef.profileRaw];   

tempProfilesAll(:,~goodRuns,:) = [];

% Show the distribution of lags and peak negative correlations 
allChCorr(:,1) = [tempProfileNoRef.magLowTheta]'; 
allChCorr(:,2) = [tempProfileNoRef.magLowAlpha]'; 
allChCorr(:,3) = [tempProfileNoRef.magLowBeta]'; 
allChCorr(:,4) = [tempProfileNoRef.magLow]'; 
allChCorr(:,5) = [tempProfileNoRef.magLowRaw]'; 
allChCorr(~goodRuns,:) = []; 

allChLagVal(:,1)= [tempProfileNoRef.lagLowTheta]'; 
allChLagVal(:,2)= [tempProfileNoRef.lagLowAlpha]'; 
allChLagVal(:,3)= [tempProfileNoRef.lagLowBeta]'; 
allChLagVal(:,4) = [tempProfileNoRef.lagLow]'; 
allChLagVal(:,5) = [tempProfileNoRef.lagLowRaw]'; 
allChLagVal(~goodRuns,:) = []; 

% Show the distribution of lags and peak negative correlations - Superficial
allChCorrSuper(:,1) = [tempProfileSuper_6Ch_NoRef.magLowTheta]'; 
allChCorrSuper(:,2) = [tempProfileSuper_6Ch_NoRef.magLowAlpha]'; 
allChCorrSuper(:,3) = [tempProfileSuper_6Ch_NoRef.magLowBeta]'; 
allChCorrSuper(:,4) = [tempProfileSuper_6Ch_NoRef.magLow]'; 
allChCorrSuper(:,5) = [tempProfileSuper_6Ch_NoRef.magLowRaw]'; 
allChCorrSuper(~goodRuns,:) = []; 

% Show the distribution of lags and peak negative correlations - Middle
allChCorrMid(:,1) = [tempProfileMid_6Ch_NoRef.magLowTheta]'; 
allChCorrMid(:,2) = [tempProfileMid_6Ch_NoRef.magLowAlpha]'; 
allChCorrMid(:,3) = [tempProfileMid_6Ch_NoRef.magLowBeta]'; 
allChCorrMid(:,4) = [tempProfileMid_6Ch_NoRef.magLow]'; 
allChCorrMid(:,5) = [tempProfileMid_6Ch_NoRef.magLowRaw]'; 
allChCorrMid(~goodRuns,:) = []; 

% Show the distribution of lags and peak negative correlations - Deep
allChCorrDeep(:,1) = [tempProfileDeep_6Ch_NoRef.magLowTheta]'; 
allChCorrDeep(:,2) = [tempProfileDeep_6Ch_NoRef.magLowAlpha]'; 
allChCorrDeep(:,3) = [tempProfileDeep_6Ch_NoRef.magLowBeta]'; 
allChCorrDeep(:,4) = [tempProfileDeep_6Ch_NoRef.magLow]'; 
allChCorrDeep(:,5) = [tempProfileDeep_6Ch_NoRef.magLowRaw]'; 
allChCorrDeep(~goodRuns,:) = []; 


% Show the distributions of lags and correlations for 10/10 and 6/6 split
% for superficial and deep channels
superCorr_10Ch_Gamma = [tempProfileSuperNoRef.magLow]'; 
deepCorr_10Ch_Gamma  = [tempProfileDeepNoRef.magLow]';
superCorr_6Ch_Gamma  = [tempProfileSuper_6Ch_NoRef.magLow]'; 
deepCorr_6Ch_Gamma   = [tempProfileDeep_6Ch_NoRef.magLow]'; 

superCorr_10Ch_Gamma(~(goodRuns & ~singleChFlag)) = []; 
deepCorr_10Ch_Gamma(~(goodRuns & ~singleChFlag))  = []; 
superCorr_6Ch_Gamma(~(goodRuns & ~singleChFlag))  = []; 
deepCorr_6Ch_Gamma(~(goodRuns & ~singleChFlag))  = [];  

superLag_10Ch_Gamma = [tempProfileSuperNoRef.lagLow]'./10; 
deepLag_10Ch_Gamma  = [tempProfileDeepNoRef.lagLow]'./10;
superLag_6Ch_Gamma  = [tempProfileSuper_6Ch_NoRef.lagLow]'./10; 
deepLag_6Ch_Gamma   = [tempProfileDeep_6Ch_NoRef.lagLow]'./10;

superLag_10Ch_Gamma(~(goodRuns & ~singleChFlag)) = []; 
deepLag_10Ch_Gamma(~(goodRuns & ~singleChFlag))  = []; 
superLag_6Ch_Gamma(~(goodRuns & ~singleChFlag))  = []; 
deepLag_6Ch_Gamma(~(goodRuns & ~singleChFlag))   = [];  

% Get the median +/- mad lag 
x = -200:200;
allChLag = [tempProfileNoRef.lagLow]'; allChLag(~(goodRuns & ~singleChFlag)) = []; 
lagFrameRange = find(x == round(median(allChLag)- mad(allChLag))) : find(x == round(median(allChLag)+ mad(allChLag)));
%169:189; % Median +/- MAD across two monkeys

smFlagFOV = smFlag; smFlagFOV(~(goodRunsSpatial & ~singleChFlag)) = [];
smFlagROI = smFlag; smFlagROI(~(goodRuns)) = [];

if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat'],'file') 
    save(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat'],...
        'smFlagFOV','tempProfilesAll','allChCorr','allChLagVal','superCorr_10Ch_Gamma','deepCorr_10Ch_Gamma',...
        'superCorr_6Ch_Gamma','deepCorr_6Ch_Gamma','superLag_10Ch_Gamma','deepLag_10Ch_Gamma',...
        'superLag_6Ch_Gamma','deepLag_6Ch_Gamma','lagFrameRange','allChLag','allChCorrSuper',...
        'allChCorrMid','allChCorrDeep','smFlagROI','-append');
end

%%  Cross correlating physiology and imaging for FOV
videoFlag      = 0; % To save videos of hybrids - change if you need to
corrFCHybrid   = NaN(size(processedDat,1),size(processedDat,2),5,5,401);
corrFCSuper    = NaN(size(processedDat,1),size(processedDat,2),5,5,401);
corrFCDeep     = NaN(size(processedDat,1),size(processedDat,2),5,5,401);
super_DeepCorr = NaN(size(processedDat,1),size(processedDat,2),5,5,401);
corrFCMid      = NaN(size(processedDat,1),size(processedDat,2),5,5,401);
super_MidCorr  = NaN(size(processedDat,1),size(processedDat,2),5,5,401);
mid_DeepCorr   = NaN(size(processedDat,1),size(processedDat,2),5,5,401);

super_DeepAvgFrames = NaN(size(processedDat,1),size(processedDat,2),5,5);
super_MidAvgFrames  = NaN(size(processedDat,1),size(processedDat,2),5,5);
deep_MidAvgFrames   = NaN(size(processedDat,1),size(processedDat,2),5,5);
peakNegValsAll      = NaN(size(processedDat,1),size(processedDat,2),5,5);
peakNegTimesAll     = NaN(size(processedDat,1),size(processedDat,2),5,5);

superHybridAllBands = NaN(size(processedDat,1),size(processedDat,2),5,5,5);
midHybridAllBands   = NaN(size(processedDat,1),size(processedDat,2),5,5,5);
deepHybridAllBands  = NaN(size(processedDat,1),size(processedDat,2),5,5,5);
bandNames           = {'Theta'; 'Alpha';'Beta';'Gamma';'Spiking'};

tic;
for iDate = 1:size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);
    for iRun = 1: size(allRuns{iDate,1},1)       
        clear runName dataDir clipMask elecMask clipMaskCortex corrMask...
            x negIdx lowIdx serverDir
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName];
        serverDir = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\kem294\Data\' ...
            monkeyName '_SqM\' hemisphere ' Hemisphere\'  expDate '\' runName ];

        clc; disp(['Analyzing data for ' monkeyName ' '  expDate ' run: ' runName]);

        if ~exist([serverDir '\crossCorrFOV_10_NoRef.mat'],'file') 
            disp('No reference - channel split 10(superficial)/10(Deep)');
            getCrossCorrFOV(monkeyName,expDate,runName,serverDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,10,'NoRef');
        end

        if ~exist([serverDir '\crossCorrFOV_6_NoRef.mat'],'file') 
            disp('No reference - channel split top 6(superficial)/bottom 6(Deep)');
            getCrossCorrFOV(monkeyName,expDate,runName,serverDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,6,'NoRef');
        end

        if ~exist([serverDir '\crossCorrFOV_6_BipolarRef.mat'],'file')
            disp('Bipolar reference -  channel split top 6(superficial)/bottom 6(Deep)');
            getCrossCorrFOV(monkeyName,expDate,runName,serverDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,6,'BipolarRef');
        end

        if ~exist([serverDir '\crossCorrFOV_6_AvgRef.mat'],'file')
            disp('Avg reference - channel split top 6(superficial)/bottom 6(Deep)');
            getCrossCorrFOV(monkeyName,expDate,runName,serverDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,6,'AvgRef');
        end

        if ~exist([serverDir '\crossCorrFOV_6_CSDRef.mat'],'file') 
            disp('CSD reference - channel split top 6(superficial)/bottom 6(Deep)');
            getCrossCorrFOV(monkeyName,expDate,runName,serverDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,6,'CSDRef');
        end
      %
        clear varInfo;
        varFlag = 0; 
        try matfile([dataDir '\processedHybridMapVars.mat']);
            varInfo = who('-file',[dataDir '\processedHybridMapVars.mat']);
            if sum(ismember(varInfo,'crossFreqCrossLayerHybridT'))==0
                varFlag = 1;
            end
        catch
            varFlag = 1;
        end 

        if varFlag
            % IMAGING: Load the appropriate masks for the imaging data
            % Load clipmask
            if exist([dataDir '\clipMask0' runName(end) '.BMP'],'file') == 0
                clipMask = imread([dataDir '\clipMask0' runName(end) '.png']);
            else
                clipMask = imread([dataDir '\clipMask0' runName(end) '.bmp']);
            end

            % This mask includes vessels and electrode(s)
            if exist([dataDir '\maskSkull0' runName(end) '.bmp'],'file') == 0
                elecMask = imread([dataDir '\maskSkull0' runName(end) '.png']);
            else
                elecMask = imread([dataDir '\maskSkull0' runName(end) '.bmp']);
            end

            clipMask       = imresize(clipMask,1/3); % Resize mask
            clipMask       = clipMask(:,:,1)>0; % Converting to 0s and 1s
            elecMask       = imresize(elecMask-100,1/3); % Binarize image
            elecMask       = elecMask(:,:,1)>0;
            clipMaskCortex = clipMask & ~elecMask; % Mask of cortex sans vessels
            imSize         = size(clipMask);

            % Mask for comparing the imaging and hybrid map (after removing
            % areas beyond lateral sulcus)
            if exist([dataDir '\corrMask0' runName(end) '.BMP'],'file')
                corrMask = imread([dataDir '\corrMask0' runName(end) '.bmp']);
                corrMaskFlag = 1;
            elseif exist([dataDir '\corrMask0' runName(end) '.PNG'],'file')
                corrMask = imread([dataDir '\corrMask0' runName(end) '.png']);
                corrMaskFlag = 1;
            else
                corrMaskFlag = 0;
            end

            if corrMaskFlag
                corrMask = imresize(corrMask,1/3); % Resize the mask
                corrMask = corrMask(:,:,1)>0;
                corrMask = corrMask & ~elecMask;
            else
                corrMask = clipMaskCortex;
            end

            corrMaskT = reshape(corrMask,[imSize(1)*imSize(2) 1]);

            % Get ROI location and FC map
            seedLocIn    = load([dataDir '\roiCenterLoc.mat']);
            seedLocProbe = seedLocIn.seedLocProbe;
            seedLocIn    = seedLocIn.seedLocIn;
            circleRad    = round(roiSize{iDate}(iRun)/(spatialBin)); % 500um radius
            greenFig     = imresize(greenIm{iDate,iRun},1/spatialBin,'OutputSize',[imSize(1) imSize(2)]);

            pDatTemp = processedDat{iDate,iRun}.tempBandPass;
            seedSigT = calculateSeedSignal(greenFig,corrMask,...
                seedLocIn,circleRad,pDatTemp); % Get Gaussian weighted seed signal

            fcMap             = plotCorrMap(seedSigT,pDatTemp,0);
            fcMap             = reshape(fcMap,[361*438 1]);
            fcMap(~corrMaskT) = NaN;

            chLen = estChInCortex{1,iDate}(iRun,2) - estChInCortex{1,iDate}(iRun,1);

            % Initialize run-related variables
            % Correlations within the same frequencies            
            corrFCHybridT = NaN(5,4,401); corrFCSuperT    = NaN(5,4,401);
            corrFCDeepT   = NaN(5,4,401); super_DeepCorrT = NaN(5,4,401);

            corrFCMidT          = NaN(5,5,401); super_MidCorrT       = NaN(5,5,401);
            mid_DeepCorrT       = NaN(5,5,401); peakNegValsAllT      = NaN(5,5);
            peakNegTimesAllT    = NaN(5,5);     super_DeepAvgFramesT = NaN(5,5);
            super_MidAvgFramesT = NaN(5,5);     deep_MidAvgFramesT   = NaN(5,5);
           
            % Correlations between frequencies    
            superHybridAllBandsT = NaN(5,5,5);  deepHybridAllBandsT = NaN(5,5,5);
            midHybridAllBandsT   = NaN(5,5,5); crossFreqCrossLayerHybridT = NaN(5,5,5,3,3);

            for iType = 1:5 % 10/10 split or 6/6 split of superficial/deep channels
                clear crossCorrFOV allXCorr superXCorr deepXCorr allLags fileName
                switch iType
                    case 1
                        crossCorrFOV = matfile([serverDir '\crossCorrFOV_10_NoRef.mat']);
                        fileName     = '10_NoRef';

                    case 2
                        crossCorrFOV = matfile([serverDir '\crossCorrFOV_6_NoRef.mat']);
                        fileName     = '6_NoRef';

                    case 3
                        crossCorrFOV = matfile([serverDir '\crossCorrFOV_6_BipolarRef.mat']);
                        fileName     = '6_BipolarRef';

                    case 4
                        crossCorrFOV = matfile([serverDir '\crossCorrFOV_6_AvgRef.mat']);
                        fileName     = '6_AvgRef';

                    case 5
                        crossCorrFOV = matfile([serverDir '\crossCorrFOV_6_CSDRef.mat']);
                        fileName     = '6_CSDRef';
                end

                allXcorr     = crossCorrFOV.spatialProfile;
                superXcorr   = crossCorrFOV.spatialProfileSuper;
                deepXcorr    = crossCorrFOV.spatialProfileDeep;

                if chLen~=0 && iType~=1
                    midXCorr     = crossCorrFOV.spatialProfileMid;
                end
                allLags      = crossCorrFOV.lagFull;

                x = allLags;
                negIdx = x<0 & x>=-150; negVals = x(negIdx);
                lowIdx = x<0 & x>= -80; xLow = x(lowIdx);

                % Show and save the full field of view for the hybrid map
                % Within frequency comparison

                % Theta
                crossCorrTheta      = reshape(allXcorr.ccFullTheta,[401 imSize(1)*imSize(2)]);   crossCorrTheta(:,~corrMaskT)      = NaN;
                crossCorrSuperTheta = reshape(superXcorr.ccFullTheta,[401 imSize(1)*imSize(2)]); crossCorrSuperTheta(:,~corrMaskT) = NaN;
                crossCorrDeepTheta  = reshape(deepXcorr.ccFullTheta,[401 imSize(1)*imSize(2)]);  crossCorrDeepTheta(:,~corrMaskT)  = NaN;
                lagLowTheta         = tempProfileNoRef(iDate,iRun).lagLowTheta;
            
                % Alpha
                crossCorrAlpha      = reshape(allXcorr.ccFullAlpha,[401 imSize(1)*imSize(2)]);   crossCorrAlpha(:,~corrMaskT)      = NaN;
                crossCorrSuperAlpha = reshape(superXcorr.ccFullAlpha,[401 imSize(1)*imSize(2)]); crossCorrSuperAlpha(:,~corrMaskT) = NaN;
                crossCorrDeepAlpha  = reshape(deepXcorr.ccFullAlpha,[401 imSize(1)*imSize(2)]);  crossCorrDeepAlpha(:,~corrMaskT)  = NaN;
                lagLowAlpha         = tempProfileNoRef(iDate,iRun).lagLowAlpha;
              
                % Beta
                crossCorrBeta      = reshape(allXcorr.ccFullBeta,[401 imSize(1)*imSize(2)]);   crossCorrBeta(:,~corrMaskT)      = NaN;
                crossCorrSuperBeta = reshape(superXcorr.ccFullBeta,[401 imSize(1)*imSize(2)]); crossCorrSuperBeta(:,~corrMaskT) = NaN;
                crossCorrDeepBeta  = reshape(deepXcorr.ccFullBeta,[401 imSize(1)*imSize(2)]);  crossCorrDeepBeta(:,~corrMaskT)  = NaN;
                lagLowBeta         = tempProfileNoRef(iDate,iRun).lagLowBeta;
               
                % Gamma
                crossCorrGamma      = reshape(allXcorr.ccFull,[401 imSize(1)*imSize(2)]);   crossCorrGamma(:,~corrMaskT)      = NaN;
                crossCorrSuperGamma = reshape(superXcorr.ccFull,[401 imSize(1)*imSize(2)]); crossCorrSuperGamma(:,~corrMaskT) = NaN;
                crossCorrDeepGamma  = reshape(deepXcorr.ccFull,[401 imSize(1)*imSize(2)]);  crossCorrDeepGamma(:,~corrMaskT)  = NaN;
                lagLowGamma         = tempProfileNoRef(iDate,iRun).lagLow;
               
                % Spiking
                crossCorrSpiking      = reshape(allXcorr.ccFullRaw,[401 imSize(1)*imSize(2)]);   crossCorrSpiking(:,~corrMaskT)      = NaN;
                crossCorrSuperSpiking = reshape(superXcorr.ccFullRaw,[401 imSize(1)*imSize(2)]); crossCorrSuperSpiking(:,~corrMaskT) = NaN;
                crossCorrDeepSpiking  = reshape(deepXcorr.ccFullRaw,[401 imSize(1)*imSize(2)]);  crossCorrDeepSpiking(:,~corrMaskT)  = NaN;
                lagLowSpiking         = tempProfileNoRef(iDate,iRun).lagLowRaw;               

                if iType >= 2 && chLen~=0
                    crossCorrMidTheta = reshape(midXCorr.ccFullTheta,[401 imSize(1)*imSize(2)]); crossCorrMidTheta(:,~corrMaskT)   = NaN;
                    crossCorrMidAlpha = reshape(midXCorr.ccFullAlpha,[401 imSize(1)*imSize(2)]); crossCorrMidAlpha(:,~corrMaskT)   = NaN;
                    crossCorrMidBeta = reshape(midXCorr.ccFullBeta,[401 imSize(1)*imSize(2)]);   crossCorrMidBeta(:,~corrMaskT)    = NaN;
                    crossCorrMidGamma = reshape(midXCorr.ccFull,[401 imSize(1)*imSize(2)]);      crossCorrMidGamma(:,~corrMaskT)   = NaN;
                    crossCorrMidSpiking = reshape(midXCorr.ccFullRaw,[401 imSize(1)*imSize(2)]); crossCorrMidSpiking(:,~corrMaskT) = NaN;
                end               
                   
                % Correlate each hybrid map(all,super,deep,mid maps)
                % with FC map and correlate super/deep/mid maps to one
                % another
                for iBand = 1:5
                    clear bandName mapsAll crossCorrSuperR crossCorrDeepR crossCorrMidR
                    switch iBand
                        case 1
                            bandName        = 'Theta';
                            mapsAll         = crossCorrTheta;
                            crossCorrSuperR = crossCorrSuperTheta;
                            crossCorrDeepR  = crossCorrDeepTheta;
                            lagLow          = lagLowTheta;

                            if iType >= 2 && chLen~=0
                                crossCorrMidR   = crossCorrMidTheta;                                
                            end

                        case 2
                            bandName        = 'Alpha';
                            mapsAll         = crossCorrAlpha;
                            crossCorrSuperR = crossCorrSuperAlpha;
                            crossCorrDeepR  = crossCorrDeepAlpha;
                            lagLow          = lagLowAlpha;

                            if iType >= 2 && chLen~=0
                                crossCorrMidR   = crossCorrMidAlpha;
                            end

                        case 3
                            bandName        = 'Beta';
                            mapsAll         = crossCorrBeta;
                            crossCorrSuperR = crossCorrSuperBeta;
                            crossCorrDeepR  = crossCorrDeepBeta;
                            lagLow          = lagLowBeta;

                            if iType >= 2 && chLen~=0
                                crossCorrMidR   = crossCorrMidBeta;
                            end

                        case 4
                            bandName        = 'Gamma';
                            mapsAll         = crossCorrGamma;
                            crossCorrSuperR = crossCorrSuperGamma;
                            crossCorrDeepR  = crossCorrDeepGamma;
                            lagLow          = lagLowGamma;

                            if iType >= 2 && chLen~=0
                                crossCorrMidR   = crossCorrMidGamma;
                            end

                        case 5
                            bandName        = 'Spiking';
                            mapsAll         = crossCorrSpiking;
                            crossCorrSuperR = crossCorrSuperSpiking;
                            crossCorrDeepR  = crossCorrDeepSpiking;
                            lagLow          = lagLowSpiking;

                            if iType >= 2 && chLen~=0
                                crossCorrMidR   = crossCorrMidSpiking;
                            end
                    end

                    for iMap = 1:size(crossCorrGamma,1)
                        corrFCHybridT(iBand,iType,iMap)   = corr(fcMap,mapsAll(iMap,:)','rows','complete');
                        corrFCSuperT(iBand,iType,iMap)    = corr(fcMap,crossCorrSuperR(iMap,:)','rows','complete');
                        corrFCDeepT(iBand,iType,iMap)     = corr(fcMap,crossCorrDeepR(iMap,:)','rows','complete');
                        super_DeepCorrT(iBand,iType,iMap) = corr(crossCorrSuperR(iMap,:)',crossCorrDeepR(iMap,:)','rows','complete');

                        if iType >= 2 && chLen~=0
                            corrFCMidT(iBand,iType,iMap)     = corr(fcMap,crossCorrMidR(iMap,:)','rows','complete');
                            super_MidCorrT(iBand,iType,iMap) = corr(crossCorrSuperR(iMap,:)',crossCorrMidR(iMap,:)','rows','complete');
                            mid_DeepCorrT(iBand,iType,iMap)  = corr(crossCorrMidR(iMap,:)',crossCorrDeepR(iMap,:)','rows','complete');
                        end
                    end

                    % Get the peak negative correlation value and the lags
                    % for all channel hybrids
                    [peakNegValsAllT(iBand,iType),peakNegTimesAllT(iBand,iType)] = min(squeeze(corrFCHybridT(iBand,iType,negIdx)));
                    peakNegTimesAllT(iBand,iType) = negVals(peakNegTimesAllT(iBand,iType))./10;

                    % Average the super/deep hybrid maps for a certain range and correlate them
                    % with one another for 10/10 and 6/6
                    clear superHybridAvg deepHybridAvg
                    superHybridAvg = mean(crossCorrSuperR(lagFrameRange,:),1,'omitnan');
                    deepHybridAvg  = mean(crossCorrDeepR(lagFrameRange,:),1,'omitnan');
                    super_DeepAvgFramesT(iBand,iType) = corr(superHybridAvg',deepHybridAvg','rows','complete');

                    if iType >= 2 && chLen~=0
                        midHybridAvg = mean(crossCorrMidR(lagFrameRange,:),1,'omitnan');
                        super_MidAvgFramesT(iBand,iType) = corr(superHybridAvg',midHybridAvg','rows','complete');
                        deep_MidAvgFramesT(iBand,iType)  = corr(deepHybridAvg',midHybridAvg','rows','complete');
                    end

                    if ~exist([dataDir '\HybridMapFOV_' fileName '_' bandName '.png'],'file') || 1
                        clear frameLow grayImFull crossCorr
                        frameLow       = (x == lagLow);
                        crossCorr      = reshape(mapsAll,[401 imSize(1) imSize(2)]);
                        figure; imagesc(squeeze(crossCorr(frameLow,:,:))); hold on; axis image off;
                        colormap(flipud(jet)); caxis([-0.5 0.5]);
                        grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));
                        h = imagesc(grayImFull); hold off; set(h,'AlphaData',~corrMask);
                        title(['Peak negative at ' num2str(lagLow/10) 's']);
                        f = gcf; exportgraphics(f,[dataDir '\HybridMapFOV_' fileName '_' bandName '.png'],'Resolution',300); close gcf;
                    end

                    if ~exist([dataDir '\HybridMapFOV_Super_' fileName '_' bandName '.png'],'file')  || 1
                        clear frameLow grayImFull crossCorrSuper
                        frameLow       = (x == lagLow);
                        crossCorrSuper = reshape(crossCorrSuperR,[401 imSize(1) imSize(2)]);
                        figure; imagesc(squeeze(crossCorrSuper(frameLow,:,:))); hold on; axis image off;
                        colormap(flipud(jet)); caxis([-0.5 0.5]);
                        grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));
                        h = imagesc(grayImFull); hold off; set(h,'AlphaData',~corrMask);
                        title(['Superficial channels-Peak negative at ' num2str(lagLow/10) 's']);
                        f = gcf; exportgraphics(f,[dataDir '\HybridMapFOV_Super_' fileName '_' bandName '.png'],'Resolution',300); close gcf;
                    end

                    if ~exist([dataDir '\HybridMapFOV_Deep_' fileName '_' bandName '.png'],'file')  || 1
                        clear frameLow grayImFull crossCorrDeep
                        frameLow      = (x == lagLow);
                        crossCorrDeep = reshape(crossCorrDeepR,[401 imSize(1) imSize(2)]);
                        figure; imagesc(squeeze(crossCorrDeep(frameLow,:,:))); hold on; axis image off;
                        colormap(flipud(jet)); caxis([-0.5 0.5]);
                        grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));
                        h = imagesc(grayImFull); hold off; set(h,'AlphaData',~corrMask);
                        title(['Deep channels-Peak negative at ' num2str(lagLow/10) 's']);
                        f = gcf; exportgraphics(f,[dataDir '\HybridMapFOV_Deep_' fileName '_' bandName '.png'],'Resolution',300); close gcf;
                    end

                    if (~exist([dataDir '\HybridMapFOV_Mid_' fileName '_' bandName '.png'],'file')) && iType >= 2 && chLen~=0 && 1
                        clear frameLow grayImFull crossCorrMid 
                        frameLow     = (x == lagLow);
                        crossCorrMid = reshape(crossCorrMidR,[401 imSize(1) imSize(2)]);
                        figure; imagesc(squeeze(crossCorrMid(frameLow,:,:))); hold on; axis image off;
                        colormap(flipud(jet)); caxis([-0.5 0.5]);
                        grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));
                        h = imagesc(grayImFull); hold off; set(h,'AlphaData',~corrMask);
                        title(['Middle channels-Peak negative at ' num2str(lagLow/10) 's']);
                        f = gcf; exportgraphics(f,[dataDir '\HybridMapFOV_Mid_' fileName '_' bandName '.png'],'Resolution',300); close gcf;
                    end
                end

                % Correlate maps across frequencies for each layer
                % compartment
                for iBand1 = 1:5
                    clear super1 deep1 mid1
                    switch iBand1
                        case 1
                            super1 = mean(crossCorrSuperTheta(lagFrameRange,:),1,'omitnan');
                            deep1  = mean(crossCorrDeepTheta(lagFrameRange,:),1,'omitnan');
                            if iType >= 2 && chLen~=0
                                mid1 = mean(crossCorrMidTheta(lagFrameRange,:),1,'omitnan');
                            end
                        
                        case 2
                            super1 = mean(crossCorrSuperAlpha(lagFrameRange,:),1,'omitnan');
                            deep1  = mean(crossCorrDeepAlpha(lagFrameRange,:),1,'omitnan');
                            if iType >= 2 && chLen~=0
                                mid1 = mean(crossCorrMidAlpha(lagFrameRange,:),1,'omitnan');
                            end
                        
                        case 3
                            super1 = mean(crossCorrSuperBeta(lagFrameRange,:),1,'omitnan');
                            deep1  = mean(crossCorrDeepBeta(lagFrameRange,:),1,'omitnan');
                            if iType >= 2 && chLen~=0
                                mid1 = mean(crossCorrMidBeta(lagFrameRange,:),1,'omitnan');
                            end
                       
                        case 4
                            super1 = mean(crossCorrSuperGamma(lagFrameRange,:),1,'omitnan');
                            deep1  = mean(crossCorrDeepGamma(lagFrameRange,:),1,'omitnan');
                            if iType >= 2 && chLen~=0
                                mid1 = mean(crossCorrMidGamma(lagFrameRange,:),1,'omitnan');
                            end
                        
                        case 5
                            super1 = mean(crossCorrSuperSpiking(lagFrameRange,:),1,'omitnan');
                            deep1  = mean(crossCorrDeepSpiking(lagFrameRange,:),1,'omitnan');
                            if iType >= 2 && chLen~=0
                                mid1 = mean(crossCorrMidSpiking(lagFrameRange,:),1,'omitnan');
                            end
                    end 
                        
                    for iBand2 = 1:5
                        clear super2 deep2 mid2
                        switch iBand2
                            case 1
                                super2 = mean(crossCorrSuperTheta(lagFrameRange,:),1,'omitnan');
                                deep2  = mean(crossCorrDeepTheta(lagFrameRange,:),1,'omitnan');
                                if iType >= 2 && chLen~=0
                                    mid2 = mean(crossCorrMidTheta(lagFrameRange,:),1,'omitnan');
                                end

                            case 2
                                super2 = mean(crossCorrSuperAlpha(lagFrameRange,:),1,'omitnan');
                                deep2  = mean(crossCorrDeepAlpha(lagFrameRange,:),1,'omitnan');
                                if iType >= 2 && chLen~=0
                                    mid2 = mean(crossCorrMidAlpha(lagFrameRange,:),1,'omitnan');
                                end

                            case 3
                                super2 = mean(crossCorrSuperBeta(lagFrameRange,:),1,'omitnan');
                                deep2  = mean(crossCorrDeepBeta(lagFrameRange,:),1,'omitnan');
                                if iType >= 2 && chLen~=0
                                    mid2 = mean(crossCorrMidBeta(lagFrameRange,:),1,'omitnan');
                                end

                            case 4
                                super2 = mean(crossCorrSuperGamma(lagFrameRange,:),1,'omitnan');
                                deep2  = mean(crossCorrDeepGamma(lagFrameRange,:),1,'omitnan');
                                if iType >= 2 && chLen~=0
                                    mid2 = mean(crossCorrMidGamma(lagFrameRange,:),1,'omitnan');
                                end

                            case 5
                                super2 = mean(crossCorrSuperSpiking(lagFrameRange,:),1,'omitnan');
                                deep2  = mean(crossCorrDeepSpiking(lagFrameRange,:),1,'omitnan');
                                if iType >= 2 && chLen~=0
                                    mid2 = mean(crossCorrMidSpiking(lagFrameRange,:),1,'omitnan');
                                end
                        end

                        superHybridAllBandsT(iType,iBand1,iBand2) = corr(super1',super2','rows','complete');
                        deepHybridAllBandsT(iType,iBand1,iBand2)  = corr(deep1',deep2','rows','complete');
                        if iType>=2 && chLen~=0
                            midHybridAllBandsT(iType,iBand1,iBand2) = corr(mid1',mid2','rows','complete');
                        end

                        if iBand1~=iBand2
                            crossFreqCrossLayerHybridT(iType,iBand1,iBand2,1,1) = corr(super1',super2','rows','complete');
                            crossFreqCrossLayerHybridT(iType,iBand1,iBand2,1,3) = corr(super1',deep2','rows','complete');
                            crossFreqCrossLayerHybridT(iType,iBand1,iBand2,3,1) = corr(deep1',super2','rows','complete');
                            crossFreqCrossLayerHybridT(iType,iBand1,iBand2,3,3) = corr(deep1',deep2','rows','complete');

                            if  iType>=2 && chLen~=0
                                crossFreqCrossLayerHybridT(iType,iBand1,iBand2,1,2) = corr(super1',mid2','rows','complete');
                                crossFreqCrossLayerHybridT(iType,iBand1,iBand2,2,1) = corr(mid1',super2','rows','complete');
                                crossFreqCrossLayerHybridT(iType,iBand1,iBand2,2,2) = corr(mid1',mid2','rows','complete');
                                crossFreqCrossLayerHybridT(iType,iBand1,iBand2,2,3) = corr(mid1',deep2','rows','complete');
                                crossFreqCrossLayerHybridT(iType,iBand1,iBand2,3,2) = corr(deep1',mid2','rows','complete');
                            end

                        end
                        

                    end
                end
            end

            save([dataDir '\processedHybridMapVars.mat'],'corrFCHybridT','corrFCSuperT',...
                'corrFCDeepT','super_DeepCorrT','corrFCMidT','super_MidCorrT','mid_DeepCorrT',...
                'peakNegValsAllT','peakNegTimesAllT','super_DeepAvgFramesT','super_MidAvgFramesT',...
                'deep_MidAvgFramesT','superHybridAllBandsT','deepHybridAllBandsT','midHybridAllBandsT',...
                'crossFreqCrossLayerHybridT','fcMap');

            corrFCHybrid(iDate,iRun,:,:,:)   = corrFCHybridT;
            corrFCSuper(iDate,iRun,:,:,:)    = corrFCSuperT;
            corrFCDeep(iDate,iRun,:,:,:)     = corrFCDeepT;
            corrFCMid(iDate,iRun,:,:,:)      = corrFCMidT;
           
            super_DeepCorr(iDate,iRun,:,:,:) = super_DeepCorrT;           
            super_MidCorr(iDate,iRun,:,:,:)  = super_MidCorrT;
            mid_DeepCorr(iDate,iRun,:,:,:)   = mid_DeepCorrT ;

            peakNegValsAll(iDate,iRun,:,:)  = peakNegValsAllT;
            peakNegTimesAll(iDate,iRun,:,:) = peakNegTimesAllT;

            super_DeepAvgFrames(iDate,iRun,:,:) = super_DeepAvgFramesT;
            super_MidAvgFrames(iDate,iRun,:,:)  = super_MidAvgFramesT;
            deep_MidAvgFrames(iDate,iRun,:,:)   = deep_MidAvgFramesT;

            superHybridAllBands(iDate,iRun,:,:,:) = superHybridAllBandsT;
            midHybridAllBands(iDate,iRun,:,:,:)   = midHybridAllBandsT;
            deepHybridAllBands(iDate,iRun,:,:,:)  = deepHybridAllBandsT;

            crossFreqCrossLayerHybrid(iDate,iRun,:,:,:,:,:) = crossFreqCrossLayerHybridT; 
        else
            clear allHybridVars
            allHybridVars                    = matfile([dataDir '\processedHybridMapVars.mat']);
            corrFCHybrid(iDate,iRun,:,:,:)   = allHybridVars.corrFCHybridT;
            corrFCSuper(iDate,iRun,:,:,:)    = allHybridVars.corrFCSuperT;
            corrFCDeep(iDate,iRun,:,:,:)     = allHybridVars.corrFCDeepT;
            super_DeepCorr(iDate,iRun,:,:,:) = allHybridVars.super_DeepCorrT;

            corrFCMid(iDate,iRun,:,:,:)     = allHybridVars.corrFCMidT;
            super_MidCorr(iDate,iRun,:,:,:) = allHybridVars.super_MidCorrT;
            mid_DeepCorr(iDate,iRun,:,:,:)  = allHybridVars.mid_DeepCorrT ;

            peakNegValsAll(iDate,iRun,:,:)  = allHybridVars.peakNegValsAllT;
            peakNegTimesAll(iDate,iRun,:,:) = allHybridVars.peakNegTimesAllT;

            super_DeepAvgFrames(iDate,iRun,:,:) = allHybridVars.super_DeepAvgFramesT;
            super_MidAvgFrames(iDate,iRun,:,:)  = allHybridVars.super_MidAvgFramesT;
            deep_MidAvgFrames(iDate,iRun,:,:)   = allHybridVars.deep_MidAvgFramesT;

            superHybridAllBands(iDate,iRun,:,:,:) = allHybridVars.superHybridAllBandsT;
            midHybridAllBands(iDate,iRun,:,:,:)   = allHybridVars.midHybridAllBandsT;
            deepHybridAllBands(iDate,iRun,:,:,:)  = allHybridVars.deepHybridAllBandsT;
           
            crossFreqCrossLayerHybrid(iDate,iRun,:,:,:,:,:) = allHybridVars.crossFreqCrossLayerHybridT;

        end
    end
end
toc;

%%
corrFCHybridT = reshape(corrFCHybrid,[size(corrFCHybrid,1)*size(corrFCHybrid,2) size(corrFCHybrid,3) size(corrFCHybrid,4) size(corrFCHybrid,5)]);
nanRow        = isnan(corrFCHybridT(:,1,1,1));
corrFCHybridT(nanRow,:,:,:) = [];

gammaProfile  = squeeze(corrFCHybridT(:,4,2,:));
x = -200:200; 
negIdx = x<0 & x>=-70; negValsT = x(negIdx);
negIdx2 = x<0 & x>=-150; negValsT2 = x(negIdx2);
for iR = 1:size(gammaProfile,1)
    [negVals(iR,1),negTimes(iR,1)] = min(gammaProfile(iR,negIdx));
    negTimes(iR,1) = negValsT(negTimes(iR,1))./10;

    [negVals2(iR,1),negTimes2(iR,1)] = min(gammaProfile(iR,negIdx2));
    negTimes2(iR,1) = negValsT2(negTimes2(iR,1))./10;
end



% figure; plot(negTimes(goodRunsSpatial),negVals(goodRunsSpatial),'.b','MarkerSize',20); hold on;
% plot(negTimes(~goodRunsSpatial),negVals(~goodRunsSpatial),'.r','MarkerSize',20); hold on;
% box off; xline(0,'LineWidth',2); yline(0,'LineWidth',2); ylim([-1 0.2]); xlim([-20 20]); 
% peakNegValsAllT = reshape(peakNegValsAll,[size(peakNegValsAll,1)*size(peakNegValsAll,2) size(peakNegValsAll,3) size(peakNegValsAll,4)]);
% peakNegValsAllT(nanRow,:,:) = [];
% peakNegTimesAllT = reshape(peakNegTimesAll,[size(peakNegTimesAll,1)*size(peakNegTimesAll,2) size(peakNegTimesAll,3) size(peakNegTimesAll,4)]);
% peakNegTimesAllT(nanRow,:,:) = [];
% plot(squeeze(peakNegTimesAllT(goodRunsSpatial,4,2)),squeeze(peakNegValsAllT(goodRunsSpatial,4,2)),'.b','MarkerSize',20); hold on;
% plot(squeeze(peakNegTimesAllT(~goodRunsSpatial,4,2)),squeeze(peakNegValsAllT(~goodRunsSpatial,4,2)),'.r','MarkerSize',20); hold on;
% box off; xline(0,'LineWidth',2); yline(0,'LineWidth',2); xlim([-20 20]); ylim([-1 0.2]);


figure;subplot(121); plot(x,gammaProfile(goodRunsSpatial,:),'b'); hold on; 
plot(x,gammaProfile(~goodRunsSpatial,:),'r'); box off; ylim([-1 1]);
box off; xline(0,'LineWidth',2); yline(0,'LineWidth',2); xticklabels(-20:10:20); xlim([-160 20]);

subplot(122); plot(negTimes(goodRunsSpatial),negVals(goodRunsSpatial),'.b','MarkerSize',25); hold on;
plot(negTimes2(goodRunsSpatial),negVals2(goodRunsSpatial),'.b','MarkerSize',25); hold on;

plot(negTimes(~goodRunsSpatial),negVals(~goodRunsSpatial),'.r','MarkerSize',25); hold on;
box off; xline(0,'LineWidth',2); yline(0,'LineWidth',2); ylim([-1 0.5]); yticks(-1:0.2:0.5); xlim([-20 20]);
ylabel('Correlations between FC and cross-modal maps');xlabel('Lag (s)');xticks(-20:5:20); xlim([-16 2]);
sgtitle([ monkeyName ' - FOV all runs']); legend ('Accepted runs','Rejected runs','Location','northeast')



% gammaTimes = squeeze(peakNegTimesAllT(:,4,2));
% gammaLags  = squeeze(peakNegValsAllT(:,4,2));
% 
% X = [gammaTimes gammaLags]; Y = single(goodRunsSpatial); Y(Y==0) =2; 
% mdl = fitcnb(X,Y,'ClassNames',[1 2]);
% colorVals = [0.85 0.325 0.098; 0 0.447 0.741];
% 
% x1Range = -20:0.01:20;
% x2Range = -1:0.01:0.2;
% 
% [xx1, xx2] = meshgrid(x1Range,x2Range);
% XGrid = [xx1(:) xx2(:)];
% 
% predVals = predict(mdl,XGrid); 
% figure; gscatter(xx1(:),xx2(:),predVals,colorVals); hold on; legend off;
% plot(squeeze(peakNegTimesAllT(goodRunsSpatial,4,2)),squeeze(peakNegValsAllT(goodRunsSpatial,4,2)),'.r','MarkerSize',20); 
% plot(squeeze(peakNegTimesAllT(~goodRunsSpatial,4,2)),squeeze(peakNegValsAllT(~goodRunsSpatial,4,2)),'.b','MarkerSize',20);
% box off; xline(0,'LineWidth',2); yline(0,'LineWidth',2); xlim([-20 20]); ylim([-1 0.2]);


%% Check run 02/20/2024 - run 07 for Whiskey -- don't run it for Charlie
expDate   = '02_20_2024';
runName   = 'run07';

serverDir = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\kem294\Data\' ...
    monkeyName '_SqM\' hemisphere ' Hemisphere\'  expDate '\' runName];
dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName];

% Get masks and greens
if exist([dataDir '\green0' runName(end) '_Edited.png'],'file')
    greenTemp = imread([dataDir '\green0' runName(end) '_Edited.png']);

elseif exist([dataDir '\green0' runName(end) '_Edited.bmp'],'file')
    greenTemp = imread([dataDir '\green0' runName(end) '_Edited.bmp']);

else
    error('Greens are not edited...');
end

greenTemp = greenTemp(:,:,1);

if exist([dataDir '\clipMask0' runName(end) '.BMP'],'file') == 0
    clipMask = imread([dataDir '\clipMask0' runName(end) '.png']);
else
    clipMask = imread([dataDir '\clipMask0' runName(end) '.bmp']);
end

% This mask includes vessels and electrode(s)
if exist([dataDir '\maskSkull0' runName(end) '.bmp'],'file') == 0
    elecMask = imread([dataDir '\maskSkull0' runName(end) '.png']);
else
    elecMask = imread([dataDir '\maskSkull0' runName(end) '.bmp']);
end

clipMask       = imresize(clipMask,1/3); % Resize mask
clipMask       = clipMask(:,:,1)>0; % Converting to 0s and 1s
elecMask       = imresize(elecMask-100,1/3); % Binarize image
elecMask       = elecMask(:,:,1)>0;
clipMaskCortex = clipMask & ~elecMask; % Mask of cortex sans vessels
imSize         = size(clipMask);

% Mask for comparing the imaging and hybrid map (after removing
% areas beyond lateral sulcus)
if exist([dataDir '\corrMask0' runName(end) '.BMP'],'file')
    corrMask = imread([dataDir '\corrMask0' runName(end) '.bmp']);
    corrMaskFlag = 1;
elseif exist([dataDir '\corrMask0' runName(end) '.PNG'],'file')
    corrMask = imread([dataDir '\corrMask0' runName(end) '.png']);
    corrMaskFlag = 1;
else
    corrMaskFlag = 0;
end

% Mask for comparing the imaging and hybrid map (after removing
% areas beyond lateral sulcus)...
if exist([dataDir '\corrMask0' runName(end) '.BMP'],'file')
    corrMask = imread([dataDir '\corrMask0' runName(end) '.bmp']);
    corrMaskFlag = 1;
elseif exist([dataDir '\corrMask0' runName(end) '.PNG'],'file')
    corrMask = imread([dataDir '\corrMask0' runName(end) '.png']);
    corrMaskFlag = 1;
else
    corrMaskFlag = 0;
end

if corrMaskFlag
    corrMask = imresize(corrMask,1/3); % Resize the mask
    corrMask = corrMask(:,:,1)>0;
    corrMask = corrMask & ~elecMask;
else
    corrMask = clipMaskCortex;
end

corrMaskT = reshape(corrMask,[imSize(1)*imSize(2) 1]);

% Get ROI location and session-wise FC map
processedDat = matfile([dataDir '\processedFrames.mat']);
seedLocIn    = load([dataDir '\roiCenterLoc.mat']);
seedLocProbe = seedLocIn.seedLocProbe;
seedLocIn    = seedLocIn.seedLocIn;
circleRad    = round(38/(spatialBin)); % 500um radius
greenFig     = imresize(greenTemp,1/spatialBin,'OutputSize',[imSize(1) imSize(2)]);

pDatTemp = processedDat.tempBandPass;
seedSigT = calculateSeedSignal(greenFig,corrMask,...
    seedLocIn,circleRad,pDatTemp); % Get Gaussian weighted seed signal

fcMap = plotCorrMap(seedSigT,pDatTemp,0);
fcMap = reshape(fcMap,[imSize(1)*imSize(2) 1]);
fcMap(~corrMaskT) = NaN; 


% Get cross-modal maps
crossCorrFOV = matfile([serverDir '\crossCorrFOV_6_NoRef.mat']);
allXCorr     = crossCorrFOV.spatialProfile;
allLags      = crossCorrFOV.lagFull;

allHybridVars = matfile([dataDir '\processedHybridMapVars.mat']);
corrFCHybrid  = allHybridVars.corrFCHybridT;

peakNegValsAllRun  = allHybridVars.peakNegValsAllT;  peakNegValsAllRun  = peakNegValsAllRun(:,2);
peakNegTimesAllRun = allHybridVars.peakNegTimesAllT; peakNegTimesAllRun = peakNegTimesAllRun(:,2);
x = -200:200;

for iX = 1:size(peakNegTimesAll,1)
    peakNegIdx(iX,1) = find(x== peakNegTimesAll(iX).*10);
end

crossCorrAlpha = reshape(allXCorr.ccFullAlpha,[401 imSize(1)*imSize(2)]); crossCorrAlpha(:,~corrMaskT)   = NaN;
crossCorrBeta  = reshape(allXCorr.ccFullBeta,[401 imSize(1)*imSize(2)]);  crossCorrBeta(:,~corrMaskT)    = NaN;
crossCorrGamma = reshape(allXCorr.ccFull,[401 imSize(1)*imSize(2)]);      crossCorrGamma(:,~corrMaskT)   = NaN;
crossCorrSpiking = reshape(allXCorr.ccFullRaw,[401 imSize(1)*imSize(2)]); crossCorrSpiking(:,~corrMaskT) = NaN;

% alphaCrossModal   = squeeze(crossCorrAlpha(peakNegIdx(2),:));
% betaCrossModal    = squeeze(crossCorrBeta(peakNegIdx(3),:));
% gammaCrossModal   = squeeze(crossCorrGamma(peakNegIdx(4),:));
% spikingCrossModal = squeeze(crossCorrGamma(peakNegIdx(5),:));

map1 = crossCorrGamma; map1(:,~corrMaskT) = NaN;
for iM = 1:401    
    corrVals(iM,1) = corr(fcMap,map1(iM,:)','rows','complete');
end


%% Compile all FOV level variables
corrFCHybridT = reshape(corrFCHybrid,[size(corrFCHybrid,1)*size(corrFCHybrid,2) size(corrFCHybrid,3) size(corrFCHybrid,4) size(corrFCHybrid,5)]);
nanRow        = isnan(corrFCHybridT(:,1,1,1));
corrFCHybridT(nanRow,:,:,:) = [];
corrFCHybridT(~goodRunsSpatial,:,:,:) = [];

% Compile correlations between layer compartments and FC map
corrFCSuperT = mean(corrFCSuper(:,:,:,:,lagFrameRange),5,'omitnan'); 
corrFCSuperT = reshape(corrFCSuperT,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]); 
corrFCSuperT(nanRow,:,:) = []; 
corrFCSuperT(~(goodRunsSpatial & ~singleChFlag),:,:) = []; 

corrFCMidT = mean(corrFCMid(:,:,:,:,lagFrameRange),5,'omitnan'); 
corrFCMidT = reshape(corrFCMidT,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]); 
corrFCMidT(nanRow,:,:) = []; 
corrFCMidT(~(goodRunsSpatial & ~singleChFlag),:,:) = []; 

corrFCDeepT = mean(corrFCDeep(:,:,:,:,lagFrameRange),5,'omitnan'); 
corrFCDeepT = reshape(corrFCDeepT,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]); 
corrFCDeepT(nanRow,:,:) = []; 
corrFCDeepT(~(goodRunsSpatial & ~singleChFlag),:,:) = []; 

% Compile the peak negative correlations and lags for good runs
peakNegValsAllT = reshape(peakNegValsAll,[size(peakNegValsAll,1)*size(peakNegValsAll,2) size(peakNegValsAll,3) size(peakNegValsAll,4)]);
peakNegValsAllT(nanRow,:,:) = [];
peakNegValsAllT(~goodRunsSpatial,:,:) = [];

% Compile the across layers hybrid map correlations ie; S/D, S/M and M/D 
super_DeepAvgFramesT = reshape(super_DeepAvgFrames,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]);
super_DeepAvgFramesT(nanRow,:,:)           = [];
super_DeepAvgFramesT(~(goodRunsSpatial & ~singleChFlag),:,:) = [];

super_MidAvgFramesT = reshape(super_MidAvgFrames,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]);
super_MidAvgFramesT(nanRow,:,:)           = [];
super_MidAvgFramesT(~(goodRunsSpatial & ~singleChFlag),:,:) = [];

deep_MidAvgFramesT = reshape(deep_MidAvgFrames,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]);
deep_MidAvgFramesT(nanRow,:,:)           = [];
deep_MidAvgFramesT(~(goodRunsSpatial & ~singleChFlag),:,:) = [];

% Cross-frequency correlations between hybrid maps for each layer 
superHybridAllBandsT = reshape(superHybridAllBands,[size(superHybridAllBands,1)*size(superHybridAllBands,2) size(superHybridAllBands,3) size(superHybridAllBands,4) size(superHybridAllBands,5)]);
nanRow = isnan(superHybridAllBandsT(:,2,1,1));
superHybridAllBandsT(nanRow,:,:,:) = [];
superHybridAllBandsT(~(goodRunsSpatial & ~singleChFlag),:,:,:)   = [];

midHybridAllBandsT = reshape(midHybridAllBands,[size(superHybridAllBands,1)*size(superHybridAllBands,2) size(superHybridAllBands,3) size(superHybridAllBands,4) size(superHybridAllBands,5)]);
midHybridAllBandsT(nanRow,:,:,:) = [];
midHybridAllBandsT(~(goodRunsSpatial & ~singleChFlag),:,:,:)   = [];

deepHybridAllBandsT = reshape(deepHybridAllBands,[size(superHybridAllBands,1)*size(superHybridAllBands,2) size(superHybridAllBands,3) size(superHybridAllBands,4) size(superHybridAllBands,5)]);
deepHybridAllBandsT(nanRow,:,:,:) = [];
deepHybridAllBandsT(~(goodRunsSpatial & ~singleChFlag),:,:,:)   = [];

clear varInfo;
varFlag = 0;
try matfile(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat']);
    varInfo = who('-file',['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat']);
    if sum(ismember(varInfo,'corrFCHybridT'))==0 || sum(ismember(varInfo,'superHybridAllBandsT'))==0
        varFlag = 1;
    end
catch
    varFlag = 1;
end
        
if 1%varFlag 
    save(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat'],...
       'corrFCHybridT','corrFCSuperT','corrFCMidT','corrFCDeepT','peakNegValsAllT',...
       'super_DeepAvgFramesT','super_MidAvgFramesT','deep_MidAvgFramesT',...
       'superHybridAllBandsT','midHybridAllBandsT','deepHybridAllBandsT','-append');
end

%% Controls - Disambiguating channel position vs channel count to solidify the results
posControlVsRunCorr = NaN(size(processedDat,1),size(processedDat,2),3,5);
superDeepPosCheck   = NaN(size(processedDat,1),size(processedDat,2),3, 401,5);
fs                  = 1e3;
for iDate = 1: size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);

    for iRun = 1:size(allRuns{iDate,1})
        clear runName dataDir
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

        clc; disp([monkeyName ' '  expDate ' run: ' runName]);

        if ~exist([dataDir '\superDeepControls.mat'],'file')
            if exist([dataDir '\clipMask0' runName(end) '.BMP'],'file') == 0
                clipMask = imread([dataDir '\clipMask0' runName(end) '.png']);
            else
                clipMask = imread([dataDir '\clipMask0' runName(end) '.bmp']);
            end

            % This mask includes vessels and electrode(s)
            if exist([dataDir '\maskSkull0' runName(end) '.bmp'],'file') == 0
                elecMask = imread([dataDir '\maskSkull0' runName(end) '.png']);
            else
                elecMask = imread([dataDir '\maskSkull0' runName(end) '.bmp']);
            end

            clipMask       = imresize(clipMask,1/3); % Resize mask
            clipMask       = clipMask(:,:,1)>0; % Converting to 0s and 1s
            elecMask       = imresize(elecMask-100,1/3); % Binarize image
            elecMask       = elecMask(:,:,1)>0;
            clipMaskCortex = clipMask & ~elecMask; % Mask of cortex sans vessels
            imSize         = size(clipMask);

            % Mask for comparing the imaging and hybrid map (after removing
            % areas beyond lateral sulcus)
            if exist([dataDir '\corrMask0' runName(end) '.BMP'],'file')
                corrMask = imread([dataDir '\corrMask0' runName(end) '.bmp']);
                corrMaskFlag = 1;
            elseif exist([dataDir '\corrMask0' runName(end) '.PNG'],'file')
                corrMask = imread([dataDir '\corrMask0' runName(end) '.png']);
                corrMaskFlag = 1;
            else
                corrMaskFlag = 0;
            end

            if corrMaskFlag
                corrMask = imresize(corrMask,1/3); % Resize the mask
                corrMask = corrMask(:,:,1)>0;
                corrMask = corrMask & ~elecMask;
            else
                corrMask = clipMaskCortex;
            end

            corrMaskT = reshape(corrMask,[imSize(1)*imSize(2) 1]);

            % Obtaining ephys data
            clear probeCh timeStamp badTimes badTimesThIm chInCortex szLFP...
                szIm szMin badTimeIm badTimes10Hz
            probeCh      = probe{iRun,iDate}.probeCh;
            timeStamp    = probe{iRun,iDate}.timeStamp;
            badTimes     = badTimesLFP{iDate,iRun};
            badTimesThIm = badTimeThresh{iDate,iRun};
            chInCortex   = estChInCortex{1,iDate}(iRun,:);

            % Removing bad channels
            probeCh(:,badCh{iDate,iRun}) = [];
            szLFP   = size(probeCh,1);

            if chInCortex(1)-chInCortex(2) ~= 0

                % Get imaging data
                clear pDatTemp imSize processedDat10
                pDatTemp = processedDat{iDate,iRun}.tempBandPass;
                imSize   = size(pDatTemp);
                pDatTemp = reshape(pDatTemp,[imSize(1)*imSize(2) imSize(3)]);

                % Upsampling imaging data to 10 Hz
                disp('Upsampling imaging data to 10 Hz... ');
                parfor iP = 1:size(pDatTemp,1)
                    processedDat10(iP,:) = interp(pDatTemp(iP,:),5);
                end

                szIm = size(processedDat10,2)*100;

                if ~(szLFP == szIm)
                    szMin          = min([szLFP, szIm]);
                    probeCh        = probeCh(1:szMin,:);
                    processedDat10 = processedDat10(:,1:floor(szMin/100));

                    badTimes(badTimes>szMin)         = [];
                    badTimesThIm(badTimesThIm>szMin) = [];
                else
                    szMin = szLFP;
                end

                timeStampSorted = timeStamp- timeStamp(1);
                badTimes10Hz    = unique(badTimesThIm./1000);
                badTimeIm       = [];

                for iT = 1: length(badTimes10Hz)
                    badTimeIm(iT) = find((floor(abs(timeStampSorted - badTimes10Hz(iT))*100)./100)<=0.05,1,'first');
                end

                badTimeIm = unique(badTimeIm);
                badTimeIm(badTimeIm>size(processedDat10,2)) = [];
                processedDat10(:,badTimeIm) = [];
                probeCh(badTimes,:) = []; % Remove bad times from LFP

                % Remove bad times determined visually from spectrogram
                [probeCh,~,processedDat10] = removeBadTimesFromSpec(monkeyName,expDate,runName,probeCh,...
                    [],processedDat10);

                clear ch
                ch = chInCortex(1):chInCortex(2);

                % 1. Pick 6 random channels 5 times and get the mean hybrid maps
                % Pick 5 iterations of groups of 6 channels - Check if position of
                % channels affect the super/deep hybrids
                clear superDeepPosCheckTemp posControlVsRunCorrTemp
                for iBand = 1:3
                    for iIter = 1:5
                        clear superCh deepCh
                        rng("shuffle"); superCh = ch(randperm(6,6));
                        rng("shuffle"); deepCh  = ch(randperm(6,6));

                        % Get infraslow ephys of gamma bands
                        clear bandVals ephysSuper ephysDeep
                        switch iBand
                            case 1
                                bandVals = [8 12];
                            case 2
                                bandVals = [13 30];
                            case 3
                                bandVals = [30 90]; 
                        end
                        [bT,aT]   = butter(3,bandVals./(1e3/2),'bandpass');  % Gamma band filtering parameters

                        ephysSuper = abs(single(filtfilt(bT,aT,double(probeCh(:,superCh)))));
                        ephysDeep  = abs(single(filtfilt(bT,aT,double(probeCh(:,deepCh)))));

                        envelopeDatSuper = envelope(ephysSuper,5);
                        envelopeDatDeep  = envelope(ephysDeep,5);

                        % Get filter parameters...
                        clear z p k sos g enSize                      

                        % Bandpass - 0.01 Hz - 0.1 Hz
                        [z,p,k] = butter(3,[0.01 0.1]./(fs/2),'bandpass');
                        [sos,g] = zp2sos(z,p,k);
                        enSize  = size(envelopeDatSuper);

                        envelopeFiltSuper = filtfilt(sos,g,double([envelopeDatSuper; envelopeDatSuper; envelopeDatSuper]));
                        envelopeFiltSuper = envelopeFiltSuper(enSize(1)+1:(end-enSize(1)),:);

                        envelopeFiltDeep = filtfilt(sos,g,double([envelopeDatDeep; envelopeDatDeep; envelopeDatDeep]));
                        envelopeFiltDeep = envelopeFiltDeep(enSize(1)+1:(end-enSize(1)),:);

                        infraEphysSuper = mean(single(downsample(envelopeFiltSuper,100)),2,'omitnan');
                        infraEphysDeep  = mean(single(downsample(envelopeFiltDeep,100)),2,'omitnan');

                        % Check size of timeseries of both modalities
                        clear szIm szLFP processedDat10Temp
                        szIm  = size(processedDat10,2);
                        szLFP = size(infraEphysSuper,1);

                        if ~(szLFP == szIm)
                            szMin               = min([szLFP, szIm]);
                            infraEphysSuper     = infraEphysSuper(1:szMin,:);
                            infraEphysDeep      = infraEphysDeep(1:szMin,:);
                            processedDat10Temp  = processedDat10(:,1:szMin);
                        else
                            processedDat10Temp  = processedDat10;
                        end

                        tic;
                        disp('Performing cross correlations...');

                        parfor iP = 1:size(processedDat10Temp,1)
                            [ccFullSuper(:,iP),~] = xcorr(infraEphysSuper',processedDat10Temp(iP,:),200,'normalized');
                            [ccFullDeep(:,iP),~]  = xcorr(infraEphysDeep',processedDat10Temp(iP,:),200,'normalized');
                        end

                        disp('Correlating superficial and deep hybrids for all lags...');

                        parfor iMap = 1:size(ccFullSuper,1)
                            superDeepPosCheckTemp(iBand,iMap,iIter) = corr(ccFullSuper(iMap,:)',ccFullDeep(iMap,:)','rows','complete');
                        end
                        toc;
                    end

                    runSuperDeepCorr(:,iBand)        = squeeze(super_DeepCorr(iDate,iRun,iBand,2,:));
                    posControlVsRunCorrTemp(iBand,:) = corr(runSuperDeepCorr(:,iBand) ,squeeze(superDeepPosCheckTemp(iBand,:,:)),'rows','complete');
                end
               
                save([dataDir '\superDeepControls.mat'],'superDeepPosCheckTemp','posControlVsRunCorrTemp');

                posControlVsRunCorr(iDate,iRun,:,:) = posControlVsRunCorrTemp;
                superDeepPosCheck(iDate,iRun,:,:,:) = superDeepPosCheckTemp;
            end
            
        else
            clear controlVars
            controlVars                       = load([dataDir '\superDeepControls.mat']);
            posControlVsRunCorr(iDate,iRun,:,:) = controlVars.posControlVsRunCorrTemp;
            superDeepPosCheck(iDate,iRun,:,:,:) = controlVars.superDeepPosCheckTemp;
            runSuperDeepCorr                    = squeeze(super_DeepCorr(iDate,iRun,:,2,:))';
        end

        if ~exist([dataDir '\Super_DeepPosControls.png'],'file') && (chInCortex(1)-chInCortex(2) ~= 0)
            figure; plot(-200:200,runSuperDeepCorr(:,3),'r','LineWidth',1); hold on;
            plot(-200:200,squeeze(superDeepPosCheck(iDate,iRun,3,:,:)),'k'); ylim([0 1.2]);
            xlabel('Lags (s)'); ylabel('Correlation between superficial and deep hybrids');
            f = gcf; exportgraphics(f,[dataDir '\Super_DeepPosControls.png'],'Resolution',300); close gcf;
        end
    end
end

clear varInfo;
varFlag = 0;
try matfile(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat']);
    varInfo = who('-file',['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat']);
    if sum(ismember(varInfo,'superDeepPosCheckT'))==0
        varFlag = 1;
    end
catch
    varFlag = 1;
end
        
if varFlag
superDeepPosCheckT = reshape(superDeepPosCheck,[size(super_DeepCorr,1)*size(super_DeepCorr,2) 3 401 5]); 
superDeepPosCheckT(nanRow,:,:,:)           = []; 
superDeepPosCheckT(~(goodRunsSpatial & ~singleChFlag),:,:,:) = []; 
superDeepPosCheckT = squeeze(median(superDeepPosCheckT,4,'omitnan'));
superDeepPosCheckT = mean(superDeepPosCheckT(:,:,lagFrameRange),3,'omitnan'); 
 
    save(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat'],...
       'superDeepPosCheckT','-append');
end

%% Correlate between superficial/middle/deep layer compartments
% Initialize variables
lfpLayerCorr = NaN(size(probe,2),size(probe,1),3,5);
infraLayerCorr = NaN(size(probe,2),size(probe,1),3,5);
powerLayerCorr = NaN(size(probe,2),size(probe,1),3,5);

superInfraAllBands = NaN(size(probe,2),size(probe,1),5,5,6,6);
deepInfraAllBands  = NaN(size(probe,2),size(probe,1),5,5,9,9);
midInfraAllBands   = NaN(size(probe,2),size(probe,1),5,5,6,6);
superPowerAllBands = NaN(size(probe,2),size(probe,1),5,5);
deepPowerAllBands  = NaN(size(probe,2),size(probe,1),5,5);
midPowerAllBands   = NaN(size(probe,2),size(probe,1),5,5);

crossFreqCrossLayerInfra = NaN(size(probe,2),size(probe,1),5,5,3,3);
crossFreqCrossLayerPower = NaN(size(probe,2),size(probe,1),5,5,3,3);

% Get filter parameters...
fs = 1e3; 
gammaBand   = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); % Gamma band filtering parameters
alphaBand   = [8 12];  [bA,aA] = butter(3,alphaBand./(fs/2),'bandpass'); % Alpha band filtering parameters
betaBand   = [13 30]; [bB,aB] = butter(3,betaBand./(fs/2),'bandpass');  % Beta band filtering parameters
thetaBand  = [6 8];   [bT,aT] = butter(3,thetaBand./(fs/2),'bandpass'); % Theta band filtering parameters

tic;
for iDate = 1:size(allDates,1)
    clear expDate
    expDate = allDates(iDate,:);

    for iRun = 1:size(allRuns{iDate,1})
        clear runName dataDir probeCh rawCh lfpBadTimes lfpBadCh chInCortex 
      
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

        clc; disp(['Analyzing data for ' monkeyName ' '  expDate ' run: ' runName]);
        
        % Get run related variables
        probeCh     = probe{iRun,iDate}.probeCh;
        rawCh       = probe{iRun,iDate}.rawCh; 
        lfpBadTimes = badTimesLFP{iDate,iRun};
        lfpBadCh    = badCh{iDate,iRun};
        chInCortex  = estChInCortex{1,iDate}(iRun,:);
        
        probeCh(:,lfpBadCh)    = []; % Remove bad channels from LFP
        probeCh(lfpBadTimes,:) = []; % Remove bad times from LFP  

        rawCh(:,lfpBadCh)    = []; % Remove bad channels from raw
        rawCh(lfpBadTimes,:) = []; % Remove bad times from raw

        % Remove bad times determined visually from spectrogram
        [probeCh,rawCh,~] = removeBadTimesFromSpec(monkeyName,expDate,runName,probeCh,rawCh,[]);
      
        if size(probeCh,2)~=1
            % Get frequency band based information
            clear thetaLFP alphaLFP betaLFP gammaLFP infraEphysTheta infraEphysAlpha...
                infraEphysBeta infraEphysGamma infraEphysSpiking powerTheta powerAlpha...
                powerBeta powerGamma powerSpiking

            thetaLFP = single(filtfilt(bT,aT,double(probeCh(:,chInCortex(1):chInCortex(2)))));
            alphaLFP = single(filtfilt(bA,aA,double(probeCh(:,chInCortex(1):chInCortex(2)))));
            betaLFP  = single(filtfilt(bB,aB,double(probeCh(:,chInCortex(1):chInCortex(2)))));
            gammaLFP = single(filtfilt(bG,aG,double(probeCh(:,chInCortex(1):chInCortex(2)))));

            powerTheta   = envelope(abs(thetaLFP),5); 
            powerAlpha   = envelope(abs(alphaLFP),5); 
            powerBeta    = envelope(abs(betaLFP),5);
            powerGamma   = envelope(abs(gammaLFP),5); 
            powerSpiking = envelope(abs(rawCh(:,chInCortex(1):chInCortex(2))),5);

            infraEphysTheta   = getInfraSlowPowerLFP(probeCh,bT,aT,chInCortex);
            infraEphysAlpha   = getInfraSlowPowerLFP(probeCh,bA,aA,chInCortex);
            infraEphysBeta    = getInfraSlowPowerLFP(probeCh,bB,aB,chInCortex);
            infraEphysGamma   = getInfraSlowPowerLFP(probeCh,bG,aG,chInCortex);
            infraEphysSpiking = getInfraSlowPowerLFP(rawCh,[],[],chInCortex);

            probeCh  = single(probeCh(:,chInCortex(1):chInCortex(2)));
            rawCh    = single(rawCh(:,chInCortex(1):chInCortex(2)));

            % Separate into superficial/middle/deep layer compartments for the
            % different frequency bands
            clear lfpSuper lfpMid lfpMid thetaSuper thetaMid thetaDeep alphaSuper...
                alphaMid alphaDeep betaSuper betaMid betaDeep gammaSuper gammaMid...
                gammaDeep spikingSuper spikingMid spikingDeep infraThetaSuper...
                infraThetaMid infraThetaDeep infraAlphaSuper infraAlphaMid...
                infraAlphaDeep infraBetaSuper infraBetaMid infraBetaDeep...
                infraGammaSuper infraGammaMid infraGammaDeep infraSpikingSuper...
                infraSpikingMid infraSpikingDeep thetaPowerSuper thetaPowerMid...
                thetaPowerDeep alphaPowerSuper alphaPowerMid alphaPowerDeep...
                betaPowerSuper betaPowerMid betaPowerDeep gammaPowerSuper...
                gammaPowerMid gammaPowerDeep spikingPowerSuper spikingPowerMid...
                spikingPowerDeep

            thetaSuper = mean(thetaLFP(:,1:6),2,'omitnan');
            thetaMid   = mean(thetaLFP(:,7:12),2,'omitnan');
            thetaDeep  = mean(thetaLFP(:,13:end),2,'omitnan');

            alphaSuper = mean(alphaLFP(:,1:6),2,'omitnan');
            alphaMid   = mean(alphaLFP(:,7:12),2,'omitnan');
            alphaDeep  = mean(alphaLFP(:,13:end),2,'omitnan');

            betaSuper = mean(betaLFP(:,1:6),2,'omitnan');
            betaMid   = mean(betaLFP(:,7:12),2,'omitnan');
            betaDeep  = mean(betaLFP(:,13:end),2,'omitnan');

            gammaSuper = mean(gammaLFP(:,1:6),2,'omitnan');
            gammaMid   = mean(gammaLFP(:,7:12),2,'omitnan');
            gammaDeep  = mean(gammaLFP(:,13:end),2,'omitnan');

            spikingSuper = mean(rawCh(:,1:6),2,'omitnan');
            spikingMid   = mean(rawCh(:,7:12),2,'omitnan');
            spikingDeep  = mean(rawCh(:,13:end),2,'omitnan');

            thetaPowerSuper = mean(powerTheta(:,1:6),2,'omitnan');
            thetaPowerMid   = mean(powerTheta(:,7:12),2,'omitnan');
            thetaPowerDeep  = mean(powerTheta(:,13:end),2,'omitnan');

            alphaPowerSuper = mean(powerAlpha(:,1:6),2,'omitnan');
            alphaPowerMid   = mean(powerAlpha(:,7:12),2,'omitnan');
            alphaPowerDeep  = mean(powerAlpha(:,13:end),2,'omitnan');

            betaPowerSuper = mean(powerBeta(:,1:6),2,'omitnan');
            betaPowerMid   = mean(powerBeta(:,7:12),2,'omitnan');
            betaPowerDeep  = mean(powerBeta(:,13:end),2,'omitnan');

            gammaPowerSuper = mean(powerGamma(:,1:6),2,'omitnan');
            gammaPowerMid   = mean(powerGamma(:,7:12),2,'omitnan');
            gammaPowerDeep  = mean(powerGamma(:,13:end),2,'omitnan');

            spikingPowerSuper = mean(powerSpiking(:,1:6),2,'omitnan');
            spikingPowerMid   = mean(powerSpiking(:,7:12),2,'omitnan');
            spikingPowerDeep  = mean(powerSpiking(:,13:end),2,'omitnan');

            infraThetaSuper = mean(infraEphysTheta(:,1:6),2,'omitnan');
            infraThetaMid   = mean(infraEphysTheta(:,7:12),2,'omitnan');
            infraThetaDeep  = mean(infraEphysTheta(:,13:end),2,'omitnan');

            infraAlphaSuper = mean(infraEphysAlpha(:,1:6),2,'omitnan');
            infraAlphaMid   = mean(infraEphysAlpha(:,7:12),2,'omitnan');
            infraAlphaDeep  = mean(infraEphysAlpha(:,13:end),2,'omitnan');

            infraBetaSuper = mean(infraEphysBeta(:,1:6),2,'omitnan');
            infraBetaMid   = mean(infraEphysBeta(:,7:12),2,'omitnan');
            infraBetaDeep  = mean(infraEphysBeta(:,13:end),2,'omitnan');

            infraGammaSuper = mean(infraEphysGamma(:,1:6),2,'omitnan');
            infraGammaMid   = mean(infraEphysGamma(:,7:12),2,'omitnan');
            infraGammaDeep  = mean(infraEphysGamma(:,13:end),2,'omitnan');

            infraSpikingSuper = mean(infraEphysSpiking(:,1:6),2,'omitnan');
            infraSpikingMid   = mean(infraEphysSpiking(:,7:12),2,'omitnan');
            infraSpikingDeep  = mean(infraEphysSpiking(:,13:end),2,'omitnan');

            infraThetaSuperAllCh = (infraEphysTheta(:,1:6));
            infraThetaMidAllCh   = (infraEphysTheta(:,7:12));
            infraThetaDeepAllCh  = (infraEphysTheta(:,13:end));

            infraAlphaSuperAllCh = (infraEphysAlpha(:,1:6));
            infraAlphaMidAllCh   = (infraEphysAlpha(:,7:12));
            infraAlphaDeepAllCh  = (infraEphysAlpha(:,13:end));

            infraBetaSuperAllCh = (infraEphysBeta(:,1:6));
            infraBetaMidAllCh   = (infraEphysBeta(:,7:12));
            infraBetaDeepAllCh  = (infraEphysBeta(:,13:end));

            infraGammaSuperAllCh = (infraEphysGamma(:,1:6));
            infraGammaMidAllCh   = (infraEphysGamma(:,7:12));
            infraGammaDeepAllCh  = (infraEphysGamma(:,13:end));

            infraSpikingSuperAllCh = (infraEphysSpiking(:,1:6));
            infraSpikingMidAllCh   = (infraEphysSpiking(:,7:12));
            infraSpikingDeepAllCh  = (infraEphysSpiking(:,13:end));


            % Correlate between superficial/middle/deep layers
            % Time-series
            lfpLayerCorr(iDate,iRun,1,1) = corr(thetaSuper,thetaDeep,'rows','complete');
            lfpLayerCorr(iDate,iRun,2,1) = corr(thetaSuper,thetaMid,'rows','complete');
            lfpLayerCorr(iDate,iRun,3,1) = corr(thetaMid,thetaDeep,'rows','complete');

            lfpLayerCorr(iDate,iRun,1,2) = corr(alphaSuper,alphaDeep,'rows','complete');
            lfpLayerCorr(iDate,iRun,2,2) = corr(alphaSuper,alphaMid,'rows','complete');
            lfpLayerCorr(iDate,iRun,3,2) = corr(alphaMid,alphaDeep,'rows','complete');

            lfpLayerCorr(iDate,iRun,1,3) = corr(betaSuper,betaDeep,'rows','complete');
            lfpLayerCorr(iDate,iRun,2,3) = corr(betaSuper,betaMid,'rows','complete');
            lfpLayerCorr(iDate,iRun,3,3) = corr(betaMid,betaDeep,'rows','complete');

            lfpLayerCorr(iDate,iRun,1,4) = corr(gammaSuper,gammaDeep,'rows','complete');
            lfpLayerCorr(iDate,iRun,2,4) = corr(gammaSuper,gammaMid,'rows','complete');
            lfpLayerCorr(iDate,iRun,3,4) = corr(gammaMid,gammaDeep,'rows','complete');

            lfpLayerCorr(iDate,iRun,1,5) = corr(spikingSuper,spikingDeep,'rows','complete');
            lfpLayerCorr(iDate,iRun,2,5) = corr(spikingSuper,spikingMid,'rows','complete');
            lfpLayerCorr(iDate,iRun,3,5) = corr(spikingMid,spikingDeep,'rows','complete');

            [gammaLayerXCorr(iDate,iRun,:,1),lags] = xcorr(detrend(gammaSuper,0),detrend(gammaDeep,0),1000,'coeff');
            [gammaLayerXCorr(iDate,iRun,:,2),~]    = xcorr(detrend(gammaSuper,0),detrend(gammaMid,0),1000,'coeff');
            [gammaLayerXCorr(iDate,iRun,:,3),~]    = xcorr(detrend(gammaMid,0),detrend(gammaDeep,0),1000,'coeff');
%              
%             clear superH midH deepH
%             superH = hilbert(gammaSuper); midH = hilbert(gammaMid); deepH = hilbert(gammaDeep); 
%             phaseRad{iDate,iRun}(1) = angle(superH./deepH);
%             phaseRad{iDate,iRun}(2) = angle(superH./midH);
%             phaseRad{iDate,iRun}(3)= angle(midH./deepH); 

            % Power
            powerLayerCorr(iDate,iRun,1,1) = corr(thetaPowerSuper,thetaPowerDeep,'rows','complete');
            powerLayerCorr(iDate,iRun,2,1) = corr(thetaPowerSuper,thetaPowerMid,'rows','complete');
            powerLayerCorr(iDate,iRun,3,1) = corr(thetaPowerMid,thetaPowerDeep,'rows','complete');

            powerLayerCorr(iDate,iRun,1,2) = corr(alphaPowerSuper,alphaPowerDeep,'rows','complete');
            powerLayerCorr(iDate,iRun,2,2) = corr(alphaPowerSuper,alphaPowerMid,'rows','complete');
            powerLayerCorr(iDate,iRun,3,2) = corr(alphaPowerMid,alphaPowerDeep,'rows','complete');

            powerLayerCorr(iDate,iRun,1,3) = corr(betaPowerSuper,betaPowerDeep,'rows','complete');
            powerLayerCorr(iDate,iRun,2,3) = corr(betaPowerSuper,betaPowerMid,'rows','complete');
            powerLayerCorr(iDate,iRun,3,3) = corr(betaPowerMid,betaPowerDeep,'rows','complete');

            powerLayerCorr(iDate,iRun,1,4) = corr(gammaPowerSuper,gammaPowerDeep,'rows','complete');
            powerLayerCorr(iDate,iRun,2,4) = corr(gammaPowerSuper,gammaPowerMid,'rows','complete');
            powerLayerCorr(iDate,iRun,3,4) = corr(gammaPowerMid,gammaPowerDeep,'rows','complete');

            powerLayerCorr(iDate,iRun,1,5) = corr(spikingPowerSuper,spikingPowerDeep,'rows','complete');
            powerLayerCorr(iDate,iRun,2,5) = corr(spikingPowerSuper,spikingPowerMid,'rows','complete');
            powerLayerCorr(iDate,iRun,3,5) = corr(spikingPowerMid,spikingPowerDeep,'rows','complete');

            [gammaPowerLayerXCorr(iDate,iRun,:,1),~] = xcorr(detrend(gammaPowerSuper,0),detrend(gammaPowerDeep,0),1000,'coeff');
            [gammaPowerLayerXCorr(iDate,iRun,:,2),~] = xcorr(detrend(gammaPowerSuper,0),detrend(gammaPowerMid,0),1000,'coeff');
            [gammaPowerLayerXCorr(iDate,iRun,:,3),~] = xcorr(detrend(gammaPowerMid,0),detrend(gammaPowerDeep,0),1000,'coeff');
            
%             clear superH midH deepH
%             superH = hilbert(gammaPowerSuper); midH = hilbert(gammaPowerMid); deepH = hilbert(gammaPowerDeep);
%             phaseRadPower(iDate,iRun,:,1) = angle(superH./deepH);
%             phaseRadPower(iDate,iRun,:,2) = angle(superH./midH);
%             phaseRadPower(iDate,iRun,:,3) = angle(midH./deepH);

            % Infraslow power
            infraLayerCorr(iDate,iRun,1,1) = corr(infraThetaSuper,infraThetaDeep,'rows','complete');
            infraLayerCorr(iDate,iRun,2,1) = corr(infraThetaSuper,infraThetaMid,'rows','complete');
            infraLayerCorr(iDate,iRun,3,1) = corr(infraThetaMid,infraThetaDeep,'rows','complete');

            infraLayerCorr(iDate,iRun,1,2) = corr(infraAlphaSuper,infraAlphaDeep,'rows','complete');
            infraLayerCorr(iDate,iRun,2,2) = corr(infraAlphaSuper,infraAlphaMid,'rows','complete');
            infraLayerCorr(iDate,iRun,3,2) = corr(infraAlphaMid,infraAlphaDeep,'rows','complete');

            infraLayerCorr(iDate,iRun,1,3) = corr(infraBetaSuper,infraBetaDeep,'rows','complete');
            infraLayerCorr(iDate,iRun,2,3) = corr(infraBetaSuper,infraBetaMid,'rows','complete');
            infraLayerCorr(iDate,iRun,3,3) = corr(infraBetaMid,infraBetaDeep,'rows','complete');

            infraLayerCorr(iDate,iRun,1,4) = corr(infraGammaSuper,infraGammaDeep,'rows','complete');
            infraLayerCorr(iDate,iRun,2,4) = corr(infraGammaSuper,infraGammaMid,'rows','complete');
            infraLayerCorr(iDate,iRun,3,4) = corr(infraGammaMid,infraGammaDeep,'rows','complete');

            infraLayerCorr(iDate,iRun,1,5) = corr(infraSpikingSuper,infraSpikingDeep,'rows','complete');
            infraLayerCorr(iDate,iRun,2,5) = corr(infraSpikingSuper,infraSpikingMid,'rows','complete');
            infraLayerCorr(iDate,iRun,3,5) = corr(infraSpikingMid,infraSpikingDeep,'rows','complete');

            [infraGammaLayerXCorr(iDate,iRun,:,1),~] = xcorr(detrend(infraGammaSuper,0),detrend(infraGammaDeep,0),1000,'coeff');
            [infraGammaLayerXCorr(iDate,iRun,:,2),~] = xcorr(detrend(infraGammaSuper,0),detrend(infraGammaMid,0),1000,'coeff');
            [infraGammaLayerXCorr(iDate,iRun,:,3),~] = xcorr(detrend(infraGammaMid,0),detrend(infraGammaDeep,0),1000,'coeff');

%             clear superH midH deepH
%             superH = hilbert(infraGammaSuper); midH = hilbert(infraGammaMid); deepH = hilbert(infraGammaDeep);
%             phaseRadInfra(iDate,iRun,:,1) = angle(superH./deepH);
%             phaseRadInfra(iDate,iRun,:,2) = angle(superH./midH);
%             phaseRadInfra(iDate,iRun,:,3) = angle(midH./deepH);

            % Correlate between frequencies given a layer compartment
            clear chLen
            chLen = estChInCortex{1,iDate}(iRun,2) - estChInCortex{1,iDate}(iRun,1);
            for iBand1 = 1:5
                clear super1 deep1 mid1 superPower1 midPower1 deepPower1
                switch iBand1
                    case 1
                        super1 = infraThetaSuperAllCh;
                        deep1  = infraThetaDeepAllCh;

                        superPower1 = thetaPowerSuper;
                        deepPower1  = thetaPowerDeep;

                        if  chLen~=0
                            mid1      = infraThetaMidAllCh;
                            midPower1 = thetaPowerMid;
                        end

                    case 2
                        super1 = infraAlphaSuperAllCh;
                        deep1  = infraAlphaDeepAllCh;

                        superPower1 = alphaPowerSuper;
                        deepPower1 = alphaPowerDeep;

                        if chLen~=0
                            mid1 = infraAlphaMidAllCh;
                            midPower1 = alphaPowerMid;
                        end

                    case 3
                        super1 = infraBetaSuperAllCh;
                        deep1  = infraBetaDeepAllCh;

                        superPower1 = betaPowerSuper;
                        deepPower1  = betaPowerDeep;

                        if chLen~=0
                            mid1      = infraBetaMidAllCh;
                            midPower1 = betaPowerMid;
                        end

                    case 4
                        super1 = infraGammaSuperAllCh;
                        deep1  = infraGammaDeepAllCh;

                        superPower1 = gammaPowerSuper;
                        deepPower1  = gammaPowerDeep;

                        if  chLen~=0
                            mid1      = infraGammaMidAllCh;
                            midPower1 = gammaPowerMid;
                        end

                    case 5
                        super1 = infraSpikingSuperAllCh;
                        deep1  = infraSpikingDeepAllCh;

                        superPower1 = spikingPowerSuper;
                        deepPower1  = spikingPowerDeep;

                        if  chLen~=0
                            mid1      = infraSpikingMidAllCh;
                            midPower1 = spikingPowerMid;
                        end
                end

                for iBand2 = 1:5
                    clear super2 deep2 mid2  superPower2 midPower2 deepPower2
                    switch iBand2
                        case 1
                            super2 = infraThetaSuperAllCh;
                            deep2  = infraThetaDeepAllCh;

                            superPower2 = thetaPowerSuper;
                            deepPower2  = thetaPowerDeep;

                            if  chLen~=0
                                mid2      = infraThetaMidAllCh;
                                midPower2 = thetaPowerMid;
                            end

                        case 2
                            super2 = infraAlphaSuperAllCh;
                            deep2  = infraAlphaDeepAllCh;

                            superPower2 = alphaPowerSuper;
                            deepPower2  = alphaPowerDeep;

                            if chLen~=0
                                mid2      = infraAlphaMidAllCh;
                                midPower2 = alphaPowerMid;
                            end

                        case 3
                            super2 = infraBetaSuperAllCh;
                            deep2  = infraBetaDeepAllCh;

                            superPower2 = betaPowerSuper;
                            deepPower2  = betaPowerDeep;

                            if chLen~=0
                                mid2      = infraBetaMidAllCh;
                                midPower2 = betaPowerMid;
                            end

                        case 4
                            super2 = infraGammaSuperAllCh;
                            deep2  = infraGammaDeepAllCh;

                            superPower2 = gammaPowerSuper;
                            deepPower2  = gammaPowerDeep;

                            if  chLen~=0
                                mid2      = infraGammaMidAllCh;
                                midPower2 = gammaPowerMid;
                            end

                        case 5
                            super2 = infraSpikingSuperAllCh;
                            deep2  = infraSpikingDeepAllCh;

                            superPower2 = spikingPowerSuper;
                            deepPower2  = spikingPowerDeep;

                            if chLen~=0
                                mid2      = infraSpikingMidAllCh;
                                midPower2 = spikingPowerMid;
                            end
                    end

                    superInfraAllBands(iDate,iRun,iBand1,iBand2,:,:) = ((corr(super1,super2,'rows','complete')));
                    deepInfraAllBands(iDate,iRun,iBand1,iBand2,:,:)  = ((corr(deep1,deep2,'rows','complete')));

                    superPowerAllBands(iDate,iRun,iBand1,iBand2) = corr(superPower1,superPower2,'rows','complete');
                    deepPowerAllBands(iDate,iRun,iBand1,iBand2)  = corr(deepPower1,deepPower2,'rows','complete');

                    if  chLen~=0
                        midInfraAllBands(iDate,iRun,iBand1,iBand2,:,:) = ((corr(mid1,mid2,'rows','complete')));
                        midPowerAllBands(iDate,iRun,iBand1,iBand2) = corr(midPower1,midPower2,'rows','complete');
                    end

%                     if iBand1~=iBand2
%                         crossFreqCrossLayerInfra(iDate,iRun,iBand1,iBand2,1,1) = corr(super1,super2,'rows','complete');                         
%                         crossFreqCrossLayerInfra(iDate,iRun,iBand1,iBand2,1,3) = corr(super1,deep2,'rows','complete');                        
%                         crossFreqCrossLayerInfra(iDate,iRun,iBand1,iBand2,3,1) = corr(deep1,super2,'rows','complete');                      
%                         crossFreqCrossLayerInfra(iDate,iRun,iBand1,iBand2,3,3) = corr(deep1,deep2,'rows','complete'); 
% 
%                         crossFreqCrossLayerPower(iDate,iRun,iBand1,iBand2,1,1) = corr(superPower1,superPower2,'rows','complete');
%                         crossFreqCrossLayerPower(iDate,iRun,iBand1,iBand2,1,3) = corr(superPower1,deepPower2,'rows','complete');
%                         crossFreqCrossLayerPower(iDate,iRun,iBand1,iBand2,3,1) = corr(deepPower1,superPower2,'rows','complete');
%                         crossFreqCrossLayerPower(iDate,iRun,iBand1,iBand2,3,3) = corr(deepPower1,deepPower2,'rows','complete');
% 
%                          if  chLen~=0
%                              crossFreqCrossLayerInfra(iDate,iRun,iBand1,iBand2,1,2) = corr(super1,mid2,'rows','complete');
%                              crossFreqCrossLayerInfra(iDate,iRun,iBand1,iBand2,2,1) = corr(mid1,super2,'rows','complete');
%                              crossFreqCrossLayerInfra(iDate,iRun,iBand1,iBand2,2,2) = corr(mid1,mid2,'rows','complete');
%                              crossFreqCrossLayerInfra(iDate,iRun,iBand1,iBand2,2,3) = corr(mid1,deep2,'rows','complete');
%                              crossFreqCrossLayerInfra(iDate,iRun,iBand1,iBand2,3,2) = corr(deep1,mid2,'rows','complete');
% 
%                              crossFreqCrossLayerPower(iDate,iRun,iBand1,iBand2,1,2) = corr(superPower1,midPower2,'rows','complete');
%                              crossFreqCrossLayerPower(iDate,iRun,iBand1,iBand2,2,1) = corr(midPower1,superPower2,'rows','complete');
%                              crossFreqCrossLayerPower(iDate,iRun,iBand1,iBand2,2,2) = corr(midPower1,midPower2,'rows','complete');
%                              crossFreqCrossLayerPower(iDate,iRun,iBand1,iBand2,2,3) = corr(midPower1,deepPower2,'rows','complete');
%                              crossFreqCrossLayerPower(iDate,iRun,iBand1,iBand2,3,2) = corr(deepPower1,midPower2,'rows','complete');
%                          end
%                     end 
                end
            end
          
          
        else
            lfpLayerCorr(iDate,iRun,:,:)   = 0;
            infraLayerCorr(iDate,iRun,:,:) = 0;
            powerLayerCorr(iDate,iRun,:,:) = 0;

            superInfraAllBands(iDate,iRun,:,:) = 0;
            deepInfraAllBands(iDate,iRun,:,:)  = 0;
            midInfraAllBands(iDate,iRun,:,:)   = 0;
            superPowerAllBands(iDate,iRun,:,:) = 0;
            deepPowerAllBands(iDate,iRun,:,:)  = 0;
            midPowerAllBands(iDate,iRun,:,:)   = 0;
% 
%             crossFreqCrossLayerInfra(iDate,iRun,:,:,:,:) = 0;
%             crossFreqCrossLayerPower(iDate,iRun,:,:,:,:) = 0;
            
        end
    end
end
toc;
%%

clear varInfo;
varFlag = 0;
try matfile(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat']);
    varInfo = who('-file',['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat']);
    if sum(ismember(varInfo,'superInfraAllBandsT'))==0
        varFlag = 1;
    end
catch
    varFlag = 1;
end
        
if varFlag
    superInfraAllBandsT = reshape(superInfraAllBands,[size(probe,2)*size(probe,1) 5 5 6 6]);
    deepInfraAllBandsT  = reshape(deepInfraAllBands,[size(probe,2)*size(probe,1) 5 5 9 9]);
    midInfraAllBandsT   = reshape(midInfraAllBands,[size(probe,2)*size(probe,1) 5 5 6 6]);

    % crossFreqCrossLayerInfraT = reshape(crossFreqCrossLayerInfra,[size(probe,2)*size(probe,1) 5 5 3 3]);
    % crossFreqCrossLayerPowerT = reshape(crossFreqCrossLayerPower,[size(probe,2)*size(probe,1) 5 5 3 3]);
    lfpLayerCorrT   = reshape(lfpLayerCorr,[size(probe,2)*size(probe,1) 3 5]);
    nanLocs = (isnan(lfpLayerCorrT(:,1,1)));

    lfpLayerCorrT(nanLocs,:,:)         = [];
    superInfraAllBandsT(nanLocs,:,:,:,:)   = [];
    deepInfraAllBandsT(nanLocs,:,:,:,:)    = [];
    midInfraAllBandsT(nanLocs,:,:,:,:)     = [];

    superInfraAllBandsT(~(goodRunsSpatial & ~singleChFlag),:,:,:,:)   = [];
    deepInfraAllBandsT(~(goodRunsSpatial & ~singleChFlag),:,:,:,:)    = [];
    midInfraAllBandsT(~(goodRunsSpatial & ~singleChFlag),:,:,:,:)     = [];

    save(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat'],...
        'superInfraAllBandsT','deepInfraAllBandsT','midInfraAllBandsT','-append');
end


superInfraAllBandsTMed = median(superInfraAllBandsT,[4 5],'omitnan');
deepInfraAllBandsTMed  = median(deepInfraAllBandsT,[4 5],'omitnan');
midInfraAllBandsTMed   = median(midInfraAllBandsT,[4 5],'omitnan');

for iBand = 1:5
    figure;
    subplot(131); boxplot(squeeze(superInfraAllBandsTMed(:,iBand,:)),bandLabels);
    ylim([-0.5 1.2]); title('Superficial - infraslow power'); box off;

    subplot(132); boxplot(squeeze(deepInfraAllBandsTMed(:,iBand,:)),bandLabels);
    ylim([-0.5 1.2]); title('Middle - infraslow power'); box off;

    subplot(133); boxplot(squeeze(midInfraAllBandsTMed(:,iBand,:)),bandLabels);
    ylim([-0.5 1.2]);title('Deep - infraslow power'); box off;      
    
end

%%
dMat = logical(diag([1 1 1 1 1])); dMat = reshape(dMat,[25 1]);

ss = squeeze(reshape(superInfraAllBandsT,[size(superInfraAllBandsT,1) 25 36]));
ss = squeeze(mean(ss,1,'omitnan'));
ss(~dMat,:) = [];
ss = reshape(ss,[5 6 6]);

mm = squeeze(reshape(midInfraAllBandsT,[size(superInfraAllBandsT,1) 25 36]));
mm = squeeze(mean(mm,1,'omitnan'));
mm(~dMat,:) = [];
mm = reshape(mm,[5 6 6]);

dd = squeeze(reshape(deepInfraAllBandsT,[size(superInfraAllBandsT,1) 25 81]));
dd = squeeze(mean(dd,1,'omitnan'));
dd(~dMat,:) = [];
dd = reshape(dd,[5 9 9]);

figure;
for iB = 1:5
    subplot(5,3,3*(iB-1)+1); imagesc(squeeze(ss(iB,:,:)));
    colormap jet; caxis([0 1]); colorbar; axis square; 
    ylabel(bandLabels{iB});
    if iB == 1; title('Superficial'); end

    subplot(5,3,3*(iB-1)+2);imagesc(squeeze(mm(iB,:,:)));
    colormap jet; caxis([0 1]); colorbar;  axis square;
     if iB == 1; title('Middle'); end


    subplot(5,3,3*(iB-1)+3);imagesc(squeeze(dd(iB,:,:)));
    colormap jet; caxis([0 1]); colorbar;  axis square;
     if iB == 1; title('Deep'); end
end

%% Calculating EEG powers
% Get filter parameters...
fs = 1e3; 
gammaBand   = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); % Gamma band filtering parameters
params.Fs       = fs;
params.fpass    = [1 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;

for iDate = 1:size(allDates,1)
    clear expDate
    expDate = allDates(iDate,:);

    for iRun = 1:size(allRuns{iDate,1})
        clear runName dataDir probeCh rawCh lfpBadTimes lfpBadCh chInCortex 
      
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

        clc; disp(['Calculating powers and correlations for  ' monkeyName ' '  expDate ' run: ' runName]);
        
        % Get run related variables
        probeCh     = probe{iRun,iDate}.probeCh; 
        eegCh       = probe{iRun,iDate}.eegCh;
        lfpBadTimes = badTimesLFP{iDate,iRun};
        probeCh(lfpBadTimes,:) = [];
        chInCortex   = estChInCortex{1,iDate}(iRun,:);
        gammaCh = filtfilt(bG,aG,double(probeCh(:,chInCortex(1):chInCortex(2))));
        gIntraProbeCorr(iDate,iRun) = median(corr(gammaCh),'all'); 
        
        eegCh(lfpBadTimes,:)   = []; % Remove bad times from LFP
        [spec,timeValsSpec,freqValsSpec] = mtspecgramc(eegCh,[5 2],params); 

        % Get gamma power
        gIdx = freqValsSpec>=30 & freqValsSpec<=90;
        gPower(iDate,iRun) = 10.*log10(median(spec(:,gIdx),'all'));
     end
end

gPower = reshape(gPower,[size(gPower,1)*size(gPower,2) 1]);
gPower(gPower==0) = [];  gPower(~goodRuns) = [];

gIntraProbeCorr = reshape(gIntraProbeCorr,[size(gIntraProbeCorr,1)*size(gIntraProbeCorr,2) 1]);
gIntraProbeCorr(gIntraProbeCorr==0) = [];  gIntraProbeCorr(~goodRuns) = [];

figure; plot(isoLevelGoodRuns,gPower,'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]);
figure; plot(isoLevelGoodRuns,gIntraProbeCorr,'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]);


%% Calculate FC maps for seeds placed <1mm near the probe - Full FOV controls
%%% and correlate with hybrid map
x = -200:200;
negIdx = (-100<=x)&(x<=0); negVals = x(negIdx);

for iDate = 1:size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);
    for iRun = 1: size(allRuns{iDate,1},1)
        clear runName dataDir clipMask elecMask clipMaskCortex corrMask lowIdx ...
            pDatTemp greenFig seedLocIn crossCorrTemp allHybridMaps mapsAll...
            peakNegHybridMap mapsAllTemp probeCh badTimes szLFP skullMask ...
            inDatSize infraEphys allCortexMask

        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
        serverDir = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\kem294\Data\' monkeyName '_SqM\' ...
            hemisphere ' Hemisphere\' expDate '\' runName];
        % IMAGING: Load the appropriate masks for the imaging data
        % Load clipmask
        if exist([dataDir '\clipMask0' runName(end) '.BMP'],'file') == 0
            clipMask = imread([dataDir '\clipMask0' runName(end) '.png']);
        else
            clipMask = imread([dataDir '\clipMask0' runName(end) '.bmp']);
        end

        % Load the cortex mask
        if exist([dataDir '\skullMask.bmp'],'file') == 0
            skullMask = imread([dataDir '\skullMask.png']);
        else
            skullMask = imread([dataDir '\skullMask.bmp']);
        end

        % This mask includes vessels and electrode(s)
        if exist([dataDir '\maskSkull0' runName(end) '.bmp'],'file') == 0
            elecMask = imread([dataDir '\maskSkull0' runName(end) '.png']);
        else
            elecMask = imread([dataDir '\maskSkull0' runName(end) '.bmp']);
        end

        clipMask       = imresize(clipMask,1/3); % Resize mask
        clipMask       = clipMask(:,:,1)>0; % Converting to 0s and 1s
        elecMask       = imresize(elecMask-100,1/3); % Binarize image
        elecMask       = elecMask(:,:,1)>0;
        skullMask      = imresize(skullMask,1/3);
        skullMask      = skullMask(:,:,1)>0; 
        allCortexMask  = skullMask & ~elecMask; % Mask of cortex with vessels
        clipMaskCortex = clipMask & ~elecMask; % Mask of cortex without vessels
        
        % Mask for comparing the imaging and hybrid map (after removing
        % areas beyond lateral sulcus)
        if exist([dataDir '\corrMask0' runName(end) '.BMP'],'file')
            corrMask = imread([dataDir '\corrMask0' runName(end) '.bmp']);
            corrMaskFlag = 1;
        elseif exist([dataDir '\corrMask0' runName(end) '.PNG'],'file')
            corrMask = imread([dataDir '\corrMask0' runName(end) '.png']);
            corrMaskFlag = 1;
        else
            corrMaskFlag = 0;
        end

        if corrMaskFlag
            corrMask = imresize(corrMask,1/3); % Resize the mask
            corrMask = corrMask(:,:,1)>0;
            corrMask = corrMask & ~elecMask;
        else
            corrMask = clipMaskCortex;
        end

        % Get imaging data for the run
        pDatTemp = processedDat{iDate,iRun}.tempBandPass;
        imSize   = size(pDatTemp);
        greenFig = imresize(greenIm{iDate,iRun},1/spatialBin,'OutputSize',[imSize(1) imSize(2)]);
        
        % Get ROI location and FC map
        seedLocIn    = load([dataDir '\roiCenterLoc.mat']);
        seedLocProbe = seedLocIn.seedLocProbe;
        seedLocIn    = seedLocIn.seedLocIn;
        circleRad    = 6;

        seedSigT     = calculateSeedSignal(greenFig,corrMask,...
            seedLocIn,circleRad,pDatTemp); % Get Gaussian weighted seed signal

        fcMap            = plotCorrMap(seedSigT,pDatTemp,0);
        corrMaskT        = reshape(corrMask,[imSize(1)*imSize(2) 1]);
        fcMap            = reshape(fcMap,[361*438 1]);
        fcMap(~corrMaskT) = NaN;

        % Get the hybrid maps and determine lag  for the run
        try
            allHybridMaps  = matfile([serverDir '\crossCorrFOV_6_NoRef.mat']);
        catch
            getCrossCorrFOV(monkeyName,expDate,runName,serverDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp);

            allHybridMaps = matfile([serverDir '\crossCorrFOV_6_NoRef.mat']);

        end

        mapsAllTemp           = allHybridMaps.spatialProfile;
        mapsAll               = mapsAllTemp.ccFull;
        mapsAll               = reshape(mapsAll,[401 361*438]);
        mapsAll(:,~corrMaskT) = NaN;
%         peakNegHybridMap      = squeeze(mapsAll((x == peakNegTimesAll(iDate,iRun,3,2).*10),:));

        clc; disp(['Analyzing data for ' monkeyName ' '  expDate ' run: ' runName]);

%         %%% 3. TEMPORAL CONTROL: Shuffle Ephys using variable windows and
%         %%% determine the correlation at peak negative
% 
%         if ~exist([dataDir '\tempControlVarsFOV.mat'],'file') 
%             % Get infraslow ephys powers and imaging data
%             clear  runWiseCorrShuffled runWiseLagShuffled gammaEphys probeCh...
%                 processedDat10 ch badChannels badTimes szLFP timeStamp...
%                 badTimeThreshTemp badTimes
% 
%             disp('Obtaining temporal controls...');
% 
%             % Upsampling imaging data to 10 Hz           
%             pDatTemp   = reshape(pDatTemp,[imSize(1)*imSize(2) imSize(3)]);
%             
%             parfor iP = 1:size(pDatTemp,1)
%                 processedDat10(iP,:) = interp(pDatTemp(iP,:),5);
%             end
% 
%             szIm = size(processedDat10,2)*100;
% 
% 
%             % Get ephys data
%             probeCh               = probe{iRun,iDate}.probeCh;
%             ch                     = estChInCortex{iDate}(iRun,:);
%             badChannels            = badCh{iDate,iRun};
%             badTimes               = badTimesLFP{iDate,iRun};
%             probeCh(:,badChannels) = [];
%             szLFP                  = size(probeCh,1);
% 
%             % Make both matrices equal...
%             badTimeThreshTemp = badTimeThresh{iDate,iRun};
% 
%             if ~(szLFP == szIm)
%                 szMin          = min([szLFP, szIm]);
%                 probeCh        = probeCh(1:szMin,:);
%                 processedDat10 = processedDat10(:,1:floor(szMin/100));
% 
%                 badTimes(badTimes>szMin)           = [];
%                 badTimeThreshTemp(badTimeThreshTemp>szMin) = [];
%             else
%                 szMin = szLFP;
%             end
% 
%             timeStamp       = probe{iRun,iDate}.timeStamp;
%             timeStampSorted = timeStamp- timeStamp(1);
%             badTimes10Hz    = unique(badTimeThreshTemp./1000);
%             badTimeIm       = [];
% 
%             % Identifying frames to be removed from RS-ISOI
%             for iT = 1: length(badTimes10Hz)
%                 badTimeIm(iT) = find((floor(abs(timeStampSorted - badTimes10Hz(iT))*100)./100)<=0.05,1,'first');
%             end
% 
%             % Remove bad times determined from LFP
%             badTimeIm = unique(badTimeIm);
%             badTimeIm(badTimeIm>size(processedDat10,2)) = [];
%             processedDat10(:,badTimeIm) = [];
%             probeCh(badTimes,:) = [];
% 
%             % Remove bad times determined visually from spectrogram
%             [probeCh,~,processedDat10] = removeBadTimesFromSpec(monkeyName,expDate,runName,probeCh,[],processedDat10);
% 
%             % Get infraslow powers
%             gammaBand  = [30 90]; [bG,aG] = butter(3,gammaBand./(1e3/2),'bandpass');% Gamma band filtering parameters
% %             gammaEphys = single(filtfilt(bG,aG,double(probeCh(:,ch(1):ch(2)))));
% 
%             % Set the time windows to perform cross correlations
%             % between ISOI and Ephys
%             tic;
%             timeLen = min([size(gammaEphys,1) size(processedDat10,2)*1e2]);
%             winLen  = [0.001 0.1 1 5 10 50 100 500 timeLen./1e3].*1e3;
% %
%             for iWin = 1:length(winLen) % shuffle time series in windows of 10s
%                 clear comb1 infraEphysS roiS ccT ccFull z p k sos...
%                     enSize envelopeFiltered gammaEphysS
%                 for iRep = 1: 5
%                     clear comb1 infraEphysS roiS ccT ccFull z p k sos...
%                         enSize envelopeFiltered gammaEphysS envelopeDat
%                     rng('shuffle');
%                     comb1 = randperm(round(timeLen/winLen(iWin)));
% 
%                     if iWin == 1
%                         gammaEphysS  = single(filtfilt(bG,aG,double(probeCh(comb1,ch(1):ch(2)))));%gammaEphys(comb1,:);
%                     elseif iWin == 9
%                         gammaEphysS  = single(filtfilt(bG,aG,double(probeCh(:,ch(1):ch(2)))));
%                     else
%                         probeEphysS = [];
%                         for iL = 1:length(comb1)
%                             clear win1
%                             win1 = (comb1(iL)-1)*winLen(iWin)+1 : (comb1(iL)-1)*winLen(iWin)+winLen(iWin);
%                             win1(win1>timeLen) = [];
%                             probeEphysS = [probeEphysS; probeCh(win1,ch(1):ch(2))];
%                         end
%                          gammaEphysS  = single(filtfilt(bG,aG,double(probeEphysS)));
%                     end
% 
%                     % Get infraslow ephys
%                     envelopeDat = envelope(abs(gammaEphysS),5);
% 
%                     % Bandpass - 0.01 Hz - 0.1 Hz
%                     [z,p,k] = butter(3,[0.01 0.1]./(1e3/2),'bandpass');
%                     [sos,g] = zp2sos(z,p,k);
%                     enSize  = size(envelopeDat);
% 
%                     envelopeFiltered = filtfilt(sos,g,double([envelopeDat; envelopeDat; envelopeDat ]));
%                     envelopeFiltered = envelopeFiltered(enSize(1)+1:(end-enSize(1)),:);
%                     infraEphysS      = mean(single(downsample(envelopeFiltered,100)),2);
% 
%                     % Check size of imaging and ephys after this operation
%                     clear processedDat10R szMin
%                     szIm  = size(processedDat10,2);
%                     szLFP = size(infraEphysS,1);
%                     if ~(szLFP == szIm)
%                         szMin        = min([  szLFP, szIm]);
%                         infraEphysS   = infraEphysS(1:szMin,:);
%                         processedDat10R = processedDat10(:,1:szMin);
%                     else
%                         processedDat10R = processedDat10;
%                     end
% 
%                     parfor iP = 1:size(processedDat10R,1)
%                         if iP == 1
%                             [ccFull(:,iP),lagFull(:,iP)]  = xcorr(infraEphysS',processedDat10R(iP,:),200,'normalized');
%                         else
%                             [ccFull(:,iP),~]      = xcorr(infraEphysS',processedDat10R(iP,:),200,'normalized');
%                         end
%                     end
% 
%                     ccFull(:,~corrMaskT) = NaN;
% 
%                     for iMap = 1:401; mapVals(iMap) = corr(fcMap,ccFull(iMap,:)','rows','complete'); end
% 
%                     [runWiseCorrShuffled(iWin,iRep),runWiseLagShuffled(iWin,iRep)] = min(mapVals(negIdx));
%                     runWiseLagShuffled(iWin,iRep) = negVals(runWiseLagShuffled(iWin,iRep))./10;
%                 end
%             end
%             
%             toc;
%             runWiseCorrAllShuffled{iDate,iRun} = runWiseCorrShuffled;
%             runWiseLagAllShuffled{iDate,iRun}  = runWiseLagShuffled;
%             corrNegShuffle(iDate,iRun,:)       = median(runWiseCorrShuffled,2,'omitnan');%runWiseCorrShuffled;% 
%             corrNegTimes(iDate,iRun,:)         = median(runWiseLagShuffled,2,'omitnan');%runWiseLagShuffled;%
% 
%             save([dataDir '\tempControlVarsFOV.mat'],'runWiseCorrShuffled','runWiseLagShuffled');
% 
%         else
%             clear allVars runWiseCorrShuffled
%             allVars = load([dataDir '\tempControlVarsFOV.mat']);
%             runWiseCorrAllShuffled{iDate,iRun} = allVars.runWiseCorrShuffled;
%             runWiseLagAllShuffled{iDate,iRun}  = allVars.runWiseLagShuffled;
%             runWiseCorrShuffled                = allVars.runWiseCorrShuffled;
%             corrNegShuffle(iDate,iRun,:)       = median(runWiseCorrShuffled,2,'omitnan');
%             corrNegTimes(iDate,iRun,:)         = median(allVars.runWiseLagShuffled,2,'omitnan');
%         end
% 
%         if ~exist([dataDir '\temporalControlFOV_v2.png'],'file')
%            winLen  = {0.001 0.1 1 5 10 50 100 500 'Unshuffled'}; 
%             figure; plot(movmean(runWiseCorrShuffled,[0 1]),'Color',[0.65 0.65 0.65],'LineWidth',1); hold on;
%             plot(movmean(squeeze(median(runWiseCorrShuffled,2,'omitnan')),[0 1]),'k','LineWidth',2);
% %             semAll  = std(runWiseCorrShuffled,0,2)./sqrt(size(runWiseCorrShuffled,2));
% %             xVar = [(1:length(winLen)) fliplr((1:length(winLen)))];
% %             yVar = [(squeeze(corrNegShuffle(iDate,iRun,:))-2.*semAll)' ...
% %                 flipud((squeeze(corrNegShuffle(iDate,iRun,:))+2.*semAll))'];
% %             patch(xVar,yVar,'blue','FaceAlpha',0.3,'EdgeColor','none');
% 
%             xticklabels(winLen);  hold on;  ylim([-1 0.5]); yticks(-1:0.2:0.5);
%             xlabel('Length of window for shuffling (s)'); ylabel('Cross correlation'); grid on;
%             title(strrep(['Full FOV temporal control for ' monkeyName ' ' expDate ' ' runName],'_','\_'));
%             f = gcf; exportgraphics(f,[dataDir '\temporalControlFOV_v2.png'],'Resolution',300); close gcf;
%         end

        %%% 4.SPATIAL CONTROL (FC MAP BASED): Correlating FC maps
        %%% obtained from seeds placed at varying distances from the
        %%% probe location with the hybrid map
        disp('Obtaining spatial controls...');
        pDatTemp = reshape(pDatTemp,[imSize(1) imSize(2) imSize(3)]);
        tic;
        if ~exist([dataDir '\spatialControlVarsFOV.mat'],'file') 
            clear  lagLow  hybridLowMap runWiseSpatialCorr runWiseSpatialTimes locAll corrHybridMap
            minShift  = round(roiSize{iDate}(iRun)/spatialBin);
            theta     = linspace(0,2*pi, round(pi*minShift)); % number of angles
            maxPoints = length(theta);
            distShift = {'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'};

            for iShift = 1:5
                clear loc locShift row col numPoints pixelLoc
                if iShift == 1
                    locShift = round(roiSize{iDate}(iRun)/spatialBin);
                else
                    locShift = round(roiSize{iDate}(iRun)*2*(iShift-1)/spatialBin);
                end

                % Get coordinates for a circle
                loc(:,1) = round(locShift * cos(theta) + seedLocProbe(1));
                loc(:,2) = round(locShift * sin(theta) + seedLocProbe(2));

                % Get the FC map for the locations and correlate with hybrid map
                for iPoint = 1:size(loc,1)
                    clear seedSigT corrMapT mapVals seedRad seed 
                    seedRad =  round(roiSize{iDate}(iRun)/(spatialBin)); % seed radius for FC map- 500um

                    if any((fliplr(loc(iPoint,:))+seedRad)>size(greenFig)) || any(loc(iPoint,:)-seedRad<= 0)
                        runWiseSpatialCorr(iShift,iPoint)  = NaN;
                        runWiseSpatialTimes(iShift,iPoint) = NaN;
                        loc(iPoint,:) = NaN;
                    
                    else
                        seed = loc(iPoint,:);
                        clipMask_seed = ~corrMask(seed(2)-seedRad:seed(2)+seedRad,seed(1)-seedRad:seed(1)+seedRad);
                        allCortex_seed = ~allCortexMask(seed(2)-seedRad:seed(2)+seedRad,seed(1)-seedRad:seed(1)+seedRad);
                        maskSize = size(clipMask_seed,1) *size(clipMask_seed,2);

                        % Check if the masks occupy more than 50% of seed
                        % or if the seeds are close to the edge (the edges
                        % of the image should at max occupy 25% of seed)
                        if ((sum(clipMask_seed,'all')/maskSize)>0.5) || ((sum(allCortex_seed,'all')/maskSize)>0.25)
                            runWiseSpatialCorr(iShift,iPoint)  = NaN;
                            runWiseSpatialTimes(iShift,iPoint) = NaN;
                            loc(iPoint,:) = NaN;
                            continue;

                        else
                            % Get Gaussian weighted seed signal
                            seedSigT = calculateSeedSignal(greenFig,clipMaskCortex,loc(iPoint,:),seedRad,pDatTemp);
                            corrMapT = reshape(plotCorrMap(seedSigT,pDatTemp,0),[imSize(1)*imSize(2) 1]);

                            % Correlate FC map with all hybrid maps. This
                            % correlation makes sure that we remain agnostic
                            % about the peak negative lag and this gives a
                            % distribution of lags and/or correlations. The
                            % variables that change are lags, distance between
                            % probe and seed locations.
                            for iMap = 1:401; mapVals(iMap) = corr(corrMapT,mapsAll(iMap,:)','rows','complete'); end

                            [runWiseSpatialCorr(iShift,iPoint),runWiseSpatialTimes(iShift,iPoint)] = min(mapVals(negIdx));
                            runWiseSpatialTimes(iShift,iPoint) = negVals(runWiseSpatialTimes(iShift,iPoint))./10;

                            % Correlate FC map with peak negative hybrid map -
                            % the lag is fixed here, the distance between the
                            % probe and the seed is the variable that is
                            % changing here.
                            corrHybridMap(iShift,iPoint) = corr(corrMapT,peakNegHybridMap','rows','complete');
                        end
                    end
                    locAll{iShift} = loc;
                end
            end

            save([dataDir '\spatialControlVarsFOV.mat'],'runWiseSpatialCorr','runWiseSpatialTimes','locAll','corrHybridMap');
            spCorrControl{iRun,iDate}      = runWiseSpatialCorr;
            spCorrControlTimes{iRun,iDate} = runWiseSpatialTimes;
            spCorrHybridMap{iRun,iDate}    = corrHybridMap;

            toc;
        else
            clear allSpatialVars
            allSpatialVars = load([dataDir '\spatialControlVarsFOV.mat']);
            spCorrControl{iRun,iDate}      = allSpatialVars.runWiseSpatialCorr;
            spCorrControlTimes{iRun,iDate} = allSpatialVars.runWiseSpatialTimes;
            spCorrHybridMap{iRun,iDate}    = allSpatialVars.corrHybridMap;
            locAll                         = allSpatialVars.locAll;
        end

        if ~exist([dataDir '\spatialControlFOV_v2.png'],'file')
            cVals = {'w','k','b','g','m'};
            figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(131); imagesc(greenFig); hold on; colormap gray; axis image off;
            plot(seedLocProbe(1),seedLocProbe(2),'Marker','pentagram','MarkerSize',15,...
                'MarkerFaceColor','r','MarkerEdgeColor','none');

            for iShift = 1:5
                if ~isempty(locAll{iShift})
                    plot(locAll{iShift}(:,1),locAll{iShift}(:,2),'.','Color',cVals{iShift},'MarkerSize',10);
                end
            end

            subplot(132);boxplot(spCorrControl{iRun,iDate}',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});
            xlabel('Distance from probe (mm)'); ylabel('Correlation between FC map and peak negative map');
            ylim([-1 0.3]);

            subplot(133);boxplot(spCorrControlTimes{iRun,iDate}',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});
            xlabel('Distance from probe (mm)'); ylabel('Lag at peak negative correlation'); ylim([-11 3]);

            sgtitle(strrep(['FC maps vs Hybrid maps (varying lag) for ' monkeyName ' ' expDate ' ' runName],'_','\_'));
            f = gcf; exportgraphics(f,[dataDir '\spatialControlFOV_v3.png'],'Resolution',300); close gcf;
        end

        if ~exist([dataDir '\spatialControl_HybridMap_v2.png'],'file')
            cVals = {'w','k','b','g','m'};
            figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(121); imagesc(greenFig); hold on; colormap gray; axis image off;
            plot(seedLocProbe(1),seedLocProbe(2),'Marker','pentagram','MarkerSize',15,...
                'MarkerFaceColor','r','MarkerEdgeColor','none');

            for iShift = 1:5
                if ~isempty(locAll{iShift})
                    plot(locAll{iShift}(:,1),locAll{iShift}(:,2),'.','Color',cVals{iShift},'MarkerSize',10);
                end
            end

            subplot(122);boxplot(spCorrHybridMap{iRun,iDate}',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});
            xlabel('Distance from probe (mm)'); ylabel('Correlation between FC map and peak negative map');
            sgtitle(strrep(['FC maps vs Hybrid map (fixed lag) for ' monkeyName ' ' expDate ' ' runName],'_','\_')); ylim([-1 0.3]);

            f = gcf; exportgraphics(f,[dataDir '\spatialControl_HybridMap_v2.png'],'Resolution',300); close gcf;
        end
    end
end

%% Grouping temporal controls

clear runWiseCorrAllShuffledT 
winLen  = {0.001 0.1 1 5 10 50 100 500 'Unshuffled'};

runWiseCorrAllShuffledT = reshape(runWiseCorrAllShuffled,[size(runWiseCorrAllShuffled,1)*size(runWiseCorrAllShuffled,2) 1]);
zeroInd   = cell2mat(cellfun(@(x) isempty(x),runWiseCorrAllShuffledT,'un',0));
runWiseCorrAllShuffledT(zeroInd) = [];
runWiseCorrAllShuffledT(~goodRunsSpatial) = [];

runWiseCorrAllShuffledT = cellfun(@(x) x./abs(min(x(9,:))),runWiseCorrAllShuffledT,'un',0)';
runWiseCorrAllShuffledT = cell2mat(cellfun(@(x) median(x,2,'omitnan'),runWiseCorrAllShuffledT,'un',0));

% figure; plot(smoothdata(runWiseCorrAllShuffledT,2,'movmean',10),'Color',[0.65 0.65 0.65],'LineWidth',1); hold on;
figure; plot(smoothdata(runWiseCorrAllShuffledT,2,'movmean',10),'Color',[0.65 0.65 0.65],'LineWidth',1); hold on;
plot(median(smoothdata(runWiseCorrAllShuffledT,2,'movmean',2),2),'k','LineWidth',2);
xticklabels(winLen);  hold on;  ylim([-1 0.5]); yticks(-1:0.1:0.5);box off; ylim([-0.8 -0.1]);
xlabel('Length of window for shuffling (s)'); ylabel('Cross correlation'); grid off;
box off; ylim([-1 -0.25])

stdAll = std(runWiseCorrAllShuffledT,[],2)/sqrt(size(runWiseCorrAllShuffledT,2));
c95 = tinv([0.025 0.975],size(runWiseAll,2)-1);
y95 = bsxfun(@times,stdAll', c95(:));

meanAll = mean(runWiseAll,2,'omitnan');
posValAll = meanAll+y95';

figure; plot(median(smoothdata(runWiseAll,2,'movmean',2),2),'k','LineWidth',2);
hold on
%  xVar = [(1:length(winLen)) fliplr((1:length(winLen)))];
%  yVar = [(squeeze(corrNegShuffle(iDate,iRun,:))-2.*semAll)' ...
%  flipud((squeeze(corrNegShuffle(iDate,iRun,:))+2.*semAll))'];
%  patch(xVar,yVar,'blue','FaceAlpha',0.3,'EdgeColor','none');
xVar = [1:size(posValAll,1) fliplr((1:size(posValAll,1)))];
patch(xVar,[posValAll(:,1)' fliplr(posValAll(:,2)')],[0.65 0.65 0.65],'FaceAlpha',0.3,'EdgeColor','none')
clear corrNegShuffleT

corrNegShuffleT = reshape(corrNegShuffle,[size(corrNegShuffle,1)*size(corrNegShuffle,2) 9]);
zeroInd =  (corrNegShuffleT(:,1)) == 0; 
corrNegShuffleT(zeroInd,:) = [];
corrNegShuffleT(~goodRunsSpatial,:) = [];

figure; plot(smoothdata(corrNegShuffleT',2,'movmean',8),'Color',[0.65 0.65 0.65],'LineWidth',1); hold on;
plot(median(smoothdata(corrNegShuffleT,2,'movmeanv',3),1),'k','LineWidth',2);
xticklabels(winLen);  hold on;  ylim([-1 0.5]); yticks(-1:0.2:0.5);
xlabel('Length of window for shuffling (s)'); ylabel('Cross correlation'); grid off;
box off; ylim([-0.8 -0.1]);


%% Group population data from controls
% Get medians (for each distance group) for each run and combine across all
% runs
clear spCorrNorm medSPCorrControl medSpCorrAll spCorrMinT
spCorrMin = cellfun(@(x) x./abs(min(x,[],'all','omitnan')),spCorrControl,'un',0);  
spCorrMin = reshape(spCorrMin,[size(spCorrControl,1)*size(spCorrControl,2) 1]);
zeroInd   = cell2mat(cellfun(@(x) isempty(x),spCorrMin,'un',0));
spCorrMin(zeroInd) = [];
spCorrMin(~goodRunsSpatial) = [];
spCorrMin = cellfun(@(x) median(x,2,'omitnan'),spCorrMin,'un',0);
spCorrMin = (cat(2,spCorrMin{:}))'; 

figure; boxplot(spCorrMin,{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'}); hold on; 
 ylim([-1 0.5]); yticks(-1:0.2:1);
xlabel('Distance from probe (mm)'); ylabel('Correlation between FC map and peak negative map');

figure; violin(spCorrMin,'facecolor','b','edgecolor','none','bw',0.1);
box off; legend off; ylim([-1.2 0.5]);yticks(-1:0.2:1); 
xticklabels([0.5 1 2 3 4]);

[pSpCorr,tblSpCorr,statsSpCorr] = anova1(spCorrMin,{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'},'off');
[rSpCorr,mSpCorr,~,gnamesSpCorr] = multcompare(statsSpCorr,"CriticalValueType","bonferroni");

tblSpCorrM = array2table(rSpCorr,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblSpCorrM.("Group") = gnamesSpCorr(tblSpCorrM.("Group"));
tblSpCorrM.("Control Group") = gnamesSpCorr(tblSpCorrM.("Control Group"));

spCorrMinT = reshape(cellfun(@(x) x./abs(min(x,[],'all','omitnan')),spCorrControl,'un',0),[size(spCorrControl,1)*size(spCorrControl,2) 1]);
spCorrMinT(zeroInd) = [];
spCorrMinT(~goodRunsSpatial) = [];
spCorrMinT = -(cat(2,spCorrMinT{:}));

figure; boxplot(spCorrMinT',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'}); hold on; 
ylim([-1 0.5]);yticks(-1:0.2:0.5);box off;
xlabel('Distance from probe (mm)'); ylabel('Correlation between FC map and peak negative map');

figure; violin(spCorrMinT','facecolor','b','edgecolor','none','bw',0.1);
box off; legend off; hold on; 

x1 = (reshape(repmat(1:5,[934 1]),[934*5 1]));
y1 = reshape(spCorrMinT',[934*5 1]);
s = swarmchart(x1,y1,5,'b','filled');
s.XJitterWidth = 0.5;
ylim([-1 1.3]);yticks(-1:0.1:1.3);box off;
xticks(1:5);xticklabels({'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});

% swarmchart(repmat(1:5,[934 1]),spCorrMinT',20,'.b');
% ylim([-1 1.3 ]);yticks(-1:0.1:1.3); hold on;

% Plot individual run data 
clear testSpCorr
testSpCorr = spCorrControl{6,2};
figure; violin(-testSpCorr','facecolor','b','edgecolor','none','bw',0.1);
box off; legend off; ylim([-1 1.3]); hold on; 

x2 = (reshape(repmat(1:5,[41 1]),[41*5 1]));
y2 = reshape(-testSpCorr',[41*5 1]);
s = swarmchart(x2,y2,5,'b','filled');
s.XJitterWidth = 0.5;
ylim([-1 1.3]);yticks(-1:0.1:1.3);box off;
xticks(1:5);xticklabels({'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});




[pSpCorrT,tblSpCorrT,statsSpCorrT] = anova1(spCorrMinT',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'},'off');
[rSpCorrT,mSpCorrT,~,gnamesSpCorrT] = multcompare(statsSpCorrT,"CriticalValueType","bonferroni","Alpha", 0.008);

tblSpCorrMT = array2table(rSpCorrT,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblSpCorrMT.("Group") = gnamesSpCorrT(tblSpCorrMT.("Group"));
tblSpCorrMT.("Control Group") = gnamesSpCorrT(tblSpCorrMT.("Control Group"));


medSpCorrControl = cellfun(@(x) median(x,2,'omitnan'),spCorrControl,'un',0);
medSpCorrControl = reshape(medSpCorrControl,[size(spCorrControl,1)*size(spCorrControl,2) 1]);
zeroInd = cell2mat(cellfun(@(x) (isempty(x)), medSpCorrControl,'un',0)); 
medSpCorrControl(zeroInd) = [];
medSpCorrControl(~goodRunsSpatial) = [];
medSpCorrControl = cat(2,medSpCorrControl{:})';

figure; boxplot(medSpCorrControl,{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'}); hold on; 
ylim([-1 0.5]); yticks(-1:0.2:0.5);
xlabel('Distance from probe (mm)'); ylabel('Correlation between FC map and peak negative map');

[pSpCorr,tblSpCorr,statsSpCorr] = anova1(medSpCorrControl,{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'},'off');
[rSpCorr,mSpCorr,~,gnamesSpCorr] = multcompare(statsSpCorr,"CriticalValueType","bonferroni");

tblSpCorrM = array2table(rSpCorr,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblSpCorrM.("Group") = gnamesSpCorr(tblSpCorrM.("Group"));
tblSpCorrM.("Control Group") = gnamesSpCorr(tblSpCorrM.("Control Group"));

% Another method
medSpControlT = reshape(spCorrControl,[size(spCorrControl,1)*size(spCorrControl,2) 1]);
medSpControlT(zeroInd) = [];
medSpControlT(~goodRunsSpatial) = [];
medSpControlT = cat(2,medSpControlT{:});

figure; boxplot(medSpControlT',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'}); hold on; 
ylim([-1 0.5]);yticks(-1:0.2:0.5);
xlabel('Distance from probe (mm)'); ylabel('Correlation between FC map and peak negative map');

[pSpCorrT,tblSpCorrT,statsSpCorrT] = anova1(medSpControlT',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'},'off');
[rSpCorrT,mSpCorrT,~,gnamesSpCorrT] = multcompare(statsSpCorrT,"CriticalValueType","bonferroni");

tblSpCorrMT = array2table(rSpCorrT,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblSpCorrMT.("Group") = gnamesSpCorrT(tblSpCorrMT.("Group"));
tblSpCorrMT.("Control Group") = gnamesSpCorrT(tblSpCorrMT.("Control Group"));





