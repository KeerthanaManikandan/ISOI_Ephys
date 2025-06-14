% This script is used to calculate edge frequencies and power (20-60 Hz)
% from spectrogram of EEG. These values are then used to classify animal
% states into light, semi-deep, and deep states.
% November 13,2024 - KM
% Codes for ISOI_Ephys Paper Figure 8 and Supplementary Figure 7

% Set paths
clc; clear;
commonDir = 'C:\Users\kem294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\ISOI_Ephys\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12\fly_track\videoIO']));

%% Load monkey data aand initalize variables
clear allMonkeyVars
hemisphere = 'Left'; 
monkeys = {'CharlieSheen';'Whiskey'};

for iM = 1:2
    allMonkeyVars(iM) =  load(['X:\Data\' monkeys{iM} '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat']); %#ok<SAGROW> 
end
bandLabels = {'Theta';'Alpha';'Beta';'Gamma';'Spiking'}; 

% Gamma band filtering parameters
gammaBand   = [30 90]; [bG,aG] = butter(3,gammaBand./(1e3/2),'bandpass'); 

% PSD parameter initialization
params.Fs       = 1e3; 
params.fpass    = [2 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;

paramsNew = params;
paramsNew.fpass = [2 120];

%% Load the isoflurane levels (%) for each run for both animals
% Whiskey
isoLevelGoodRunsW = [1;1.2;0.8;0.9;0.7;1.1;1;1;1.25;1;1.75;1.1;1.75;1;1.1;1.1;1;1.3;1.1;1;1.3;1.1;1.3];
isoLevelSpatialW  = [1;1.2;0.8;0.9;0.7;1.1;1;1;1.25;1;1.75;1.1;1.75;1;1.1;1.1;1;1.3;1.1;1;1.3;1.1;1];

% Charlie Sheen
isoLevelGoodRunsC = [0.75;0.75;0.9;0.75;0.75;0.75;0.9;0.8;0.9;0.8;0.9;0.9];
isoLevelSpatialC  = [0.75;0.75;0.9;0.75;0.75;0.75;0.9;0.8;0.9;0.8;0.9;0.9];

% Combined
combinedIsoGoodRuns = [isoLevelGoodRunsC; isoLevelGoodRunsW];
combinedIsoSpatial  = [isoLevelSpatialC; isoLevelSpatialW];

% Obtain spatial and temporal correlations 
fovCorr = [allMonkeyVars(1).peakNegValsAllT(:,:,2) ; allMonkeyVars(2).peakNegValsAllT(:,:,2)];     
roiCorr = [allMonkeyVars(1).allChCorr ; allMonkeyVars(2).allChCorr];

% Plot isoflurane vs spatial correlations for all frequencies
figure;
for iBand = 1:5
    clear coeff xFit yFit mdl
    subplot(2,3,iBand);
    plot(combinedIsoSpatial,fovCorr(:,iBand),'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on;

    % Line fitting
    coeff = polyfit(combinedIsoSpatial,fovCorr(:,iBand),1);
    xFit  = linspace(min(combinedIsoSpatial),max(combinedIsoSpatial),1000);
    yFit  = polyval(coeff,xFit); mdl = fitlm(combinedIsoSpatial,fovCorr(:,iBand));
    plot(xFit,yFit,'-k','LineWidth',1);

    text(1.5, 0.2,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
    text(1.5,0.15,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);

    xlabel('iso %'); ylabel('Correlations'); xlim([0.6 1.8]); xticks(0.6:0.1:1.8);
    title(bandLabels{iBand}); box off; ylim([-1 0.3]);sgtitle('FOV');
end

% Plot isoflurane vs temporal correlations for all frequencies
figure;
for iBand = 1:5
    clear coeff xFit yFit mdl
    subplot(2,3,iBand);
    plot(combinedIsoGoodRuns,roiCorr(:,iBand),'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on;
    
    % Line fitting
    coeff = polyfit(combinedIsoGoodRuns,roiCorr(:,iBand),1);
    xFit  = linspace(min(combinedIsoGoodRuns),max(combinedIsoGoodRuns),1000);
    yFit  = polyval(coeff,xFit); mdl = fitlm(combinedIsoGoodRuns,roiCorr(:,iBand));
    plot(xFit,yFit,'-k','LineWidth',1);

    text(1.5, 0,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
    text(1.5,-0.05,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);

    xlabel('iso %'); ylabel('Correlations'); xlim([0.6 1.8]); xticks(0.6:0.1:1.8);
    title(bandLabels{iBand}); box off; ylim([-0.7 0.1]); sgtitle('ROI');
end

%% Calculate EEG power spectra for Whiskey - 10/17/22 - runs 5-8
% Initialize variabless
expDate    = '10_17_2022';
fileNum    = 5:8;
monkeyName = 'Whiskey'; 
hemisphere = 'Left';
saveFolder = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Processed Data'];

eegPSD = NaN(120,length(fileNum)); % Pre-allocating variables
animalStateIso = [1.2 2.25 3.25 1]; % Isoflurane levels for procedure on 10_17_2022,runs 5-8

% Extract/load EEG for the specified runs
for iFile = 1:length(fileNum)
    clear eegCh eegPSDT spec powMeansS badTimeThresh badTimeIndOld badTimeInd
    var =  load([saveFolder '\datafile000'  num2str(fileNum(iFile)) '_lfp.mat']);
    eegCh = var.eegCh;
    
    % Obtain the spectrogram
    [spec,timeValsSpec,~] = mtspecgramc(eegCh,[5 3],params); 
    
    % Remove bad time segments
    powMeanS = squeeze(sum(spec,2)); % Calculate cumulative powers
    badTimeThresh   = (median(powMeanS,1)+4.5*mad(powMeanS,1,1)); % Determine threshold bounds
    badTimeIndOld = floor(timeValsSpec(powMeanS>badTimeThresh))*1e3; % Identify threshold crossings
    badTimes = [];
    
    % Extract 1 s interval before and after each threshold crossing
    if ~isempty(badTimeIndOld)
        badTimes = [];
        badTimeInd =[(badTimeIndOld-1e3)'  (badTimeIndOld+1e3)']; 

        for iL = 1:size(badTimeInd,1)
            badTimes = [badTimes badTimeInd(iL,1): badTimeInd(iL,2)];
        end
         badTimes = unique(badTimes);
    end
    
    % Remove bad time segments
    eegCh(badTimes,:) = []; 

    % Calculate edge frequencies (frequency where cumulative power>90%)
    [spec,timeValsSpec,freqValsSpec] = mtspecgramc(eegCh,[2 1],paramsNew);
    cumPowFreq = cumsum(spec,2);
    mat = cumPowFreq>= 0.9.*cumPowFreq(:,end);
    [~,highPowerIdx] = max(mat,[],2);
    meanFreqW(iFile) = median(freqValsSpec(highPowerIdx)); % Edge frequency  

    % Get EEG power for frequencies from 20-60 Hz
    idx = freqValsSpec>=20 & freqValsSpec<=60;
    animalStateCheckW(iFile) = median(10.*log10(spec(:,idx)),'all','omitnan');
end

% Plot edge frequency and EEG powers 
figure; scatter(animalStateIso,meanFreqW,45,[1 2 3 1],'filled'); xlabel('Iso %'); ylabel('Edge frequency');
figure; scatter(animalStateIso,animalStateCheckW,'filled'); xlabel('Iso %'); ylabel('EEG powers');

%%  Calculate EEG power spectra for Bordeaux and Rambo 
fs = 1e3;

% Get EEG powers and edge frequencies for the two animals 
for iM = 1:2
    clear monkeyName allDates allRuns rippleName meanFreq isoLevels eegPSD specAll animalStateCheck powEnvelope    
    
    % Get monkey info
    switch iM
        case 1 % Bordeaux
            monkeyName = '15-18_Bordeaux';
            allDates   = ['01_06_2020'; '12_03_2019';'07_09_2019'];
            allRuns    = {['run01';'run02';'run03';'run04'];'run00';['run00';'run01';'run02']};
            rippleName  = {[];[];['0001';'0002';'0003']};
            isoLevels  = [1.1 2.25 1.25 1.75; 0.6 NaN NaN NaN; 1 2 1 NaN];

        case 2 % Rambo
            monkeyName = '16-18_Rambo';
            allDates   = ['07_15_2019'; '11_19_2019'];
            allRuns    = {['run00'; 'run01';'run02';'run03'];['run00'; 'run01';'run02';'run03']};
            rippleName  = {['0001';'0002';'0003';'0004'];['0001';'0002';'0003';'0004']};
            isoLevels  = [0.9 2 2 1; 1.5 2 1.4 1.9];
    end
    serverDataPath = '\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\Euthanized Animal Data\';
    
    % Get EEG data for all recordings
    for iDate = 1:size(allDates,1)
        clear eegChT eegCh expDate
        expDate  = allDates(iDate,:);       

        for iRun = 1:size(allRuns{iDate},1)
            runName = allRuns{iDate}(iRun,:);
            
            % Find if data has been stored via Ripple or NI DAQ
            if strcmp(expDate,'01_06_2020')
                dat = load([serverDataPath monkeyName '_SqM\Left Hemisphere\' expDate '\' runName ...
                    '\Bordeaux_01062020_' runName '_EEG_EKG_data.mat']); 
                eegChT = dat.storedData(:,1).*1e3;        

            elseif strcmp(expDate,'12_03_2019')
                dat = load([serverDataPath monkeyName '_SqM\Left Hemisphere\' expDate '\' runName ...
                    '\Bordeaux_1232019_' runName '_EEG_EKG_data.mat']);
                eegChT = dat.storedData(:,1).*1e3;
           
            else % Extract data from Ripple
                 datName = ([serverDataPath monkeyName '_SqM\Left Hemisphere\' expDate  ...
                    '\RippleFiles\' rippleName{iDate}(iRun,:)]);
                      
                 % Load Ripple File
                 [nsResult,hFile] = ns_OpenFile(datName);
                 [nsResult2, fileInfoVar] = ns_GetFileInfo(hFile);
                
                 % Get entity info
                 clear entityInfo
                 for iEntity = 1: fileInfoVar.EntityCount
                     [~,entityInfo(iEntity,1)] = ns_GetEntityInfo(hFile,iEntity);
                 end

                 % Get the label for all channels
                 clear elecList elecID
                 elecList  = find([entityInfo.EntityType] == 2);
                 elecLabel = {entityInfo(elecList).EntityLabel}; % Gets the label of all the channels in lfpList
                 
                 % Identify the label for EEG
                 elecIdx   = cell2mat(cellfun(@(x) strcmp(x(1:end-2),'analog'),elecLabel,'un',0));
                 elecList(~elecIdx)  = [];
                 elecLabel(~elecIdx) = [];
                 
                 clear listVal
                 if length(elecList)==2
                     listVal = 2;
                 else
                     listVal = 1;
                 end
       
                 % Extract EEG data
                 [~, analogInfo2] = ns_GetAnalogInfo(hFile, elecList(listVal));
                 elecCount        = entityInfo(elecList(listVal)).ItemCount;
                 fs_ns5           = analogInfo2.SampleRate;
                 [~, ~, eegChT]    = ns_GetAnalogData(hFile,elecList(listVal),1,elecCount);
                 
                 if fs_ns5 == 30e3 % downsample to 1 kHz
                     eegChT = downsample(eegChT,fs_ns5/1e3);                
                 end
                 
            end
            
            % Remove bad time segments 
            clear spec powMeansS badTimeThresh badTimeIndOld badTimeInd
            [spec,timeValsSpec,freqValsSpec] = mtspecgramc(eegChT,[5 2],params);

            powMeanS = squeeze(sum(spec,2)); % Calculate cumulative powers
            badTimeThresh   = (median(powMeanS,1)+4.5*mad(powMeanS,1,1)); % Calculate threshold bounds
            badTimeIndOld = floor(timeValsSpec(powMeanS>badTimeThresh))*1e3; % Identify threshold crossings
            badTimes = [];
            
            % Extract 1 s interval before and after each threshold crossing
            if ~isempty(badTimeIndOld)
                badTimes = [];
                badTimeInd =[(badTimeIndOld-1e3)'  (badTimeIndOld+1e3)']; % Taking one second before and one second after bad time segments

                for iL = 1:size(badTimeInd,1)
                    badTimes = [badTimes badTimeInd(iL,1): badTimeInd(iL,2)];
                end
                badTimes = unique(badTimes);
            end
            
            % Remove bad times
            eegChT(badTimes,:) = [];
            eegCh(:,iRun) = eegChT(1:500e3);
        end           

        eegChNorm = eegCh;
        
        % Calculate edge frequencies and powers
        for iRun = 1:size(allRuns{iDate},1)
            clear eegChT 
            eegChT = (eegChNorm(:,iRun));

            [eegPSDT,psdFreq]   = mtspectrumc(buffer(eegChT,1000),params); % Power spectrum
            eegPSD(:,iDate,iRun) = median(eegPSDT,2,'omitnan');

            % Calculate edge frequency (frequency where cumulative
            % power>90% of total power)
            [spec,timeValsSpec,freqValsSpec] = mtspecgramc(eegChT,[2 1],paramsNew);
            cumPowFreq = cumsum(spec,2);
            mat = cumPowFreq>= 0.9.*cumPowFreq(:,end);
            [~,highPowerIdx] = max(mat,[],2);

            meanFreq(iDate,iRun) = median(freqValsSpec(highPowerIdx));
            specAll(:,:,iDate,iRun)= spec; % Store spectrogram for visualization

            % Get power for frequencies from 20-60 Hz
            idx = freqValsSpec>=20 & freqValsSpec<=60;
            animalStateCheck(iDate,iRun) = median(10.*log10(spec(:,idx)),'all','omitnan');        
        end
    end

    switch iM % Assign variables to the corresponding animals
        case 1
            eegB  = eegPSD; % Power spectrum
            specB = specAll; % Spectrogram
            isoB  = [1.1 2.25 1.25 1.75; 0.6 NaN NaN NaN; 1 2 1 NaN];%reshape([1.1 2.25 1.25 1.75; 0.6 NaN NaN NaN; 1 2 1 NaN],[size(eegPSD,2)*size(eegPSD,3) 1]);
            animalB = animalStateCheck; % Average power between 20-60 Hz
            meanFreqB = meanFreq; % Edge frequency
        case 2
            eegR  = eegPSD; 
            specR = specAll; 
            isoR  = [0.9 2 2 1; 1.5 2 1.4 1.9];
            animalR = animalStateCheck;
            meanFreqR = meanFreq; 
    end
end 

% Bordeaux
bordeauxLight = [squeeze(eegB(:,1,[1 3])) squeeze(eegB(:,2,1)) squeeze(eegB(:,3,[1 3]))]; % Light
bordeauxSemi  = [squeeze(eegB(:,1,4)) squeeze(eegB(:,3,2))]; % Semi-deep
bordeauxDeep  = [squeeze(eegB(:,1,2))]; % Deep

medBordeauxLight = 10.*log10(median(bordeauxLight(:,4:5),2,'omitnan')); % Median PSD for light recordings
medBordeauxSemi  = 10.*log10(bordeauxSemi(:,2)); % PSD for semi-deep recordings

bordeauxLightSpec = cat(3,squeeze(specB(:,:,1,[1 3])),squeeze(specB(:,:,2,1)),squeeze(specB(:,:,3,[1 3]))); % Light spectrogram
bordeauxSemiSpec  = cat(3, squeeze(specB(:,:,1,4)),squeeze(specB(:,:,3,2))); % Semi-deep spectrogram
bordeauxDeepSpec  = squeeze(specB(:,:,1,2)); % Deep spectrogram

medBordeauxLightSpec = 10.*log10(median(bordeauxLightSpec(:,:,1:3),3,'omitnan')); % Average spectrogram
medBordeauxSemiSpec  = 10.*log10(bordeauxSemiSpec(:,:,1)); 
medBordeauxDeepSpec  = 10.*log10(bordeauxDeepSpec); 

% Rambo 
ramboLight = [squeeze(eegR(:,1,[1 4])) squeeze(eegR(:,2,[1 3]))] ; % Light
ramboSemi  = [squeeze(eegR(:,1,2)) squeeze(eegR(:,2,4))] ; % Semi-deep
ramboDeep  = [squeeze(eegR(:,1,3)) squeeze(eegR(:,2,2))] ; % Deep

medRamboLight = 10.*log10(median(ramboLight,2,'omitnan')); % Median PSD 
medRamboDeep  = 10.*log10(median(ramboDeep,2,'omitnan'));
dR = medRamboLight - medRamboDeep; 

ramboLightSpec = cat(3,squeeze(specR(:,:,1,[1 4])),squeeze(specR(:,:,2,[1 3]))); % Light spectrogram 
ramboSemiSpec  = cat(3,squeeze(specR(:,:,1,2)),squeeze(specR(:,:,2,4))) ; % Semi-deep spectrogram
ramboDeepSpec  = cat(3,squeeze(specR(:,:,1,3)),squeeze(specR(:,:,2,2))) ; % Deep spectrogram

medRamboLightSpec = 10.*log10(median(ramboLightSpec,3,'omitnan')); % Average spectrogram
medRamboSemiSpec  = 10.*log10(median(ramboSemiSpec,3,'omitnan'));
medRamboDeepSpec  = 10.*log10(median(ramboDeepSpec,3,'omitnan')); 

% Set colors for plotting
colorR = [1 2 3 1; 1 3 1 2];
colorB = [1 3 1 2; 1 NaN NaN NaN; 1 2 1 NaN];

colorB = reshape(colorB,[12 1]);
colorR = reshape(colorR,[8 1]);

isoR    = reshape(isoR,[8 1]);
isoB    = reshape(isoB,[12 1]);
nanVals = isnan(isoB);

meanFreqB = reshape(meanFreqB,[12 1]); % Compile edge frequency 
meanFreqR = reshape(meanFreqR,[8 1]); 

freqIdx = freqValsSpec>=20 & freqValsSpec<=60; % Compile spectral powers
specBPow = reshape(10.*log10(squeeze(median(specB(:,freqIdx,:,:),[1 2],'omitnan'))),[12 1]);
specRPow = reshape(10.*log10(squeeze(median(specR(:,freqIdx,:,:),[1 2],'omitnan'))),[8 1]);

isoB(nanVals)         = [];
animalB(nanVals)      = [];
colorB(nanVals)       = [];
specBPow(nanVals)     = [];
meanFreqB(nanVals)    = [];

figure; scatter(isoB,specBPow,35,colorB,'filled'); ylim([-20 70]);
figure; scatter(isoR,specRPow,35,colorR,'filled');ylim([-20 70]);

figure; scatter(isoB,meanFreqB,35,colorB,'filled'); %ylim([-20 70]);
figure; scatter(isoR,meanFreqR,35,colorR,'filled');%ylim([-20 70]);

%% Combining all three monkeys - Whiskey, Bordeaux, and Rambo

isoAllMonkeys   = [isoB ;isoR; animalStateIso']; % Isoflurane
colorAllMonkeys = [colorB; colorR; [1;2;3;1]]; % Plot colors
specPowAll      = [specBPow; specRPow; animalStateCheckW']; % Powers between 20-60 Hz
edgeFreqAll     = [meanFreqB; meanFreqR; meanFreqW']; % Edge Frequencies

% Plot the powers as a function of isoflurane
figure; scatter(isoAllMonkeys,specPowAll,60,colorAllMonkeys,'filled');
ylim([-25 45]);ylabel('Power (dB)'); xlabel('Iso (%)');xlim([0.5 3.5]);

% Plot the edge frequencies as a function of isoflurane
figure; scatter(isoAllMonkeys,edgeFreqAll,60,colorAllMonkeys,'filled');
ylim([0 20]);ylabel('Edge Frequency (Hz)'); xlabel('Iso (%)');xlim([0.5 3.5]);

% K means clustering to see if the edge frequencies can be separable
[idx,c] = kmeans(edgeFreqAll,3);
figure; plot(isoAllMonkeys(idx==1),edgeFreqAll(idx==1),'r.','MarkerSize',25); hold on; 
plot(isoAllMonkeys(idx==2),edgeFreqAll(idx==2),'m.','MarkerSize',25);
plot(isoAllMonkeys(idx==3),edgeFreqAll(idx==3),'y.','MarkerSize',25);

%% Load Ephys+imaging data and get EEG powers from 20-60 Hz for simultaneous ISOI and Ephys experiments
hemisphere = 'Left'; spatialBin = 3;
for iM = 1:2
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

    isoLevel = reshape(isoLevel,[size(isoLevel,1)*size(isoLevel,2) 1]);
    isoLevel(isnan(isoLevel)) = [];
    isoLevelGoodRuns = isoLevel; isoLevelGoodRuns(~goodRuns) = [];

    % Get monkey data and parameters
    [allDates,allRuns, refDate, refDir,lensCombo, roiSize, ephysFileNameAll, serverPath,probeLabel,...
        chInCortexNotes, greenMapRef] = getMonkeyParams_Imaging_Ephys(monkeyName,commonDir, hemisphere);

    % Get monkey data....
    [~,~,probe,badCh,badTimesLFP,~,estChInCortex] = ...
        getAllData_Imaging_Ephys(monkeyName,hemisphere,commonDir,serverPath,allDates,allRuns,...
        ephysFileNameAll,greenMapRef,chInCortexNotes,probeLabel,spatialBin);

    clc; disp(['All physiology and imaging data for ' monkeyName ' loaded']);

    % Calculating EEG powers, edge frequencies, pairwise correlations 
    % Get filter parameters...
    clear fovCorr roiCorr gIntraProbeCorr gPowerCorr gInfraSlowCorr eegMaxPow meanFreq
    % Initialize variables
    fs = 1e3;
    gIntraProbeCorr = NaN(size(probe,2), size(probe,1));
    gPowerCorr      = NaN(size(probe,2), size(probe,1));
    gInfraSlowCorr  = NaN(size(probe,2), size(probe,1));
    eegMaxPow       = NaN(size(probe,2), size(probe,1));
    meanFreq        = NaN(size(probe,2), size(probe,1));
    
    clear fovCorr roiCorr

    fovCorr = [allMonkeyVars(iM).peakNegValsAllT(:,:,2)]; % Spatial correlations
    roiCorr = [allMonkeyVars(iM).allChCorr]; % Temporal correlations

    % Get powers and edge frequencies for all recordings. 
    for iDate = 1:size(allDates,1)
        clear expDate eegCh
        expDate = allDates(iDate,:);

        for iRun = 1:size(allRuns{iDate,1})
            clear runName dataDir probeCh rawCh lfpBadTimes lfpBadCh chInCortex

            runName = allRuns{iDate,1}(iRun,:);
            dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

            clc; disp(['Calculating powers and correlations for  ' monkeyName ' '  expDate ' run: ' runName]);

            % Get run related variables
            probeCh     = probe{iRun,iDate}.probeCh; % LFP
            eegChT       = probe{iRun,iDate}.eegCh; % EEG 
            lfpBadTimes = badTimesLFP{iDate,iRun}; % Bad times
            badChDat    = badCh{iDate,iRun};

            % Get gamma band time series, power, and infraslow powers
            probeCh(:,badChDat)    = [];  
            probeCh(lfpBadTimes,:) = [];

            chInCortex = estChInCortex{1,iDate}(iRun,:);
            gammaCh    = filtfilt(bG,aG,double(probeCh(:,chInCortex(1):chInCortex(2)))); % Time series
            gPower     = envelope(gammaCh,5); % Power
            gInfraSlow = getInfraSlowPowerLFP(probeCh,bG,aG,chInCortex); % Infraslow powers

            gIntraProbeCorr(iDate,iRun) = median(corr(gammaCh),'all','omitnan'); % Time series correlations
            gPowerCorr(iDate,iRun)      = median(corr(gPower),'all','omitnan'); % Power correlations
            gInfraSlowCorr(iDate,iRun)  = median(corr(gInfraSlow),'all','omitnan'); % Infraslow power correlations

            eegChT(lfpBadTimes,:) = []; % Remove bad times from LFP
            paramsNew             = params;%paramsNew.fpass = [2 120];
            [spec,~,freqValsSpec] = mtspecgramc(eegChT,[2 1],paramsNew); % Get spectrogram
            cumPowFreq            = cumsum(spec,2); % Cumulative powers
            mat                   = (cumPowFreq>= 0.9.*cumPowFreq(:,end)); % Find edge frequencies
            [~,highPowerIdx]      = max(mat,[],2); 
            meanFreq(iDate,iRun)  = median(freqValsSpec(highPowerIdx));
            
            % Get EEG powers between 20-60 Hz
            clear spec timeValsSpec freqValsSpec
            [spec,timeValsSpec,freqValsSpec] = mtspecgramc(eegChT,[5 2],params);          
            freqIdx = freqValsSpec>=20 & freqValsSpec<=60;
            eegMaxPow(iDate,iRun) = 10.*log10(squeeze(median(spec(:,freqIdx),[1 2],'omitnan')));
        end 
    end
    
    % Compile across all recordings
    eegMaxPow = reshape(eegMaxPow,[size(eegMaxPow,1)*size(eegMaxPow,2) 1]); % Powers between 20-60 Hz
    eegMaxPow(isnan(eegMaxPow)) = [];  eegMaxPow(~goodRuns) = [];

    gIntraProbeCorr = reshape(gIntraProbeCorr,[size(gIntraProbeCorr,1)*size(gIntraProbeCorr,2) 1]); % Timeseries correlations
    gIntraProbeCorr(isnan(gIntraProbeCorr)) = [];  gIntraProbeCorr(~goodRuns) = [];

    gPowerCorr = reshape(gPowerCorr,[size(gPowerCorr,1)*size(gPowerCorr,2) 1]); % Power correlations
    gPowerCorr(isnan(gPowerCorr)) = [];  gPowerCorr(~goodRuns) = [];

    gInfraSlowCorr = reshape(gInfraSlowCorr,[size(gInfraSlowCorr,1)*size(gInfraSlowCorr,2) 1]); % Infraslow power correlations
    gInfraSlowCorr(isnan(gInfraSlowCorr)) = [];  gInfraSlowCorr(~goodRuns) = [];

    meanFreq = reshape(meanFreq,[size(meanFreq,1)*size(meanFreq,2) 1]); % Edge frequencies
    meanFreq(isnan(meanFreq)) = [];  meanFreq(~goodRuns) = [];

    % Plotting pairwise correlations vs isoflurane % for individual monkeys
    figure; subplot(131); showLinearFit(isoLevelGoodRuns,gIntraProbeCorr,1.5,0.8,0.7); % Iso vs Ephys params
    xlabel('Iso %'); ylabel('Correlations');
    axis square;title('iso vs gamma time series'); xlim([0.6 1.8]); ylim([-1 1]);

    subplot(132); showLinearFit(isoLevelGoodRuns,gPowerCorr,1.5,0.8,0.7); % Power
    xlabel('Iso %'); ylabel('Correlations');
    axis square; title('iso vs gamma power'); xlim([0.6 1.8]);ylim([-1 1]);

    subplot(133);showLinearFit(isoLevelGoodRuns,gInfraSlowCorr,1.5,0.8,0.7); % Infraslow 
    xlabel('Iso %'); ylabel('Correlations');
    axis square;title('iso vs gamma infraslow power');xlim([0.6 1.8]);ylim([-1 1]);

    sgtitle(['iso % vs ephys parameters for ' monkeyName]);

    % Edge Frequency vs Ephys params
    maxLim = round(max(meanFreq)); minLim = floor(min(meanFreq)); 
    figure; subplot(131); showLinearFit(meanFreq,gIntraProbeCorr,maxLim-10,0.9,0.8); 
    xlabel('EEG metric'); ylabel('Correlations');
    axis square; title('EEG power vs gamma time series');xlim([minLim-10 maxLim+10]);ylim([-1 1]);

    subplot(132); showLinearFit(meanFreq,gPowerCorr,maxLim-10,0.9,0.8);
    xlabel('EEG metric'); ylabel('Correlations');
    axis square; title('EEG power vs gamma power');xlim([minLim-10 maxLim+10]);ylim([-1 1]);
    
    subplot(133); showLinearFit(meanFreq,gInfraSlowCorr,maxLim-10,0.9,0.8);
    axis square;xlabel('EEG metric'); ylabel('Correlations'); 
    title('EEG power vs gamma infraslow power');xlim([minLim-10 maxLim+10]);ylim([-1 1]);
    
    sgtitle(['EEG Power (au) vs ephys parameters for ' monkeyName]);


    % EEG power vs Ephys params
    maxLim = round(max(eegMaxPow)); minLim = floor(min(eegMaxPow)); 
    figure; subplot(131); showLinearFit(eegMaxPow,gIntraProbeCorr,maxLim-10,0.9,0.8); 
    xlabel('Edge Frequency'); ylabel('Correlations');
    axis square; title('Edge Frequency vs gamma time series');xlim([minLim-10 maxLim+10]);ylim([-1 1]);

    subplot(132); showLinearFit(eegMaxPow,gPowerCorr,maxLim-10,0.9,0.8);
    xlabel('Edge Frequency'); ylabel('Correlations');
    axis square; title('Edge Frequency vs gamma power');xlim([minLim-10 maxLim+10]);ylim([-1 1]);
    
    subplot(133); showLinearFit(eegMaxPow,gInfraSlowCorr,maxLim-10,0.9,0.8);
    axis square;xlabel('EEG metric'); ylabel('Correlations'); 
    title('Edge Frequency vs gamma infraslow power');xlim([minLim-10 maxLim+10]);ylim([-1 1]);
    
    sgtitle(['EEG Power (au) vs ephys parameters for ' monkeyName]);

    figure; % Edge frequencies vs Spatial correlations
    for iBand = 1:5
        clear coeff xFit yFit mdl
        subplot(2,3,iBand);
        showLinearFit(meanFreq,fovCorr(:,iBand),maxLim-10,0.2,0.1);
        xlabel('Edge Frequency'); ylabel('Correlations'); xlim([minLim-10 maxLim+10]);
        title(bandLabels{iBand}); box off; ylim([-1 0.3]);
    end
     sgtitle(['FOV - ' monkeyName]);

    figure; % Edge Frequencies vs Temporal correlations
    for iBand = 1:5
        clear coeff xFit yFit mdl
        subplot(2,3,iBand);
        showLinearFit(meanFreq,roiCorr(:,iBand),maxLim-10,0,-0.1);
        xlabel('Edge Frequency'); ylabel('Correlations'); xlim([minLim-10 maxLim+10]);
        title(bandLabels{iBand}); box off; ylim([-0.7 0.1]); 
    end
  sgtitle(['ROI - ' monkeyName]);

    figure; % EEG Power vs Spatial correlations
    for iBand = 1:5
        clear coeff xFit yFit mdl
        subplot(2,3,iBand);
        showLinearFit(eegMaxPow,fovCorr(:,iBand),maxLim-10,0.2,0.1);
        xlabel('EEG metric'); ylabel('Correlations'); xlim([minLim-10 maxLim+10]);
        title(bandLabels{iBand}); box off; ylim([-1 0.3]);
    end
     sgtitle(['FOV - ' monkeyName]);

    figure; % EEG Power vs Temporal correlations
    for iBand = 1:5
        clear coeff xFit yFit mdl
        subplot(2,3,iBand);
        showLinearFit(eegMaxPow,roiCorr(:,iBand),maxLim-10,0,-0.1);
        xlabel('EEG metric'); ylabel('Correlations'); xlim([minLim-10 maxLim+10]);
        title(bandLabels{iBand}); box off; ylim([-0.7 0.1]); 
    end
  sgtitle(['ROI - ' monkeyName]);

    figure; % Iso vs Spatial correlations
    for iBand = 1:5
        clear coeff xFit yFit mdl
        subplot(2,3,iBand);
        showLinearFit(isoLevelGoodRuns,fovCorr(:,iBand),1.5,0.2,0.1);
        xlabel('Iso %'); ylabel('Correlations'); xlim([0.6 1.8]);
        title(bandLabels{iBand}); box off; ylim([-1 0.3]);
    end
   sgtitle(['FOV - ' monkeyName]);

    figure; % Iso vs Temporal correlations
    for iBand = 1:5
        clear coeff xFit yFit mdl
        subplot(2,3,iBand);
        showLinearFit(isoLevelGoodRuns,roiCorr(:,iBand),1.5,0,-0.1);
        xlabel('Iso %'); ylabel('Correlations'); xlim([0.6 1.8]);
        title(bandLabels{iBand}); box off; ylim([-0.7 0.1]);
    end
    sgtitle(['ROI - ' monkeyName]);
  
    if iM==1 % Temporarily assign the EEG metrics to a variable 
        eegMaxPowC       = eegMaxPow;
        gIntraProbeCorrC = gIntraProbeCorr;
        gInfraSlowCorrC = gInfraSlowCorr;
        gPowerCorrC      = gPowerCorr;
        isoC             = isoLevelGoodRuns;
        meanFreqC        = meanFreq;
    else
        eegMaxAll          = [eegMaxPowC;eegMaxPow]; % Powers between 20-60 Hz
        gIntraProbeCorrAll = [gIntraProbeCorrC; gIntraProbeCorr]; % Timeseries correlations
        gInfraSlowCorrAll = [gInfraSlowCorrC; gInfraSlowCorr]; % Infraslow correlations
        gPowerCorrAll      = [gPowerCorrC; gPowerCorr]; % Power correlations
        isoAll             = [isoC; isoLevelGoodRuns]; % Isoflurane %
        meanFreqAll        = [meanFreqC; meanFreq]; % Edge frequencies

    end
end

%% Plot data from all monkeys - Charlie, Whiskey, Bordeaux and Rambo

isoHigh   = isoAll>=2; % to ensure that all data points are included
iso       = [isoAll(~isoHigh); isoAllMonkeys]; % Isoflurane %
mFreq     = [meanFreqAll(~isoHigh); edgeFreqAll]; % Edge frequencies
eegPowAll = [eegMaxAll(~isoHigh); specPowAll]; % EEG Powers

fovCorr = [allMonkeyVars(1).peakNegValsAllT(:,:,2) ; allMonkeyVars(2).peakNegValsAllT(:,:,2)]; % Spatial correlations     
roiCorr = [allMonkeyVars(1).allChCorr ; allMonkeyVars(2).allChCorr]; % Temporal correlations

oneIdx         = (gInfraSlowCorrAll== 1); 
infraPow       = gInfraSlowCorrAll; infraPow(oneIdx) =[]; % Infraslow correlations
meanFreqAllNew = meanFreqAll; meanFreqAllNew(oneIdx) = []; % Edge Frequencies
isoAllNew      = isoAll; isoAllNew(oneIdx) = []; % Isoflurane
eegPowNew      = eegMaxAll; eegPowNew(oneIdx) = []; % EEG powers


%% Figure 8 panels
figure; subplot(141); % Isoflurane vs edge frequencies
scatter(iso,mFreq,70,[ones(size(meanFreqAll(~isoHigh))); colorAllMonkeys],'filled');hold on;
ylabel('Edge Frequency (Hz)'); xlabel('Iso (%)');
xlim([0.5 3.5]); xticks(0.5:0.5:3.5); ylim([2 22]); yticks(2:2:22);

% Fit an exponential function
modelfun = @(b,x) b(1) * exp(-b(2).*x);  
beta0 = [10 2]; 
mdl = fitnlm(iso,mFreq, modelfun, beta0);
X = 0.5:0.1:3.5;
coefficients = mdl.Coefficients{:, 'Estimate'};

yFitted = coefficients(1) * exp(-coefficients(2).*X);
plot(X, yFitted, 'r-', 'LineWidth', 2);
text(2.5, 22,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
text(2.5, 20,['p-val: ' num2str(mdl.Coefficients.pValue(2))]); axis square;

subplot(142); % Plot edge frequencies vs pairwise correlations
lowID = infraPow<0.1; 
showLinearFit(meanFreqAllNew(~lowID),infraPow(~lowID),20,0.4,0.35);
xlabel('Edge Frequency'); ylabel('Correlations');
axis square;title('Edge freq vs infraslow power correlations'); 
xlim([5 22]); ylim([0.4 1]); xticks(0:2:22); yticks(0:0.1:1);

subplot(143); % Plot edge frequencies vs temporal correlations
showLinearFit(meanFreqAll, roiCorr(:,4),20,0,-0.1)
xlabel('Edge Frequency'); ylabel('Correlations');
axis square;title('Edge freq vs ROI correlations - gamma');
xlim([5 22]); ylim([-0.5 0]); xticks(0:2:22); yticks(-1:0.1:0);

subplot(144); % Plot edge frequencies vs spatial correlations
showLinearFit(meanFreqAll,fovCorr(:,4),20,0,-0.1)
xlabel('Edge Frequency'); ylabel('Correlations');
axis square;title('Edge freq vs FOV correlations - gamma'); 
xlim([5 22]); ylim([-1 -0.2]); xticks(0:2:22); yticks(-1:0.1:0);

%% Supplementary Figure 7 panels
figure; subplot(141); % Isoflurane vs EEG Power
scatter(iso,eegPowAll,70,[ones(size(eegMaxAll(~isoHigh))); colorAllMonkeys],'filled'); hold on;
ylabel('EEG Metric'); xlabel('Iso (%)'); 
xlim([0.5 3.5]); xticks(0.5:0.5:3.5); ylim([-20 40]); yticks(-20:5:40);

% Fit exponential function
modelfun = @(b,x) b(1) * exp(-b(2).*x);  
beta0 = [10 2]; 
mdl = fitnlm(iso,eegPowAll, modelfun, beta0);
X = 0.5:0.1:3.5;
coefficients = mdl.Coefficients{:, 'Estimate'};

yFitted = coefficients(1) * exp(-coefficients(2).*X) ;
plot(X, yFitted, 'r-', 'LineWidth', 2);
text(2.5, 22,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
text(2.5, 20,['p-val: ' num2str(mdl.Coefficients.pValue(2))]); axis square;

subplot(142); % Plot EEG Powers vs pairwise correlations
lowID = infraPow<0.1; 
showLinearFit(eegPowNew(~lowID),infraPow(~lowID),35,0.4,0.35);
xlabel('EEG Metric'); ylabel('Correlations');
axis square; xlim([20 40]); xticks(-20:2:38); ylim([0.4 1]);

subplot(143); % Plot EEG Powers vs temporal correlations
showLinearFit(eegMaxAll, roiCorr(:,4),35,0,-0.1)
xlabel('EEG Metric'); ylabel('Correlations');
axis square; xlim([20 38]);xticks(-20:2:38);ylim([-0.6 0]);

subplot(144); % Plot EEG Powers vs spatial correlations
showLinearFit(eegMaxAll,fovCorr(:,4),35,0,-0.1)
xlabel('EEG Metric'); ylabel('Correlations');
axis square; xlim([20 38]);xticks(-20:2:38);ylim([-1 -0.1]);

%% Show pairwise correlations as a function of isoflurane
figure; subplot(131); showLinearFit(isoAll,gIntraProbeCorrAll,1.5,0.8,0.7); % Iso vs Time series
xlabel('Iso %'); ylabel('Correlations');
axis square;title('iso vs gamma time series'); xlim([0.6 1.8]); ylim([-1 1]);

subplot(132); showLinearFit(isoAll,gPowerCorrAll,1.5,0.8,0.7); % Iso vs Power
xlabel('Iso %'); ylabel('Correlations');
axis square; title('iso vs gamma power'); xlim([0.6 1.8]);ylim([-1 1]);

subplot(133);showLinearFit(isoAll,gInfraSlowCorrAll,1.5,0.8,0.7); % Iso vs infraslow powers
xlabel('Iso %'); ylabel('Correlations');
axis square;title('iso vs gamma infraslow power');xlim([0.6 1.8]);ylim([-1 1]);

sgtitle('iso % vs ephys parameters for all monkeys');

%%  Function to fit a line
function showLinearFit(xVal,yVal,textLocX,textLocY1,textLocY2)
    plot(xVal,yVal,'o','MarkerSize',7,'MarkerFaceColor',[0 0.4470 0.7410]); hold on; box off; 
    coeff = polyfit(xVal,yVal,1);
    xFit = linspace(min(xVal),max(xVal),1000); 
    yFit = polyval(coeff,xFit); mdl = fitlm(xVal,yVal);
    plot(xFit,yFit,'-k','LineWidth',1);
    text(textLocX, textLocY1,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
    text(textLocX, textLocY2,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);
end
