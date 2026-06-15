clc; clear;
commonDir = 'C:\Users\kem294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\ISOI_Ephys\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12\fly_track\videoIO']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\spectral_analysis\continuous\dupes']));
rmpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12\spectral_analysis\continuous\dupes']));

%% Initialize variables and get monkey data
hemisphere = 'Left'; spatialBin = 3;
iM = 2; % 1 - Charlie Sheen, 2 - Whiskey

% Get good runs, channel info, location of electorde in a date x run format
switch iM
    case 1 
        monkeyName = 'CharlieSheen';
        goodRuns  = ([1 1 NaN NaN NaN;... % Good runs
            1 1 NaN NaN NaN; ...
            1 1 1 0 1; ...
            1 1 1 1 NaN]);

        singleChFlag  = ([0 0 NaN NaN NaN;... % Channel type
            0 0 NaN NaN NaN; ...
            0 0 0 0 0; ...
            0 0 0 0 NaN]);

        goodRunsSpatial  = ([1 1 NaN NaN NaN;... % Spatial correlations
            1 1 NaN NaN NaN; ...
            1 1 1 0 1; ...
            1 1 1 1 NaN]);

        smFlag  = (['S' 'M' '#' '#' '#';... % Electrode location
            'S' 'M' '#' '#' '#'; ...
            'S' 'S' 'M' 'M' 'M'; ...
            'S' 'M' 'M' 'M' '#']); 

        isoLevel = ([0.75 0.75 NaN NaN NaN;... % Isoflurane levels (%)
            0.75 0.75 NaN NaN NaN; ...
            0.9 0.9 0.9 0.9 0.9; ...
            0.75 0.8 0.8 0.9 NaN]); 

    case 2
        monkeyName = 'Whiskey';
        goodRuns   = ([1 1 1 NaN NaN NaN NaN; ... % Good runs
            1 1 1 1 1 1 1 ; ...
            1 0 1 1 NaN NaN NaN; ...
            1 1 0 1 1 1 0; ...
            1 0 0 1 1 1 1]);

        singleChFlag = ([1 1 1 NaN NaN NaN NaN; ... % Spatial correlations
            0 0 0 0 0 0 0 ; ...
            0 0 0 0 NaN NaN NaN; ...
            0 0 0 0 0 0 0; ...
            0 0 0 0 0 0 0]);

        goodRunsSpatial = ([1 1 1 NaN NaN NaN NaN; ... % Date x run
            1 1 1 1 1 1 1 ; ...
            1 0 1 1 NaN NaN NaN; ...
            1 1 0 1 1 1 1; ...
            1 0 0 1 1 1 0]);

         smFlag = (['M' 'M' 'M' '#' '#' '#' '#'; ... % Electrode location
            'M' 'M' 'M' 'M' 'S' 'S' 'S' ; ...
            'S' 'M' 'M' 'S' '#' '#' '#'; ...
            'S' 'M' 'M' 'S' 'M' 'M' 'M'; ...
            'S' 'S' 'M' 'M' 'M' 'S' 'M']);
         
         isoLevel = ([1 1.1 1.25 NaN NaN NaN NaN; ... % Isoflurane levels (%)
            1.2 1 1 1.1 1.1 1.1 1.1 ; ...
            0.8 1.05 1.75 1.75 NaN NaN NaN; ...
            0.9 1 1 1 1 1 1; ...
            0.7 0.8 1 1.1 1.3 1.3 1.3]);
end

% Compile and reshape variables 
goodRuns = reshape(goodRuns,[size(goodRuns,1)*size(goodRuns,2) 1]); % Good run flag
goodRuns(isnan(goodRuns)) = []; goodRuns = logical(goodRuns);

singleChFlag = reshape(singleChFlag,[size(singleChFlag,1)*size(singleChFlag,2) 1]); % Channel type flag
singleChFlag(isnan(singleChFlag)) = []; singleChFlag = logical(singleChFlag); 

goodRunsSpatial = reshape(goodRunsSpatial,[size(goodRunsSpatial,1)*size(goodRunsSpatial,2) 1]); % Good spatial run flag
goodRunsSpatial(isnan(goodRunsSpatial)) = []; goodRunsSpatial = logical(goodRunsSpatial);

% Sort electrode locations into sensory and motor areas
smFlag = reshape(smFlag,[size(smFlag,1)*size(smFlag,2) 1]);
smFlag((smFlag == '#')) = [];

sensoryGoodRuns = (smFlag =='S') & goodRuns; % Sensory location flag
sensoryGoodSpatialRuns = (smFlag == 'S') & goodRunsSpatial;

motorGoodRuns = (smFlag == 'M') & goodRuns; % Motor location flag
motorGoodSpatialRuns = (smFlag == 'M') & goodRuns; 

% Reshape isoflurane levels and remove the values for the bad runs
isoLevel = reshape(isoLevel,[size(isoLevel,1)*size(isoLevel,2) 1]);
isoLevel(isnan(isoLevel)) = []; 
isoLevelGoodRuns = isoLevel; isoLevelGoodRuns(~goodRuns) = [];
isoLevelSpatial = isoLevel; isoLevelSpatial(~goodRunsSpatial) = [];

% Get monkey data and parameters
[allDates,allRuns, refDate, refDir,lensCombo, roiSize, ephysFileNameAll, serverPath,probeLabel,...
    chInCortexNotes, greenMapRef] = getMonkeyParams_Imaging_Ephys(monkeyName,commonDir, hemisphere);

% Get monkey data....
[processedDat,greenIm,probe,badCh,badTimesLFP,~,estChInCortex] = ...
    getAllData_Imaging_Ephys(monkeyName,hemisphere,commonDir,serverPath,allDates,allRuns,...
    ephysFileNameAll,greenMapRef,chInCortexNotes,probeLabel,spatialBin);

clc; disp(['All physiology and imaging data for ' monkeyName ' loaded']);

if iM == 1
    monkeyID   = '302-19'; % Charlie Sheen
    daqOnly    = ['08_07_2023';'01_11_2022';'11_20_2023'];
    daqOnlyExp = ['08072023';'01112022';'11202023'];
    rippleOnly = ['11_29_2021';'01_11_2022';'08_07_2023'];
else
    monkeyID   = '303-19'; % Whiskey
    daqOnly    = ['08_14_2023';'10_16_2023'];
    daqOnlyExp = ['08142023';'10162023'];
    rippleOnly = ['12_04_2023';'02_20_2024';'04_29_2024'];
end

fs = 1e3;
[b,a] = butter(3,[1 250]./(fs/2),'bandpass'); % Bandpass filtering parameters across 1-250 Hz
[bS,aS] = butter(3,[57 62]./(fs/2),'stop'); % Bandstop filtering between 57-62 Hz

params.Fs       = 1e3; 
params.fpass    = [2 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;

%% Get data from NI-DAQ only
varHours     = NaN(5,12,'single');
varMovAvg    = NaN(5,12,'single');
varDiffEnv   = NaN(5,12,'single');
meanMovAvg   = NaN(5,12,'single');
meanEdgeFreq = NaN(5,12,'single');

% DAQ 
for iDate = 1:size(daqOnly,1)
    clear edgeFreqStart movAvg diffRange
  
    for iR = 1:10

        clear l eegStart spec timeValsSpec powMeanS badTimeThresh badTimeIndOld badTimes

        if strcmp(daqOnly(iDate,:),'08_07_2023') && iR == 2
            break;  % Include only the first recording for this experiment from NI-DAQ
        end 
        
        if strcmp(daqOnly(iDate,:),'11_20_2023') && (iR == 2 || iR == 4)
            continue; 
        end

         if strcmp(daqOnly(iDate,:),'01_11_2022') && iR>4
            break; 
        end

        try
            l = load(['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\' monkeyID '_' monkeyName '_SqM' ...
                '\Left Hemisphere\' daqOnly(iDate,:) '\run0' num2str(iR-1) '\' monkeyName '_' ...
                daqOnlyExp(iDate,:) '_run0' num2str(iR-1) '_EEG_EKG_data.mat'],'storedData');

            eegStart = l.storedData(:,1);
            eegStart = filtfilt(b,a,eegStart);
            eegStart    = filtfilt(bS,aS,eegStart);
            
            [spec,timeValsSpec,~] = mtspecgramc(eegStart,[5 3],params);

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
            eegStart(badTimes,:) = [];

             if iM==2
                if iDate ==1 && iR == 7
                    eegStart = eegStart(1:400e3);

                elseif iDate ==1 && iR == 2
                    eegStart = eegStart([150e3:400e3 600e3:850e3]);
                
                elseif iDate==1 && iR ==3
                    eegStart = eegStart([1:300e3 450e3:700e3]);
                
                elseif iDate == 1 && iR == 4
                    eegStart = eegStart([1:400e3 600e3:end]);
                end
            end

            % Calculate edge frequencies (frequency where cumulative power>90% of total power)
            [spec,~,freqValsSpec] = mtspecgramc(eegStart,[2 1],params);
            cumPowFreq = cumsum(spec,2);
            mat = cumPowFreq>= 0.9.*cumPowFreq(:,end);
            [~,highPowerIdx] = max(mat,[],2);

           % %%
           %  figure;imagesc(timeValsSpec,freqValsSpec,10.*log10(spec)')
           %  set(gca,'YDir','normal')
           %  colorbar; hold on;
           %  plot(timeValsSpec,freqValsSpec(highPowerIdx),'w'); 
           %  var(freqValsSpec(highPowerIdx))
           %  figure; plot(eegStart)
           %  %%
            edgeFreqStart{iR} = freqValsSpec(highPowerIdx);
            movAvg{iR} = smoothdata(movmean(edgeFreqStart{iR},60));

            [up,lo] = envelope(edgeFreqStart{iR},60,'peak');
            diffRange{iR}= up-lo;
            

        catch
            edgeFreqStart{iR} = NaN;
            movAvg{iR} = NaN;
            diffRange{iR} = NaN;
            continue;
        end
    end

    numElements                       = numel(cell2mat(cellfun(@(x) var(x,'omitnan'),edgeFreqStart,'un',0)));
    varHours(iDate,1:numElements)     = cell2mat(cellfun(@(x) var(x,'omitnan'),edgeFreqStart,'un',0));
    varMovAvg(iDate,1:numElements)    = cell2mat(cellfun(@(x) var(x,'omitnan'),movAvg,'un',0));
    varDiffEnv(iDate,1:numElements)   = cell2mat(cellfun(@(x) var(x,'omitnan'),diffRange,'un',0));
    meanMovAvg(iDate,1:numElements)   = cell2mat(cellfun(@(x) median(x,'omitnan'),movAvg,'un',0)); 
    meanEdgeFreq(iDate,1:numElements) = cell2mat(cellfun(@(x) median(x,'omitnan'),edgeFreqStart,'un',0)); 

end

%% Get EEG data from Ripple
clear movAvg 
for iDate = 1:size(rippleOnly,1)
    clear edgeFreqStart diffRange movAvg
    disp(rippleOnly(iDate,:)); 

    for iR = 1:10
        try
            datLoc = [serverPath rippleOnly(iDate,:) '\Electrophysiology\run0' num2str(iR-1)];
            targetFile = dir(fullfile(datLoc, 'datafile*.ns2'));

                 try
                     datName = [datLoc '\' targetFile.name(1:end-4)];
                     [edgeFreqStart{iR},diffRange{iR},movAvg{iR},eegDat] = getEEGDataRipple(datName);
                 
                 catch
                     datLoc = [serverPath rippleOnly(iDate,:) '\Electrophysiology'];
                     targetFile = dir(fullfile(datLoc, 'datafile*.ns2'));

                     if isempty(targetFile)
                         targetFile = dir(fullfile(datLoc, 'run0*.ns2'));
                         if size(targetFile,1)>1
                             fileNames = {targetFile.name}; idx = 1;
                             for iRec = 1:size(targetFile,1)
                                 fileName = fileNames{iRec};
                                 if size(fileName,2)>=15
                                     continue;
                                 else
                                     datName = [datLoc '\' fileName(1:end-4)];
                                     [edgeFreqStart{idx},diffRange{idx},movAvg{idx},eegDat] = getEEGDataRipple(datName);
                                     
                                     if strcmp(rippleOnly(iDate,:),'01_11_2022') && iRec == 2
                                         edgeFreqStart{idx} = edgeFreqStart{idx}(1:900); 
                                         diffRange{idx}     = diffRange{idx}(1:900); 
                                         movAvg{idx}        = movAvg{idx}(1:900); 
                                     end
                                     idx = idx+1;
                                 end
                             end
                             break;
                         end
                     end
                 end            
        catch
            edgeFreqStart{iR} = NaN;
            movAvg{iR} = NaN;
            diffRange{iR} = NaN;
        end
        if iM==2
            if iDate == 1 && (iR ==3 ||iR==4)
                edgeFreqStart{iR} = NaN;
                movAvg{iR} = NaN;
                diffRange{iR} = NaN;
            end
        end
    end

    numElements  = numel(cell2mat(cellfun(@(x) var(x,'omitnan'),edgeFreqStart,'un',0)));

    if strcmp(rippleOnly(iDate,:),'08_07_2023')
        edgeFreqStart = edgeFreqStart(2:end); % Remove the extra entry for run00
        movAvg        = movAvg(2:end);
        diffRange     = diffRange(2:end);

        numElements  = numel(cell2mat(cellfun(@(x) var(x,'omitnan'),edgeFreqStart,'un',0)));
        varHours(1,2:numElements+1)     = cell2mat(cellfun(@(x) var(x,'omitnan'),edgeFreqStart,'un',0));
        varMovAvg(1,2:numElements+1)    = cell2mat(cellfun(@(x) var(x,'omitnan'),movAvg,'un',0));
        varDiffEnv(1,2:numElements+1)   = cell2mat(cellfun(@(x) var(x,'omitnan'),diffRange,'un',0));
        meanMovAvg(1,2:numElements+1)   = cell2mat(cellfun(@(x) median(x,'omitnan'),movAvg,'un',0));
        meanEdgeFreq(1,2:numElements+1) = cell2mat(cellfun(@(x) median(x,'omitnan'),edgeFreqStart,'un',0));

    elseif strcmp(rippleOnly(iDate,:),'01_11_2022')
        varHours(2,4:numElements+3)     = cell2mat(cellfun(@(x) var(x,'omitnan'),edgeFreqStart,'un',0));
        varMovAvg(2,4:numElements+3)    = cell2mat(cellfun(@(x) var(x,'omitnan'),movAvg,'un',0));
        varDiffEnv(2,4:numElements+3)   = cell2mat(cellfun(@(x) var(x,'omitnan'),diffRange,'un',0));
        meanMovAvg(2,4:numElements+3)   = cell2mat(cellfun(@(x) median(x,'omitnan'),movAvg,'un',0));
        meanEdgeFreq(2,4:numElements+3) = cell2mat(cellfun(@(x) median(x,'omitnan'),edgeFreqStart,'un',0));

    else
        varHours(iDate+size(daqOnly,1),1:numElements)     = cell2mat(cellfun(@(x) var(x,'omitnan'),edgeFreqStart,'un',0));
        varMovAvg(iDate+size(daqOnly,1),1:numElements)    = cell2mat(cellfun(@(x) var(x,'omitnan'),movAvg,'un',0));
        varDiffEnv(iDate+size(daqOnly,1),1:numElements)   = cell2mat(cellfun(@(x) var(x,'omitnan'),diffRange,'un',0));
        meanMovAvg(iDate+size(daqOnly,1),1:numElements)   = cell2mat(cellfun(@(x) median(x,'omitnan'),movAvg,'un',0));
        meanEdgeFreq(iDate+size(daqOnly,1),1:numElements) = cell2mat(cellfun(@(x) median(x,'omitnan'),edgeFreqStart,'un',0));
    end
      
end

% Group data across both animals
animalData(iM).varHours     = varHours;
animalData(iM).varMovAvg    = varMovAvg;
animalData(iM).varDiffEnv   = varDiffEnv;
animalData(iM).meanMovAvg   = meanMovAvg;
animalData(iM).meanEdgeFreq = meanEdgeFreq;


%% Plotting data animal by animal

for iM = 1:3
    if iM~=3
            edgeVar      = animalData(iM).varHours;
            movAvgVar    = animalData(iM).varMovAvg;
            envelopeVar  = animalData(iM).varDiffEnv;
            edgeMedian   = animalData(iM).meanEdgeFreq;
            movAvgMedian = animalData(iM).meanMovAvg;
    else
        edgeVar      = [animalData(1).varHours; animalData(2).varHours];
        movAvgVar    = [animalData(1).varMovAvg; animalData(2).varMovAvg];
        envelopeVar  = [animalData(1).varDiffEnv; animalData(2).varDiffEnv];
        edgeMedian   = [animalData(1).meanEdgeFreq; animalData(2).meanEdgeFreq];
        movAvgMedian = [animalData(1).meanMovAvg; animalData(2).meanMovAvg];
    end

    nanRow =  all(isnan(edgeVar), 2);
    if sum(nanRow)~=0
        edgeVar(nanRow,:)      = [];
        movAvgVar(nanRow,:)    = [];
        envelopeVar(nanRow,:)  = [];
        edgeMedian(nanRow,:)   = [];
        movAvgMedian(nanRow,:) = [];
    end
% 
    figure; 
    subplot(321); plotData(edgeVar,'Time (hours)', 'Variance of edge frequency');
    if iM~=3; ylim([0 50]); end
    subplot(322); plotData(movAvgVar,'Time (hours)', 'Variance of moving average edge frequency');
     if iM~=3; ylim([-5 30]); end
    subplot(323); plotData(envelopeVar,'Time (hours)', 'Variance of range (upper-lower envelope)');
    if iM~=3; ylim([0 70]); end
    subplot(324); plotData(edgeMedian,'Time (hours)', 'Median of edge frequency');
    if iM~=3; ylim([0 30]); end
    subplot(325); plotData(movAvgMedian,'Time (hours)', 'Median of moving average edge frequency');
    if iM~=3; ylim([0 30]); end


    if iM == 1
        sgtitle('Charlie Sheen');
    elseif iM == 2
        sgtitle('Whiskey');
    else
        sgtitle('Combined');
    end
end

%% Plotting all data for both animals
 cols = 1:2:12;

for iVar = 1:5
    switch iVar
        case 1
            var = edgeVar; 
            titleName ='Variance of edge frequencies';
        case 2
            var = movAvgVar;
            titleName = 'Variance of edge frequencies (moving average)';
        case 3
            var = envelopeVar;
            titleName = 'Variance of upper-lower envelope';
        case 4
            var = edgeMedian;
            titleName = 'Median of edge frequencies';
        case 5
            var = movAvgMedian; 
            titleName = 'Median of edge frequencies (moving average)';
    end 


    data = cell2mat(arrayfun(@(c) var(:,[c c+1]), cols, 'UniformOutput', false));
    groupLabels = repelem(1:length(cols), 9*2);

    figure;
    boxplot(data(:), groupLabels);
    xticklabels(arrayfun(@(c) sprintf('%d-%d',c,c+1), cols, 'UniformOutput', false));
    hold on;

    y1 = data(:);
    x1 = groupLabels(:);  

    validIdx = ~isnan(y1);
    x1 = x1(validIdx);
    y1 = y1(validIdx);

    s = swarmchart(x1, y1, 25, 'o', 'filled');
    s.XJitterWidth = 0.5;
    box off;
    title(titleName);
    xlabel('Hours'); ylabel('Edge frequency (Hz)');
end

%% Calculate the Mann-Kendall Tau statistical test to identify trend 
for iVar = 1: 5
    switch iVar
        case 1
            var = edgeVar;
        case 2
            var = movAvgVar;
        case 3
            var = envelopeVar;
        case 4
            var = edgeMedian;
        case 5
            var = movAvgMedian;
    end  

    for iExp = 1:size(movAvgMedian,1)
        yVal = var(iExp,:);
        yVal = yVal(~isnan(yVal));
        xVal = 1:numel(yVal);
        [taub(iExp,iVar), tau, h(iExp,iVar), sig(iExp,iVar), Z, S, sen ]   = ktaub([yVal; xVal]', 0.001, 0);
    end
end


%% Extract EEG data from Ripple 
function [edgeFreq,diffRange,movAvg,eegDat] = getEEGDataRipple(datName)
    % edgeFreq: edge frequency
    % diffRange: difference between upper and lower envelope
    % movAvg: smoothed version of edgeFreq
    % eegDat: EEG data for the recording of interest 

    fs = 1e3;
    [b,a] = butter(3,[1 250]./(fs/2),'bandpass'); % Bandpass filtering parameters across 1-250 Hz
    [bS,aS] = butter(3,[57 62]./(fs/2),'stop'); % Bandstop filtering between 57-62 Hz
    
    params.Fs       = 1e3; 
    params.fpass    = [2 120];
    params.pad      = -1;
    params.tapers   = [3 5];
    params.trialave = 0;

    [~,hFile] = ns_OpenFile(datName);
    [~, fileInfo] = ns_GetFileInfo(hFile);

    for iEntity = 1: fileInfo.EntityCount
        [~,entityInfo(iEntity,1)] = ns_GetEntityInfo(hFile,iEntity); 
    end

    
    elecList   = find([entityInfo.EntityType] == 2);
    
    [~, analogInfo2] = ns_GetAnalogInfo(hFile, elecList(1));
    fs_ns5  = analogInfo2.SampleRate;
    
    elecLabel  = {entityInfo(elecList).EntityLabel}; % Gets the label of all the channels in lfpList
    elecIdx    = cell2mat(cellfun(@(x) strcmp(x(1:end-2),'analog'),elecLabel,'un',0));
    elecList(~elecIdx) = [];
    elecLabel(~elecIdx) = [];
    elecCount = entityInfo(elecList(1)).ItemCount;
    
    if ~isempty(elecList) && strcmp(elecLabel{1},'analog 1')
    
        [~, ~, eegDat] = ns_GetAnalogData(hFile,elecList(1),1,elecCount);
        if fs_ns5 == 30e3; eegDat = downsample(eegDat,fs_ns5/fs); end
    
        eegDat = filtfilt(b,a,eegDat);
        eegDat = single(filtfilt(bS,aS,eegDat));
    end
    [spec,timeValsSpec,~] = mtspecgramc(eegDat,[5 3],params);
    
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
    eegDat(badTimes,:) = [];
    
    % Calculate edge frequencies (frequency where cumulative power>90%)
    [spec,~,freqValsSpec] = mtspecgramc(eegDat,[2 1],params);
    cumPowFreq = cumsum(spec,2);
    mat = cumPowFreq>= 0.9.*cumPowFreq(:,end);
    [~,highPowerIdx] = max(mat,[],2);
    edgeFreq= freqValsSpec(highPowerIdx);
    movAvg = smoothdata(movmean(edgeFreq,60));
    
    [up,lo] = envelope(edgeFreq,60,'peak'); % Envelope
    diffRange = up-lo;  % Difference in upper and lower envelope

end

%% Boxplot function with data points
function plotData(var,xLabel, yLabel)

    boxplot(var); hold on;
    % Show the data points
    pointSize = size(var);
    x1 = reshape(repmat(1:pointSize(2), [pointSize(1) 1]), [pointSize(1)*pointSize(2) 1]);
    y1 = reshape(var, [pointSize(1)*pointSize(2) 1]);
    validIdx = ~isnan(y1);
    x1 = x1(validIdx);
    y1 = y1(validIdx);
    
    s = swarmchart(x1, y1, 25, 'o', 'filled');
    s.XJitterWidth = 0.5;
    box off;
    xlabel(xLabel);
    ylabel(yLabel);

end

%% Check data from NI-DAQ with LFP data
iDate = 2; 
iR = 4; iRun = 1;

l = load(['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\' monkeyID '_' monkeyName '_SqM' ...
    '\Left Hemisphere\' daqOnly(iDate,:) '\run0' num2str(iR-1) '\' monkeyName '_' ...
    daqOnlyExp(iDate,:) '_run0' num2str(iR-1) '_EEG_EKG_data.mat'],'storedData');

 eegDAQ = l.storedData(:,1);

 [spec,timeValsSpec,~] = mtspecgramc(eegDAQ,[5 3],params);

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
 eegDAQ(badTimes,:) = [];

 eegRipple = probe{iRun,iDate}.eegCh; 