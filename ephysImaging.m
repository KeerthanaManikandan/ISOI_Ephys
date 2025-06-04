% ephysImaging.m
% This code analyses electrophysiological and imaging data recorded simultaneously
% June 23,2023 - Keerthana Manikandan
% Modified: October 1, 2023
clear;clc
commonDir = 'C:\Users\KEM294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\Imaging']));
addpath(genpath([commonDir '\Codes\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\fly_track\videoIO']));

%% Initializing all the relevant variables
% Imaging parameters
monkeyName   = 'Whiskey';
spatialBin   = 3;
hemisphere   = 'Left';
saveFlag     = 0;
checkBadRuns = 0;
fs           = 1e3; % LFP 

% Ephys parameters
gammaBand      = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); % Gamma band filtering parameters
alphaBand      = [8 12];  [bA,aA] = butter(3,alphaBand./(fs/2),'bandpass'); % Alpha band filtering parameters
betaBand       = [13 30]; [bB,aB] = butter(3,betaBand./(fs/2),'bandpass');
probeBList     = [1:13 15:21 23:32];
chOutCortex    = 1:3;
chDeep         = 30:32;

% Load all the runs, experiment dates, green reference images, reference directory for the corresponding monkey
if strcmp(monkeyName,'CharlieSheen')
    allDates     = ['11_29_2021'; '01_11_2022'; '08_07_2023';'11_20_2023'];
    allRuns      = {['run00'; 'run01']; ['run03'; 'run04']; ['run02'; 'run05';'run06'; 'run07';'run08']...
                    ['run02'; 'run05'; 'run06'; 'run07' ]};

    % Get imaging parameters
    refDate      = '08_31_2021';
    refDir       = [commonDir '\' monkeyName '_SqM\' hemisphere ' Hemisphere\' refDate '\Master Green Images\'];
    refImageName = 'Charlie Sheen Combined Green 08_31_2021';

    % Get ephys parameters
    ephysFileNameAll = {['run00'; 'run01'];['run0003';'run0004'];['datafile0002'; 'datafile0005';'datafile0006' ;...
        'datafile0007';'datafile0008']; ['datafile0002'; 'datafile0005', 'datafile0006'; 'datafile0007']};
    serverPath       = '\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\302-19_CharlieSheen_SqM\Left Hemisphere\';
    probeLabel       = {[];[];['CDE1'; 'CDE1';'CDE1' ;'CDE1';'CDE1']; ['B25A'; 'B25A', 'B25A'; 'B25A']};

elseif strcmp(monkeyName,'Whiskey')
    allDates = ['08_14_2023'; '10_16_2023'];
    allRuns  = {[ 'run02'; 'run04'; 'run05'; 'run06'];[ 'run03'; 'run04'; 'run05'; ...
        'run06'; 'run07'; 'run08'; 'run09']};

    % Get imaging parameters
    refDate      = '05_09_2022';
    refDir       = [commonDir '\' monkeyName '_SqM\' hemisphere ' Hemisphere\' refDate '\Master Green Images\'];
    refImageName = 'Combined Green';

    % Get ephys parameters
    ephysFileNameAll = {['datafile0002' ; 'datafile0004'; 'datafile0005'; 'datafile0006']; ['datafile0003' ;...
        'datafile0004'; 'datafile0005'; 'datafile0006'; 'datafile0007'; 'datafile0008' ;  'datafile0009']};
    serverPath       = '\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\303-19_Whiskey_SqM\Left Hemisphere\';
    probeLabel       = {[]; ['CDE1';'B25A'; 'B25A'; 'B25A'; 'B25A'; 'B25A' ;'B25A']};

end

%% Store/Retrieve imaging data
% Get the master green image
if exist([refDir refImageName '.bmp'],'file') == 0 % Make sure to get what the reference run is outside the function
    greenMapRef = imread([refDir refImageName '.png']);
else
    greenMapRef = imread([refDir refImageName '.bmp']);
end

greenMapRef  = greenMapRef(:,:,1);
greenMapSize = size(greenMapRef);
greenMap     = double(greenMapRef);
greenMap_RGB = ind2rgb(greenMap,gray(256));

% Store/retrieve imaging data
for iDate =2% 1:size(allDates,1)
    for iRun = 4%1:7%size(allRuns{iDate,1},1)
        
        % Load all run data 
        expDate    = allDates(iDate,:);
        runName    = allRuns{iDate,1}(iRun,:);
        saveFolder = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Processed Data'];
        datFileNum = ephysFileNameAll{iDate,1};
        fileNum    = str2double(datFileNum(iRun,end));

        % Get the directory and filenames
        dataDir     = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
        if ~exist(dataDir,'dir'); [~,~] = mkdir(dataDir); end 

        numFiles    = length(dir([dataDir '\Spatial Downsample SS3' '/*.mat'])); % loading the downsampled data only
        datName     = 'Data_RS_10Hz_SS3_';
        templateDir = ['X:\Data\' monkeyName '_SqM\Left Hemisphere\'];

        fileInfo = dir(dataDir); 
        if isempty(find(strcmp({fileInfo.name},['green' runName(end-1:end) '.png']), 1)) && isempty(find(strcmp({fileInfo.name},['green' runName(end-1:end) '.bmp']), 1))
            try
                copyfile([serverPath '\' expDate '\' runName '\green' runName(end-1:end) '.png'], dataDir);
            catch
                copyfile([serverPath '\' expDate '\' runName '\green' runName(end-1:end) '.bmp'], dataDir);
            end

        end
          

        if exist([templateDir  expDate '\' allRuns{iDate,1}(iRun,:) '\green0' allRuns{iDate,1}(iRun,5) '_Edited.png'],'file')
            greenTemp = imread([dataDir '\green0' allRuns{iDate,1}(iRun,5) '_Edited.png']);
        elseif exist([templateDir  expDate '\' allRuns{iDate,1}(iRun,:) '\green0' allRuns{iDate,1}(iRun,5) '_Edited.bmp'],'file')
            greenTemp = imread([dataDir '\green0' allRuns{iDate,1}(iRun,5) '_Edited.bmp']);
        else
            error('Greens are not edited...');
            
        end
     
        greenIm{iDate,iRun}  = greenTemp(:,:,1);

        % Store/Retrieve processed imaging data
        if strcmp(expDate,'06_22_2021')
            serverDataPath    = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\302-19_CharlieSheen_SqM\Left Hemisphere\06_22_2021 Craniotomy\' runName];
        else
            serverDataPath = [serverPath expDate '\' runName];
        end

        % Store processed data 
        
        if ~exist([dataDir '\processedFrames.mat'],'file')
            clear tempBandPass
            disp(['Processing imaging data for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)]);
            [~,~,~,~,tempBandPass] = getPreProcessedDataRestingState(serverDataPath,dataDir,runName,numFiles,spatialBin,datName,0);
            disp(['Storing imaging data for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)]);
            save([dataDir '\processedFrames.mat'],'tempBandPass'); 
            processedDat{iDate,iRun} = matfile([dataDir '\processedFrames.mat']); clear tempBandPass;
            
        else
            % Retrieve processed data 
            disp(['Loading imaging data for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)]);
            processedDat{iDate,iRun} = matfile([dataDir '\processedFrames.mat']);
        end

        cd(commonDir);
    end
end

%% Store/Retrieve LFP data - Need to make this into a function... 
for iDate = 2%1:size(allDates,1)
    for iRun = 4%1:7%1:size(allRuns{iDate,1},1)
        clear entityInfo datFileName timeStamp

        % Load all run information 
        expDate    = allDates(iDate,:);
        runName    = allRuns{iDate,1}(iRun,:);
        saveFolder = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Processed Data'];
        dataDir    = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
        datFileNum = ephysFileNameAll{iDate,1};
        fileNum    = str2double(datFileNum(iRun,end));

        % Get the name of stored file
        datFileName = ephysFileNameAll{iDate,1}(iRun,:);
        datFileName = datFileName(1:end-1);

        if (fileNum>=10)
            datFileName = datFileName(1:end-1);
        end

        % Check if LFP is already stored and if camera timestamps are also
        % stored... 
        try load([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'])
            vInfo = who('-file',[saveFolder '\' datFileName num2str(fileNum) '_lfp.mat']);
        catch
            vInfo = [];
        end
        
        % Preprocess and store data... 
        if  ~exist([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'],'file') || ~ismember('timeStamp',vInfo)
            disp(['Preprocessing LFP data for : ' runName]);

            if exist([serverPath expDate '\Electrophysiology\' runName '\' datFileName num2str(fileNum) '.nev'],'file')
                datName = [serverPath expDate '\Electrophysiology\' runName '\' datFileName num2str(fileNum)];

            elseif exist([serverPath expDate '\Electrophysiology\' datFileName num2str(fileNum) '.nev'],'file')
                datName = [serverPath expDate '\Electrophysiology\' datFileName num2str(fileNum)];
            end

            try % Open the data name
                try
                    [nsResult,hFile] = ns_OpenFile(datName);
                catch
                    disp(['Data did not load for : ' runName]);
                    continue;
                end

                if ~strcmp(nsResult,'ns_OK')
                    disp('Data file did not open! - going to the next datafile');
                    continue;
                end

                % Get file information

                [nsResult2, fileInfo] = ns_GetFileInfo(hFile);
                if ~strcmp(nsResult2,'ns_OK')
                    disp('Data file information did not load!');
                    continue;
                end

                % Get entity information
                for iEntity = 1: fileInfo.EntityCount
                    [~,entityInfo(iEntity,1)] = ns_GetEntityInfo(hFile,iEntity);
                end

                % Sort the entities whether they contain events, neural data, segment data
                lfpList   = find([entityInfo.EntityType] == 2);
                lfpLabel  = {entityInfo(lfpList).EntityLabel}; % Gets the label of all the channels in lfpList

                % Remove the raw signal channels and get the indices of the LFP only
                if strcmp(lfpLabel{1}(1:3),'lfp')
                    lfpIdx = cell2mat(cellfun(@(x) strcmp(x(1:3),'lfp'),lfpLabel,'un',0)); % 32-channel electrode

                elseif strcmp(lfpLabel{1}(end-2:end),'lfp')
                    lfpIdx = cell2mat(cellfun(@(x) strcmp(x(end-2:end),'lfp'),lfpLabel,'un',0)); % 32-channel electrode

                elseif strcmp(lfpLabel{1}(1:6),'analog')
                    lfpIdx = cell2mat(cellfun(@(x) strcmp(x,'analog 2'),lfpLabel,'un',0)); % Single electrode
                end

                lfpList(~lfpIdx)  = [];
                lfpLabel(~lfpIdx) = [];

                if strcmp(lfpLabel{1},lfpLabel{2}) % To remove 30kHz sampled data
                    lfpList(2) = [];
                    lfpLabel(2) = [];
                end

                if strcmp(lfpLabel{1}(1:3),'lfp')
                    chNum = str2num(cell2mat(cellfun(@(x) x(end-1:end),lfpLabel','un',0))); %#ok<ST2NM>
                else
                    chNum = str2num(cell2mat(cellfun(@(x) x(2:3),lfpLabel','un',0))); %#ok<ST2NM>
                end
                [~,elecID] = sort(chNum);

                % Get the LFP for all the channels in sorted order
                clear b a bS aS probeCh
                [b,a]   = butter(3,[1 250]./(fs/2),'bandpass'); % Bandpass filtering parameters across 1-250 Hz
                [bS,aS] = butter(3,[57 62]./(fs/2),'stop'); % Bandstop filtering between 57-62 Hz

                disp('Getting the LFP for the single probe and filtering them ... ')

                for iElec = 1:length(lfpLabel) % Get LFP 
                    clear elecEntityID lfpEntityID lfpCount

                    if ~isempty(elecID)
                        lfpEntityID   = lfpList(elecID(iElec));
                    else
                        lfpEntityID   = lfpList(iElec);
                    end

                    lfpCount      = entityInfo(lfpEntityID).ItemCount;

                    if ~exist('analogInfo','var')
                        [~, analogInfo] = ns_GetAnalogInfo(hFile, lfpEntityID);
                        fs  = analogInfo.SampleRate;
                    end

                    % Get LFP data
                    [~, ~, probeCh(:,iElec)] = ns_GetAnalogData(hFile,lfpEntityID,1,lfpCount);
                end

                probeCh = filtfilt(b,a,probeCh); % Bandpass filtering across 1-250 Hz
                probeCh = filtfilt(bS,aS,probeCh); % Bandstop filtering between 57-62 Hz

                % Save LFP that is recorded simultaneously with imaging
                % Obtain camera frame information
                eventList   = find([entityInfo.EntityType] == 1);

                % Grab the times where the camera frame occurred
                timeStamp = [];
                for iT = 1:entityInfo(eventList(2)).ItemCount
                    [nsFrame, timeStamp(iT), ~, ~] = ns_GetEventData(hFile, eventList(2),iT);
                end

                if ~strcmp(nsFrame,'ns_OK')
                    disp('Camera frames not recorded...');
                end
                
                % Keep data between first and last camera timestamp
                t1k = unique(floor(timeStamp.*1000));
                probeCh = single(probeCh(t1k(1):t1k(end),:));

                if size(probeCh,1)>=905000;  probeCh = probeCh(1:905000,:); end

                % Getting the EEG
                clear elecList elecID
                elecList  = find([entityInfo.EntityType] == 2);
                elecLabel = {entityInfo(elecList).EntityLabel}; % Gets the label of all the channels in lfpList
                elecIdx   = cell2mat(cellfun(@(x) strcmp(x(1:end-2),'analog'),elecLabel,'un',0));
                elecList(~elecIdx) = [];
                elecLabel(~elecIdx) = [];

                if ~isempty(elecList)
                    if strcmp(elecLabel{1},'analog 1') % analog 1 has EEG
                        clear probe2Temp
                        [~, analogInfo2] = ns_GetAnalogInfo(hFile, elecList(1));
                        elecCount      = entityInfo(elecList(1)).ItemCount;
                        fs_ns5           = analogInfo2.SampleRate;
                        [~, ~, eegCh]    = ns_GetAnalogData(hFile,elecList(1),1,elecCount);
                        if fs_ns5 == 30e3; eegCh = downsample(eegCh,fs_ns5/fs); end
                        eegCh = filtfilt(b,a,eegCh);
                        eegCh = eegCh(t1k(1):t1k(end));
                    end
                end

                % Save the LFP
                if ~exist('saveFolder','dir'); [~,~] = mkdir(saveFolder); end
                disp('Storing data... ' );

                if exist('eegCh','var')
                    save([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'] ,'probeCh','eegCh','timeStamp');
                else
                    save([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'] ,'probeCh','timeStamp');
                end

                % Filter the LFP to remove 1-5 Hz.
                if ~exist('fs','var'); fs = 1e3; end
                [bL,aL] = butter(3,([6 250]./(fs/2)),'bandpass');
                probe{iRun,iDate} = single(filtfilt(bL,aL,double(probeCh)));
                eeg{iRun,iDate}   = filtfilt(bL,aL,eegCh);
                cameraTimeStamp{iRun,iDate} = timeStamp;
                
            catch
                continue;
            end

        else
            % Retrieve LFP Data
            disp(['Retrieving electrophysiology data for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)]);
            load([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat']);
            probe{iRun,iDate} = probeCh;

            % Filter the LFP to remove 1-5 Hz.
            fs = 1e3;
            [bL,aL] = butter(3,([6 250]./(fs/2)),'bandpass');
            probe{iRun,iDate} = single(filtfilt(bL,aL,double(probe{iRun,iDate})));

            if exist('eegCh','var')
                eeg{iRun,iDate}   = filtfilt(bL,aL,eegCh);
            else
                eeg{iRun,iDate} = [];
            end
            cameraTimeStamp{iRun,iDate} = timeStamp; % Camera frame acquisition timestamps

        end
    end
end

%% Identify noisy data from LFP and remove them from both LFP and imaging data
fs              = 1000;
params.Fs       = fs;
params.fpass    = [1 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;

for iDate = 2%:size(allDates,1)
    for iRun = 4%1:7%:size(allRuns{iDate,1},1)
        clear probeCh processedDatR clipMask allCortexMask elecMask vesselMask clipMaskCortex channels spec...
            powMeanS badTimeThreshHigh processedDat10 badTimeIm

        badTimeInd = []; badTimes = [];
        runName    = allRuns{iDate,1}(iRun,:);
        dataDir    = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

        disp(['Identifying noisy data from LFP for:  ' allDates(iDate,:) ' ' runName] );
        probeCh = probe{iRun,iDate};
        channels  = 1:size(probeCh,2);

        % Determine the bad time segments and bad channels in LFP data
        % Remove bad channels...
        badChVal = [];

        % Remove bad channels that are already known from probe datasheet
        if ~isempty(probeLabel{iDate,1}(iRun,:)) 
            if strcmp(probeLabel{iDate,1}(iRun,:),'BD29')
                channels(ismember(channels,8)) = [];
                badChVal = [badChVal 8];

            elseif strcmp(probeLabel{iDate,1}(iRun,:),'CDE1')
                channels(ismember(channels,9)) = [];
                badChVal = [badChVal 9];

            elseif strcmp(probeLabel{iDate,1}(iRun,:),'D553')
                channels(ismember(channels,5)) = [];
                badChVal = [badChVal 5];

            elseif strcmp(probeLabel{iDate,1}(iRun,:),'B25A')
                 channels(ismember(channels,27)) = [];
                badChVal = [badChVal 27];
            end
        end

        [spec,timeValsSpec,freqValsSpec] = mtspecgramc(probeCh(:,channels),[5 2],params);
        powTimeBin = squeeze(sum(10.*log10(abs(spec)),2));

        if length(channels)~= 1              
            badElecThreshHigh = (median(powTimeBin,2)+5*mad(powTimeBin,1,2));
            badElecThreshLow  = (median(powTimeBin,2)-5*mad(powTimeBin,1,2));

            badChVal = [badChVal channels(sum((powTimeBin>badElecThreshHigh),1) >= floor(0.75*size(powTimeBin,1)) |(sum((powTimeBin<badElecThreshLow),1) >= floor(0.75*size(powTimeBin,1))))];
            if ~isempty(badChVal); channels(ismember(channels,badChVal)) = []; end
        end
                
        % Remove bad time segments
        if length(channels)~= 1 % Linear array
            badTimeThreshHigh = mean(probeCh(:,channels(channels>=15)),'all','omitnan') + 5*std(probeCh(:,channels(channels>=15)),[],[1,2]);
            badTimeThreshLow  = mean(probeCh(:,channels(channels>=15)),'all','omitnan') - 5*std(probeCh(:,channels(channels>=15)),[],[1,2]);
        else % Single probe 
            badTimeThreshHigh = mean(probeCh,'omitnan') + 5*std(probeCh,[]);
            badTimeThreshLow= mean(probeCh,'omitnan')- 5*std(probeCh,[]);
        end

        % Determine threshold crossings...
        badTimeIndOld = find(mean(probeCh(:,channels(channels>=15)),2,'omitnan')>badTimeThreshHigh | mean(probeCh(:,channels(channels>=15)),2,'omitnan')<badTimeThreshLow);

        % Determine the bad time segments 
        badTimes = [];
        if ~isempty(badTimeIndOld)
            badTimeInd = [(badTimeIndOld-fs/20)  (badTimeIndOld+fs/20)]; % Taking 50 ms before and after each threshold crossing

            for iL = 1:size(badTimeInd,1)
                badTimes = [badTimes badTimeInd(iL,1): badTimeInd(iL,2)];
            end

        end
        badTimes = unique(badTimes);

        % Getting spectrograms to check if the LFP is good before and after
        % removing bad time segments 
        if ~exist([dataDir '\AverageLFPSpec.png'],'file') % Plotting and saving the spectrogram 
            clear spec timeValsSpec freqValsSpec

            if length(channels)~= 1 % Linear array 
                [spec,timeValsSpec,freqValsSpec] = mtspecgramc(probeCh(:,channels(channels>=15)),[5 2],params);
                figure; subplot(211);
                imagesc(timeValsSpec,freqValsSpec, 10.*log10(squeeze(mean(spec,3)')));
                caxis([-20 20]); set(gca,'YDir','normal'); colormap jet; colorbar;
                xlabel('Time (s)'); ylabel('Power (dB)'); title('Before removing bad time segments');
                               
                if ~isempty(badTimes)                    
                    probeCh(badTimes,:) = []; % Removing bad Time segments... 

                    clear spec timeValsSpec freqValsSpec 
                    [spec,timeValsSpec,freqValsSpec] = mtspecgramc(probeCh(:,channels(channels>=15)),[5 2],params);
                    subplot(212); imagesc(timeValsSpec,freqValsSpec, 10.*log10(squeeze(mean(spec,3)')));
                    caxis([-20 20]); set(gca,'YDir','normal');colormap jet; colorbar;
                    xlabel('Time (s)'); ylabel('Power (dB)'); title('After removing bad time segments');
                else
                    subplot(212); imagesc(timeValsSpec,freqValsSpec, 10.*log10(squeeze(mean(spec(:,:,15:end),3)')));
                    caxis([-20 20]); set(gca,'YDir','normal'); title('No bad time segments detected');colormap jet; colorbar;
                end

            else % Single probe 
                [spec,timeValsSpec,freqValsSpec] = mtspecgramc(probeCh,[5 2],params);
                figure; subplot(211); imagesc(timeValsSpec,freqValsSpec, 10.*log10(spec'));
                caxis([-20 20]); set(gca,'YDir','normal'); xlabel('Time (s)'); ylabel('Power (dB)'); colormap jet; colorbar;

                if ~isempty(badTimes)                    
                    probeCh(badTimes,:) = []; % Removing bad Time segments...

                    clear spec timeValsSpec freqValsSpec
                    [spec,timeValsSpec,freqValsSpec] = mtspecgramc(probeCh,[5 2],params);
                    subplot(212);imagesc(timeValsSpec,freqValsSpec, 10.*log10(spec'));
                    caxis([-20 20]); set(gca,'YDir','normal'); colormap jet; 
                    xlabel('Time (s)'); ylabel('Power (dB)'); title('After removing bad time segments'); colorbar;
                else
                    subplot(212); imagesc(timeValsSpec,freqValsSpec, 10.*log10(squeeze(mean(spec(:,:,15:end),3)')));
                    caxis([-20 20]); set(gca,'YDir','normal'); title('No bad time segments detected'); colormap jet; colorbar;
                end
            end
            sgtitle(strrep([monkeyName ' Date: ' expDate ' - Run: ' runName(end)],'_','\_'));
            f = gcf; exportgraphics(f,[dataDir '\AverageLFPSpec.png'],'Resolution',300); close gcf;

        else % Removing bad Time segments...
            if ~isempty(badTimes)
                probeCh(badTimes,:)=[]; 
            end
        end

        probe{iRun,iDate} = probeCh;  % Replacing LFP with its noise-free version

        % Remove bad segments from imaging data...
        % Load the appropriate masks 
        if exist([dataDir '\clipMask0' runName(end) '.BMP'],'file') == 0 % Clipmask 
            clipMask = imread([dataDir '\clipMask0' runName(end) '.png']);
        else
            clipMask = imread([dataDir '\clipMask0' runName(end) '.bmp']);
        end

        if exist([dataDir '\skullMask.bmp'],'file') == 0
            allCortexMask = imread([dataDir '\skullMask.png']); % Includes vessels in this mask
        else
            allCortexMask = imread([dataDir '\skullMask.bmp']); % Includes vessels in this mask
        end

        if exist([dataDir '\maskSkull0' runName(end) '.bmp'],'file') == 0
            elecMask = imread([dataDir '\maskSkull0' runName(end) '.png']); % Includes vessels in this mask
        else
            elecMask = imread([dataDir '\maskSkull0' runName(end) '.bmp']); % Includes vessels in this mask
        end

        clipMask       = imresize(clipMask,1/3); % Resize mask
        clipMask       = clipMask(:,:,1)>0; % Converting to 0s and 1s
        allCortexMask  = imresize(allCortexMask,1/3); % Resizing cortex mask with vessels
        allCortexMask  = allCortexMask(:,:,1)>0;
        elecMask       = imresize(elecMask-100,1/3); % so as to binarize image
        elecMask       = elecMask(:,:,1)>0;
        clipMaskCortex = clipMask & ~elecMask;

        processedDatR = reshape(processedDat{iDate,iRun},[361*438 size(processedDat{iDate,iRun},3)]);

        % Upsample the imaging dataset
        parfor iP = 1:size(processedDatR,1)
            processedDat10(iP,:) = interp(processedDatR(iP,:),5);
        end

        % Identify and remove bad times from ephys and imaging dataset
        timeStamp       = cameraTimeStamp{iRun,iDate};
        timeStampSorted = timeStamp- timeStamp(1);
        badTimes10Hz    = unique(badTimeIndOld./1000);
        badTimeIm       = [];

        % Identifying frames to be removed from RS-ISOI 
        for iT = 1: length(badTimes10Hz)
            badTimeIm(iT) = find((floor(abs(timeStampSorted - badTimes10Hz(iT))*100)./100)<=0.05,1,'first');
        end
        badTimeIm = unique(badTimeIm);


        % Check if FC maps from a particular run are good
        if ~exist([dataDir '\FCMaps.png'],'file')
            figure('units','normalized','outerposition',[0 0 1 1]); 
            imagesc(greenIm{iDate,iRun}); colormap gray; axis image off; title('Pick a seed to get a FC map'); hold on;
            seedLoc = ginput(1);
            seedLoc = fliplr(round(seedLoc));
            plot(seedLoc(2),seedLoc(1),'w.','MarkerSize',15); close gcf;
           
            seedSig = calculateSeedSignal(imresize(greenIm{iDate,iRun},1/spatialBin), clipMask,fliplr(round(seedLoc./spatialBin)),12,reshape(processedDat10,[361 438 size(processedDat10,2)])); % Get Gaussian weighted seed signal
            corrMapBefore = plotCorrMap(seedSig,reshape(processedDat10,[361 438 size(processedDat10,2)]),0);
          
            % Plotting FC maps... 
            greenImRGB = ind2rgb(greenIm{iDate,iRun},gray(256));
            frameSize = size(greenImRGB);

            figure('units','normalized','outerposition',[0 0 1 1]); subplot(121); imagesc(greenImRGB); hold on;
            corrMapR = imresize(corrMapBefore,spatialBin,'OutputSize',[frameSize(1) frameSize(2)]);
            imagesc(corrMapR,'alphadata',corrMapR.*0.8);
            colormap jet; colorbar; caxis([0 1]); axis image off;
            plot(seedLoc(2),seedLoc(1),'w.','MarkerSize',15); title('FC map before removing bad frames'); 
            
            if ~isempty(badTimeIm)
                processedDat10(:,badTimeIm) = [];
                seedSig = calculateSeedSignal(imresize(greenIm{iDate,iRun},1/spatialBin), clipMask,fliplr(round(seedLoc./spatialBin)),12,reshape(processedDat10,[361 438 size(processedDat10,2)])); % Get Gaussian weighted seed signal
                corrMapAfter = plotCorrMap(seedSig,reshape(processedDat10,[361 438 size(processedDat10,2)])); 

                subplot(122); imagesc(greenImRGB); hold on;
                corrMapR = imresize(corrMapAfter,spatialBin,'OutputSize',[frameSize(1) frameSize(2)]);
                imagesc(corrMapR,'alphadata',corrMapR.*0.8);
                colormap jet; colorbar; caxis([0 1]); axis image off;
                plot(seedLoc(2),seedLoc(1),'w.','MarkerSize',15); title('FC map after removing bad frames'); 
            else
                subplot(122); imagesc(greenImRGB); hold on;
                corrMapR = imresize(corrMapBefore,spatialBin,'OutputSize',[frameSize(1) frameSize(2)]);
                imagesc(corrMapR,'alphadata',corrMapR.*0.8);
                colormap jet; colorbar; caxis([0 1]); axis image off;
                plot(seedLoc(2),seedLoc(1),'w.','MarkerSize',15); title('No bad frames detected ');
            end

            sgtitle(strrep([monkeyName ' Date: ' expDate ' - Run: ' runName(end)],'_','\_'));
            f = gcf; exportgraphics(f,[dataDir '\FCMaps.png.png'],'Resolution',300); close gcf;

        else % Removing the bad frames
            if ~isempty(badTimeIm)
                processedDat10(:,badTimeIm) = [];
            end
        end

        % Remove vessels, skull and non-cortex from the imaging data
        clipMaskCortexR = reshape(clipMaskCortex,[size(clipMaskCortex,1)*size(clipMaskCortex,2) 1]);
        cortexOnlyMat   = processedDat10;

        cortexOnlyMat(~clipMaskCortexR,:) = NaN;        
        imData10{iDate,iRun}              = reshape(cortexOnlyMat,[361 438 size(cortexOnlyMat,2)]);

        % Determine the transition channels for LFP data 
        if size(probeCh,2)~=1

            % Get mean intra probe marginals from wideband range for Probe
            marginalVal = mean(corr(probeCh,'Rows','complete'),2,'omitnan');

            % Get the slope of the marginals for Probe A
            fx = abs(movmean(gradient(marginalVal),2,'omitnan'));

            % Find the channel with maximum slope
            if strcmp(expDate,'11_01_2021') || strcmp(expDate,'01_01_2021')
                estChInCortex{iDate}(iRun,1) = 1;
            else
                estChInCortex{iDate}(iRun,1) = find(fx == max(fx));
            end

            if estChInCortex{iDate}(iRun,1)==20
                estChInCortex{iDate}(iRun,1) = (estChInCortex{iDate}(iRun,1) - 19);
            elseif estChInCortex{iDate}(iRun,1)>20
                estChInCortex{iDate}(iRun,1) = (estChInCortex{iDate}(iRun,1) - 20);
            end

            estChInCortex{iDate}(iRun,2) = size(probeCh,2);
        else
            continue;
        end
        %{
        % Identify noisy frames
        % creating a mask of ones of size 100x100 pixels
        clear mask
        mask = zeros(size(clipMask));
        mask(150:250, 200:300) = 1;
        maskR = reshape(mask,[361*438 1]);
        processedDatTest = processedDat10;
%         processedDatTest(~clipMaskCortexR,:) = NaN;
        processedDatTest(~maskR,:) = [];

        medDatTime = squeeze(median(processedDatTest,1,'omitnan'));
        highT = mean(medDatTime) + 3*std(medDatTime);
        lowT  = mean(medDatTime) - 3*std(medDatTime);
        crossing =  find(medDatTime>=highT | medDatTime<=lowT);

        
        

        % Find bad frames
        % Use only pixels in cortex
        medDatTime = squeeze(median(cortexOnlyMat,1,'omitnan'));
        highT = mean(medDatTime) + 3*std(medDatTime);
        lowT  = mean(medDatTime) - 3*std(medDatTime);
        crossing =  find(medDatTime>=highT | medDatTime<=lowT);

        % Remove bad frames identified from both ephys and imaging dataset
        % To try doing the opposite - imaging first followed by ephys
        processedDat10(:,crossing) =[];
        cortexOnlyMat(:,crossing) = [];

        crossingEphys = crossing.*100;
        diffDuration = [1 diff(crossingEphys)];
        diffIdx = find(diffDuration~=100);
        ephysBadTimes = [];

        % Taking 50 ms before and after bad frame
        for iL = 1:length(crossingEphys)
            ephysBadTimes = [ephysBadTimes crossingEphys(iL)-50: crossingEphys(iL)+50];
        end

        ephysBadTimes = unique(ephysBadTimes);
        probeCh(ephysBadTimes,:) =[];

%         % Finding the median frame
%         medianFrame = median(processedDat10,2,'omitnan');
%         madFrame = mad(processedDat10,1,2);
% 
%         threshFrame = medianFrame+3.*madFrame;
%         threshFrameLow = medianFrame - 3.*madFrame;
%         numPixels = numel(threshFrame);
%         % threshFrameR = reshape(threshFrame,[numPixels 1]);
% 
%         badFrameIdx = zeros(size(processedDat10,2),1);
%         for iL = 1: size(processedDat10,2)
%             if (numel(find(processedDat10(:,iL)>threshFrame))> numPixels/10)|| (numel(find(processedDat10(:,iL)<threshFrameLow))> numPixels/10)
%                 badFrameIdx(iL) = 1;
%             end
%         end
%         badFrameVals = find(badFrameIdx);
        %}
    end
end
%% Check for the quality of imaging data by obtaining FC maps
iDate = 2; iRun = 4;
expDate    = allDates(iDate,:);
runName    = allRuns{iDate,1}(iRun,:);
spatialBin = 3;
figure; imagesc(greenIm{iDate,iRun}); colormap gray; axis image off; hold on;
disp('Pick a seed pixel for getting a FC map ');
seedLoc = ginput(1);
seedLoc = fliplr(round(seedLoc));
plot(seedLoc(2),seedLoc(1),'w.','MarkerSize',15); close gcf;

if exist([templateDir '\' allDates(iDate,:) '\' allRuns{iDate,1}(iRun,:) '\clipMask0' allRuns{iDate,1}(iRun,end) '.BMP'],'file') == 0
    clipMask = imread([templateDir  allDates(iDate,:) '\' allRuns{iDate,1}(iRun,:) '\clipMask0' allRuns{iDate,1}(iRun,end) '.png']);
else
    clipMask = imread([templateDir allDates(iDate,:) '\' allRuns{iDate,1}(iRun,:) '\clipMask0' allRuns{iDate,1}(iRun,end) '.bmp']);
end
clipMask   = imresize(clipMask(:,:,1),1/spatialBin);

seedSig = calculateSeedSignal(imresize(greenIm{iDate,iRun},1/spatialBin), clipMask,fliplr(round(seedLoc./spatialBin)),12,processedDat{iDate,iRun}); % Get Gaussian weighted seed signal
corrMap = plotCorrMap(seedSig,processedDat{iDate,iRun},0);

% Overlay correlation map on top of greens

greenImRGB = ind2rgb(greenIm{iDate,iRun},gray(256));
frameSize = size(greenImRGB);
figure; imagesc(greenImRGB);
hold on;
corrMapR = imresize(corrMap,spatialBin,'OutputSize',[frameSize(1) frameSize(2)]);
imagesc(corrMapR,'alphadata',corrMapR.*0.8);
colormap jet; colorbar; caxis([0 1]);
axis image off;
plot(seedLoc(2),seedLoc(1),'w.','MarkerSize',15);

%% Cross correlate infra slow changes in bandlimited LFP for the entire FOV
[z,p,k] = cheby1(2,1,[0.01 0.1]./(fs/2),'bandpass');
[sos,g] = zp2sos(z,p,k);

tic;
for iDate = 2%1:size(allDates,1)
    clear expDate datFileNum
    expDate    = allDates(iDate,:);
    datFileNum = ephysFileNameAll{iDate,1};

    for iRun = 4%1:7%:3%:size(allRuns{iDate,1},1)
        clear probeCh probeB infraA infraB envelopeBandLimited runName clipMask cc lags

        clc; disp(['Obtaining slow changes in band limited LFP for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(iRun)]);

        runName = allRuns{iDate,1}(iRun,:);
        probeCh = probe{iRun,iDate};
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

        if isempty(probeCh); continue; end

        % Get channels inside cortex...
        if size(probeCh,2)~=1
            ch = estChInCortex{iDate}(iRun,:);
            if ch(1)== 0
                ch= 1:32;
            end
            probeCh = probeCh(:,ch(1):ch(end));
        end

        % Obtain infraslow component of gamma band power
        probeBandLimited = single(filtfilt(bG,aG,double(probeCh)));
       
        % Rectify the signal
        probeBandLimited = abs(probeBandLimited);

        % Get the envelope of the signal
        envelopeBandLimited = envelope(probeBandLimited);

        % Filter envelope from 0.01-0.1 Hz and downsample it...

        envelopeBandLimited = [envelopeBandLimited(1:1e3,:) ;envelopeBandLimited ;envelopeBandLimited(1:1e3,:) ];
        envelopeFiltered = filtfilt(sos,g,double(envelopeBandLimited));
        envelopeFiltered = envelopeFiltered(1e3+1:(end-1e3),:);
        infraEphys = single(downsample(envelopeFiltered,100));
        clear envelopeFiltered envelopeBandLimited

        % Get imaging data
        processedDat10 =  reshape(imData10{iDate,iRun},[361*438 size(imData10{iDate,iRun},3)]);

        % Make sure the size of the array is the same...
        sz = min([size(processedDat10,2) size(infraEphys,1)]);
        infraEphys = infraEphys(1:sz,:);
        infraEphys = mean(infraEphys(1:sz,:),2,'omitnan');
        processedDat10 = processedDat10(:,1:sz);

        % Determine correlation between imaging and electrophysiology
        % for different lags...
        tic;
        clear cc lags

        lP = size(infraEphys,2);
        lIm = size(processedDat10,1);
        lT  = size(processedDat10,2);

        %         if size(infraEphys,2)== 1
        
        parfor iP = 1:size(processedDat10,1)
            [cc(:,iP),lags(:,iP)] = xcorr(infraEphys',processedDat10(iP,:),300,'normalized');
        end

        
        %         else
        %         parfor iP = 1:lIm
        %             for iCh = 1:lP
        %                 [cc(:,iCh,iP),lags(:,iCh,iP)] = xcorr(infraEphys(:,iCh)',processedDat10(iP,:),150,'normalized');
        %             end

        %         end
        %
        %         end
        lags = single(lags(:,1));
        toc;

        % Store the variables
        crossCorr{iDate,iRun} = cc;
        lagVals{iDate,iRun}   = lags;
        crossCorr{iDate,iRun} = reshape(crossCorr{iDate,iRun},[size(cc,1) 361 438]);

        %         % Save the variables
%         save([dataDir '\crossCorrMean15_32Avg.mat'],'cc','lags','-v7.3');
% 
%         % Save the map of averaged correlations for different lags
%         v = VideoWriter([dataDir '\corrVideoMean15_32Avg.mp4']);
%         open(v);
%         figure; 
%         for iLag = 1:size(lagVals{iDate,iRun},1)
%             imagesc(squeeze(mean(crossCorr{iDate,iRun}(iLag,:,:,:),2,'omitnan'))); axis image off;
%             colormap jet; colorbar; caxis([-0.5 0.5]);
%             title(['Cross correlation at  lag ' num2str(lagVals{iDate,iRun}(iLag,21,150,100)./10) 's']);
%             pause(0.05);
%             frame = getframe(gcf);
%             writeVideo(v,frame);
%         end
%         close(v);
%         close gcf;
    end
end
toc;

%%
v = VideoWriter([dataDir '\corrVideoMean_15mins.mp4']);
open(v);
figure('units','normalized','outerposition',[0 0 1 1]);

for iLag = 1:size(lags,1)
    imagesc(squeeze((crossCorr{iDate,iRun}(iLag,:,:)))); axis image off;
    colormap jet; colorbar; caxis([-0.5 0.5]);
    title(['Cross correlation at  lag ' num2str(lagVals{iDate,iRun}(iLag,150,100)./10) 's']);

%     subplot(122); imagesc(squeeze(cc_R(iLag,:,:))); axis image off;
%     colormap jet; colorbar; caxis([-0.5 0.5]);
%     title(['Average of infraslow signals: Cross correlation at  lag ' num2str(lagVals{iDate,iRun}(iLag,21,150,100)./10) 's']);
%       pause(0.05);
      frame = getframe(gcf);
            writeVideo(v,frame);
end
close(v);
close gcf;

%% Get the lag where highest correlation occurs 
% Get an ROI around the probe

figure; imagesc(greenIm{iDate,iRun}); colormap gray; axis image off; 
roiLag = drawrectangle; 
roiPos = roiLag.Position; close gcf;
figure; imagesc(greenIm{iDate,iRun}); colormap gray; hold on; axis image off; rectangle('Position',roiPos,'EdgeColor','w','LineWidth',3);
roiPos = round(roiPos./3); 

% Crop the ROI and average the time series in the pixels
pDat = squeeze(mean(imData10{iDate,iRun}(roiPos(2):roiPos(2)+ roiPos(4),roiPos(1):roiPos(1)+ roiPos(3),:),[1 2],'omitnan'));

% Cross - Correlate ISOI and LFP 
[ccT,lagsT] = xcorr(infraEphys',pDat','normalized');

figure; plot(lagsT,ccT);
xticks(lagsT(1): 1000:lagsT(end))
xline(0,'LineWidth',2);
xticklabels(lagsT(1)/10: 100:lagsT(end)/10);

% Truncating the correlation maps to within -30 to 30 sec
figure; plot(lagsT,ccT); xlim([-302 302]); xline(0,'LineWidth',2);
xticks(-300:30:300); xticklabels(-30:3:30);



% Correlate truncated corr map with FC map
cropMap = corrMap(roiPos(2):roiPos(2)+ roiPos(4),roiPos(1):roiPos(1)+ roiPos(3));
% greenCrop = 

% corrMapT =  reshape(corrMap(roiPos(2):roiPos(2)+ roiPos(4),roiPos(1):roiPos(1)+ roiPos(3)),[23*25 1]);
% crop10T = reshape(crop10,[601 25*23]);

for iL = 1:size(cropROI,1)
    fcCorr(iL) = corr(cropROI(iL,:)',corrMapT,'rows','complete' ); 
end 

%%

v = VideoWriter([dataDir '\corrVideoMean_Probe.mp4']);
open(v);
figure('units','normalized','outerposition',[0 0 1 1]);

for iLag = 1:size(crop10,1)
    imagesc(squeeze((crop10(iLag,:,:)))); axis image off;
    colormap jet; colorbar; caxis([-0.5 0.5]);
    title(['Cross correlation at  lag ' num2str(lagVals{iDate,iRun}(iLag+8725-1,150,100)./10) 's']);

%     subplot(122); imagesc(squeeze(cc_R(iLag,:,:))); axis image off;
%     colormap jet; colorbar; caxis([-0.5 0.5]);
%     title(['Average of infraslow signals: Cross correlation at  lag ' num2str(lagVals{iDate,iRun}(iLag,21,150,100)./10) 's']);
%       pause(0.05);
      frame = getframe(gcf);
            writeVideo(v,frame);
end
close(v);
close gcf;


%% Compare FC map - ISOI with LFP-ISOI map
meanCorrR = reshape(meanCorr,[size(meanCorr,1) 361*438]);
ccR = reshape(cc_R,[size(cc_R,1) 361*438]);
c = diag(corr(meanCorrR',ccR','rows','complete'));

c1 =corr(ccR','rows','complete'); % Ephys - imaging map with ISOI FC map for that run 


% Crop the master green map 
figure; imagesc(greenMapRef); 
roi = drawrectangle; 

fcMap = imcrop(rsConnMatrix,roi.Position); 
greenRefCropped = imcrop(greenMapRef, roi.Position); 

% Register the greens 
greenRun = greenIm{2,5}; 

[movPointsTemp,fixedPointsTemp] = cpselect(greenRun,greenRefCropped,'Wait',true);
[OTemp,spacingTemp,~] = point_registration(size(greenRefCropped),fliplr(fixedPointsTemp),fliplr(movPointsTemp));
temp = zeros(size(greenRefCropped));
temp(1:size(greenRun,1),1:size(greenRun,2)) = squeeze(clipMaskCortex);
cortexMatT = bspline_transform(OTemp,temp,spacingTemp,3);

for iL = 1: size(ccR,1)
    clear temp; 
    temp = zeros(size(greenRefCropped));
    temp(1:size(greenRun,1),1:size(greenRun,2)) = squeeze(cc_R(iL,:,:));
    [imRTemp(iL,:,:),tTemp] = bspline_transform(OTemp,temp,spacingTemp,3); % bicubic spline interpolation
end 


% Find the lag where cross correlations near the probe is high

        % Compare it with FC Map

%% Infra slow changes in power of bandlimited LFP for a specific ROI

[z,p,k] = cheby1(2,1,[0.01 0.1]./(fs/2),'bandpass');
[sos,g] = zp2sos(z,p,k);

% n = fir1(30,[0.01 0.1]./(fs/2));
rowIdx = 1;
 tic;
for iDate =2%1:size(allDates,1)
    clear expDate datFileNum
    expDate    = allDates(iDate,:);
    datFileNum = ephysFileNameAll{iDate,1};

    for iRun =1:3%:size(allRuns{iDate,1},1)
        clear probeCh probeB infraA infraB envelopeBandLimited envelopeB runName
        clc; disp(['Obtaining slow changes in band limited LFP for ' monkeyName ': Date: ' allDates(iDate,:) ' ; File: ' num2str(iRun)]);

        runName = allRuns{iDate,1}(iRun,:);
        probeCh = probe{iRun,iDate};
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

        if isempty(probeCh); rowIdx = rowIdx+1; continue; end

        % Get channels inside cortex...
        if size(probeCh,2)~=1
            ch = estChInCortex{iDate}(iRun,:);

            if ch(1)== 0
                rowIdx = rowIdx+1;
                ch= 1:32;
            end

            probeCh = probeCh(:,ch(1):ch(end));
        end

        % Obtain infraslow component of band power
        for iBand = 3%1:3
            clear probeBandLimited  envelopeBandLimited
            switch iBand
                case 1 % Alpha band
                    probeBandLimited = filtfilt(bA,aA,probeCh);
                case 2 % Beta
                    probeBandLimited = filtfilt(bB,aB,probeCh);
                case 3 % Gamma band
                    probeBandLimited = single(filtfilt(bG,aG,double(probeCh)));
            end

            % Rectify the signal
            probeBandLimited = abs(probeBandLimited);

            % Get the envelope of the signal
            envelopeBandLimited = envelope(probeBandLimited);

            % Filter envelope from 0.01-0.1 Hz and downsample it...
            envelopeBandLimited = [envelopeBandLimited(1:1e3,:); envelopeBandLimited ;envelopeBandLimited(1:1e3,:) ];
            envelopeFiltered = filtfilt(sos,g,double(envelopeBandLimited)); 
            envelopeFiltered = envelopeFiltered(1e3+1:(end-1e3),:);
            infraEphys = single(downsample(envelopeFiltered,100)); clear envelopeFiltered 

            processedDat10 = imData10{iDate,iRun}; %reshape(processedDat{iDate,iRun},[361*438 size(processedDat{iDate,iRun},3)]);

            % Load the appropriate masks
            if exist([dataDir '\clipMask0' runName(end) '.BMP'],'file') == 0
                clipMask = imread([dataDir '\clipMask0' runName(end) '.png']);
            else
                clipMask = imread([dataDir '\clipMask0' runName(end) '.bmp']);
            end

            if exist([dataDir '\skullMask.bmp'],'file') == 0
                allCortexMask = imread([dataDir '\skullMask.png']); % Has vessels in this mask
            else
                allCortexMask = imread([dataDir '\skullMask.bmp']); % Has vessels in this mask
            end

            if exist([dataDir '\maskSkull0' runName(end) '.bmp'],'file') == 0
                elecMask = imread([dataDir '\maskSkull0' runName(end) '.png']); % Has vessels in this mask
            else
                elecMask = imread([dataDir '\maskSkull0' runName(end) '.bmp']); % Has vessels in this mask
            end

            clipMask       = imresize(clipMask,1/3); % Resize mask
            clipMask       = clipMask(:,:,1)>0; % Converting to 0s and 1s
            allCortexMask  = imresize(allCortexMask,1/3); % Resizing cortex mask with vessels
            allCortexMask  = allCortexMask(:,:,1)>0;
            elecMask       = imresize(elecMask,1/3);
            elecMask       = elecMask(:,:,1)>0;
            vesselMask     = ~clipMask & allCortexMask; % Vessel mask
            clipMaskCortex = clipMask & ~elecMask;

            clipMaskCortexR = reshape(clipMaskCortex,[size(clipMaskCortex,1)*size(clipMaskCortex,2) 1]);
            processedDat10(~clipMaskCortexR,:) = NaN;

            figure; imagesc(greenIm{iDate,iRun}); colormap gray; axis image off; hold on; seedRad = 6;
            title('Pick a seed pixel that is inside a patch ');
            seedLocIn = ginput(1);
            seedLocIn = fliplr(round(seedLocIn./3));
            processedDat10Im = reshape(processedDat10,[361 438 size(processedDat10,2)]);
            inDat = processedDat10Im(seedLocIn(2)-seedRad:seedLocIn(2)+seedRad,seedLocIn(1)-seedRad:seedLocIn(1)+seedRad,:);
            inDatR = median(reshape(inDat,[size(inDat,1)*size(inDat,2) size(inDat,3)]),1,'omitnan');

            title('Pick a seed pixel that is outside a patch');
            seedLocOut = ginput(1); close gcf;
            seedLocOut = fliplr(round(seedLocOut./3));
            outDat = processedDat10Im(seedLocOut(2)-seedRad:seedLocOut(2)+seedRad,seedLocOut(1)-seedRad:seedLocOut(1)+seedRad,:);
            outDatR = median(reshape(outDat,[size(outDat,1)*size(outDat,2) size(outDat,3)]),1,'omitnan');

            % Plotting the location of the ROIs
            greenImRGB = ind2rgb(greenIm{iDate,iRun},gray(256));
            figure; imagesc(imresize(greenImRGB,[361 438]));
            hold on; axis image;

            in = [seedLocIn(2)-seedRad:seedLocIn(2)+seedRad; seedLocIn(1)-seedRad:seedLocIn(1)+seedRad];
            out = [seedLocOut(2)-seedRad:seedLocOut(2)+seedRad; seedLocOut(1)-seedRad:seedLocOut(1)+seedRad];

            plot([in(1,1) in(1,1) in(1,1)+(seedRad*2) in(1,1)+(seedRad*2) in(1,1)],[in(2,1) in(2,1)+(seedRad*2) in(2,1)+(seedRad*2) in(2,1) in(2,1)],'LineWidth',4);
            plot([out(1,1) out(1,1) out(1,1)+(seedRad*2) out(1,1)+(seedRad*2) out(1,1)],[out(2,1) out(2,1)+(seedRad*2) out(2,1)+(seedRad*2) out(2,1) out(2,1)],'LineWidth',4);

            % Determine correlation between imaging and electrophysiology
            % for different lags...
            clear cc lags ccOut lagsOut
            sz = length(inDatR);

            for iCh = 1:size(probeCh,2)
                [cc(:,iCh),lags(:,iCh)]       = xcorr(infraEphys(1:sz,iCh)',inDatR(1:sz),10,'normalized');
                [ccOut(:,iCh),lagsOut(:,iCh)] = xcorr(infraEphys(1:sz,iCh)',outDatR(1:sz),10,'normalized');
            end

            % Find the maximum lag
            [maxC,ind] = max(cc,[],1,'omitnan');
            maxLag = lags(ind);

            [maxCOut,ind] = max(ccOut,[],1,'omitnan');
            maxLagOut = lagsOut(ind);

            figure; stem(maxLag, maxC, 'filled'); hold on;
            stem(maxLagOut, maxCOut, 'filled');
            xlabel('Lag (s)'); ylabel('Correlation');
        end
    end
end
% Video code for specific ROI 

%         for iPlot = 1%:2
%             clear plotVar titleLabel
%             switch iPlot
%                 case 1
%                     plotVar = crossCorr{iDate,iRun};
%                     titleLabel = 'Cross Correlation';
%                 case 2
%                     plotVar = lagVals{iDate,iRun};
%                     titleLabel = 'Lags';
%             end
% 
%             imagesc(squeeze(plotVar(1,15,:,:))); axis image off; 
%             roi = drawrectangle; 
%             posVals = floor(roi.Position); close gcf; 
%             v = VideoWriter([dataDir '\corrVideoElec50.mp4']);
%             open(v); 
% 
%             figure; iFig = 1;
%             for iLag = 1:size(lagVals{iDate,iRun},1)
%                 imagesc(squeeze(plotVar(iLag,15,posVals(3):posVals(3)+posVals(4), posVals(1):posVals(1)+ posVals(2)))); axis image off;
%                 colormap jet;colorbar; %caxis([-101 101]);pause(0.1); colorbar
%                 caxis([-0.5 0.5]); title([titleLabel ' at  lag ' num2str(lagVals{iDate,iRun}(iLag,100,100)./10) 's']); %lagVals{iDate,iRun}(iLag,100,100)
%                 pause(0.5);
%                 iFig = iFig+1;
%                 frame = getframe(gcf);
%                 writeVideo(v,frame);
%             end
%             close(v);
% %%
%           
%         end

%         Find maximum lag and correlation corresponding to the maximum lag
%         if size(probeCh,2) == 1
%              cc(:,~clipMaskCortexR) = NaN;
%             [maxCorr{iDate,iRun},ind] = max(cc,[],1,'omitnan');
%             maxLag{iDate,iRun} = lags(ind);
% 
%             figure; ccR = reshape(maxCorr{iDate,iRun},[361 438]);
%             greenR = ind2rgb(imresize(greenIm{iDate,iRun},1/3),gray(256));
%             imagesc(ccR); hold on; axis image off; colormap jet; colorbar; caxis([-0.5 0.5]); title([ runName ': Maximum correlation']) ;
% 
%             %             %,'alphadata',ccR.*5); colormap jet;
%             %             imagesc(ccR,'alphadata',ccR.*5); colormap jet; colorbar; caxis([-1 1]);
% 
% 
%             figure; lagR = reshape(maxLag{iDate,iRun},[361 438]);
%             imagesc(greenR); hold on;  axis image off;
%             imagesc(lagR);colormap jet; colorbar; caxis([-10 10]);title([runName ': Lag at which maximum correlation occured']);  %,'alphadata',lagR.*5);
%         else
%              cc(:,:,~clipMaskCortexR) = NaN;
%         end

toc;

%%
%   clear corrinfra10
%             corrInfra10 = corr(imDataClean{iDate,iRun}', infraBandLimited{iDate,iRun,iBand}(1:size(imData10,2),:));
%             corrInfra10Im{iDate,iRun,iBand} = reshape(corrInfra10,[361 438 size(corrInfra10,2)]);
% meanCorrInfra10 = cellfun(@(x)mean(x,3),corrInfra10Im, 'un',0);
% bandLabels = {'Alpha'; 'Beta'; 'Gamma'};
% for iDate = 1%:size(allDates,1)
%     expDate = allDates(iDate,:);
%     for iRun = 1%:size(allRuns{iDate,1},1)
%
%         for iBand = 1:3
%             figure; greenImRGB = ind2rgb(greenIm{iDate,iRun},gray(256));
%             imagesc(imresize(greenImRGB,[361 438]));  hold on; axis image;
%             imagesc(imgaussfilt(meanCorrInfra10{iDate,iRun,iBand},3),'alphadata',imgaussfilt(meanCorrInfra10{iDate,iRun,iBand}.*4,1)); colormap jet; caxis([0 0.5]); colorbar;
%             title(strrep([expDate ' Run: ' num2str(iRun) ' Band: ' bandLabels{iBand}],'_','\_'));
%         end
%     end
% end
%% Check for the quality of imaging data by obtaining FC maps
clear downSampledDat corrMap 

spatialBin = 3; iDate = 2; iRun =1; 

figure; imagesc(greenIm{iDate,iRun}); colormap gray; axis image off; hold on;
disp('Pick a seed pixel for getting a FC map ');
seedLoc = ginput(1);
seedLoc = fliplr(round(seedLoc));
plot(seedLoc(2),seedLoc(1),'w.','MarkerSize',15); close gcf;

if exist([templateDir '\' allDates(iDate,:) '\' allRuns{iDate,1}(iRun,:) '\clipMask0' allRuns{iDate,1}(iRun,end) '.BMP'],'file') == 0
    clipMask = imread([templateDir  allDates(iDate,:) '\' allRuns{iDate,1}(iRun,:) '\clipMask0' allRuns{iDate,1}(iRun,end) '.png']);
else
    clipMask = imread([templateDir allDates(iDate,:) '\' allRuns{iDate,1}(iRun,:) '\clipMask0' allRuns{iDate,1}(iRun,end) '.bmp']);
end
clipMask   = imresize(clipMask(:,:,1),1/spatialBin);
%

% downSampledDat =  imData10{iDate,iRun};
% downSampledDat = reshape(downSampledDat,[361 438 size(downSampledDat,2)]); 
%%
downSampledDat =  processedDat{iDate,iRun};
seedSig = calculateSeedSignal(imresize(greenIm{iDate,iRun},1/spatialBin), clipMask,fliplr(round(seedLoc./spatialBin)),12,downSampledDat); % Get Gaussian weighted seed signal
corrMap = plotCorrMap(seedSig,downSampledDat,0);

%% Overlay correlation map on top of greens

greenImRGB = ind2rgb(greenIm{iDate,iRun},gray(256));
frameSize = size(greenImRGB);
figure; imagesc(greenImRGB);
hold on;
corrMapR = imresize(corrMap,spatialBin,'OutputSize',[frameSize(1) frameSize(2)]);
imagesc(corrMapR,'alphadata',corrMapR.*0.8);
colormap jet; colorbar; caxis([0 1]);
axis image off;
plot(seedLoc(2),seedLoc(1),'w.','MarkerSize',15);


%% Select pixels in and out of patches
% allBT = allBadTimes{1,1};
% allBT = unique(round(allBT./100)); % 100 fold downsampling of bad time points

figure; imagesc(greenIm{1,1}); colormap gray; axis image off; hold on;
title('Pick a seed pixel that is inside a patch ');
seedLocIn = ginput(1);
% [seedLocInMaster,~] = cpselect(greenMapRef,greenIm{iDate,1},'Wait',true);
seedLocIn = fliplr(round(seedLocIn./3));
% seedLocInMaster = fliplr(round(seedLocInMaster));

seedRad = 12;
inDat = processedDat(seedLocIn(2)-seedRad:seedLocIn(2)+seedRad,seedLocIn(1)-seedRad:seedLocIn(1)+seedRad,:,1,1);

title('Pick a seed pixel that is outside a patch');
% [seedLocOutMaster,~] = cpselect(greenMapRef,greenIm{iDate,1},'Wait',true);
% seedLocOutMaster = fliplr(round(seedLocOutMaster));
seedLocOut = ginput(1); close gcf;
seedLocOut = fliplr(round(seedLocOut./3));
outDat = processedDat(seedLocOut(2)-seedRad:seedLocOut(2)+seedRad,seedLocOut(1)-seedRad:seedLocOut(1)+seedRad,:,1,1);

% % Upsample imaging data from the superpixels
% % Upsample imaging data and obtain median time course to 10 Hz
% inDat10 = interp(median(reshape(inDat,[numel(inDat(:,:,1)) 1810])',2,'omitnan'),5);
% outDat10 = interp(median(reshape(outDat,[numel(inDat(:,:,1)) 1810])',2,'omitnan'),5);
%
% % Upsample median time course to 1kHz
% inDat1k = interp(inDat10,100);
% outDat1k = interp(outDat10,100);
%
% badTimeSample = allBadTimes{1,1};
% inDat1k(badTimeSample(badTimeSample<=length(inDat1k))) = [];
% outDat1k(badTimeSample(badTimeSample<=length(outDat1k))) = [];
%
% % Correlate upsampled imaging data with neural data
%
% gEnvelopeIn = corr(envelopeBandLimited(1:length(inDat1k),:),inDat1k);
% gEnvelopeOut = corr(envelopeBandLimited(1:length(outDat1k),:),outDat1k);
%
% gInfraIn  = corr(infraBandLimited(1:length(inDat1k),:),inDat1k);
% gInfraOut = corr(infraBandLimited(1:length(outDat1k),:),outDat1k);

% Upsample imaging data and obtain median time course
inDat10 = interp(median(reshape(inDat,[numel(inDat(:,:,1)) 1810]),1,'omitnan'),5);
outDat10 = interp(median(reshape(outDat,[numel(inDat(:,:,1)) 1810]),1,'omitnan'),5);

% Remove bad frames (downsampled bad time segments)
inDat10(allBT(allBT<=length(inDat10)))= [];
outDat10(allBT(allBT<=length(outDat10))) = [];


% Correlate envelope and infra slow fluctuations with imaging data
infra10Hz = infraBandLimited{1,1,3};
infra10Hz = infra10Hz(1:length(inDat10),:);
gEnvelopeIn = corr(infra10Hz(1:length(inDat10),:),inDat10');
gEnvelopeOut = corr(infra10Hz(1:length(inDat10),:),outDat10');
%%
clear cc lags ccOut lagsOut
for iCh = 1:32
    [cc(:,iCh),lags(:,iCh)] = xcorr(infraGammaEphys(:,iCh),inDat10(:),'normalized');
    [ccOut(:,iCh),lagsOut(:,iCh)] = xcorr(infraGammaEphys(:,iCh),outDat10(:),'normalized');
end
[maxC,ind] = max(cc,[],1);
maxLag = lags(ind);

[maxCOut,ind] = max(ccOut,[],1);
maxLagOut = lagsOut(ind);


%%

corrInfraSuper = imgaussfilt(median(corrInfra10R(:,:,1:10),3,'omitnan'),2);
corrInfraDeep  = imgaussfilt(mean(corrInfra10R(:,:,14:20),3,'omitnan'),2);
corrInfra4     = imgaussfilt(mean(corrInfra10R(:,:,11:13),3,'omitnan'),2);

greenImRGB = ind2rgb(greenIm{iDate,iRun},gray(256));
figure; imagesc(imresize(greenImRGB,[361 438]));
hold on; axis image;
imagesc(corrInfraSuper,'alphadata',corrInfraSuper.*4); colormap jet; caxis([0 0.4]); colorbar; title ('Superficial channels');

figure; imagesc(imresize(greenImRGB,[361 438]));
hold on; axis image;
imagesc(corrInfraDeep,'alphadata',corrInfraDeep.*4); colormap jet; caxis([0 0.4]); colorbar;title ('Deep channels');

figure; imagesc(imresize(greenImRGB,[361 438]));
hold on; axis image;
imagesc(corrInfra4,'alphadata',corrInfra4.*10); colormap jet; caxis([-0.2 0.2]); colorbar;title ('Layer 4 channels');

%% Moving average correlation
seedRad = 12;
if exist([templateDir '\' allDates(iDate,:) '\' allRuns{iDate,1}(iRun,:) '\clipMask0' allRuns{iDate,1}(iRun,end) '.BMP'],'file') == 0
    clipMask = imread([templateDir  allDates(iDate,:) '\' allRuns{iDate,1}(iRun,:) '\clipMask0' allRuns{iDate,1}(iRun,end) '.png']);
else
    clipMask = imread([templateDir allDates(iDate,:) '\' allRuns{iDate,1}(iRun,:) '\clipMask0' allRuns{iDate,1}(iRun,end) '.bmp']);
end
clipMask   = imresize(clipMask(:,:,1),1/spatialBin);
corrValsAll = zeros(361,438,31);
imDat10_R = reshape(imData10,[361 438 8722]);
greenR = imresize(greenIm{1,1},[361 438]);
infra10HzR = single(infra10Hz(1:size(imDat10_R,3),:));

%%
L= size(imDat10_R,2);
gaussianDisk_R = fspecial('gaussian',2*seedRad+1,seedRad);
tic;
%  seedSig = zeros(size(imDat10_R,3),1);
parfor iR = 1:size(imDat10_R,1)
    for iC = 1:L
        seed = ([iC,iR]);
        clc; disp(['Seed ' num2str(seed)]);
        if iR<= seedRad || iC<= seedRad;   continue; end
        if seed(2)+seedRad> size(imDat10_R,1); continue; end
        if seed(1)+seedRad> size(imDat10_R,2); continue; end
        seedSig = calculateSeedSignal(greenR, clipMask, seed, seedRad,imDat10_R);
        corrValsAll(iR,iC,:) = corr(seedSig, infra10HzR);
    end
end
toc;

% inDat = processedDat(seedLocIn(2)-seedRad:seedLocIn(2)+seedRad,seedLocIn(1)-seedRad:seedLocIn(1)+seedRad,:,1,1);



% gInfraIn  = corr(infraBandLimited(1:length(inDat10),:),inDat10);
% gInfraOut = corr(infraBandLimited(1:length(inDat10),:),outDat10);

%%
% clear movePointsTemp fixedPointsTemp
% [movPointsTemp,fixedPointsTemp] = cpselect(greenIm{iDate,1},greenMapRef,'Wait',true);
% movPoints{iDate,1} = movPointsTemp;
% fixedPoints{iDate,1} = fixedPointsTemp;
%
% temp = zeros(size(greenMapRef));
% temp(1:size(greenIm{iDate,1},1),1:size(greenIm{iDate,1},2)) = greenIm{iDate,1};
%
% [OTemp,spacingTemp,~] = point_registration(size(greenMapRef),fliplr(fixedPoints{iDate,1}),fliplr(movPoints{iDate,1}));
% [imRTemp{iDate,1},tTemp] = bspline_transform(OTemp,temp,spacingTemp,3); % bicubic spline interpolation
%
% O{iDate,1} = OTemp;
% spacing{iDate,1} = spacingTemp;
% figure; imshowpair(greenMapRef,imRTemp{iDate,1},'blend'); %colormap gray;
% % title(['Green map from Date: ' allDates(iDate,:)]);
% %             t{iDate,1} = tTemp;
% save([templateDir expDate '\imageRegisPoints_RS_Ephys.mat'],'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp');

%% Obtain FC map for location where electrode is placed
figure; imagesc(greenMapRef);axis image off; colormap gray; hold on;
title('Select the reference site');
refSeed(1,:,:) = ginput(1);

[rsConnMatrix,corrMapFinal] = getRSConnectivityMaps(squeeze(refSeed(1,:,:))',monkeyName);

%%  Plot greens with superpixels
iDate = 1; iRun = 1;
greenImRGB = ind2rgb(greenIm{iDate,iRun},gray(256));
figure; imagesc(imresize(greenImRGB,[361 438]));
hold on; axis image;

in = [seedLocIn(2)-seedRad:seedLocIn(2)+seedRad; seedLocIn(1)-seedRad:seedLocIn(1)+seedRad];
out = [seedLocOut(2)-seedRad:seedLocOut(2)+seedRad; seedLocOut(1)-seedRad:seedLocOut(1)+seedRad];
% plot([in(2,1) in(2,1)+(seedRad*2) in(2,1)+(seedRad*2) in(2,1) in(2,1)],[in(1,1) in(1,1) in(1,1)+(seedRad*2) in(1,1)+(seedRad*2) in(1,1)],'LineWidth',4);
% plot( [out(2,1) out(2,1)+(seedRad*2) out(2,1)+(seedRad*2) out(2,1) out(2,1)],[out(1,1) out(1,1) out(1,1)+(seedRad*2) out(1,1)+(seedRad*2) out(1,1)],'LineWidth',4);

plot([in(1,1) in(1,1) in(1,1)+(seedRad*2) in(1,1)+(seedRad*2) in(1,1)],[in(2,1) in(2,1)+(seedRad*2) in(2,1)+(seedRad*2) in(2,1) in(2,1)],'LineWidth',4);
plot( [out(1,1) out(1,1) out(1,1)+(seedRad*2) out(1,1)+(seedRad*2) out(1,1)],[out(2,1) out(2,1)+(seedRad*2) out(2,1)+(seedRad*2) out(2,1) out(2,1)],'LineWidth',4);

%%  Plot FC map with superpixels
figure; imagesc(ind2rgb(greenMapRef,gray(256))); %imagesc(imresize(greenMapRef,[361 438]));
hold on; axis image;
imagesc(corrMapFinal,'AlphaData',corrMapFinal.*1.2); colormap jet;
seedRadR = seedRad*3;
in = [seedLocInMaster(2)-seedRadR:seedLocInMaster(2)+seedRadR; seedLocInMaster(1)-seedRadR:seedLocInMaster(1)+seedRadR];
out = [seedLocOutMaster(2)-seedRadR:seedLocOutMaster(2)+seedRadR; seedLocOutMaster(1)-seedRadR:seedLocOutMaster(1)+seedRadR];
% plot([in(2,1) in(2,1)+(seedRad*2) in(2,1)+(seedRad*2) in(2,1) in(2,1)],[in(1,1) in(1,1) in(1,1)+(seedRad*2) in(1,1)+(seedRad*2) in(1,1)],'LineWidth',4);
% plot( [out(2,1) out(2,1)+(seedRad*2) out(2,1)+(seedRad*2) out(2,1) out(2,1)],[out(1,1) out(1,1) out(1,1)+(seedRad*2) out(1,1)+(seedRad*2) out(1,1)],'LineWidth',4);

plot([in(1,1) in(1,1) in(1,1)+(seedRadR*2) in(1,1)+(seedRadR*2) in(1,1)],[in(2,1) in(2,1)+(seedRadR*2) in(2,1)+(seedRadR*2) in(2,1) in(2,1)],'LineWidth',2);
plot( [out(1,1) out(1,1) out(1,1)+(seedRadR*2) out(1,1)+(seedRadR*2) out(1,1)],[out(2,1) out(2,1)+(seedRadR*2) out(2,1)+(seedRadR*2) out(2,1) out(2,1)],'LineWidth',2);


%%
figure; subplot(121); imagesc(imgaussfilt(gEnvelopeIn,1)); colormap jet; title('Test pixels');
colorbar; caxis([-0.2 0.2]);
subplot(122); imagesc(imgaussfilt(gEnvelopeOut,1)); colormap jet; colorbar; caxis([-0.2 0.2]); title('Control pixels');
sgtitle('Gamma power vs RS hemodynamics');

figure; subplot(121); imagesc(imgaussfilt((gInfraIn),1)); colormap jet;
colorbar;caxis([-0.2 0.2]); title('Test pixels');
subplot(122); imagesc(imgaussfilt((gInfraOut),1)); colormap jet; colorbar;caxis([-0.2 0.2]); title('Control pixels');
sgtitle('Infra slow Gamma fluctuations vs RS hemodynamics');

%%
figure;  plot(movmean(gEnvelopeIn,3),1:length(gEnvelopeIn),'LineWidth',2);   hold on;
% colorbar; caxis([-0.02 0.02]);
plot(movmean(gEnvelopeOut,3),1:length(gEnvelopeOut),'LineWidth',2); set(gca,'YDir','reverse');
sgtitle('Infra slowGamma fluctuations vs RS hemodynamics'); xlim([-0.08 0.08]);
%%

figure;  plot(movmean(gInfraIn,2),1:length(gInfraIn),'LineWidth',2); hold on;
% colorbar; caxis([-0.02 0.02]);
plot(movmean(gInfraOut,2),1:length(gInfraIn),'LineWidth',2); set(gca,'YDir','reverse');
sgtitle('Infra slowGamma fluctuations vs RS hemodynamics');xlim([-0.08 0.08]);








%% Calculating the averaged seed signal over a few pixels
function seedSig = calculateSeedSignal(green, clipMask, seed, seedRad,frames)
% choose seed location (gaussian average of a radius around chosen point)
% seedRad = 6; % 6 = 250um radius, 12 = 500um radius, 24 = 1000um radius
clear gaussianDisk green_Seed clipMask_Seed
if exist('seedRad','var') == 0; seedRad = 6; end
seedSig = zeros(size(frames,3),1);

gaussianDisk = fspecial('gaussian',2*seedRad+1,seedRad);
green_seed = green(seed(2)-seedRad:seed(2)+seedRad,seed(1)-seedRad:seed(1)+seedRad,:);
clipMask_seed = clipMask(seed(2)-seedRad:seed(2)+seedRad,seed(1)-seedRad:seed(1)+seedRad);
gaussianDisk = gaussianDisk .* double(clipMask_seed);
gaussianDisk = gaussianDisk ./ sum(gaussianDisk(:));

for x = 1:2*seedRad+1     % rows are observations, columns are variables
    for y = 1:2*seedRad+1
        seedSig = seedSig + gaussianDisk(x,y) * squeeze(frames(seed(2)+x-seedRad-1,seed(1)+y-seedRad-1,:));

    end
end
end

