% ephysImaging_v2.m
% This script analyzes electrophysiological and imaging data recorded
% simultaneously for one monkey
% December 13, 2023 - KM
% See ephysImaging.m for previous versions
clear;clc
commonDir = 'C:\Users\KEM294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\Imaging']));
addpath(genpath([commonDir '\Codes\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\fly_track\videoIO']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\spectral_analysis\continuous\dupes']));

%% Initializing all the relevant variables
monkeyName   = 'Whiskey';
spatialBin   = 3;
hemisphere   = 'Left';
saveFlag     = 0;
checkBadRuns = 0;
fs           = 1e3; % LFP

% Ephys filtering parameters
gammaBand   = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); % Gamma band filtering parameters
alphaBand   = [8 12];  [bA,aA] = butter(3,alphaBand./(fs/2),'bandpass'); % Alpha band filtering parameters
betaBand    = [13 30]; [bB,aB] = butter(3,betaBand./(fs/2),'bandpass');
chOutCortex = 1:3;
chDeep      = 30:32;
[z,p,k]     = cheby1(2,1,[0.01 0.1]./(fs/2),'bandpass');
[sos,g]     = zp2sos(z,p,k);

% Load all the runs, experiment dates, green reference images, reference directory for the corresponding monkey

if strcmp(monkeyName,'CharlieSheen') % Charlie Sheen
    allDates     = ['11_29_2021'; '01_11_2022'; '08_07_2023';'11_20_2023'];
    allRuns      = {['run00'; 'run01']; ['run03'; 'run04']; ['run02'; 'run05';'run06'; 'run07'; 'run08'];...
        ['run02'; 'run05'; 'run06'; 'run07' ]};

    % Get imaging parameters
    refDate      = '08_31_2021';
    refDir       = [commonDir '\' monkeyName '_SqM\' hemisphere ' Hemisphere\' refDate '\Master Green Images\'];
    refImageName = 'Charlie Sheen Combined Green 08_31_2021';
    lensCombo        = {['50/50 MF'; '50/50 MF'];['20/50';'20/50'];['50/50 MF';'50/50 MF';'50/50 MF'; '50/50 MF'; '50/50 MF'];...
        ['50/50 MF';'50/80 MF'; '50/80 MF';'50/80 MF']};
    roiSize          = {[64 64];[24 24];[64 64 64 64 64];[38 38 38 38];}; %  500um x 500um

    % Get ephys parameters
    ephysFileNameAll = {['run00'; 'run01'];['run0003';'run0004'];['datafile0002'; 'datafile0005';'datafile0006';...
        'datafile0007';'datafile0008']; ['datafile0002';'datafile0005'; 'datafile0006'; 'datafile0007']};
    serverPath       = '\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\302-19_CharlieSheen_SqM\Left Hemisphere\';
    probeLabel       = {[]; []; ['CDE1';'CDE1';'CDE1';'CDE1';'CDE1';'CDE1']; ['B25A';'B25A'; 'B25A'; 'B25A']};
    chInCortexNotes  = {[1 1];[1 1];[7 10 8 6 4 7];[3 4 5 7]}; % transition determined during experiment


elseif strcmp(monkeyName,'Whiskey') % Whiskey
    allDates = ['08_14_2023'; '10_16_2023'; '12_04_2023';'02_20_2024'];
    allRuns  = {['run02'; 'run04'; 'run05';];[ 'run03'; 'run04'; 'run05'; ...
        'run06'; 'run07'; 'run08'; 'run09'];['run04'; 'run05';'run06'; 'run07']; ...
        ['run01'; 'run02';'run03'; 'run04'; 'run05'; 'run06'; 'run07']};

    % Get imaging parameters
    refDate      = '05_09_2022';
    refDir       = [commonDir '\' monkeyName '_SqM\' hemisphere ' Hemisphere\' refDate '\Master Green Images\'];
    refImageName = 'Combined Green';
    lensCombo        = {['28/50'; '28/50'; '28/50'];[ '50/80 MF'; '50/80 MF'; '50/80 MF'; ...
        '50/80 MF'; '50/80 MF'; '50/80 MF'; '50/80 MF'];['50/80 MF'; '50/80 MF';'50/80 MF'; '50/80 MF'];...
        ['50/80 MF'; '50/80 MF';'50/80 MF'; '50/80 MF';'50/80 MF' ; '50/80 MF'; '50/80 MF']};
    roiSize          = {[36 36 36];[38 38 38 38 38 38 38];[38 38 38 38];[38 38 38 38 38 38 38]}; % for 500 um x 500 um

    % Get ephys parameters
    ephysFileNameAll = {['datafile0002'; 'datafile0004'; 'datafile0005']; ['datafile0003' ;...
        'datafile0004'; 'datafile0005'; 'datafile0006'; 'datafile0007'; 'datafile0008' ;  'datafile0009']; ...
        ['datafile0004'; 'datafile0005'; 'datafile0006'; 'datafile0007']; ['datafile0002'; 'datafile0003';...
        'datafile0004'; 'datafile0005'; 'datafile0006'; 'datafile0007';'datafile0008' ]};
    serverPath       = '\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\303-19_Whiskey_SqM\Left Hemisphere\';
    probeLabel       = {[]; ['CDE1';'B25A'; 'B25A'; 'B25A'; 'B25A'; 'B25A' ;'B25A'];['B25A'; 'B25A'; '0763'; '0763'];...
        ['0763'; '0763'; '0763'; '0763'; '0763'; '0763'; '0763']};
    chInCortexNotes  = {[1 1 1]; [5 5 7 4 5 5 6];[6 6 8 5]; [7 6 6 8 7 7 9]};  % transition determined during experiment

end

% Get the master green image
if exist([refDir refImageName '.bmp'],'file') == 0 % Make sure to get what the reference run is outside the function
    greenMapRef = imread([refDir refImageName '.png']);
else
    greenMapRef = imread([refDir refImageName '.bmp']);
end

greenMapRef  = greenMapRef(:,:,1);

%% Store/Retrieve imaging data
for iDate = 1:size(allDates,1)
    clear expDate
    expDate    = allDates(iDate,:);

    for iRun = 1:size(allRuns{iDate,1},1)
        % Load all the information
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
        fileInfo    = dir(dataDir);

        % Check if greens exist in the data folder
        if isempty(find(strcmp({fileInfo.name},['green' runName(end-1:end) '.png']), 1)) && isempty(find(strcmp({fileInfo.name},['green' runName(end-1:end) '.bmp']), 1))
            try
                copyfile([serverPath '\' expDate '\' runName '\green' runName(end-1:end) '.png'], dataDir);
            catch
                copyfile([serverPath '\' expDate '\' runName '\green' runName(end-1:end) '.bmp'], dataDir);
            end
        end

        % Check if greens are edited
        if exist([templateDir  expDate '\' allRuns{iDate,1}(iRun,:) '\green0' allRuns{iDate,1}(iRun,5) '_Edited.png'],'file')
            greenTemp = imread([dataDir '\green0' allRuns{iDate,1}(iRun,5) '_Edited.png']);

        elseif exist([templateDir  expDate '\' allRuns{iDate,1}(iRun,:) '\green0' allRuns{iDate,1}(iRun,5) '_Edited.bmp'],'file')
            greenTemp = imread([dataDir '\green0' allRuns{iDate,1}(iRun,5) '_Edited.bmp']);

        else
            error('Greens are not edited...');
        end

        greenIm{iDate,iRun}  = greenTemp(:,:,1); % Getting greens

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
            [~,~,~,spSmooth,tempBandPass] = getPreProcessedDataRestingState(serverDataPath,dataDir,runName,numFiles,spatialBin,datName);
            disp(['Storing imaging data for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)]);
            save([dataDir '\processedFrames.mat'],'tempBandPass');
            processedDat{iDate,iRun} = matfile([dataDir '\processedFrames.mat']); clear tempBandPass;

        else
            % Retrieve processed data
            disp(['Loading imaging data for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)]);
            processedDat{iDate,iRun} = matfile([dataDir '\processedFrames.mat']);
        end

        cd(commonDir);

        % Get the functional connectivity map for the run to assess the
        % quality of imaging...
        clear clipMask pDatTemp seedLoc corrMap corrMapR
        if exist([dataDir '\clipMask0' runName(end) '.BMP'],'file') == 0 % Clipmask
            clipMask = imread([dataDir '\clipMask0' runName(end) '.png']);
        else
            clipMask = imread([dataDir '\clipMask0' runName(end) '.bmp']);
        end

        clipMask = imresize(clipMask,1/3); % Resize mask
        clipMask = clipMask(:,:,1)>0; % Converting to 0s and 1s

        if ~exist([dataDir '\FCMap.png'],'file')|| ~exist([dataDir '\FCMapQC.png'],'file')
            clear pDatTemp
            pDatTemp = processedDat{iDate,iRun}.tempBandPass;
        end

        if ~exist([dataDir '\FCMap.png'],'file') % Near the electrode
            figure('units','normalized','outerposition',[0 0 1 1]);
            imagesc(greenIm{iDate,iRun}); colormap gray; axis image off;
            title(strrep(['Pick a seed to get a FC map for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)],'_','\_')); hold on;
            seedLoc = ginput(1);
            seedLoc = fliplr(round(seedLoc));
            plot(seedLoc(2),seedLoc(1),'w.','MarkerSize',15); close gcf;

            seedSig = calculateSeedSignal(imresize(greenIm{iDate,iRun},1/spatialBin), clipMask,fliplr(round(seedLoc./spatialBin)),12,pDatTemp); % Get Gaussian weighted seed signal
            corrMap = plotCorrMap(seedSig,pDatTemp,0);

            % Plotting FC maps...
            greenImRGB = ind2rgb(greenIm{iDate,iRun},gray(256));
            frameSize = size(greenImRGB);

            figure('units','normalized','outerposition',[0 0 1 1]); imagesc(greenImRGB); hold on;
            corrMapR = imresize(corrMap,spatialBin,'OutputSize',[frameSize(1) frameSize(2)]);
            imagesc(corrMapR,'alphadata',corrMapR.*0.8);
            colormap jet; colorbar; caxis([0 1]); axis image off;
            plot(seedLoc(2),seedLoc(1),'w.','MarkerSize',35);

            sgtitle(strrep([monkeyName ' Date: ' expDate ' - Run: ' runName(end)],'_','\_'));
            f = gcf; exportgraphics(f,[dataDir '\FCMap.png'],'Resolution',300); close gcf;
        end

        if ~exist([dataDir '\FCMapQC.png'],'file') % Imaging quality check
            figure('units','normalized','outerposition',[0 0 1 1]);
            imagesc(greenIm{iDate,iRun}); colormap gray; axis image off;
            title(strrep(['Pick a seed to get a FC map to check imaging quality for ' monkeyName ],'_','\_')); hold on;
            seedLoc = ginput(1);
            seedLoc = fliplr(round(seedLoc));
            plot(seedLoc(2),seedLoc(1),'w.','MarkerSize',15); close gcf;

            seedSig = calculateSeedSignal(imresize(greenIm{iDate,iRun},1/spatialBin), clipMask,fliplr(round(seedLoc./spatialBin)),12,pDatTemp); % Get Gaussian weighted seed signal
            corrMap = plotCorrMap(seedSig,pDatTemp,0);

            % Plotting FC maps...
            greenImRGB = ind2rgb(greenIm{iDate,iRun},gray(256));
            frameSize = size(greenImRGB);

            figure('units','normalized','outerposition',[0 0 1 1]); imagesc(greenImRGB); hold on;
            corrMapR = imresize(corrMap,spatialBin,'OutputSize',[frameSize(1) frameSize(2)]);
            imagesc(corrMapR,'alphadata',corrMapR.*0.8);
            colormap jet; colorbar; caxis([0 1]); axis image off;
            plot(seedLoc(2),seedLoc(1),'w.','MarkerSize',35);

            sgtitle(strrep([monkeyName ' Date: ' expDate ' - Run: ' runName(end)],'_','\_'));
            f = gcf; exportgraphics(f,[dataDir '\FCMapQC.png'],'Resolution',300); close gcf;
        end
        clear pDatTemp

        % Get the functional connectivity map averaged across multiple
        % resting state runs
        if ~exist([dataDir '\AvgRS_FCMap.png'],'file')
            figure('units','normalized','outerposition',[0 0 1 1]);
            imagesc(greenMapRef);axis image off; colormap gray; hold on;
            title(strrep(['Select the location where the probe is located for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)],'_','\_'));
            refSeed(1,:,:) = ginput(1); close gcf;
            [rsConnMatrix,corrMapFinal] = getRSConnectivityMaps(squeeze(refSeed(1,:,:))',monkeyName);

            figure('units','normalized','outerposition',[0 0 1 1]);
            imagesc(ind2rgb(greenMapRef,gray(256))); hold on;
            imagesc(corrMapFinal,'alphaData',corrMapFinal.*0.9); axis image off; hold on;
            plot(refSeed(1,:,1),refSeed(1,:,2),'w.','MarkerSize',35);
            colormap jet; caxis([0 1]); f = gcf;
            sgtitle(strrep([monkeyName ' Date: ' expDate ' - Run: ' runName(end)],'_','\_'));
            exportgraphics(f,[dataDir '\AvgRS_FCMap.png'],'Resolution',300); close gcf;
        end
    end
end

%% Store/Retrieve LFP data
for iDate = 1:size(allDates,1)
    clear expDate
    expDate    = allDates(iDate,:);
    for iRun = 1:size(allRuns{iDate,1},1)
        clear entityInfo datFileName timeStamp
        % Load all run information
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

            [probe{iRun,iDate}, eeg{iRun,iDate}, cameraTimeStamp{iRun,iDate}] = saveLFPSingleProbe(runName,datName,...
                saveFolder,datFileName,fileNum,fs);
        else
            % Retrieve LFP Data
            disp(['Retrieving electrophysiology data for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)]);
            load([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat']); % Use matfile here?
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

%% Determining bad time segments, channels and transition channels for LFP time series
fs              = 1000;
params.Fs       = fs;
params.fpass    = [1 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;

for iDate =1:size(allDates,1)
    clear expDate
    expDate    = allDates(iDate,:);

    for iRun =  1:size(allRuns{iDate,1},1)
        clear  runName dataDir probeCh transitionCh probeTemp
        expDate    = allDates(iDate,:);
        badTimeInd = []; badTimes = [];
        runName    = allRuns{iDate,1}(iRun,:);
        dataDir    = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

        % Determine the bad time segments and bad channels in LFP data
        disp(['Identifying noisy data from LFP for:  ' allDates(iDate,:) ' ' runName]);
        probeCh      = probe{iRun,iDate};
        channels     = 1:size(probeCh,2);
        transitionCh = chInCortexNotes{iDate}(iRun);
        badChVal     = [];

        % Remove bad channels that are already known from probe datasheet
        if ~isempty(probeLabel{iDate,1})
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

        % Obtain the spectrogram and determine bad channels
        [spec,~,~] = mtspecgramc(probeCh(:,channels),[5 2],params); % Get the spectrogram
        powTimeBin = squeeze(sum(10.*log10(abs(spec)),2));          % Sum of powers

        if length(channels)~= 1       % Determine the minimum and maximum threshold to determine bad channels
            badElecThreshHigh = (median(powTimeBin,2)+5*mad(powTimeBin,1,2));
            badElecThreshLow  = (median(powTimeBin,2)-5*mad(powTimeBin,1,2));

            badChVal = [badChVal channels(sum((powTimeBin>badElecThreshHigh),1) >= floor(0.75*size(powTimeBin,1)) |(sum((powTimeBin<badElecThreshLow),1) >= floor(0.75*size(powTimeBin,1))))];
            if ~isempty(badChVal); channels(ismember(channels,badChVal)) = []; end
        end

        % Remove bad channels from the probe
        badCh{iDate,iRun} = badChVal;
        probeCh(:,badChVal) = [];

        % Determine the bad time segments
        clear probeBL
        probeBL = single(filtfilt(bG,aG,double(probeCh))); % Filtering the gamma band

        % Bad time segment threshold bounds
        if length(channels)~= 1 % Linear array
            badTimeThreshHigh = mean(probeBL(:,15:end),'all','omitnan') + 5*std(probeBL(:,15:end),[],[1,2]);
            badTimeThreshLow  = mean(probeBL(:,15:end),'all','omitnan') - 5*std(probeBL(:,15:end),[],[1,2]);
        else % Single probe
            badTimeThreshHigh = mean(probeBL,'omitnan') + 5*std(probeBL,[]);
            badTimeThreshLow= mean(probeBL,'omitnan')- 5*std(probeBL,[]);
        end

        % Determine threshold crossings...
        badTimeIndOld = find(mean(probeBL(:,15:end),2,'omitnan')>badTimeThreshHigh | mean(probeBL(:,15:end),2,'omitnan')<badTimeThreshLow);

        % Determine the bad time segments
        badTimes = [];

        if ~isempty(badTimeIndOld)
            badTimeInd = [(badTimeIndOld-fs/10)  (badTimeIndOld+fs/10)]; % Taking 100 ms before and after each threshold crossing

            for iL = 1:size(badTimeInd,1)
                badTimes = [badTimes badTimeInd(iL,1): badTimeInd(iL,2)];
            end
        end

        badTimes = unique(badTimes);
        badTimesLFP{iDate,iRun} = badTimes;

        % *Checkpoint*: Get the spectrograms before and after bad time
        % segment removal to check if noise has been removed.
        if ~exist([dataDir '\AverageLFPSpec.png'],'file')% Plotting and saving the spectrogram
            clear spec timeValsSpec freqValsSpec

            if length(channels)~= 1 % Linear array
                [spec,timeValsSpec,freqValsSpec] = mtspecgramc(probeCh(:,15:end),[5 2],params);
                figure; subplot(211);
                imagesc(timeValsSpec,freqValsSpec, 10.*log10(squeeze(mean(spec,3)')));
                caxis([-20 20]); set(gca,'YDir','normal'); colormap jet; colorbar;
                xlabel('Time (s)'); ylabel('Power (dB)'); title('Before removing bad time segments');

                if ~isempty(badTimes)
                    probeChTemp = probeCh;
                    probeChTemp(badTimes,:) = []; % Removing bad Time segments...

                    clear spec timeValsSpec freqValsSpec
                    [spec,timeValsSpec,freqValsSpec] = mtspecgramc(probeChTemp(:,15:end),[5 2],params);
                    subplot(212); imagesc(timeValsSpec,freqValsSpec, 10.*log10(squeeze(mean(spec,3)')));
                    caxis([-20 20]); set(gca,'YDir','normal');colormap jet; colorbar;
                    xlabel('Time (s)'); ylabel('Power (dB)'); title('After removing bad time segments');
                else
                    subplot(212); imagesc(timeValsSpec,freqValsSpec, 10.*log10(squeeze(mean(spec(:,:,15:end),3)')));
                    caxis([-20 20]); set(gca,'YDir','normal'); title('No bad time segments detected');colormap jet; colorbar;
                end

            else % Single probe
                probeChTemp = probeCh;
                [spec,timeValsSpec,freqValsSpec] = mtspecgramc(probeChTemp,[5 2],params);
                figure; subplot(211); imagesc(timeValsSpec,freqValsSpec, 10.*log10(spec'));
                caxis([-20 20]); set(gca,'YDir','normal'); xlabel('Time (s)'); ylabel('Power (dB)'); colormap jet; colorbar;
                if strcmp(expDate,'08_14_2023'); caxis([-60 60]); end

                if ~isempty(badTimes)
                    probeChTemp(badTimes,:) = []; % Removing bad Time segments...

                    clear spec timeValsSpec freqValsSpec
                    [spec,timeValsSpec,freqValsSpec] = mtspecgramc(probeChTemp,[5 2],params);
                    subplot(212);imagesc(timeValsSpec,freqValsSpec, 10.*log10(spec'));
                    caxis([-20 20]); set(gca,'YDir','normal'); colormap jet;
                    xlabel('Time (s)'); ylabel('Power (dB)'); title('After removing bad time segments'); colorbar;
                    if strcmp(expDate,'08_14_2023'); caxis([-60 60]); end
                else
                    subplot(212); imagesc(timeValsSpec,freqValsSpec, 10.*log10(spec'));
                    caxis([-20 20]); set(gca,'YDir','normal'); title('No bad time segments detected'); colormap jet; colorbar;
                    if strcmp(expDate,'08_14_2023'); caxis([-60 60]); end
                end
            end
            sgtitle(strrep([monkeyName ' Date: ' expDate ' - Run: ' runName(end)],'_','\_'));
            f = gcf; exportgraphics(f,[dataDir '\AverageLFPSpec.png'],'Resolution',300); close gcf;
        end

        probeTemp = probeCh;
        probeTemp(badTimes,:) = [];

        % Determine the transition channels for LFP data
        if size(probeTemp,2)~=1

            % Get mean intra probe marginals from wideband range for Probe
            marginalVal = mean(imgaussfilt(corr(probeTemp,'Rows','complete'),1),2,'omitnan');

            if ~exist(fullfile(dataDir,'pairCorr.png'),'file')
                t = tiledlayout(1,3,'TileSpacing','Compact');
                nexttile([1 2]); imagesc(imgaussfilt(corr(probeTemp,'Rows','complete'),1));
                colormap jet; colorbar; axis image tight;
                xticks(1:2:size(probeTemp,2));  yticks(1:2:size(probeTemp,2));
                nexttile; plot(marginalVal,1:length(marginalVal));set(gca,'YDir','reverse');
                axis padded;xlim([0 round((max(marginalVal)+0.1).*10)/10]);
                xticks(0:0.2:round((max(marginalVal)+0.1).*10)/10);

                sgtitle(strrep(['Within probe correlations for: ' monkeyName ' Date: ' expDate ' - Run: ' runName(end)],'_','\_'));
                f = gcf; exportgraphics(f,[dataDir '\pairCorr.png'],'Resolution',300); close gcf;
            end

            % Get the slope of the marginals for Probe A
            fx = abs(movmean(gradient(marginalVal),2,'omitnan'));

            % Find the channel with maximum slope
            clear tempCh
            tempCh = find(fx == max(fx));

            if tempCh>=20 % Set transition to 1 if max slope is after channel 20
                tempCh = 1;
            end

            if abs(tempCh - transitionCh)>=10 % Set transition to the value determined from notes if difference between estimated and observed value exceeds 10
                estChInCortex{iDate}(iRun,1) = transitionCh;

            elseif abs(tempCh - transitionCh)<= 3 % Set transition to estimated value if difference between estimated and observed value is less than or equal to 3
                estChInCortex{iDate}(iRun,1) = tempCh;
            else
                estChInCortex{iDate}(iRun,1) = floor((tempCh + transitionCh)/2); % Set transition to mean between estimated and observed value otherwise
            end

            if estChInCortex{iDate}(iRun,1)+ 20 > size(probeTemp,2) % Limit to the first 20 channels only
                estChInCortex{iDate}(iRun,2) = size(probeTemp,2);
            else
                estChInCortex{iDate}(iRun,2) = estChInCortex{iDate}(iRun,1)+ 20;
            end

        else
            estChInCortex{iDate}(iRun,:) = [1 1];
        end

        % Animal state: plot spectrogram of EEG
        clear eegGood spec timeValsSpec freqValsSpec
        if ~exist(fullfile(dataDir,'eegSpec.png'),'file')
            segLen = 500; winSize = 250;
            eegGood             = eeg{iRun,iDate};
            eegGood(badTimes,:) = [];

            [spec,timeValsSpec,freqValsSpec] = mtspecgramc(eegGood,[5 2],params); % Spectrogram

            figure; imagesc(timeValsSpec,freqValsSpec, 10.*log10((spec)'));
            caxis([-20 40]); set(gca,'YDir','normal');colormap jet; colorbar;
            xlabel('Time (s)'); ylabel('Frequency (Hz)');
            title(strrep([monkeyName ' Date: ' expDate ' Run: ' runName(end) ' - Spectrogram of EEG'],'_','\_'));
            f = gcf; exportgraphics(f,[dataDir '\eegSpec.png'],'Resolution',300); close gcf;
        end
    end
end

%% Computing cross correlations
for iDate = 2%1:size(allDates,1)
    clear expDate
    expDate    = allDates(iDate,:);
    for iRun =5% 1:size(allRuns{iDate,1},1)

        clc; clear  runName dataDir probeCh badTimes badChVal crossCorrAll processedDat10...
            processedDatR probeTemp clipMaskROI grayIm clearMaskFlag varInfo probeLocFlag
        
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
        disp([ monkeyName ' ' allDates(iDate,:) ' ' runName]);
       
        varInfo = who('-file', fullfile(dataDir,'crossCorrROI.mat'));
        if find(ismember(varInfo,'clipMaskROI')); clipMaskFlag = 1; else; clipMaskFlag = 0; end 
       
        clear varInfo; varInfo = who('-file', fullfile(dataDir,'roiCenterLoc.mat'));
        if find(ismember(varInfo,'seedLocProbe')); probeLocFlag = 1; else; probeLocFlag = 0; end

        if ~exist(fullfile(dataDir,'crossCorrROI.mat'),'file') && ~exist(fullfile(dataDir,'crossCorrFOV.mat'),'file') && ~clipMaskFlag || ~probeLocFlag % Check if cross correlations have been computed
            disp('Loading ISOI and LFP data for this run...');
            % EPHYS: Retrieve LFP for the run
            probeCh  = probe{iRun,iDate};
            channels = 1:size(probeCh,2);
            badChVal = badCh{iDate,iRun};

            % EPHYS: Remove bad channels from LFP
            badTimes  = badTimesLFP{iDate,iRun};
            probeCh(:,badChVal)   = [];

            % IMAGING: Load the appropriate masks for the imaging data
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
            processedDatR  = reshape(processedDat{iDate,iRun}.tempBandPass,[361*438 size(processedDat{iDate,iRun},'tempBandPass',3)]);

            % IMAGING: Upsample the imaging dataset
            disp('Upsampling imaging data to 10 Hz...');
            parfor iP = 1:size(processedDatR,1)
                processedDat10(iP,:) = interp(processedDatR(iP,:),5);
            end

            % IMAGING: Remove non-cortical components from imaging data
            clipMaskCortexR = reshape(clipMaskCortex,[size(clipMaskCortex,1)*size(clipMaskCortex,2) 1]);
            processedDat10(~clipMaskCortexR,:) = NaN;
            processedDat10R = reshape(processedDat10,[361 438 size(processedDat10,2)]);

            % IMAGING: Pick the area around the electrode
            seedRad = roiSize{iDate}(iRun);%floor(roiSize{iDate}(iRun)/spatialBin) ;%
            greenFig = imresize(greenIm{iDate,iRun},1/spatialBin,'OutputSize',[361 438]);            
           
        
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

                if ~probeLocFlag
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

            % IMAGING: Show the ROI
            if ~exist(fullfile(dataDir,'ROI.png'),'file') || 1
                figure; imagesc(greenFig); colormap gray; axis image off; hold on;
                rectangle('Position',[ seedLocIn(1)-seedRad/2,seedLocIn(2)-seedRad/2,seedRad,seedRad],'EdgeColor','r','LineWidth',2);
                title('ROI near the electrode');
                f = gcf; exportgraphics(f,[dataDir '\ROI.png'],'Resolution',300); close gcf;
            end

            % IMAGING: Isolate the ROI
            disp('Isolating ROI and removing bad frames...');
            clear inDat inDatSize
            inDat = processedDat10R(seedLocIn(2)-seedRad/2:seedLocIn(2)+seedRad/2,seedLocIn(1)-seedRad/2:seedLocIn(1)+seedRad/2,:);
            inDatSize = size(inDat);
            inDat = reshape(inDat,[inDatSize(1)*inDatSize(2) inDatSize(3)]);

            clipMaskROI = clipMaskCortex(seedLocIn(2)-seedRad/2:seedLocIn(2)+seedRad/2,seedLocIn(1)-seedRad/2:seedLocIn(1)+seedRad/2);

            % Check if the # samples are equal before removing bad time
            % segments from both LFP and imaging data...
            clear szLFP szIm szMin probeTemp

            probeTemp = probeCh;
            szLFP     = size(probeTemp,1);
            szIm      = size(inDat,2)*100;

            if ~(szLFP == szIm) % Make both matrices equal...
                szMin     = min([szLFP, szIm]);
                probeCh   = probeCh(1:szMin,:);
                badTimes(badTimes>szMin) = [];
                inDatOrig = inDat(:,1:round(szMin/100));
            else
                inDatOrig = inDat;
            end

            probeTemp(badTimes,:) = [];

            % IMAGING: Upsample ROI (x100 times) and remove concurrent bad
            % frames
            clear inDatTemp inDatUp
            tic;
            parfor iP = 1:size(inDat,1)
                inDatTemp = interp(inDat(iP,:),100);
                if ~(szLFP == szIm); inDatTemp = inDatTemp(:,1:szMin); end
                inDatTemp(:,badTimes) = []; % Remove bad frames
                inDatUp(iP,:) = downsample(inDatTemp,100); % Downsample imaging data
            end
            toc;

            % IMAGING: Upsampling full field of view and removing bad frames...
            clear varInfo
            varInfo = who('-file', fullfile(dataDir,'crossCorrFOV.mat'));
            if find(ismember(varInfo,'ccFull'));  ccFullFlag = 1;  else; ccFullFlag = 0;  end

            if ~exist(fullfile(dataDir,'crossCorrFOV.mat'),'file') || ~ccFullFlag %|| goodDatFlag
                disp('Upsampling full FOV and removing bad frames...');
                fovFlag = 1;
                clear inDatTemp goodDat ccFull ccFullAlpha
                
                tic;
                parfor iP = 1:size(processedDat10,1)
                    inDatTemp = interp(processedDat10(iP,:),100);
                    if ~(szLFP == szIm); inDatTemp = inDatTemp(:,1:szMin); end
                    inDatTemp(:,badTimes) = []; % Remove bad frames
                    goodDat(iP,:) = downsample(inDatTemp,100); % Downsample imaging data
                end
                toc;        

            else
                % Load the cross correlations for the full FOV
                clear crossCorrFull 
                crossCorrFull                 = load([dataDir '\crossCorrFOV.mat'],'ccFull','ccFullAlpha');
                crossCorrFOV{iRun,iDate}      = crossCorrFull.ccFull;
                crossCorrAlphaFOV{iRun,iDate} = crossCorrFull.ccFullAlpha;
                fovFlag                       = 0;
            end

            % Remove bad time segments determined from spectrogram
            if strcmp(monkeyName,'CharlieSheen')
                if strcmp(expDate,'01_11_2022') && strcmp(runName,'run03') % Remove last 150 s of data
                    probeTemp(end-150e3+1:end,:) = [];
                    processedDat10(:,end-1500+1:end) = [];
                    if fovFlag == 1; goodDat(:,end-1500+1:end) = [];end

                elseif strcmp(expDate,'01_11_2022') && strcmp(runName,'run04') % Remove 380-600 s data
                    probeTemp(380e3+1:600e3,:) = [];
                    processedDat10(:,3800+1:6000) = [];
                    if fovFlag == 1; goodDat(:,3800+1:6000) = [];end
                end

            elseif strcmp(monkeyName,'Whiskey')
                if strcmp(expDate,'08_14_2023') && strcmp(runName,'run04') % Remove 320 - 410 s of data
                    probeTemp(320e3+1:410e3,:) = [];
                    processedDat10(:,3200+1:4100) = [];
                    if fovFlag == 1; goodDat(:,3200+1:4100) = [];end

                elseif strcmp(expDate,'10_16_2023') && strcmp(runName,'run06') % Remove 430 - 510 s of data
                    probeTemp(430e3+1:510e3,:) = [];
                    processedDat10(:,4300+1:5100) = [];
                    if fovFlag == 1; goodDat(:,4300+1:5100) = [];end

                elseif strcmp(expDate,'12_04_2023') && strcmp(runName,'run04') % Remove 510 - 630 s of data
                    probeTemp(510e3+1:630e3,:) = [];
                    processedDat10(:,5100+1:6300) = [];
                    if fovFlag == 1; goodDat(:,5100+1:6300) = [];end

                elseif strcmp(expDate,'12_04_2023') && strcmp(runName,'run05') % Remove 750 - 800 s of data
                    probeTemp(750e3+1:800e3,:) = [];
                    processedDat10(:,7500+1:8000) = [];
                    if fovFlag == 1; goodDat(:,7500+1:8000) = [];end


                elseif strcmp(expDate,'02_20_2024') && strcmp(runName,'run01') % Remove last 100 s of data
                    probeTemp(end-100e3+1:end,:) = [];
                    processedDat10(:,end-1000+1:end) = [];
                    if fovFlag == 1; goodDat(:,end-1000+1:end)= [];end

                elseif strcmp(expDate,'02_20_2024') && strcmp(runName,'run03') % Remove 1-100; 580-650 s of data
                    probeTemp([1:100e3 580e3+1:650e3],:) = [];
                    processedDat10(:,[1:1000 5800+1:6500]) = [];
                    if fovFlag == 1; goodDat(:,[1:1000 5800+1:6500])  = [];end

                elseif strcmp(expDate,'02_20_2024') && strcmp(runName,'run05') % Remove 500 - 700 s of data
                    probeTemp(500e3+1:700e3,:) = [];
                    processedDat10(:,5000+1:7000) = [];
                    if fovFlag == 1; goodDat(:,5000+1:7000) = [];end
                end
            end

            % EPHYS: Get infraslow gamma power fluctuations after
            % removing the bad time segments
            disp('Getting infraslow power fluctuations...');
            clear envelopeBL envelopeFiltered probeBL ch
            ch  =  estChInCortex{iDate}(iRun,:);

            % % Average referencing
            %  probeTemp = probeTemp - mean(probeTemp(:,ch(1):ch(2)),2,'omitnan');

            infraEphys = getInfraSlowPowerLFP(probeTemp,bG,aG,ch);
            infraEphys = mean(infraEphys,2,'omitnan');

            % EPHYS: Get infraslow gamma power fluctuations before
            % removing bad time segments...
            clear envelopeBL sz szInfraGammaOrig szInDatOrig
            infraEphysOrig = getInfraSlowPowerLFP(probeCh,bG,aG,ch);
            infraEphysOrig = mean(infraEphysOrig,2,'omitnan');

            szInfraGammaOrig = size(infraEphysOrig,1); szInDatOrig = size(inDatOrig,2);
            if ~(szInfraGammaOrig == szInDatOrig)
                clear minSize
                minSize = min([szInfraGammaOrig szInDatOrig]);
                infraEphysOrig = infraEphysOrig(1:minSize);
                inDatOrig    = inDatOrig(:,1:minSize);
            end

            % EPHYS: Get the infraslow alpha power fluctuations after
            % removing bad time segments...
            clear infraEphysAlpha sz
            infraEphysAlpha = getInfraSlowPowerLFP(probeTemp,bA,aA,ch);
            infraEphysAlpha = mean(infraEphysAlpha,2,'omitnan');

            % Check size of vectors before performing cross correlations
            clear szInfraGamma szInDat szInfraAlpha
            szInfraGamma = size(infraEphys,1); szInDat = size(inDatUp,2); szInfraAlpha = size(infraEphysAlpha,1);
            if ~(szInfraGamma == szInDat) || ~(szInfraAlpha == szInDat) || ~(szInfraAlpha == szInfraGamma)
                clear minSize
                minSize = min([szInfraGamma szInDat szInfraAlpha ]);
                infraEphys = infraEphys(1:minSize);
                infraEphysAlpha = infraEphysAlpha(1:minSize);
                inDatUp    = inDatUp(:,1:minSize);
            end

            % EPHYS:  Scramble the data
            comb        = randperm(size(infraEphys,1));
            infraEphysS = infraEphys(comb);
            inDatS      = inDatUp(:,comb);
            infraAlphaS = infraEphysAlpha(comb);

            % XCORR: Get the cross correlation between LFP and ISOI for the area
            % around the electrode
            if ~exist(fullfile(dataDir,'crossCorrROI.mat'),'file')|| 1
                tic;
                clear cc ccS ccA ccAS ccOrig lags
                disp('Performing cross correlations...');
                parfor iP = 1:size(inDatUp,1)
                    if iP == 1
                        [cc(:,iP),lags(:,iP)]  = xcorr(infraEphys',inDatUp(iP,:),200,'normalized');
                    else
                        [cc(:,iP),~]  = xcorr(infraEphys',inDatUp(iP,:),200,'normalized');
                    end

                    [ccOrig(:,iP),~] =  xcorr(infraEphysOrig',inDatOrig(iP,:),200,'normalized');
                    [ccS(:,iP),~]    = xcorr(infraEphysS',inDatS(iP,:),200,'normalized');
                    [ccA(:,iP),~]    = xcorr(infraEphysAlpha',inDatUp(iP,:),200,'normalized');
                    [ccAS(:,iP),~]   = xcorr(infraAlphaS',inDatS(iP,:),200,'normalized');
                end

                cc  = reshape(cc,[401 inDatSize(1) inDatSize(2)]);  ccOrig = reshape(ccOrig,[401 inDatSize(1) inDatSize(2)]);
                ccS = reshape(ccS,[401 inDatSize(1) inDatSize(2)]); ccA    = reshape(ccA,[401 inDatSize(1) inDatSize(2)]);
                ccAS = reshape(ccAS,[401 inDatSize(1) inDatSize(2)]);

                crossCorr{iRun,iDate}          = cc;
                crossCorrOrig{iRun,iDate}      = ccOrig;
                crossCorrScrambled{iRun,iDate} = ccS;
                crossCorrAlpha{iRun,iDate}     = ccA;
                crossCorrAlphaS{iRun,iDate}    = ccAS;
                allLags{iRun,iDate}            = lags;
                roiMask{iRun,iDate}            = clipMaskROI; 
                toc;
                disp('Saving cross correlations around the electrode...');
                save([dataDir '\crossCorrROI.mat'],'cc','ccOrig','ccS','ccA','ccAS','lags','clipMaskROI');
            end

            % XCORR: Get the cross correlation between LFP and ISOI for the
            % full field of view...
            if ~exist(fullfile(dataDir,'crossCorrFOV.mat'),'file')
                % Upsampling full field of view and removing bad frames...            
                % Cross correlate LFP and ISOI for the entire FOV
                clear ccFull ccFullAlpha
                tic;
                parfor iP = 1:size(goodDat,1)
                    if iP == 1
                        [ccFull(:,iP),lagFull(:,iP)]  = xcorr(infraEphys',goodDat(iP,:),200,'normalized');
                    else
                        [ccFull(:,iP),~]  = xcorr(infraEphys',goodDat(iP,:),200,'normalized');
                        [ccFullAlpha(:,iP),~]  = xcorr(infraEphysAlpha',goodDat(iP,:),200,'normalized');
                    end
                end

                ccFull      = reshape(ccFull,[401 361 438]);
                ccFullAlpha = reshape(ccFullAlpha,[401 361 438]);

                crossCorrFOV{iRun,iDate}      = ccFull;
                crossCorrAlphaFOV{iRun,iDate} = ccFullAlpha;

                disp('Saving cross correlations for the entire FOV...');
                save([dataDir '\crossCorrFOV.mat'],'ccFull','ccFullAlpha','lagFull');
                toc;
            else
                % Load the cross correlations for ROI
                crossCorrAll = load([dataDir '\crossCorrROI.mat'],'cc','lags');
                crossCorr{iRun,iDate}          = crossCorrAll.cc;
                allLags{iRun,iDate}            = crossCorrAll.lags;
            end

            clear x xNew negIdx vidNew xLow minMedcorrInd maxMedcorrIndLow frameNumMedLow 
            x = allLags{iRun,iDate};
            negIdx = x<0 & x>=-150; xNew = x(negIdx);
            lowIdx = x<0 & x>= -80; xLow = x(lowIdx);
%             imSize = size(crossCorr{iRun,iDate},[2 3]);
 
            [~,minMedcorrInd] = min(median(crossCorr{iRun,iDate}(negIdx,:,:),[2,3],'omitnan'));
            lagLow          = xNew(minMedcorrInd);
            frameNumMedLow  = (x == lagLow);

           %% SPATIAL CONTROL 1 : Shift ROI 1mm by 1mm and cross correlate 
           if ~exist([dataDir '\spatialControl.png'],'file')
               clear  pDatTemp
               pDatTemp = processedDat{iDate,iRun}.tempBandPass;
               tic;
               cVals = {'w','k','b','g','m'};
               figure('units','normalized','outerposition',[0 0 1 1]);
               subplot(121); imagesc(greenFig); hold on; colormap gray; axis image off;
               plot(seedLocProbe(1),seedLocProbe(2),'Marker','pentagram','MarkerSize',15,'MarkerFaceColor','r','MarkerEdgeColor','none');

               for iShift = 1:5
                   clear locShift
                   if iShift == 1
                       locShift = round(roiSize{iDate}(iRun)/spatialBin);
                   else
                       locShift = round(roiSize{iDate}(iRun)*2*(iShift-1)/spatialBin);
                   end

                   for iDir = 1:8 % 8 directions from location of probe
                       clear loc roiControlT roiControlSize roiControlT
                       switch iDir
                           case 1
                               loc = [seedLocProbe(1)+ locShift seedLocProbe(2)]; % x+d,y
                           case 2
                               loc = [seedLocProbe(1)+ locShift seedLocProbe(2)+ locShift]; % x+d,y+d
                           case 3
                               loc = [seedLocProbe(1) seedLocProbe(2)+ locShift];% x,y+d
                           case 4
                               loc = [seedLocProbe(1)- locShift seedLocProbe(2)+ locShift];% x-d,y+d
                           case 5
                               loc = [seedLocProbe(1)- locShift seedLocProbe(2)];% x-d,y
                           case 6
                               loc = [seedLocProbe(1)- locShift seedLocProbe(2)- locShift];% x-d,y-d
                           case 7
                               loc = [seedLocProbe(1) seedLocProbe(2)- locShift];% x,y-d
                           case 8
                               loc = [seedLocProbe(1)+ locShift seedLocProbe(2)- locShift];% x+d,y-d
                       end



                       if sum(fliplr(loc)+12>size(greenFig)) || sum(loc-12<= 0)
                           spCorrControl(iRun,iDate,iShift,iDir) = NaN;
                           continue;
                       else
                           plot(loc(1),loc(2),'.','Color',cVals{iShift},'MarkerSize',30);
                           seedSigT = calculateSeedSignal(imresize(greenIm{iDate,iRun},1/spatialBin),clipMaskCortex,loc,12,pDatTemp); % Get Gaussian weighted seed signal
                           corrMapT = reshape(plotCorrMap(seedSigT,pDatTemp,0),[361*438 1]);

                           crossCorrMap = reshape(squeeze(crossCorrFOV{iRun,iDate}(frameNumMedLow,:,:)),[361*438 1]);
                           spCorrControl(iRun,iDate,iShift,iDir) =  corr(crossCorrMap,corrMapT,'rows','complete');
                       end
                   end
               end
               subplot(122);boxplot(squeeze(spCorrControl(iRun,iDate,:,:))',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});
               xlabel('Distance from probe (mm)'); ylabel('Correlation between FC map and peak negative map');
               f = gcf; exportgraphics(f,[dataDir '\spatialControl.png'],'Resolution',300); close gcf;
           end

            % TEMPORAL CONTROL
            if ~exist([dataDir '\temporalControl.png'],'file')
                timeLen = length(infraEphys);
                winLen  = [1 5 10 50 100 500 1000 5000 timeLen];
                for iSh = 1:length(winLen) % shuffle time series in windows of 10s
                    clear comb1 comb2 infraEphysS roiS  ccT

                    if iSh == 1
                        comb1 = randperm(winLen(iSh));
                        infraEphysS = infraEphys(comb1);
                        comb2 = randperm(winLen(iSh));
                        roiS        = inDatUp(:,comb2);

                    else
                        infraEphysS= [];  roiS = [];
                        comb1 = randperm(round(timeLen/winLen(iSh)));
                        comb2 = randperm(round(timeLen/winLen(iSh)));
                        for iL = 1:length(comb1)
                            clear win1 win2
                            win1 = (comb1(iL)-1)*winLen(iSh)+1 : (comb1(iL)-1)*winLen(iSh)+winLen(iSh);
                            win2 = (comb1(iL)-1)*winLen(iSh)+1 : (comb1(iL)-1)*winLen(iSh)+winLen(iSh);
                            win1(win1>timeLen) = []; win2(win2>timeLen) = [];

                            infraEphysS = [infraEphysS; infraEphys(win1)];
                            roiS        = [roiS inDatUp(:,win2)];
                        end
                    end

                    parfor iP = 1:size(roiS,1)
                        [ccT(:,iP),~]  = xcorr(infraEphysS',roiS(iP,:),200,'normalized');
                    end

                    ccT  = reshape(ccT,[401 inDatSize(1) inDatSize(2)]);
                    ccWinROI(iRun,iDate,iSh) = median(squeeze(ccT(frameNumMedLow,:,:)),'all','omitnan');

                end

                figure; plot(movmean(squeeze(ccWinROI(iRun,iDate,:)),3)); xticklabels(winLen); ylim([min(squeeze(ccWinROI(iRun,iDate,:)))-0.1 0.1]);
                xlabel('Window Length (samples)'); ylabel('Median correlation at peak negative'); grid on;
                f = gcf; exportgraphics(f,[dataDir '\temporalControl.png'],'Resolution',300); close gcf;
            end

        else
            disp('Cross correlations already computed, proceeding with the next steps...');

            % Load the cross correlations for ROI
            crossCorrAll = load([dataDir '\crossCorrROI.mat'],'cc','ccOrig','ccS','ccA',...
                'ccAS','lags','clipMaskROI');
            crossCorr{iRun,iDate}          = crossCorrAll.cc;
            crossCorrOrig{iRun,iDate}      = crossCorrAll.ccOrig;
            crossCorrScrambled{iRun,iDate} = crossCorrAll.ccS;
            crossCorrAlpha{iRun,iDate}     = crossCorrAll.ccA;
            crossCorrAlphaS{iRun,iDate}    = crossCorrAll.ccAS;
            allLags{iRun,iDate}            = crossCorrAll.lags;
            roiMask{iRun,iDate}            = crossCorrAll.clipMaskROI; 

            % Load the cross correlations for the full FOV
            crossCorrFull = load([dataDir '\crossCorrFOV.mat'],'ccFull','ccFullAlpha','lagFull');
            crossCorrFOV{iRun,iDate} = crossCorrFull.ccFull;
            crossCorrAlphaFOV{iRun,iDate} = crossCorrFull.ccFullAlpha;

        end

        %% Spatial (ROI) and temporal lag profiles for gamma band
        [medCorrLagGamma(iRun,iDate), medCorrLagGammaLow(iRun,iDate)] = plotLagProfiles(dataDir,'Gamma',...
            crossCorr{iRun,iDate},allLags{iRun,iDate},monkeyName,expDate,runName,roiMask{iRun,iDate});

        % Spatial (ROI) and temporal profiles for alpha band
        [medCorrLagAlpha(iRun,iDate), medCorrLagAlphaLow(iRun,iDate)] = plotLagProfiles(dataDir,'Alpha',...
            crossCorrAlpha{iRun,iDate},allLags{iRun,iDate},monkeyName,expDate,runName,roiMask{iRun,iDate});

        %% Spatial profiles for the FOV at peak positive and negative...
        % Gamma band
        showPeakXcorrFOV(dataDir,'Gamma',monkeyName,expDate,runName,greenFig,crossCorrFOV{iRun,iDate},...
            clipMaskCortex,medCorrLagGamma(iRun,iDate),medCorrLagGammaLow(iRun,iDate),allLags{iRun,iDate});
%%
        % Alpha band
         showPeakXcorrFOV(dataDir,'Alpha',monkeyName,expDate,runName,greenFig,crossCorrAlphaFOV{iRun,iDate},...
            clipMaskCortex,medCorrLagAlpha(iRun,iDate),medCorrLagAlphaLow(iRun,iDate),allLags{iRun,iDate});

        % Save the spatial profiles for all lags (negative) 
        saveVideoFOV(dataDir,'Gamma',monkeyName,expDate,runName,crossCorrFOV{iRun,iDate},...
            clipMaskCortex,allLags{iRun,iDate});

        
        %         % Infra slow - gamma
        %         % XCORR: Get the spatial profile of cross correlations for all lags
        %         if ~exist(fullfile(dataDir,'\corrVideoMean.mp4'),'file')
        %             disp('Plotting cross correlation profiles...');
        %             v = VideoWriter([dataDir '\corrVideoMean'],'MPEG-4');
        %             open(v);  figure('units','normalized','outerposition',[0 0 1 1]);
        %
        %             for iLag = 1:size(xNew,1)
        %                 subplot(131); imagesc(squeeze(ccOrig(iLag,:,:))); axis image off;
        %                 colormap jet; colorbar; caxis([-0.5 0.5]);
        %                 title(['No removal of bad time segments - Lag: ' num2str(lags(iLag,1)./10) 's']);
        %
        %
        %                 subplot(132); imagesc(squeeze(cc(iLag,:,:))); axis image off; hold on;
        %                 h = imagesc(grayIm); hold off;
        %                 set(h,'AlphaData',~clipMaskROI);
        %                 colormap jet; colorbar; caxis([-0.5 0.5]);
        %                 title(['Concurrent removal of bad time segments - Lag: ' num2str(lags(iLag,1)./10) 's']);
        %
        %                 subplot(133); imagesc(squeeze(ccS(iLag,:,:))); axis image off; hold on;
        %                 h = imagesc(grayIm); hold off;
        %                 set(h,'AlphaData',~clipMaskROI);
        %                 colormap jet; colorbar; caxis([-0.5 0.5]);
        %                 title(['Scrambled infraslow LFP & imaging ' num2str(lags(iLag,1)./10) 's']);
        %                 sgtitle(strrep([monkeyName ' Date: ' expDate ' - Run: ' runName(end)],'_','\_'));
        %
        %                 pause(0.01);
        %                 frame = getframe(gcf);
        %                 writeVideo(v,frame);
        %             end
        %             close(v); close gcf;
        %         end

        % Get moving window cross correlogram to pick a time series
        % Shift time series to get the maximum xcorr
        %         if iRun == 5
        %             padVal = 95;frameNum = 107;
        %         elseif iRun == 6
        %             padVal = 86;frameNum = 592;
        %         elseif iRun == 7
        %             padVal = 69;frameNum = 592;
        %         end
        %         inDatPadded = [zeros(size(inDatUp,1),padVal) inDatUp ];
        %         % Get moving window xcorr and find the time window where max xcorr
        %         % occurs
        %         framePicked = inDatPadded(frameNum,:);
        %         infraEphysPad = [infraEphys; zeros(padVal,1)];
        %         infraAlphaPad = [infraEphysAlpha; zeros(padVal,1)];
        %         winLen = 10; sliding = 1; timeLen = length(framePicked);

        %         % Infra slow  - alpha
        %         % XCORR: Get the spatial profile of cross correlations for all lags
        %         if ~exist(fullfile(dataDir,'\corrVideoMeanCheck_allCh_Alpha.mp4'),'file')
        %             v = VideoWriter([dataDir '\corrVideoMeanCheck_allCh_Alpha'],'MPEG-4');
        %             open(v);  figure;
        %
        %             for iLag = 1:size(xNew,1)
        %                 imagesc(squeeze(ccA(iLag,:,:))); axis image off;
        %                 colormap jet; colorbar; caxis([-0.5 0.5]);
        %                 title(strrep([ 'Lag: ' num2str(lags(iLag,1)./10) 's ' monkeyName ' Date: ' expDate ...
        %                     ' - Run: ' runName(end)],'_','\_'));
        %
        %                 pause(0.01);
        %                 frame = getframe(gcf);
        %                 writeVideo(v,frame);
        %             end
        %             close(v); close gcf;
        %         end

    end
end

%% Compiling multiple runs across all experiments for gamma and alpha bands
for iBand = 1:2 
    clear meanXCorr bandName gMed
    switch iBand
        case 1
            meanXCorr = cellfun(@(x)median(x,[2 3],'omitnan'),crossCorr,'un',0);

            bandName  = 'Gamma band';
        case 2
            meanXCorr = cellfun(@(x)median(x,[2 3],'omitnan'),crossCorrAlpha,'un',0);
            bandName  = 'Alpha band';
    end
    clear x xNew;
    meanXCorr = horzcat(meanXCorr{:});
    x = allLags{1,1};
    negIdx = x<-20 & x>=-200; xNew = x(negIdx);

%     [~,maxAllInd] = max(meanXCorr(negIdx,:));
%     maxCorrTime    = xNew(maxAllInd)./10;

    figure;plot(x,meanXCorr,'Color',[0.65 0.65 0.65]); hold on;
    plot(x,movmean(median(meanXCorr,2,'omitnan'),10),'k','LineWidth',2);
    grid off; xticklabels(-20:5:20);ylim([-0.4 0.4]);box off;xlim([-150 0]);
    title(['Temporal Profile for all recordings for ' bandName]); 
    xlabel('Lag (s)'); ylabel('Median cross correlation');

    gMed = median(meanXCorr,2,'omitnan');
%     [~,maxIdx] = (max(gMed));

    % Plot sem as patches around median
    sem  = std(meanXCorr,0,2)./size(meanXCorr,2);
    figure; plot(x,movmean(gMed,1),'k','LineWidth',2); hold on;
    patch([(x) ;flipud((x))],[(gMed-sem);  flipud((gMed+sem))],'b','FaceAlpha',0.5,'EdgeColor','none');
    xticklabels(-20:5:20);xlim([-150 0]);ylim([-0.4 0.4]);
    title(['Median temporal Profile across all recordings for ' bandName]);
    xlabel('Lag (s)'); ylabel('Median cross correlation');
end

%% Comparison of spatial profile with functional connectivity map 
% 1. Correlate spatial map with FC map
% 2. Weighted correlation between spatial map and FC map 
% 3. Combine information from microelectrode map and compare the spatial
% organization. - can be done with Charlie 

iDate = 2; iRun = 6;
expDate    = allDates(iDate,:);
runName    = allRuns{iDate,1}(iRun,:);
dataDir    = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

%  Load the appropriate masks for the imaging data
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

% Seed to get the FC map for the run 
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(greenIm{iDate,iRun}); colormap gray; axis image off;
title(strrep(['Pick a seed to get a FC map for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)],'_','\_')); hold on;
seedLoc = ginput(1);
seedLoc = fliplr(round(seedLoc));
plot(seedLoc(2),seedLoc(1),'w.','MarkerSize',15);pause(5); close gcf;

pDatTemp = processedDat{iDate,iRun}.tempBandPass;
seedSig  = calculateSeedSignal(imresize(greenIm{iDate,iRun},1/spatialBin), clipMask,fliplr(round(seedLoc./spatialBin)),12,pDatTemp); % Get Gaussian weighted seed signal
corrMap  = plotCorrMap(seedSig,pDatTemp,0);

corrMap(~clipMaskCortex) = NaN;
imSize                   = size(corrMap);
grayIm                   = cat(3, 0.25.*ones((imSize)),0.25.*ones((imSize)), 0.25.*ones((imSize)));

% Show the FC map
figure; imagesc(corrMap); hold on; 
h = imagesc(grayIm); hold off;
set(h,'AlphaData',~clipMaskCortex); axis image off; 
colorbar; colormap jet; caxis([0 1]);

% Show the map at peak negative cross correlation 
l = allLags{iRun,iDate};
frameIdx = find(l == -21);
peakNegFrame = squeeze(crossCorrFOV{iRun,iDate}(frameIdx,:,:));
peakNegFrame(~clipMaskCortex) = NaN;
figure; imagesc(peakNegFrame); hold on;
h = imagesc(grayIm); hold off;
set(h,'AlphaData',~clipMaskCortex); axis image off;
colorbar; colormap(flipud(jet)); caxis([-0.5 0.4]);

%% 1. Correlating FC map with spatial profile at peak negative
corrMapR    = reshape(corrMap,[imSize(1)*imSize(2) 1]);
peakNegR    = reshape(peakNegFrame,[imSize(1)*imSize(2) 1]); 
spatialCorr = corr(corrMapR, peakNegR,'rows','complete'); 

%% 2. Spatial correlations after masking out FC values<0.2 and also calculating mean separation
weightMat = zeros(size(corrMapR)); 
weightMat(corrMapR>=0.2) = 1; 

corrMapR(~weightMat) = NaN; 
peakNegR(~weightMat) = NaN; 
mCorr  = corr(corrMapR, peakNegR,'rows','complete');

%% 3. Controls to determine how accurate are the calculated values
% Control within a session - get different FC map 
pDatTemp    = processedDat{iDate,iRun}.tempBandPass;
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(greenIm{iDate,iRun}); colormap gray; axis image off;
title(strrep(['Pick a seed to get a FC map for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)],'_','\_')); hold on;
seedLoc1 = ginput(1);
seedLoc1 = fliplr(round(seedLoc1)); close gcf;
seedSig     = calculateSeedSignal(imresize(greenIm{iDate,iRun},1/spatialBin), clipMask,fliplr(round(seedLoc1./spatialBin)),12,pDatTemp); % Get Gaussian weighted seed signal
corrMapNew  = plotCorrMap(seedSig,pDatTemp,0);
corrMapNew(~clipMaskCortex) = NaN;

peakNegRNew    = reshape(peakNegFrame,[imSize(1)*imSize(2) 1]); 

corrMapNewR = reshape(corrMapNew,[imSize(1)*imSize(2) 1]); 
spatialCorrNew = corr(corrMapNewR,peakNegRNew,'rows','complete'); 
% eDistNew = sqrt(mean((corrMapNewR - peakNegRNew).^2,'omitnan')); 

weightMatNew = zeros(size(corrMapNewR)); 
weightMatNew(corrMapNewR>=0.2) = 1; 
corrMapNewR(~weightMatNew) = NaN;
peakNegRNew(~weightMatNew) = NaN;
spatialCorrNewM  = corr(corrMapNewR,peakNegRNew,'rows','complete'); 

%% Mismatch across sessions
clear pDatTemp 
pDatTemp    = processedDat{iDate,iRun-5}.tempBandPass;

figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(greenIm{iDate,iRun}); colormap gray; axis image off;
title(strrep(['Pick a seed to get a FC map for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)],'_','\_')); hold on;
seedLoc1 = ginput(1);
seedLoc1 = fliplr(round(seedLoc1)); close gcf;

seedSig     = calculateSeedSignal(imresize(greenIm{iDate,iRun-5},1/spatialBin), clipMask,fliplr(round(seedLoc./spatialBin)),12,pDatTemp); % Get Gaussian weighted seed signal
corrMapNewT  = plotCorrMap(seedSig,pDatTemp,0);
corrMapNewT(~clipMaskCortex) = NaN; 

peakNegRNew    = reshape(peakNegFrame,[imSize(1)*imSize(2) 1]); 

corrMapNewR = reshape(corrMapNewT,[imSize(1)*imSize(2) 1]); 
spatialCorrNewS = corr(corrMapNewR,peakNegRNew,'rows','complete'); 
% eDistNew = sqrt(mean((corrMapNewR - peakNegRNew).^2,'omitnan')); 

corrMapNewR(~weightMat) = NaN;
peakNegRNew(~weightMat) = NaN;
% eDistNewM = sqrt(mean(corrMapNewR-peakNegRNew).^2,'omitnan');






%% Spatial correspondence based on microelectrode mapping - Applicable to Charlie... 
zones = {'Hand','Arm','Face','Trunk','Leg'};
areas = {'Motor','Somatosensory'}; 

% Cosine similarity - refer to Nick's Prog Neuro paper Fig 4D. 
% Assign zero values to NaN
fcMapNoNaN   = corrMapR; fcMapNoNaN(isnan(fcMapNoNaN))     = 0; 
peakNegNoNaN = peakNegR; peakNegNoNaN(isnan(peakNegNoNaN)) = 0; 

d = 1-pdist2(fcMapNoNaN',peakNegNoNaN','cosine');
dSelf = 1-pdist2(fcMapNoNaN',fcMapNoNaN','cosine');

% Mean distance between pixels 
corrMapT = corrMap;      corrMapT(isnan(corrMap)) = 0;
peakNegT = peakNegFrame; peakNegT(isnan(peakNegFrame)) = 0; 

mD = mean(pdist2(corrMapT,peakNegT),'all','omitnan');
mDSelf = mean(pdist2(corrMapT,corrMapT),'all','omitnan');



%% Plot FC map 
greenImRGB = ind2rgb(greenIm{iDate,iRun},gray(256));
frameSize = size(greenImRGB);


figure('units','normalized','outerposition',[0 0 1 1]); imagesc(greenImRGB); hold on;
corrMapR = imresize(corrMap,spatialBin,'OutputSize',[frameSize(1) frameSize(2)]);
imagesc(corrMapR,'alphadata',corrMapR.*0.8);
colormap jet; colorbar; caxis([0 1]); axis image off;
plot(seedLoc(2),seedLoc(1),'w.','MarkerSize',35);

%% Functions that are used in this script

function [lagHigh, lagLow] = plotLagProfiles(dataDir,bandName,crossCorr,allLags,monkeyName,expDate,runName,clipMaskROI)
% Functions to plot the following -
% 1.  XCORR: Median xcorr vs lag
% 2.  XCORR: Find the frame with peak positive and peak negative correlations

clear x xNew negIdx vidNew xLow;
x = allLags;
negIdx = x<0 & x>=-150; xNew = x(negIdx);
lowIdx = x<0 & x>= -80; xLow = x(lowIdx);
imSize = size(crossCorr,[2 3]);
grayIm = cat(3, 0.25.*ones((imSize)),0.25.*ones((imSize)), 0.25.*ones((imSize)));

% vidIdx = x<=0 & x>=-150; vidNew = x(vidIdx);

% 1. XCORR: Median xcorr vs lag
if ~exist([dataDir '\xcorrVsLag' bandName '.png'],'file') || 1
    figure('Position',[400 400 475 600]);
    plot(x,movmean(median(crossCorr,[2,3],'omitnan'),3));
    xlabel('Lag (s)'); ylabel('cross correlation');   ylim([-0.4 0.4]);
    legend('Median xcorr','Location','northeast','Autoupdate','off');
    xline(0); xticks(-200:40:200); xticklabels(-20:4:20); grid on; xlim([-160 0]);
    title(strrep([ bandName ' xcorr vs lag for ' monkeyName ' Date: ' expDate ' - Run: ' runName(end)],'_','\_'));
    f = gcf; exportgraphics(f,[dataDir,'\xcorrVsLag' bandName '.png'],'Resolution',300); close gcf;
end

% 2. XCORR: Find the frame with peak positive and peak negative correlations
clear maxMedcorrInd maxMedcorrIndLow frameNumMed frameNumMedLow
[~,maxMedcorrInd] = max(median(crossCorr(negIdx,:,:),[2,3],'omitnan'));
lagHigh           = xNew(maxMedcorrInd);
frameNumMed       = (x == lagHigh);

[~,maxMedcorrIndLow] = min(median(crossCorr(lowIdx,:,:),[2,3],'omitnan'));
lagLow               = xLow(maxMedcorrIndLow);
frameNumMedLow       = (x == lagLow);

% Save the frame where the highest correlation occurred
if ~exist([dataDir '\highestXCorrFrames_' bandName '.png'],'file')
    figure('units','normalized','outerposition',[0 0 1 1]); tiledlayout(1,2);
    clear ax1; ax1 = nexttile;
    imagesc(squeeze(crossCorr(frameNumMed,:,:))); axis image off; hold on;
    h = imagesc(grayIm); hold off;
    set(h,'AlphaData',~clipMaskROI);
    colormap jet; colorbar; caxis([-0.5 0.5]);
    title(['Peak Positive median correlation at Lag: ' num2str(xNew(maxMedcorrInd)./10) 's ']);

    nexttile; imagesc(squeeze(crossCorr(frameNumMedLow,:,:))); axis image off; hold on;
    h = imagesc(grayIm); hold off; set(h,'AlphaData',~clipMaskROI);
    colormap(flipud(jet)); colorbar; caxis([-0.5 0.5]);
    title(['Peak Negative median correlation at Lag: ' num2str(xLow(maxMedcorrIndLow)./10) 's ']);
    colormap(ax1,'jet');

    sgtitle(strrep([ bandName ' xcorr ' monkeyName ' Date: ' expDate ...
        ' - Run: ' runName(end)],'_','\_'));

    f = gcf; exportgraphics(f,[dataDir,'\highestXCorrFrames_' bandName '.png'],'Resolution',300); close gcf;
end

end

%%
function showPeakXcorrFOV(dataDir,bandName,monkeyName,expDate,runName,greenFig,crossCorrFOV,clipMaskCortex,frameNumHighTime,frameNumLowTime,allLags)
% This function plots the following -
% 1. Correlation map at peak positive and negative
% 2. Significant pixels of correlation map on blood vessel map

imSize = size(crossCorrFOV,[2,3]);
grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));

frameNumHigh = find(allLags == frameNumHighTime);
frameNumLow  = find(allLags == frameNumLowTime);

% 1. Correlation map at peak positive and negative
if ~exist([dataDir '\XCorrFOV_' bandName '.png'],'file') || 1
    figure('units','normalized','outerposition',[0 0 1 1]); tiledlayout(1,2);
    clear ax1; ax1 = nexttile;
    imagesc(squeeze(crossCorrFOV(frameNumHigh,:,:))); axis image off; hold on;
    h = imagesc(grayImFull); hold off;
    set(h,'AlphaData',~clipMaskCortex);
    colormap(ax1,'jet'); colorbar; caxis([-0.5 0.5]);
    title(strrep([' Peak Positive: ' num2str(frameNumHigh/10)],'_','\_'));

    ax2 = nexttile; imagesc(squeeze(crossCorrFOV(frameNumLow,:,:))); axis image off; hold on;
    h = imagesc(grayImFull); hold off;
    set(h,'AlphaData',~clipMaskCortex);
    colormap(ax2,flipud(jet)); colorbar; caxis([-0.5 0.5]);
    title(strrep([' Peak Negative: ' num2str(allLags(frameNumLow)/10)],'_','\_'));

    colormap(ax1,'jet');
    sgtitle(strrep([ bandName ' xcorr ' monkeyName ' Date: ' expDate ...
        ' - Run: ' runName(end)],'_','\_'));

    f = gcf; exportgraphics(f,[dataDir,'\XCorrFOV_' bandName '.png'],'Resolution',300); close gcf;
end

% 2. Significant pixels of correlation map on blood vessel map
if ~exist([dataDir '\greenCorrFOV_' bandName '.png'],'file')
    clear figFrameHigh figFrameLow
    figFrameHigh = squeeze(crossCorrFOV(frameNumHigh,:,:));
    figFrameLow  = squeeze(crossCorrFOV(frameNumLow,:,:));

    figure('units','normalized','outerposition',[0 0 1 1]); tiledlayout(1,2);
    clear ax1; ax1 = nexttile;
    imagesc(ind2rgb(greenFig,gray(256))); hold on; axis image off;
    imagesc(figFrameHigh,'alphaData',figFrameHigh.*5);
    colormap(ax1,'jet'); caxis([-0.5 0.5]); colorbar;
    title(strrep([' Peak Positive: ' num2str(allLags(frameNumLow)/10)],'_','\_'));

    ax2 = nexttile;  imagesc(ind2rgb(greenFig,gray(256))); hold on; axis image off;
    imagesc(figFrameLow,'alphaData',figFrameLow.*-5);
    colormap(ax2,flipud(jet)); colorbar; caxis([-0.5 0.5]);
    title(strrep([' Peak Negative: ' num2str(frameNumLow/10)],'_','\_'));
    colormap(ax1,'jet');

    sgtitle(strrep([ bandName ' xcorr ' monkeyName ' Date: ' expDate ...
        ' - Run: ' runName(end)],'_','\_'));

    f = gcf; exportgraphics(f,[dataDir,'\greenCorrFOV_' bandName '.png'],'Resolution',300); close gcf;
end

end

%% 
function saveVideoFOV(dataDir,bandName,monkeyName,expDate,runName,crossCorrFOV,clipMaskCortex,allLags)

if ~exist(fullfile([dataDir,'\XCorrFOV_' bandName 'Video.avi']),'file')
    imSize = size(crossCorrFOV,[2,3]);
    grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));

    v = VideoWriter([dataDir,'\XCorrFOV_' bandName 'Video']);
    open(v);

    lagInd = find(allLags>=-150 & allLags<=0);
    lagVals = allLags(allLags>=-150 & allLags<=0); 

    for iLag = 1:length(lagVals)
        imagesc(squeeze(crossCorrFOV(lagInd(iLag),:,:)));
        hold on; h = imagesc(grayImFull); hold off; set(h,'AlphaData',~clipMaskCortex);
        axis image off; colormap(flipud(jet)); caxis([-0.5 0.5]); colorbar;

        title(strrep([ bandName ' xcorr ' monkeyName ' Date: ' expDate ...
            ' - Run: ' runName(end) ' Lag at: ' num2str(lagVals(iLag)./10) 's'],'_','\_'));

        pause(0.2); frame = getframe(gcf); writeVideo(v,frame);
    end
    close(v); close gcf;
end
end
%% 
function wCorr = weightedCorr(weightMat,xVar, yVar)
% This function performs weighted correlations
% Refer: https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#Weighted_correlation_coefficient

% Weighted means
xWMean = sum(weightMat.*xVar,'omitnan')/sum(weightMat); 
yWMean = sum(weightMat.*yVar,'omitnan')/sum(weightMat);

% Weighted covariances
covX  = sum(weightMat.*(xVar - xWMean).^2,'omitnan')/sum(weightMat);
covY  = sum(weightMat.*(yVar - yWMean).^2,'omitnan')/sum(weightMat); 
covXY = sum(weightMat.*(xVar - xWMean).*(yVar - yWMean),'omitnan')./sum(weightMat); 

% Weighted correlation
wCorr = covXY/(sqrt(covX*covY));  

end
