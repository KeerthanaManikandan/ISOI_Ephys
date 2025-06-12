function [processedDat,greenIm,probe,badCh,badTimesLFP,badTimeThresh,estChInCortex] = getAllData_Imaging_Ephys(monkeyName,hemisphere,commonDir,...
    serverPath,allDates,allRuns,ephysFileNameAll,greenMapRef,chInCortexNotes,probeLabel,spatialBin)

% This function stores/retrieves imaging and physiology data for
% simultaneous imaging and physiology experiments. This function also
% identifies bad time segments and channels and transition channels for
% physiology.
% April 03, 2024 - KM
% Set parameters for spectrogram
fs              = 1000;
params.Fs       = fs;
params.fpass    = [1 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;
gammaBand   = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); % Gamma band filtering parameters

% Store/Retrieve imaging data for all experiments and all sessions...
for iDate = 1:size(allDates,1)
    clear expDate
    expDate    = allDates(iDate,:);

    for iRun = 1:size(allRuns{iDate,1},1)
        % Load all the information
        runName    = allRuns{iDate,1}(iRun,:);
        datFileNum = ephysFileNameAll{iDate,1};
        fileNum    = str2double(datFileNum(iRun,end));

        % Get the directory and filenames
        dataDir     = ['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
        if ~exist(dataDir,'dir'); [~,~] = mkdir(dataDir); end

        numFiles    = length(dir([dataDir '\Spatial Downsample SS3' '/*.mat'])); % loading the downsampled data only
        datName     = 'Data_RS_10Hz_SS3_';
        templateDir = ['D:\Data\' monkeyName '_SqM\Left Hemisphere\'];
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
            [~,~,~,~,tempBandPass] = getPreProcessedDataRestingState(serverDataPath,dataDir,runName,numFiles,spatialBin,datName);
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
        
        % Get FC map near the electrode 
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
            colormap jet; colorbar; clim([0 1]); axis image off;
            plot(seedLoc(2),seedLoc(1),'w.','MarkerSize',35);

            sgtitle(strrep([monkeyName ' Date: ' expDate ' - Run: ' runName(end)],'_','\_'));
            f = gcf; exportgraphics(f,[dataDir '\FCMap.png'],'Resolution',300); close gcf;
        end
        
        % Check for quality of imaging by generating FC map from a random
        % seed with known connectivity patterns.
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
            colormap jet; colorbar; clim([0 1]); axis image off;
            plot(seedLoc(2),seedLoc(1),'w.','MarkerSize',35);

            sgtitle(strrep([monkeyName ' Date: ' expDate ' - Run: ' runName(end)],'_','\_'));
            f = gcf; exportgraphics(f,[dataDir '\FCMapQC.png'],'Resolution',300); close gcf;
        end
        clear pDatTemp

        % Get average FC map from multiple RS runs 
        if ~exist([dataDir '\AvgRS_FCMap.png'],'file')
            tic;
            figure('units','normalized','outerposition',[0 0 1 1]);
            imagesc(greenMapRef);axis image off; colormap gray; hold on;
            title(strrep(['Select the location where the probe is located for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)],'_','\_'));
            refSeed(1,:,:) = ginput(1); close gcf;
            [~,corrMapFinal] = getRSConnectivityMaps(squeeze(refSeed(1,:,:))',monkeyName); 

            figure('units','normalized','outerposition',[0 0 1 1]);
            imagesc(ind2rgb(greenMapRef,gray(256))); hold on;
            imagesc(corrMapFinal,'alphaData',corrMapFinal.*0.9); axis image off; hold on;
            plot(refSeed(1,:,1),refSeed(1,:,2),'w.','MarkerSize',35);
            colormap jet; clim([0 1]); f = gcf; toc;
            sgtitle(strrep([monkeyName ' Date: ' expDate ' - Run: ' runName(end)],'_','\_'));
            exportgraphics(f,[dataDir '\AvgRS_FCMap.png'],'Resolution',300); close gcf;
             
        end
    end
end

% Store/Retrieve electrophysiology data...
for iDate = 1:size(allDates,1)
    clear expDate
    expDate    = allDates(iDate,:);
    for iRun = 1:size(allRuns{iDate,1},1)
        clear entityInfo datFileName timeStamp

        % Load all run information
        runName    = allRuns{iDate,1}(iRun,:);
        saveFolder = ['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Processed Data'];
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
        if ~exist([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'],'file')
             vInfo = [];          
        else
             vInfo = who('-file',[saveFolder '\' datFileName num2str(fileNum) '_lfp.mat']);    
        end

        % Preprocess and store data...
        if  ~exist([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'],'file') || any(~ismember('timeStamp',vInfo))
            disp(['Preprocessing LFP data for : ' runName]);

            if exist([serverPath expDate '\Electrophysiology\' runName '\' datFileName num2str(fileNum) '.nev'],'file')
                datName = [serverPath expDate '\Electrophysiology\' runName '\' datFileName num2str(fileNum)];

            elseif exist([serverPath expDate '\Electrophysiology\' datFileName num2str(fileNum) '.nev'],'file')
                datName = [serverPath expDate '\Electrophysiology\' datFileName num2str(fileNum)];
            end

            [probeCh,rawCh,eegCh,timeStamp] = saveLFPSingleProbe(runName,datName,fs);

            % Save the LFP
            if ~exist('saveFolder','dir'); [~,~] = mkdir(saveFolder); end
            disp('Storing data... ' );

            if exist('eegCh','var')
                save([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'] ,'probeCh','rawCh','eegCh','timeStamp');
            else
                save([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'] ,'probeCh','rawCh','timeStamp');
            end

            probe{iRun,iDate} = matfile([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat']); clear probeCh eegCh timeStamp rawCh;        
     
        else
            % Retrieve LFP and spiking data
            disp(['Retrieving electrophysiology data for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)]);
            probe{iRun,iDate} = matfile([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat']);
        end
    end
end

% Determining bad time segments, channels and transition channels for LFP time series...
for iDate = 1:size(allDates,1)
    clear expDate
    expDate    = allDates(iDate,:);

    for iRun =  1:size(allRuns{iDate,1},1)
        clear  runName dataDir probeCh transitionCh probeTemp
        runName    = allRuns{iDate,1}(iRun,:);
        dataDir    = ['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

        % Determine the bad time segments and bad channels in LFP data
        disp(['Identifying noisy data from LFP for:  ' allDates(iDate,:) ' ' runName]);
        probeCh      = probe{iRun,iDate}.probeCh;
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
            elseif strcmp(probeLabel{iDate,1}(iRun,:),'0763')
                channels(ismember(channels,[2 4])) = [];
                badChVal = [badChVal 2 4];
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
        badTimeThresh{iDate,iRun} = badTimeIndOld; % Threshold crossings
        
        % Determine the bad time segments
        badTimeInd = []; badTimes = [];
        
        % Remove 50 ms before and after each threshold crossing
        if ~isempty(badTimeIndOld)
            badTimeInd = [(badTimeIndOld-fs/20)  (badTimeIndOld+fs/20)]; 

            for iL = 1:size(badTimeInd,1)
                badTimes = [badTimes badTimeInd(iL,1): badTimeInd(iL,2)];
            end
        end
        
        badTimes = unique(badTimes);
        badTimes(badTimes>size(probeCh,1))=[];
        badTimesLFP{iDate,iRun} = badTimes;

        % *Checkpoint*: Get the spectrograms before and after bad time
        % segment removal to check if noise has been removed.
        if ~exist([dataDir '\AverageLFPSpec.png'],'file')% Plotting and saving the spectrogram
            clear spec timeValsSpec freqValsSpec

            if length(channels)~= 1 % Linear array
                [spec,timeValsSpec,freqValsSpec] = mtspecgramc(probeCh(:,15:end),[5 2],params);
                figure; subplot(211);
                imagesc(timeValsSpec,freqValsSpec, 10.*log10(squeeze(mean(spec,3)'))); % Before removal of bad time segments
                clim([-20 20]); set(gca,'YDir','normal'); colormap jet; colorbar;
                xlabel('Time (s)'); ylabel('Power (dB)'); title('Before removing bad time segments');

                if ~isempty(badTimes)
                    probeChTemp = probeCh;
                    probeChTemp(badTimes,:) = []; % Removing bad Time segments...

                    clear spec timeValsSpec freqValsSpec
                    [spec,timeValsSpec,freqValsSpec] = mtspecgramc(probeChTemp(:,15:end),[5 2],params);
                    subplot(212); imagesc(timeValsSpec,freqValsSpec, 10.*log10(squeeze(mean(spec,3)'))); % After removal of bad time segments
                    clim([-20 20]); set(gca,'YDir','normal');colormap jet; colorbar;
                    xlabel('Time (s)'); ylabel('Power (dB)'); title('After removing bad time segments');
                else
                    subplot(212); imagesc(timeValsSpec,freqValsSpec, 10.*log10(squeeze(mean(spec(:,:,15:end),3)')));
                    clim([-20 20]); set(gca,'YDir','normal'); title('No bad time segments detected');colormap jet; colorbar;
                end

            else % Single probe
                probeChTemp = probeCh;
                [spec,timeValsSpec,freqValsSpec] = mtspecgramc(probeChTemp,[5 2],params);
                figure; subplot(211); imagesc(timeValsSpec,freqValsSpec, 10.*log10(spec'));% After removal of bad time segments
                clim([-20 20]); set(gca,'YDir','normal'); xlabel('Time (s)'); ylabel('Power (dB)'); colormap jet; colorbar;
                if strcmp(expDate,'08_14_2023'); clim([-60 60]); end

                if ~isempty(badTimes)
                    probeChTemp(badTimes,:) = []; % Removing bad Time segments...

                    clear spec timeValsSpec freqValsSpec
                    [spec,timeValsSpec,freqValsSpec] = mtspecgramc(probeChTemp,[5 2],params);
                    subplot(212);imagesc(timeValsSpec,freqValsSpec, 10.*log10(spec'));
                    clim([-20 20]); set(gca,'YDir','normal'); colormap jet;
                    xlabel('Time (s)'); ylabel('Power (dB)'); title('After removing bad time segments'); colorbar;
                    if strcmp(expDate,'08_14_2023'); clim([-60 60]); end
                else
                    subplot(212); imagesc(timeValsSpec,freqValsSpec, 10.*log10(spec'));
                    clim([-20 20]); set(gca,'YDir','normal'); title('No bad time segments detected'); colormap jet; colorbar;
                    if strcmp(expDate,'08_14_2023'); clim([-60 60]); end
                end
            end
            sgtitle(strrep([monkeyName ' Date: ' expDate ' - Run: ' runName(end)],'_','\_'));
            f = gcf; exportgraphics(f,[dataDir '\AverageLFPSpec.png'],'Resolution',300); close gcf;
        end

        probeTemp = probeCh;
        probeTemp(badTimes,:) = [];

        % Determine the channels that are within cortex
        if size(probeTemp,2)~=1

            % Get mean within probe marginals from wideband range for the
            % electrode
            marginalVal = mean(imgaussfilt(corr(probeTemp,'Rows','complete'),1),2,'omitnan');

            if ~exist(fullfile(dataDir,'pairCorr.png'),'file') % Within probe correlations 
                tiledlayout(1,3,'TileSpacing','Compact');
                nexttile([1 2]); imagesc(imgaussfilt(corr(probeTemp,'Rows','complete'),1));
                colormap jet; colorbar; axis image tight;
                xticks(1:2:size(probeTemp,2));  yticks(1:2:size(probeTemp,2));
                
                nexttile; plot(marginalVal,1:length(marginalVal));set(gca,'YDir','reverse');
                axis padded;xlim([0 round((max(marginalVal)+0.1).*10)/10]);
                xticks(0:0.2:round((max(marginalVal)+0.1).*10)/10);

                sgtitle(strrep(['Within probe correlations for: ' monkeyName ' Date: ' expDate ' - Run: ' runName(end)],'_','\_'));
                f = gcf; exportgraphics(f,[dataDir '\pairCorr.png'],'Resolution',300); close gcf;
            end

            % Get the slope of the marginals for the electrode 
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

        % Animal state: Plot spectrogram of the EEG
        clear eegGood spec timeValsSpec freqValsSpec
        if ~exist(fullfile(dataDir,'eegSpec.png'),'file')
            eegGood             = probe{iRun,iDate}.eegCh;
            eegGood(badTimes,:) = [];

            [spec,timeValsSpec,freqValsSpec] = mtspecgramc(eegGood,[5 2],params); % Spectrogram

            figure; imagesc(timeValsSpec,freqValsSpec, 10.*log10((spec)'));
            clim([-20 40]); set(gca,'YDir','normal');colormap jet; colorbar;
            xlabel('Time (s)'); ylabel('Frequency (Hz)');
            title(strrep([monkeyName ' Date: ' expDate ' Run: ' runName(end) ' - Spectrogram of EEG'],'_','\_'));
            f = gcf; exportgraphics(f,[dataDir '\eegSpec.png'],'Resolution',300); close gcf;
        end
    end
end

end