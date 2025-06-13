% ephysImaging_v3.m
% This script analyzes electrophysiological and imaging data recorded
% simultaneously for one monkey
% April 3, 2024 - KM
% See ephysImaging.m and ephysImaging.m_v2 for previous versions
% Set paths
clc; clear -except spCorrAll spCorr
commonDir = 'C:\Users\KEM294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\Imaging']));
addpath(genpath([commonDir '\Codes\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\fly_track\videoIO']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\spectral_analysis\continuous\dupes']));

%% Storing/Retrieving monkey data
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

% Get monkey parameters
[allDates,allRuns,refDate,refDir,lensCombo,roiSize,ephysFileNameAll,serverPath,probeLabel,chInCortexNotes,greenMapRef] ...
    = getMonkeyParams_Imaging_Ephys(monkeyName, commonDir, hemisphere);

%% Store/Retrieve imaging and physiology data
[processedDat,greenIm,probe,badCh,badTimesLFP,estChInCortex] = getAllData_Imaging_Ephys(monkeyName,hemisphere,commonDir,serverPath,allDates,allRuns,...
    ephysFileNameAll,greenMapRef,chInCortexNotes,probeLabel,spatialBin);   % Change the logical check at line 67 of this function
clc; disp(['Stored/retrieved imaging and physiology data for ' monkeyName]);

%% Computing cross correlations

for iDate = 3%1:size(allDates,1)
    clear expDate
    expDate    = allDates(iDate,:);

    for iRun =1%:7% 1:size(allRuns{iDate,1},1)
        clc; clear  runName dataDir probeCh badTimes badChVal crossCorrAll processedDat10...
            processedDatR probeTemp clipMaskROI grayIm clearMaskFlag varInfo probeLocFlag

        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
        disp([ monkeyName ' ' allDates(iDate,:) ' ' runName]);

        varInfo = who('-file', fullfile(dataDir,'crossCorrROI.mat'));
        if find(ismember(varInfo,'clipMaskROI')); clipMaskFlag = 1; else; clipMaskFlag = 0; end

        clear varInfo; varInfo = who('-file', fullfile(dataDir,'roiCenterLoc.mat'));
        if find(ismember(varInfo,'seedLocProbe')); probeLocFlag = 1; else; probeLocFlag = 0; end
        %          if ~exist([dataDir '\temporalControl.png'],'file') end
        %           if ~exist([dataDir '\spatialControl.png'],'file')   end

        %%% Add flags to check if you got spatial/temporal controls or not

        imSize = size(imresize(greenIm{iDate,iRun},1/spatialBin));

        %         if ~exist(fullfile(dataDir,'crossCorrROI.mat'),'file') && ~exist(fullfile(dataDir,'crossCorrFOV.mat'),'file') && ~clipMaskFlag && ~probeLocFlag % Check if cross correlations have been computed
        disp('Loading ISOI and LFP data for this run...');
        % EPHYS: Retrieve LFP for the run
        probeCh  = probe{iRun,iDate}.probeCh;
        channels = 1:size(probeCh,2);
        badChVal = badCh{iDate,iRun};

        % EPHYS: Remove bad channels from LFP
        badTimes  = badTimesLFP{iDate,iRun};
        probeCh(:,badChVal) = [];

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

        if exist([dataDir '\corrMask0' runName(end) '.BMP'],'file')% Mask for comparing the imaging and hybrid map
            corrMask = imread([dataDir '\corrMask0' runName(end) '.bmp']); corrMaskFlag= 1;
        elseif exist([dataDir '\corrMask0' runName(end) '.PNG'],'file') 
            corrMask = imread([dataDir '\corrMask0' runName(end) '.png']);  corrMaskFlag= 1;
        else
            corrMaskFlag = 0;
        end

        if corrMaskFlag
            corrMask = imresize(corrMask,1/3);
            corrMask = corrMask(:,:,1)>0;
            corrMask = corrMask & ~elecMask;
        else
            corrMask = clipMaskCortex;
        end

        % IMAGING: Get imaging data
        processedDatR  = reshape(processedDat{iDate,iRun}.tempBandPass,[imSize(1)*imSize(2) size(processedDat{iDate,iRun},'tempBandPass',3)]);

        % IMAGING: Upsample the imaging dataset
        disp('Upsampling imaging data to 10 Hz...');
        parfor iP = 1:size(processedDatR,1)
            processedDat10(iP,:) = interp(processedDatR(iP,:),5);
        end

        % IMAGING: Remove non-cortical components from imaging data
        clipMaskCortexR = reshape(corrMask,[size(clipMaskCortex,1)*size(clipMaskCortex,2) 1]);
        processedDat10(~clipMaskCortexR,:) = NaN;
        processedDat10R = reshape(processedDat10,[imSize(1) imSize(2) size(processedDat10,2)]);

        % IMAGING: Pick the area around the electrode
        seedRad = roiSize{iDate}(iRun);%floor(roiSize{iDate}(iRun)/spatialBin) ;%
        greenFig = imresize(greenIm{iDate,iRun},1/spatialBin,'OutputSize',[imSize(1) imSize(2)]);

        % IMAGING:  Check if the ROI center and the probe location has been picked
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
        if ~exist(fullfile(dataDir,'ROI.png'),'file')
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
            probeTemp = probeTemp(1:szMin,:);

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

%         % IMAGING: Upsampling full field of view and removing bad frames...
%         clear varInfo
%         if ~exist(fullfile(dataDir,'crossCorrFOV.mat'),'file') || 1
%             ccFullFlag = 0; fovFlag = 0;
%         else
%             varInfo = who('-file', fullfile(dataDir,'crossCorrFOV.mat'));
%             if find(ismember(varInfo,'ccFull'));  ccFullFlag = 1;  else; ccFullFlag = 0;  end
            fovFlag = 0;
%         end
% 
%         if ~exist(fullfile(dataDir,'crossCorrFOV.mat'),'file') || 1
%             disp('Upsampling full FOV and removing bad frames...');
%             fovFlag = 1; tic;
% 
%             clear inDatTemp goodDat ccFull ccFullAlpha
%             parfor iP = 1:size(processedDat10,1)
%                 inDatTemp = interp(processedDat10(iP,:),100);
%                 if ~(szLFP == szIm); inDatTemp = inDatTemp(:,1:szMin); end
%                 inDatTemp(:,badTimes) = []; % Remove bad frames
%                 goodDat(iP,:) = downsample(inDatTemp,100); % Downsample imaging data
%             end
%             toc;
%         end

        % Remove bad time segments from imaging and physiology as determined from spectrogram
        if strcmp(monkeyName,'CharlieSheen')
            if strcmp(expDate,'01_11_2022') && strcmp(runName,'run03') % Remove last 150 s of data
                probeTemp(end-150e3+1:end,:) = [];
                processedDat10(:,end-1500+1:end) = [];
                if fovFlag == 1; goodDat(:,end-1500+1:end) = []; end

            elseif strcmp(expDate,'01_11_2022') && strcmp(runName,'run04') % Remove 380-600 s data
                probeTemp(380e3+1:600e3,:) = [];
                processedDat10(:,3800+1:6000) = [];
                if fovFlag == 1; goodDat(:,3800+1:6000) = []; end
            end

        elseif strcmp(monkeyName,'Whiskey')
            if strcmp(expDate,'08_14_2023') && strcmp(runName,'run04') % Remove 320 - 410 s of data
                probeTemp(320e3+1:410e3,:) = [];
                processedDat10(:,3200+1:4100) = [];
                if fovFlag == 1; goodDat(:,3200+1:4100) = []; end

            elseif strcmp(expDate,'10_16_2023') && strcmp(runName,'run06') % Remove 430 - 510 s of data
                probeTemp(430e3+1:510e3,:) = [];
                processedDat10(:,4300+1:5100) = [];
                if fovFlag == 1; goodDat(:,4300+1:5100) = []; end

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

        % EPHYS:  Shuffle the data
        comb        = randperm(size(infraEphys,1));
        infraEphysS = infraEphys(comb);
        inDatS      = inDatUp(:,comb);
        infraAlphaS = infraEphysAlpha(comb);

        % XCORR: Get the cross correlation between LFP and ISOI for the area
        % around the electrode
        if ~exist(fullfile(dataDir,'crossCorrROI.mat'),'file') || 1
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
        else
            % Load the cross correlations for ROI
            crossCorrAll = load([dataDir '\crossCorrROI.mat'],'cc','ccOrig','ccS','ccA','ccAS','lags','clipMaskROI');
            crossCorr{iRun,iDate}          = crossCorrAll.cc;
            allLags{iRun,iDate}            = crossCorrAll.lags;
        end
%%
        % XCORR: Get the cross correlation between LFP and ISOI for the
        % full field of view...
        if ~exist(fullfile(dataDir,'crossCorrFOV.mat'),'file') || ~ccFullFlag || 1
            % Upsampling full field of view and removing bad frames...
            % Cross correlate LFP and ISOI for the entire FOV
            disp('Saving cross correlations for the entire field of view...')
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

            ccFull      = reshape(ccFull,[401 imSize(1) imSize(2)]);
            ccFullAlpha = reshape(ccFullAlpha,[401 imSize(1) imSize(2)]);

            disp('Saving cross correlations for the entire FOV...');
            save([dataDir '\crossCorrFOV.mat'],'ccFull','ccFullAlpha','lagFull','-v7.3');
            crossCorrFOV{iRun,iDate}      = matfile([dataDir '\crossCorrFOV.mat']);

            toc;
        else
            % Load the cross correlations for ROI
            crossCorrFOV{iRun,iDate} = matfile([dataDir '\crossCorrFOV.mat']);
        end

        clear x xNew negIdx vidNew xLow minMedcorrInd maxMedcorrIndLow frameNumMedLow
        x = allLags{iRun,iDate};
        negIdx = x<0 & x>=-150; xNew = x(negIdx);
        lowIdx = x<0 & x>= -80; xLow = x(lowIdx);

        % Find the frame where cross correlations are minimum
        [~,minMedcorrInd] = min(median(crossCorr{iRun,iDate}(negIdx,:,:),[2,3],'omitnan'));
        lagLow          = xNew(minMedcorrInd);
        frameNumMedLow  = (x == lagLow);

        % Find the correlation between FC map and map at peak negative
        clear pDatTemp seedSig ccFOV
        pDatTemp = processedDat{iDate,iRun}.tempBandPass;
        ccFOV    = crossCorrFOV{iRun,iDate}.ccFull;
        ccFOV = reshape(squeeze(ccFOV(frameNumMedLow,:,:)),[imSize(1)*imSize(2) 1]);

        % Correlate the FC map with map at peak negative...
        clear rad col row circlePixels pixelLoc spCorrAllT
        rad          = round(roiSize{iDate}(iRun)*2/spatialBin);
        [col, row]   = meshgrid(1:imSize(2), 1:imSize(1));
        circlePixels = (row - seedLocProbe(2)).^2 + (col - seedLocProbe(1)).^2 <= rad.^2;

        [pixelLoc(:,2), pixelLoc(:,1)] = find(circlePixels);
        sz1 = imSize(1); sz2 = imSize(2);  
        tic;
        parfor iMap = 1:size(pixelLoc,1) % Pixels inside 1mm radius
            seedSigT           = calculateSeedSignal(greenFig,corrMask,pixelLoc(iMap,:),12,pDatTemp); % Get Gaussian weighted seed signal
            corrMapT           = reshape(plotCorrMap(seedSigT,pDatTemp,1),[sz1*sz2 1]);
            spCorrAllT(iMap,1) = corr(corrMapT,ccFOV,'rows','complete');
        end
        toc;
        spCorrAll{iRun,iDate} = spCorrAllT;
        spCorr(iRun,iDate) = median(spCorrAll{iRun,iDate},'all','omitnan');

        % SPATIAL CONTROL 1 : Shift ROI 1mm by 1mm and cross correlate
        if ~exist([dataDir '\spatialControl.png'],'file')  || 1
            tic;
            cVals = {'w','k','b','g','m'};
            figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(131); imagesc(greenFig); hold on; colormap gray; axis image off;
            plot(seedLocProbe(1),seedLocProbe(2),'Marker','pentagram','MarkerSize',15,'MarkerFaceColor','r','MarkerEdgeColor','none');

            minShift  = round(roiSize{iDate}(iRun)/spatialBin);
            theta     = linspace(0,2*pi, round(pi*minShift)); % number of angles
            maxPoints = length(theta);

            spCorrControl{iRun,iDate} = NaN(5,maxPoints);
            spCorrDiff{iRun,iDate}    = NaN(5,maxPoints);

            for iShift = 1:5
                clear locShift loc
                if iShift == 1
                    locShift = round(roiSize{iDate}(iRun)/spatialBin);
                else
                    locShift = round(roiSize{iDate}(iRun)*2*(iShift-1)/spatialBin);
                end

                % Get coordinates for a circle
                loc(:,1) = round(locShift * cos(theta) + seedLocProbe(1));
                loc(:,2) = round(locShift * sin(theta) + seedLocProbe(2));

                % Plot the points sampled on the green blood vessel map
                plot(loc(:,1),loc(:,2),'.','Color',cVals{iShift},'MarkerSize',10);

                for iPoint = 1:size(loc,1)
                    clear seedSigT corrMapT
                    if sum(fliplr(loc(iPoint,:))+12>size(greenFig)) || (sum(loc(iPoint,:)-12<= 0))
                        spCorrControl{iRun,iDate}(iShift,iPoint) = NaN;
                    else
                        seedSigT = calculateSeedSignal(greenFig,clipMaskCortex,loc(iPoint,:),12,pDatTemp); % Get Gaussian weighted seed signal
                        corrMapT = reshape(plotCorrMap(seedSigT,pDatTemp,0),[imSize(1)*imSize(2) 1]);

                        spCorrControl{iRun,iDate}(iShift,iPoint) =  corr(corrMapT,ccFOV,'rows','complete');
                        spCorrDiff{iRun,iDate}(iShift,iPoint)    = spCorr(iRun,iDate) - corr(corrMapT,ccFOV,'rows','complete');
                    end
                end
            end

            subplot(132);boxplot(spCorrControl{iRun,iDate}',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});
            xlabel('Distance from probe (mm)'); ylabel('Correlation between FC map and peak negative map');

            subplot(133);boxplot(spCorrDiff{iRun,iDate}',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});
            xlabel('Distance from probe (mm)'); ylabel('Deviation from the actual correlation between FC map and peak negative map');

            f = gcf; exportgraphics(f,[dataDir '\spatialControl.png'],'Resolution',300); close gcf;
        end

        % TEMPORAL CONTROL: Shuffle time series in varying windows and cross correlate
        if ~exist([dataDir '\temporalControl.png'],'file') || 1
            timeLen = length(infraEphys);
            winLen  = [1 5 10 50 100 500 1000 5000 timeLen];

            for iWin = 1:length(winLen) % shuffle time series in windows of 10s
                clear comb1 comb2
                comb1 = randperm(round(timeLen/winLen(iWin)));
                comb2 = randperm(round(timeLen/winLen(iWin)));

                for iType = 1:3 % Shuffle type
                    clear infraEphysS roiS ccT

                    switch iType
                        case 1 % Shuffle ISOI
                            if iWin == 1
                                roiS = inDatUp(:,comb1);
                            elseif iWin == 9
                                roiS  = inDatUp;
                            else
                                roiS = [];
                                for iL = 1:length(comb1)
                                    clear win1 win2
                                    win1 = (comb1(iL)-1)*winLen(iWin)+1 : (comb1(iL)-1)*winLen(iWin)+winLen(iWin);
                                    win1(win1>timeLen) = [];
                                    roiS        = [roiS inDatUp(:,win1)];
                                end
                            end
                            infraEphysS = infraEphys;

                        case 2 % Shuffle Ephys
                            if iWin == 1
                                infraEphysS  = infraEphys(comb1);
                            elseif iWin == 9
                                infraEphysS  = infraEphys;
                            else
                                infraEphysS = [];
                                for iL = 1:length(comb1)
                                    clear win1 win2
                                    win1 = (comb1(iL)-1)*winLen(iWin)+1 : (comb1(iL)-1)*winLen(iWin)+winLen(iWin);
                                    win1(win1>timeLen) = [];
                                    infraEphysS = [infraEphysS; infraEphys(win1)];
                                end
                            end
                            roiS = inDatUp;

                        case 3 % Shuffle both
                            if iWin == 1
                                infraEphysS = infraEphys(comb1);
                                roiS        = inDatUp(:,comb2);

                            elseif iWin == 9
                                infraEphysS  = infraEphys;
                                roiS = inDatUp;

                            else
                                infraEphysS= [];  roiS = [];
                                for iL = 1:length(comb1)
                                    clear win1 win2
                                    win1 = (comb1(iL)-1)*winLen(iWin)+1 : (comb1(iL)-1)*winLen(iWin)+winLen(iWin);
                                    win2 = (comb1(iL)-1)*winLen(iWin)+1 : (comb1(iL)-1)*winLen(iWin)+winLen(iWin);
                                    win1(win1>timeLen) = []; win2(win2>timeLen) = [];

                                    infraEphysS = [infraEphysS; infraEphys(win1)];
                                    roiS        = [roiS inDatUp(:,win2)];
                                end
                            end
                    end

                    % Check size of time series
                    szEphys = size(infraEphysS,1); szROI = size(roiS,2);
                    if ~isequal(szEphys,szROI)
                        minSize = min(szEphys,szROI);
                        infraEphysS = infraEphysS(1:minSize);
                        roiS = roiS(:,1:minSize);
                    end

                    parfor iP = 1:size(roiS,1)
                        [ccT(:,iP),~]  = xcorr(infraEphysS',roiS(iP,:),200,'normalized');
                    end

                    ccT  = reshape(ccT,[401 inDatSize(1) inDatSize(2)]);
                    ccWinROI(iRun,iDate,iWin,iType) = median(squeeze(ccT(frameNumMedLow,:,:)),'all','omitnan');
                end
            end

            figure;
            for iType = 1:3
                plot(movmean(squeeze(ccWinROI(iRun,iDate,:,iType)),[0 1]),'LineWidth',1); xticklabels(winLen);  hold on;
            end
            xlabel('Window Length (samples)'); ylabel('Median correlation at peak negative');
            legend('ISOI Shuffled','Ephys Shuffled','Both Shuffled','Location','southwest'); grid on;
            ylim([min(squeeze(ccWinROI(iRun,iDate,:)))-0.1 0.15]);
            f = gcf; exportgraphics(f,[dataDir '\temporalControl.png'],'Resolution',300); close gcf;

        end

        %         else
        %             disp('Cross correlations already computed, proceeding with the next steps...');
        %
        %             % Load the cross correlations for ROI
        %             crossCorrAll = load([dataDir '\crossCorrROI.mat'],'cc','ccOrig','ccS','ccA',...
        %                 'ccAS','lags','clipMaskROI');
        %             crossCorr{iRun,iDate}          = crossCorrAll.cc;
        %             crossCorrOrig{iRun,iDate}      = crossCorrAll.ccOrig;
        %             crossCorrScrambled{iRun,iDate} = crossCorrAll.ccS;
        %             crossCorrAlpha{iRun,iDate}     = crossCorrAll.ccA;
        %             crossCorrAlphaS{iRun,iDate}    = crossCorrAll.ccAS;
        %             allLags{iRun,iDate}            = crossCorrAll.lags;
        %             roiMask{iRun,iDate}            = crossCorrAll.clipMaskROI;
        %
        %             % Load the cross correlations for the full FOV
        %             crossCorrFOV{iRun,iDate} = matfile([dataDir '\crossCorrFOV.mat']);

        %         end

      %% Spatial (ROI) and temporal lag profiles for gamma band
        [medCorrLagGamma(iRun,iDate), medCorrLagGammaLow(iRun,iDate)] = plotLagProfiles(dataDir,'Gamma',...
            crossCorr{iRun,iDate},allLags{iRun,iDate},monkeyName,expDate,runName,roiMask{iRun,iDate});

        % Spatial (ROI) and temporal profiles for alpha band
        [medCorrLagAlpha(iRun,iDate), medCorrLagAlphaLow(iRun,iDate)] = plotLagProfiles(dataDir,'Alpha',...
            crossCorrAlpha{iRun,iDate},allLags{iRun,iDate},monkeyName,expDate,runName,roiMask{iRun,iDate});

       %% Spatial profiles for the FOV at peak positive and negative...
        % Gamma band
        showPeakXcorrFOV(dataDir,'Gamma',monkeyName,expDate,runName,greenFig,crossCorrFOV{iRun,iDate}.ccFull,...
            corrMask,medCorrLagGamma(iRun,iDate),medCorrLagGammaLow(iRun,iDate),allLags{iRun,iDate});
        % Alpha band
        showPeakXcorrFOV(dataDir,'Alpha',monkeyName,expDate,runName,greenFig,crossCorrFOV{iRun,iDate}.ccFullAlpha,...
            corrMask,medCorrLagAlpha(iRun,iDate),medCorrLagAlphaLow(iRun,iDate),allLags{iRun,iDate});

        %         % Save the spatial profiles for all lags (negative)
        %         saveVideoFOV(dataDir,'Gamma',monkeyName,expDate,runName,crossCorrFOV{iRun,iDate},...
        %             clipMaskCortex,allLags{iRun,iDate});

    end
end


%% Functions that are used in this script
