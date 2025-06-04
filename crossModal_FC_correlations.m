%% Comparing cross-modal maps with RS-averaged FC maps
% This code compares the recording-wise FC map for ISOI+ephys experiments
% with the average FC maps
% February 28,2025 - KM
% Set paths
clear;clc
commonDir = 'C:\Users\kem294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\Imaging']));
addpath(genpath([commonDir '\Codes\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\fly_track\videoIO']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\spectral_analysis\continuous\dupes']));
clc;

%% Get cross-modal map for certain example runs

% dateVals    = ['10_16_2023'; '10_16_2023'; '10_16_2023'; '12_04_2023'; '12_04_2023'; '04_29_2024' ; '04_29_2024'; '04_29_2024']; 
% runNamesAll = ['run04'; 'run05'; 'run06'; 'run04'; 'run05'; 'run02' ; 'run04'; 'run08']; 
% roiSize     = [38 38 38 38 38 38 38 38]; 
% monkeyName  = 'Whiskey';

dateVals = ['10_16_2023'; '08_07_2023'];
runNamesAll = ['run06';'run08'];
roiSize = [38 38];
monkeyName = {'Whiskey', 'CharlieSheen'};

%% Initialize variables and get monkey data
% iM = 2;
% switch iM
%     case 1
%         monkeyName = 'CharlieSheen';
%         goodRuns  = ([1 1 NaN NaN NaN;... % Date x run
%             1 1 NaN NaN NaN; ...
%             1 1 1 0 1; ...
%             1 1 1 1 NaN]);
% 
%         singleChFlag  = ([0 0 NaN NaN NaN;... % Date x run
%             0 0 NaN NaN NaN; ...
%             0 0 0 0 0; ...
%             0 0 0 0 NaN]);
% 
%         goodRunsSpatial  = ([1 1 NaN NaN NaN;... % Date x run
%             1 1 NaN NaN NaN; ...
%             1 1 1 0 1; ...
%             1 1 1 1 NaN]);
% 
%         smFlag  = (['S' 'M' '#' '#' '#';... % Date x run
%             'S' 'M' '#' '#' '#'; ...
%             'S' 'S' 'M' 'M' 'M'; ...
%             'S' 'M' 'M' 'M' '#']);
% 
%         isoLevel = ([0.75 0.75 NaN NaN NaN;... % Date x run
%             0.75 0.75 NaN NaN NaN; ...
%             0.9 0.9 0.9 0.9 0.9; ...
%             0.75 0.8 0.8 0.9 NaN]);
% 
%         refDate      = '08_31_2021';
%         refDir       = [commonDir '\' monkeyName '_SqM\' hemisphere ' Hemisphere\' refDate '\Master Green Images\'];
%         refImageName = 'Charlie Sheen Combined Green 08_31_2021';
% 
% 
%     case 2
%         monkeyName = 'Whiskey';
%         goodRuns   = ([1 1 1 NaN NaN NaN NaN; ... % Date x run
%             1 1 1 1 1 1 1 ; ...
%             1 0 1 1 NaN NaN NaN; ...
%             1 1 0 1 1 1 0; ...
%             1 0 0 1 1 1 1]);
% 
%         singleChFlag = ([1 1 1 NaN NaN NaN NaN; ... % Date x run
%             0 0 0 0 0 0 0 ; ...
%             0 0 0 0 NaN NaN NaN; ...
%             0 0 0 0 0 0 0; ...
%             0 0 0 0 0 0 0]);
% 
%         goodRunsSpatial = ([1 1 1 NaN NaN NaN NaN; ... % Date x run
%             1 1 1 1 1 1 1 ; ...
%             1 0 1 1 NaN NaN NaN; ...
%             1 1 0 1 1 1 1; ...
%             1 0 0 1 1 1 0]);
% 
%          smFlag = (['M' 'M' 'M' '#' '#' '#' '#'; ... % Date x run; # - in place of NaNs
%             'M' 'M' 'M' 'M' 'S' 'S' 'S' ; ...
%             'S' 'M' 'M' 'S' '#' '#' '#'; ...
%             'S' 'M' 'M' 'S' 'M' 'M' 'M'; ...
%             'S' 'S' 'M' 'M' 'M' 'S' 'M']);
% 
%          isoLevel = ([1 1.1 1.25 NaN NaN NaN NaN; ... % Date x run
%             1.2 1 1 1.1 1.1 1.1 1.1 ; ...
%             0.8 1.05 1.75 1.75 NaN NaN NaN; ...
%             0.9 1 1 1 1 1 1; ...
%             0.7 0.8 1 1.1 1.3 1.3 1.3]);
% 
%          refDate      = '05_09_2022';
%          refDir       = [commonDir '\' monkeyName '_SqM\' hemisphere ' Hemisphere\' refDate '\Master Green Images\'];
%          refImageName = 'Combined Green';
% 
% end
% 
% % Get index of good runs
% goodRuns = reshape(goodRuns,[size(goodRuns,1)*size(goodRuns,2) 1]);
% goodRuns(isnan(goodRuns)) = []; goodRuns = logical(goodRuns);
% 
% % Check if there are single channel electrode recordings
% singleChFlag = reshape(singleChFlag,[size(singleChFlag,1)*size(singleChFlag,2) 1]);
% singleChFlag(isnan(singleChFlag)) = []; singleChFlag = logical(singleChFlag);
% 
% % Flag bad runs
% goodRunsSpatial = reshape(goodRunsSpatial,[size(goodRunsSpatial,1)*size(goodRunsSpatial,2) 1]);
% goodRunsSpatial(isnan(goodRunsSpatial)) = []; goodRunsSpatial = logical(goodRunsSpatial);
% 
% % Get electrode penetration location (somatosensory or motor)
% smFlag = reshape(smFlag,[size(smFlag,1)*size(smFlag,2) 1]);
% smFlag((smFlag == '#')) = [];
% sensoryGoodRuns = (smFlag =='S') & goodRuns;
% sensoryGoodSpatialRuns = (smFlag == 'S') & goodRunsSpatial;
% motorGoodRuns = (smFlag == 'M') & goodRuns;
% motorGoodSpatialRuns = (smFlag == 'M') & goodRuns;
% 
% % Get isoflurane levels
% isoLevel = reshape(isoLevel,[size(isoLevel,1)*size(isoLevel,2) 1]);
% isoLevel(isnan(isoLevel)) = [];
% isoLevelGoodRuns = isoLevel; isoLevelGoodRuns(~goodRuns) = [];
% isoLevelSpatial = isoLevel; isoLevelSpatial(~goodRunsSpatial) = [];

%% Get recording-wise and average FC maps

hemisphere = 'Left'; spatialBin = 3;

for iSession = 2%:size(dateVals,1)
    clear expDate serverDir runName crossCorrFOV
    expDate   = dateVals(iSession,:);
    runName   = runNamesAll(iSession,:);
    serverDir = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\kem294\Data\' ...
        monkeyName{iSession} '_SqM\' hemisphere ' Hemisphere\'  expDate '\' runName];
    dataDir = ['D:\Data\' monkeyName{iSession} '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName];
    
    if strcmp(monkeyName{iSession},'Whiskey')
        refDate      = '05_09_2022';
        refDir       = [commonDir '\' monkeyName{iSession} '_SqM\' hemisphere ' Hemisphere\' refDate '\Master Green Images\'];
        refImageName = 'Combined Green';

    elseif strcmp(monkeyName{iSession},'CharlieSheen')
        refDate      = '08_31_2021';
        refDir       = [commonDir '\' monkeyName{iSession} '_SqM\' hemisphere ' Hemisphere\' refDate '\Master Green Images\'];
        refImageName = 'Charlie Sheen Combined Green 08_31_2021';
    end 

    % Get the master green image
    if exist([refDir refImageName '.bmp'],'file') == 0 % Make sure to get what the reference run is outside the function
        greenMapRef = imread([refDir refImageName '.png']);
    else
        greenMapRef = imread([refDir refImageName '.bmp']);
    end
    greenMapRef  = greenMapRef(:,:,1);

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

%     % Get cross-modal maps
%     crossCorrFOV = matfile([serverDir '\crossCorrFOV_6_NoRef.mat']);
%     allXCorr     = crossCorrFOV.spatialProfile;
%     allLags      = crossCorrFOV.lagFull;
% 
%     allHybridVars = matfile([dataDir '\processedHybridMapVars.mat']);
%     corrFCHybrid  = allHybridVars.corrFCHybridT; 
% 
%     peakNegValsAll  = allHybridVars.peakNegValsAllT;  peakNegValsAll  = peakNegValsAll(:,2);
%     peakNegTimesAll = allHybridVars.peakNegTimesAllT; peakNegTimesAll = peakNegTimesAll(:,2);
%     x = -200:200;
% 
%     for iX = 1:size(peakNegTimesAll,1)
%         peakNegIdx(iX,1) = find(x== peakNegTimesAll(iX).*10);
%     end
% 
%     crossCorrAlpha = reshape(allXCorr.ccFullAlpha,[401 imSize(1)*imSize(2)]); crossCorrAlpha(:,~corrMaskT)   = NaN;
%     crossCorrBeta  = reshape(allXCorr.ccFullBeta,[401 imSize(1)*imSize(2)]);  crossCorrBeta(:,~corrMaskT)    = NaN;
%     crossCorrGamma = reshape(allXCorr.ccFull,[401 imSize(1)*imSize(2)]);      crossCorrGamma(:,~corrMaskT)   = NaN;
%     crossCorrSpiking = reshape(allXCorr.ccFullRaw,[401 imSize(1)*imSize(2)]); crossCorrSpiking(:,~corrMaskT) = NaN;
%     
%     alphaCrossModal   = squeeze(crossCorrAlpha(peakNegIdx(2),:));
%     betaCrossModal    = squeeze(crossCorrBeta(peakNegIdx(3),:));
%     gammaCrossModal   = squeeze(crossCorrGamma(peakNegIdx(4),:));
%     spikingCrossModal = squeeze(crossCorrGamma(peakNegIdx(5),:));

    % Get ROI location and session-wise FC map
    processedDat = matfile([dataDir '\processedFrames.mat']);
    seedLocIn    = load([dataDir '\roiCenterLoc.mat']);
    seedLocProbe = seedLocIn.seedLocProbe;
    seedLocIn    = seedLocIn.seedLocIn;
    circleRad    = round(roiSize(iSession)/(spatialBin)); % 500um radius
    greenFig     = imresize(greenTemp,1/spatialBin,'OutputSize',[imSize(1) imSize(2)]);

    pDatTemp = processedDat.tempBandPass;
    seedSigT = calculateSeedSignal(greenFig,corrMask,...
        seedLocIn,circleRad,pDatTemp); % Get Gaussian weighted seed signal

    fcMap             = plotCorrMap(seedSigT,pDatTemp,0);

%     [optimizer,metric] = imregconfig('multimodal');
%     tformRunWise = imregtform(imresize(greenIm{iDate,iRun},1/spatialBin),imresize(greenIm{iDate,1},1/spatialBin),'affine',optimizer,metric);

    % Register the green and the mask to the master green
    if ~exist([dataDir '\imRegisMaster.mat'],'file')
        clear movPointsTemp fixedPointsTemp OTemp spacingTemp corrMapFinal
        
        % Get FC map averaged over multiple RS runs
        figure('units','normalized','outerposition',[0 0 1 1]);
        imagesc(greenMapRef);axis image off; colormap gray; hold on;
        title(strrep(['Select the location where the probe is located for ' monkeyName{iSession} ' ' expDate ' ' runName],'_','\_'));
        refSeed(1,:,:) = ginput(1); close gcf;
        [~,corrMapFinal] = getRSConnectivityMaps(squeeze(refSeed(1,:,:))',monkeyName{iSession});

        figure('units','normalized','outerposition',[0 0 1 1]);
        imagesc(ind2rgb(greenMapRef,gray(256))); hold on;
        imagesc(corrMapFinal,'alphaData',corrMapFinal.*0.9); axis image off; hold on;
        plot(refSeed(1,:,1),refSeed(1,:,2),'w.','MarkerSize',35);
        colormap jet; clim([0 1]);
        
        % Nonlinear registration to master green
        [movPointsTemp,fixedPointsTemp] = cpselect(greenTemp,greenMapRef,'Wait',true);
        [OTemp,spacingTemp,~] = point_registration(size(greenMapRef),fliplr(fixedPointsTemp),fliplr(movPointsTemp));
        save([dataDir '\imRegisMaster.mat'],'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp','corrMapFinal');
 
    else
        regisPoints = load([dataDir '\imRegisMaster.mat'],'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp','corrMapFinal');
        OTemp = regisPoints.OTemp;
        spacingTemp = regisPoints.spacingTemp; 
        corrMapFinal = regisPoints.corrMapFinal;
    end

    greenSize = size(greenMapRef);
    tempMask  = zeros(size(greenMapRef));
    tempMask(1:size(greenTemp,1),1:size(greenTemp,2)) = imresize(corrMask,[1082 1312]);

    tempFCMap = zeros(size(greenMapRef));
    tempFCMap(1:size(greenTemp,1),1:size(greenTemp,2)) = imresize(fcMap,[1082 1312]);
    
    % Need to recheck how I do the correlations...
    [imMaskTemp,tMaskTemp]   = bspline_transform(OTemp,tempMask,spacingTemp,3); % bicubic spline interpolation
    [imFCMapTemp,tFCMapTemp] = bspline_transform(OTemp,tempFCMap,spacingTemp,3);

    imMaskTempR = logical(reshape(imMaskTemp,[greenSize(1)*greenSize(2) 1]));
    imFCMapTempR = reshape(imFCMapTemp,[greenSize(1)*greenSize(2) 1]);
    imFCMapTempR(~imMaskTempR) = NaN;
    corrMapFinalR = reshape(corrMapFinal,[greenSize(1)*greenSize(2) 1]);
    corrMapFinalR(~imMaskTempR) = NaN;

    corrFC_Avg_SessionWise(iSession) = corr(imFCMapTempR,corrMapFinalR,'rows','complete');

end


%% Plot an example average FC and recording-wise FC - get the variables for a specific run from the previous section
% Average FC map
figure; subplot(121);imagesc(ind2rgb(greenMapRef,gray(256))); hold on;
corrMapAvg = reshape(corrMapFinalR,[greenSize(1) greenSize(2)]);
imagesc(corrMapFinal,'AlphaData',corrMapFinal.*1.0); colormap jet; clim([0 1]);
axis image off; colorbar;title('Average FC map');

% Recording-wise FC map
fcMapWarped = reshape(imFCMapTempR,[greenSize(1) greenSize(2)]);
subplot(122);imagesc(ind2rgb(greenMapRef,gray(256))); hold on;
imagesc(fcMapWarped,'AlphaData',fcMapWarped.*1); colormap jet; clim([0 1]);
axis image off; colorbar;title('Recording-wise FC map');






