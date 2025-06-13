%% Comparing cross-modal maps with RS-averaged FC maps
% This code compares the recording-wise FC map for ISOI+ephys experiments
% with the average FC maps
% February 28,2025 - KM
% Set paths
clc; clear;
commonDir = 'C:\Users\kem294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\ISOI_Ephys\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12\fly_track\videoIO']));

%% Get cross-modal map for certain example runs (one accepted, one rejected)

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
% Initialize variables
hemisphere = 'Left'; spatialBin = 3;

for iSession = 1:size(dateVals,1)
    clear expDate serverDir runName crossCorrFOV
    expDate   = dateVals(iSession,:);
    runName   = runNamesAll(iSession,:);

    serverDir = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\kem294\Data\' ...
        monkeyName{iSession} '_SqM\' hemisphere ' Hemisphere\'  expDate '\' runName];

    dataDir = ['D:\Data\' monkeyName{iSession} '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName];
    
    % Get location of master blood vessel maps 
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
    
    % Get masks
    [clipMaskCortex, corrMask] = getMasks(dataDir,runName);
    corrMaskT                  = reshape(corrMask,[imSize(1)*imSize(2) 1]);

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

    fcMap = plotCorrMap(seedSigT,pDatTemp,0); % Get FC map

    % Register the green and the mask to the master green
    if ~exist([dataDir '\imRegisMaster.mat'],'file') % Check if image registration has been done already
        clear movPointsTemp fixedPointsTemp OTemp spacingTemp corrMapFinal
        
        % Get FC map averaged over multiple RS runs
        figure('units','normalized','outerposition',[0 0 1 1]);
        imagesc(greenMapRef);axis image off; colormap gray; hold on;
        title(strrep(['Select the location where the probe is located for ' monkeyName{iSession} ' ' expDate ' ' runName],'_','\_'));
        refSeed(1,:,:) = ginput(1); close gcf;
        [~,corrMapFinal] = getRSConnectivityMaps(squeeze(refSeed(1,:,:))',monkeyName{iSession}); % Get average FC map

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
    
    % Get masks 
    greenSize = size(greenMapRef);
    tempMask  = zeros(size(greenMapRef));
    tempMask(1:size(greenTemp,1),1:size(greenTemp,2)) = imresize(corrMask,[1082 1312]);

    tempFCMap = zeros(size(greenMapRef));
    tempFCMap(1:size(greenTemp,1),1:size(greenTemp,2)) = imresize(fcMap,[1082 1312]);
    
    % Warp the recording-wise map to the average
    [imMaskTemp,tMaskTemp]   = bspline_transform(OTemp,tempMask,spacingTemp,3); % bicubic spline interpolation
    [imFCMapTemp,tFCMapTemp] = bspline_transform(OTemp,tempFCMap,spacingTemp,3);

    imMaskTempR  = logical(reshape(imMaskTemp,[greenSize(1)*greenSize(2) 1]));
    imFCMapTempR = reshape(imFCMapTemp,[greenSize(1)*greenSize(2) 1]);
    imFCMapTempR(~imMaskTempR) = NaN;

    corrMapFinalR = reshape(corrMapFinal,[greenSize(1)*greenSize(2) 1]);
    corrMapFinalR(~imMaskTempR) = NaN; % Using the same masks to include the common pixels for correlationg 

    % Correlate recording-wise map to the average 
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