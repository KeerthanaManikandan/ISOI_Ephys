% This code correlates average FC map with session-wise FC maps from the
% first procedure to the last. The purpose of this analysis is to compare
% the quality of imaging data and recording stability over time. 
% Edited on June 6, 2025 - KM
% Codes for ISOI_Ephys Paper Supplementary Figure 6

% Set paths
clc; clear;
commonDir = 'C:\Users\kem294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\ISOI_Ephys\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12\fly_track\videoIO']));

%% Get monkey data
monkeyName = 'Whiskey'; hemisphere = 'Left';

[allDates,allRuns, refDate, refDir,lensCombo, roiSize, ephysFileNameAll, serverPath,probeLabel,...
    chInCortexNotes, greenMapRef] = getMonkeyParams_Imaging_Ephys(monkeyName,commonDir, hemisphere);

dateVals   = [1 3 4];
checkDates = ['04_25_2022';'08_08_2022'; allDates(dateVals,:)];
runNames00 = ['run00'; 'run00';'run00'; 'run00'; 'run00'];
spatialBin = 3;

%% Get a grid of 100 points from the last procedure: 02_20_2024
% Load the green from 02_20_2024
clear expDate dataDir templateDir
expDate    = checkDates(end,:);

% Load all the information
runName    = runNames00(end,:);
fileNum    = 0;

% Get the directory and filenames
dataDir  = ['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
templateDir = ['D:\Data\' monkeyName '_SqM\Left Hemisphere\'];
cDriveFlag  = 0;
if isempty(dir(dataDir)); cDriveFlag = 1; end

% Get greens and masks
[greenTemp,elecMask,clipMaskCortex,corrMask]= getGreensAndMasks(serverPath,dataDir,templateDir,expDate,runName,runNames00);

% Get image registration parameters
if exist([dataDir '\imRegisMaster.mat'],'file') || exist(['C:\Users\kem294\Documents\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate  '\imageRegisPoints_RS.mat'],'file')
    if ~cDriveFlag
        regisPoints = load([dataDir '\imRegisMaster.mat'],'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp');
    else
        regisPoints = load(['C:\Users\kem294\Documents\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate  '\imageRegisPoints_RS.mat'],...
            'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp');
    end

    movPointsTemp   = regisPoints.movPointsTemp;
    fixedPointsTemp = regisPoints.fixedPointsTemp; 

else % Map points to perform registration
    clear movPointsTemp fixedPointsTemp OTemp spacingTemp
    [movPointsTemp,fixedPointsTemp] = cpselect(greenTemp,greenMapRef,'Wait',true);
    [OTemp,spacingTemp,~] = point_registration(size(greenMapRef),fliplr(fixedPointsTemp),fliplr(movPointsTemp));

    if ~cDriveFlag
        save([dataDir '\imRegisMaster.mat'],'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp');
    else
        save(['C:\Users\kem294\Documents\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate  '\imageRegisPoints_RS.mat'],...
            'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp');
    end
end

% Perform registration from one session to average map
[OTemp,spacingTemp,~] = point_registration(size(greenMapRef),fliplr(movPointsTemp),fliplr(fixedPointsTemp));

% Find points to place the grid
imSize = size(greenTemp,[1,2]); % Get image size

% Get a 10x10 grid
x = linspace(1,imSize(1),10);
y = linspace(1,imSize(2),10);
[x,y] = meshgrid(x,y);

corrMaskTemp = imresize(corrMask,imSize);
seedRad =  round(roiSize{5}(1)); % seed radius for FC map- 500um

% check if any of the points are on a vessel/sulci or are out of the image
clear point seedLocNonLinear
iPoint = 1;
for iRow = 1:10
    for iCol = 1:10
        clear pointTemp
        pointTemp = [x(iRow,iCol),y(iRow,iCol)];
        if ~(any(((pointTemp(1,:))+seedRad)>size(greenTemp,[1,2])) || any(pointTemp(1,:)-seedRad<= 0)|| corrMaskTemp(round(pointTemp(1)),round(pointTemp(2)))==0)
            seedLocNonLinear(iPoint,:) = bspline_trans_points_double(OTemp,spacingTemp,[pointTemp(1) pointTemp(2)]);
            point(iPoint,:) = pointTemp;
            iPoint = iPoint+1;
        end
    end
end

% Get average FC maps for all points
clc;
if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\avgCorrMapControls_grid.mat'],'file')
    for iSeed = 1:size(seedLocNonLinear,1)
        disp(['Obtaining average FC for Seed: ' num2str(iSeed)]);
        [~,corrMapFinal(:,:,iSeed)] = getRSConnectivityMaps(seedLocNonLinear(iSeed,:)',monkeyName);

    end
        save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\avgCorrMapControls_grid.mat'],'corrMapFinal','seedLocNonLinear');

else % Load the saved average FC map and seed locations
    corrMapFinal = load(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\avgCorrMapControls_grid.mat']);
    seedLocNonLinear  = corrMapFinal.seedLocNonLinear;
    corrMapFinal = corrMapFinal.corrMapFinal;
end

%% Obtain the session-wise FC map for the chosen sessions
for iDate = 1:size(checkDates,1)
    clear expDate dataDir templateDir
    expDate    = checkDates(iDate,:);

    % Load all the information
    runName    = runNames00(iDate,:);
    fileNum    = 0;

    % Get the directory and filenames
    dataDir  = ['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
    templateDir = ['D:\Data\' monkeyName '_SqM\Left Hemisphere\'];
    cDriveFlag  = 0;

    if isempty(dir(dataDir))
        if (exist(['C:\Users\kem294\Documents\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ],'dir'))
            dataDir     = ['C:\Users\kem294\Documents\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
            templateDir = ['C:\Users\kem294\Documents\Data\' monkeyName '_SqM\Left Hemisphere\'];
            cDriveFlag  = 1;
        else
            [~,~] = mkdir(dataDir);
        end
    end

    % Loading the downsampled data
    numFiles    = length(dir([dataDir '\Spatial Downsample SS3' '/*.mat']));
    datName     = 'Data_RS_10Hz_SS3_';

    % Get greens and masks
    [greenTemp,elecMask,clipMaskCortex,corrMask]= getGreensAndMasks(serverPath,dataDir,templateDir,expDate,runName,runNames00);
   
    greenIm{iDate,1}  = greenTemp(:,:,1);

    % Store/Retrieve processed rs-ISOI data
    if strcmp(expDate,'06_22_2021')
        serverDataPath    = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\302-19_CharlieSheen_SqM\Left Hemisphere\06_22_2021 Craniotomy\' runName];
    else
        serverDataPath = [serverPath expDate '\' runName];
    end

    % Store rs-ISOI data if not saved already
    if ~exist([dataDir '\processedFrames.mat'],'file')
        clear tempBandPass
        disp(['Processing imaging data for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)]);
        [~,~,~,~,tempBandPass] = getPreProcessedDataRestingState(serverDataPath,dataDir,runName,numFiles,spatialBin,datName);
        disp(['Storing imaging data for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)]);
        save([dataDir '\processedFrames.mat'],'tempBandPass');
        processedDat{iDate,1} = matfile([dataDir '\processedFrames.mat']); clear tempBandPass;

    else % Retrieve processed data        
        disp(['Loading imaging data for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)]);
        processedDat{iDate,1} = matfile([dataDir '\processedFrames.mat']);
    end

    cd(commonDir);

    clear pDatTemp fcMap seedLoc seedSig OTemp spacingTemp

    % Get rs-ISOI data
    pDatTemp = processedDat{iDate,1}.tempBandPass; 

    % Perform nonlinear registration
    clear regisPoints

    % Check if registration has been performed already
    if exist([dataDir '\imRegisMaster.mat'],'file') || exist(['C:\Users\kem294\Documents\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate  '\imageRegisPoints_RS.mat'],'file')
        if ~cDriveFlag
            regisPoints = load([dataDir '\imRegisMaster.mat'],'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp');
        else
            regisPoints = load(['C:\Users\kem294\Documents\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate  '\imageRegisPoints_RS.mat'],...
                'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp');
        end

        OTemp = regisPoints.OTemp;
        spacingTemp = regisPoints.spacingTemp;

    else % Map points to perform registration
        clear movPointsTemp fixedPointsTemp OTemp spacingTemp
        [movPointsTemp,fixedPointsTemp] = cpselect(greenTemp,greenMapRef,'Wait',true);
        [OTemp,spacingTemp,~] = point_registration(size(greenMapRef),fliplr(fixedPointsTemp),fliplr(movPointsTemp));

        if ~cDriveFlag
            save([dataDir '\imRegisMaster.mat'],'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp');
        else
            save(['C:\Users\kem294\Documents\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate  '\imageRegisPoints_RS.mat'],...
                'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp');
        end
    end

    % Perform warping to register images
    greenSize = size(greenMapRef);
    tempMask  = zeros(size(greenMapRef));
    tempMask(1:size(greenTemp,1),1:size(greenTemp,2)) = imresize(corrMask,[1082 1312]);
    [imMaskTemp,tMaskTemp]   = bspline_transform(OTemp,tempMask,spacingTemp,3); % bicubic spline interpolation

    for iSeed = 1:size(seedLocNonLinear,1)
        seedLoc = round(bspline_trans_points_double(OTemp,spacingTemp,[seedLocNonLinear(iSeed,:)])); % Non-linear seed transformation
        try
            seedSig = calculateSeedSignal(imresize(greenIm{iDate,1},1/spatialBin), clipMaskCortex,fliplr(round(seedLoc./spatialBin)),12,pDatTemp); % Get Gaussian weighted seed signal
            fcMap(:,:,iSeed) = plotCorrMap(seedSig,pDatTemp,0); % Obtain FC map

            tempFCMap = NaN(size(greenMapRef));
            tempFCMap(1:size(greenTemp,1),1:size(greenTemp,2)) = imresize(squeeze(fcMap(:,:,iSeed)),[1082 1312]);

            [imFCMapTemp,~] = bspline_transform(OTemp,tempFCMap,spacingTemp,3); % Warp images

            imMaskTempR = logical(reshape(imMaskTemp,[greenSize(1)*greenSize(2) 1]));
            imFCMapTempR = reshape(imFCMapTemp,[greenSize(1)*greenSize(2) 1]);
            imFCMapTempR(~imMaskTempR) = NaN;

            corrMapFinalR = reshape(squeeze(corrMapFinal(:,:,iSeed)),[greenSize(1)*greenSize(2) 1]);
            corrMapFinalR(~imMaskTempR) = NaN;

            % Correlate average FC with session-wise FC
            corrFC_Avg_SessionWise(iSeed,iDate) = corr(imFCMapTempR,corrMapFinalR,'rows','complete');
        catch
            corrFC_Avg_SessionWise(iSeed,iDate) = NaN;
        end
    end
end

% Remove seeds #44 and #51 as they are on the lateral sulcus
corrFC_Avg_SessionWise([44 51],:) = NaN;

%% Plot the correlations between session-wise and average FC maps
% Plot the seeds chosen on the average green vessel map
figure;imagesc(greenMapRef); hold on;plot(seedLocNonLinear(:,2),seedLocNonLinear(:,1),'.w','MarkerSize',20)
colormap gray;axis image off;

% Plot the correlations 
numSeeds = size(corrFC_Avg_SessionWise,1);

figure; % Show the distribution
x1 = (reshape(repmat(1:5,[numSeeds 1]),[numSeeds*5 1]));
y1 = reshape(corrFC_Avg_SessionWise,[numSeeds*5 1]);
s = swarmchart(x1,y1,20,[0 0.4470 0.7410]); 
s.XJitterWidth = 0.5; hold on;

% Show the boxplot
boxplot(corrFC_Avg_SessionWise,'Labels',{0, 4, 16, 20, 23});
xlabel('Months'); ylabel('Correlations between session and average FC');
box off;ylim([0 1]); 

% ANOVA to compare the distributions
[pSpCorrT,tblSpCorrT,statsSpCorrT] = anova1(corrFC_Avg_SessionWise,{'0', '4', '16', '20', '23'},'off');
[rSpCorrT,mSpCorrT,~,gnamesSpCorrT] = multcompare(statsSpCorrT,"CriticalValueType","bonferroni","Alpha", 0.001);

tblSpCorrMT = array2table(rSpCorrT,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblSpCorrMT.("Group") = gnamesSpCorrT(tblSpCorrMT.("Group"));
tblSpCorrMT.("Control Group") = gnamesSpCorrT(tblSpCorrMT.("Control Group"));

%% Get greens and masks
function [greenTemp,elecMask,clipMaskCortex,corrMask]= getGreensAndMasks(serverPath,dataDir,templateDir,expDate,runName,runNames00)

fileInfo    = dir(dataDir);

% Get greens - Check if greens exist in the data folder
if isempty(find(strcmp({fileInfo.name},['green' runName(end-1:end) '.png']), 1)) && isempty(find(strcmp({fileInfo.name},['green' runName(end-1:end) '.bmp']), 1))
    try
        copyfile([serverPath '\' expDate '\' runName '\green' runName(end-1:end) '.png'], dataDir);
    catch
        copyfile([serverPath '\' expDate '\' runName '\green' runName(end-1:end) '.bmp'], dataDir);
    end
end

% Check if greens are edited and load greens
if exist([templateDir  expDate '\' runNames00(end,:) '\green0' runNames00(end,5) '_Edited.png'],'file')
    greenTemp = imread([dataDir '\green0' runNames00(end,5) '_Edited.png']);

elseif exist([templateDir  expDate '\' runNames00(end,:) '\green0' runNames00(end,5) '_Edited.bmp'],'file')
    greenTemp = imread([dataDir '\green0' runNames00(end,5) '_Edited.bmp']);

else
    error('Greens are not edited...');
end

% Get masks
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

% Mask for correlating pixels
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
end