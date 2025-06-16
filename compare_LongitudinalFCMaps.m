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
addpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12\fly_track\videoIO']));

%% Get monkey data
monkeyName = 'Whiskey'; hemisphere = 'Left';

[allDates,allRuns, refDate, refDir,lensCombo, roiSize, ephysFileNameAll, serverPath,probeLabel,...
    chInCortexNotes, greenMapRef] = getMonkeyParams_Imaging_Ephys(monkeyName,commonDir, hemisphere);

dateVals   = [1 3 4];
checkDates = ['04_25_2022';'04_25_2022';'04_25_2022';...
    '05_09_2022';'05_09_2022';'05_09_2022';...
    '06_28_2022';'06_28_2022';'06_28_2022';...
    '08_08_2022'; allDates(dateVals,:)];

runNames00 = ['run00'; 'run01'; 'run02';...
    'run00';'run01';'run02';...
    'run00';'run01';'run02';...
    'run00';'run00';'run00'; 'run00'; 'run00'];
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
        [corrMapFinal(:,:,iSeed)] = getAvgFCMaps(seedLocNonLinear(iSeed,:)',monkeyName);

    end
        save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\avgCorrMapControls_grid.mat'],'corrMapFinal','seedLocNonLinear');

else % Load the saved average FC map and seed locations
    corrMapFinal = load(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\avgCorrMapControls_grid.mat']);
    seedLocNonLinear  = corrMapFinal.seedLocNonLinear;
    corrMapFinal = corrMapFinal.corrMapFinal;
end

% %% Correlate session-wise FC map within the baseline to the average map
% % Check if the session-wise correlations are already saved
% varInfo = who('-file',['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\avgCorrMapControls_grid.mat']);
% 
% if ~ismember('corrFCAvgWithin',varInfo)
%     % Get rs-ISOI data
%     clear allDatesRS allRunsRS
%     [allDatesRS,allRunsRS,~,~,refImageName, serverPath,greenMapRef] = getMonkeyParamsRS(monkeyName,commonDir,hemisphere);
% 
%     % Load the run-wise FC map for all sessions and correlate with average
%     tic;
%     for iSeed = 1: size(seedLocNonLinear,1)
%         disp(['Seed: ' num2str(iSeed)]);
% 
%         % Get run-wise map for all sessions
%         [~,corrMaps] = getAvgFCMaps(seedLocNonLinear(iSeed,:)',monkeyName);
% 
%         % Reshape maps to the size of the master green
%         corrMaps   = imresize(corrMaps,[size(corrMapFinal,1) size(corrMapFinal,2)]);
%         corrMapBsl = reshape(squeeze(corrMapFinal(:,:,iSeed)),[size(corrMapFinal,1)*size(corrMapFinal,2) 1]); % Reshape average map
% 
%         % Correlate average FC with
%         for iD = 1: size(corrMaps,3)
%             for iR = 1:size(corrMaps,4)
%                 corrMapTemp = reshape(squeeze(corrMaps(:,:,iD,iR)),[size(corrMaps,1)*size(corrMaps,2) 1]);
%                 corrFCAvgWithin(iSeed,iD,iR) = corr(corrMapTemp,corrMapBsl);
%             end
%         end
%     end
% 
%     toc;
%     save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\avgCorrMapControls_grid.mat'],'corrFCAvgWithin','-append'); % Save the variable
% 
% else % Load the saved variable
%     corrFCAvgWithin = load(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\avgCorrMapControls_grid.mat'],'corrFCAvgWithin');
%     corrFCAvgWithin = corrFCAvgWithin.corrFCAvgWithin;
% end
% 
% %% Show the distributions of correlations
% figure;
% corrTempVals = reshape(permute(corrFCAvgWithin,[1 3 2 ]),[size(corrFCAvgWithin,1) 3]);
% boxplot(corrTempVals); hold on;
% numSeeds = size(corrTempVals,1);
% 
% figure; % Show the distribution
% x1 = (reshape(repmat(1:3,[numSeeds 1]),[numSeeds*5 1]));
% y1 = reshape(corrTempVals,[numSeeds*3 1]);
% s = swarmchart(x1,y1,20,[0 0.4470 0.7410]); 
% s.XJitterWidth = 0.5; hold on;


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

%% Function to get average rs-FC maps from 3 sessions instead of all sessions/.
% Modified from getRSConnectivityMaps.m to remove the last experiment from
% the baseline/average FC map.

function [avgCorrMapFinalR,corrMap] =  getAvgFCMaps(refSeed,monkeyName)
% Get the RS connectivity matrix before thresholding
% Initializing all the relevant variables
commonDir = 'C:\Users\kem294\Documents\Data';
spatialBin   = 3;
hemisphere   = 'Left';
[allDates,allRuns,~,refDir,refImageName, serverPath,greenMapRef] = getMonkeyParamsRS(monkeyName,commonDir,hemisphere);

% Get the necessary variables needed to compute the connectivity maps
clear seedSig
greenMapSize = size(greenMapRef);
greenMap     = double(greenMapRef);
greenMap_RGB = ind2rgb(greenMap,gray(256)); % Convert from grayscale to RGB

% Get the master green mask 
if exist([refDir refImageName '_mask.bmp'],'file') == 0 % Make sure to get what the reference run is outside the function
    greenMapMask = imread([refDir refImageName '_mask.png']);
else
    greenMapMask  = imread([refDir refImageName '_mask.bmp']);
end
greenMapMask = double(greenMapMask(:,:,1) == 255);

% Get the reference seed
refSeed = fliplr(round(refSeed));

% Initializing required variables
numCols = max(cellfun(@(x)size(x,1),allRuns));
tformWithinRuns = cell(size(allDates,1),numCols); 

% Getting maps for all runs individually
for iDate = 1: size(allDates,1)-1 % Remove last experiment from the average
    for iRun = 1:size(allRuns{iDate,1},1)

        % Get run details
        clear runName expDate serverDataPath dataDir numFiles templateDir datName tempBandPass
        runName       = allRuns{iDate,1}(iRun,:);
        expDate       = allDates(iDate,:);
        if strcmp(expDate,'06_22_2021')
            serverDataPath    = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\302-19_CharlieSheen_SqM\Left Hemisphere\06_22_2021 Craniotomy\' runName];
        else
            serverDataPath = [serverPath expDate '\' runName];
        end

        dataDir  = [commonDir '\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
        numFiles = length(dir([dataDir '\Spatial Downsample SS3' '/*.mat'])); % loading the downsampled data only
        datName  = 'Data_RS_10Hz_SS3_';
        templateDir = [commonDir '\' monkeyName '_SqM\Left Hemisphere\'];

                
        % Get the processed data or compute the processed data if not saved
        if ~exist([dataDir '\processedFrames.mat'],'file')
            [~,~,~,~,tempBandPass] = getPreProcessedDataRestingState(serverDataPath,dataDir,runName,numFiles,spatialBin,datName);
            save([dataDir '\processedFrames.mat'],'tempBandPass');
            processedDat = tempBandPass; 
        else
            clear processedDatTemp processedDat
            processedDatTemp = load([dataDir '\processedFrames.mat']);
            processedDatTemp = processedDatTemp.tempBandPass;
            processedDat     = processedDatTemp;
        end

        % Load the greens
        if exist([templateDir  expDate '\' allRuns{iDate,1}(iRun,:) '\green0' allRuns{iDate,1}(iRun,5) '_Edited.bmp'],'file') == 0
            greenTemp = imread([dataDir '\green0' allRuns{iDate,1}(iRun,5) '_Edited.png']);
        else
            greenTemp = imread([dataDir '\green0' allRuns{iDate,1}(iRun,5) '_Edited.bmp']);
        end
        greenIm{iDate,iRun}  = greenTemp(:,:,1);

        % Load the clipmask 
         if exist([dataDir '\clipMask0' runName(end) '.BMP'],'file') == 0
            clipMask = imread([dataDir '\clipMask0' runName(end) '.png']);
        else
            clipMask = imread([dataDir '\clipMask0' runName(end) '.bmp']);
         end
         clipMask = imresize(clipMask,1/3); % Resize mask
         clipMask = clipMask(:,:,1)>0; % Converting to 0s and 1s
        
        % Transformation of the seed point to the run-wise maps so as to
        % get individual run-wise correlation maps
        % Linear transformation
        clear tformRunWise
        if exist([dataDir '\tform_RS.mat'],'file') == 0
            [optimizer,metric] = imregconfig('multimodal');
            tformRunWise = imregtform(imresize(greenIm{iDate,iRun},1/spatialBin),imresize(greenIm{iDate,1},1/spatialBin),'affine',optimizer,metric);
            tformWithinRuns{iDate,iRun} = tformRunWise;
            save([dataDir '\tform_RS.mat'],'tformRunWise');
        else
            
            if isempty(tformWithinRuns{iDate,iRun}) 
                load([dataDir '\tform_RS.mat'],'tformRunWise');
                tformWithinRuns{iDate,iRun} = tformRunWise;
            end
        end
        
        % Non-linear transformation 
         if iRun == 1
            if (exist([templateDir expDate '\imageRegisPoints_RS.mat'],'file') == 0) 
                selectPoints = 1; addPoints = 0; checkRegisImage = 1;

                % Get the fixed and moving points and store them in a matfile
                clear temp
                if selectPoints
                    clear movePointsTemp fixedPointsTemp
                    [movPointsTemp,fixedPointsTemp] = cpselect(greenIm{iDate,1},greenMapRef,'Wait',true);
                    movPoints{iDate,1} = movPointsTemp;
                    fixedPoints{iDate,1} = fixedPointsTemp;
                end
                
                if addPoints
                    [movPointsTemp,fixedPointsTemp] = cpselect(greenIm{iDate,1},greenMapRef,movPoints{iDate,1},fixedPoints{iDate,1},'Wait',true);
                    movPoints{iDate,1} = movPointsTemp;
                    fixedPoints{iDate,1} = fixedPointsTemp;
                end

                % In case you want to save again especially when you want to add points or to debug
                if ~(selectPoints || addPoints) 
                    try
                        load([templateDir expDate '\imageRegisPoints_RS.mat']);
                        movPoints{iDate,1} = movPointsTemp;
                        fixedPoints{iDate,1} = fixedPointsTemp;
                    catch
                        disp('Image registration not done, beginning image registration now..');
                        clear movePointsTemp fixedPointsTemp
                        [movPointsTemp,fixedPointsTemp] = cpselect(greenIm{iDate,1},greenMapRef,'Wait',true);
                        movPoints{iDate,1} = movPointsTemp;
                        fixedPoints{iDate,1} = fixedPointsTemp;
                    end
                end

                clear temp OTemp spacingTemp tTemp
                fixedPointsTemp = fixedPoints{iDate,1}; movPointsTemp = movPoints{iDate,1};
                temp = zeros(size(greenMapRef));
                temp(1:size(greenIm{iDate,1},1),1:size(greenIm{iDate,1},2)) = greenIm{iDate,1};
                
                % Perform image registration
                [OTemp,spacingTemp,~] = point_registration(size(greenMapRef),fliplr(fixedPoints{iDate,1}),fliplr(movPoints{iDate,1}));
                [imRTemp{iDate,1},~] = bspline_transform(OTemp,temp,spacingTemp,3); % bicubic spline interpolation
                O{iDate,1} = OTemp;
                spacing{iDate,1} = spacingTemp;

                save([templateDir expDate '\imageRegisPoints_RS.mat'],'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp');
                
                if checkRegisImage % Check the green images
                    figure; imshowpair(greenMapRef,imRTemp{iDate,1},'blend'); %colormap gray;
                    title(['Green map from Date: ' allDates(iDate,:)]);
                end
            else % Load the registered image 
                load([templateDir expDate '\imageRegisPoints_RS.mat'],'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp');
                movPoints{iDate,1} = movPointsTemp;
                fixedPoints{iDate,1} = fixedPointsTemp;
                O{iDate,1} = OTemp;
                spacing{iDate,1} = spacingTemp;
            end
         end

         % Get connectivity maps
         % Get run-wise correlation maps
         Ox = O{iDate,1}(:,:,1); Oy = O{iDate,1}(:,:,2); dx = spacing{iDate,1}(1); dy = spacing{iDate,1}(2);
         seedLocNonLinear(iDate,:) = bspline_trans_points_double(O{iDate,1},spacing{iDate,1},[refSeed(1) refSeed(2)]); % Non-linear transformation
         [seedLocR(2),seedLocR(1)] = transformPointsInverse(tformWithinRuns{iDate,iRun},seedLocNonLinear(iDate,2),seedLocNonLinear(iDate,1)); % Linear transformation

         if sum(seedLocR<=0)>0 || sum(squeeze(round(seedLocR./spatialBin))<=1)>=1; noSeedSig = false; continue; end % if seed value is negative, incrementing the pixel values
       
         try
             seedSig = calculateSeedSignal(imresize(greenIm{iDate,iRun},1/spatialBin), clipMask,fliplr(squeeze(round(seedLocR./spatialBin))),3,processedDat); % Get Gaussian weighted seed signal
             noSeedSig = false;
         catch
             noSeedSig = true;
             continue;
         end

         corrMapTemp(:,:,iDate,iRun) = plotCorrMap(seedSig,processedDat,0); % Get FC map
        
        
         % Register the correlation map back to the master green 
         clear corrMapLinear corrMapNonLinear
         scaleRatio = round(greenMapSize(1)/size(processedDat,1),2);
         corrMapLinear = imresize(imwarp(corrMapTemp(:,:,iDate,iRun),tformWithinRuns{iDate,iRun},'OutputView',imref2d(size(corrMapTemp(:,:,iDate,iRun)))),spatialBin,'OutputSize',[greenMapSize(1) greenMapSize(2)]); % Linear transformation
         corrMapNonLinear = bspline_transform_2d_double(Ox,Oy,corrMapLinear,dx,dy,3); % Nonlinear transformation

         % Mutliply the correlation map with master green mask so as to
         % remove the corners outside the cortex
         corrMapNonLinear         = corrMapNonLinear.*greenMapMask;
         corrMap(:,:,iDate,iRun)  = imresize(corrMapNonLinear,1/scaleRatio,'OutputSize',[size(processedDat,1) size(processedDat,2)]); % Resize the image

         % % Keep a copy and mask it for correlating with the average
         % corrMapMask = reshape(corrMapTemp(:,:,iDate,iRun),[size(corrMapTemp,1)*size(corrMapTemp,2) 1]);
         % clipMaskR   = reshape(clipMask,[size(corrMapTemp,1)*size(corrMapTemp,2) 1]);
         % 
         % corrMapMask(~clipMaskR,:) = NaN; corrMapMask = reshape(corrMapMask,[size(corrMapTemp,1) size(corrMapTemp,2)]);
         % 
         % corrMaskLinear = imresize(imwarp(corrMapMask(:,:,iDate,iRun),tformWithinRuns{iDate,iRun},'OutputView',imref2d(size(corrMapTemp(:,:,iDate,iRun)))),spatialBin,'OutputSize',[greenMapSize(1) greenMapSize(2)]); % Linear transformation
         % corrMaskNonLinear = (bspline_transform_2d_double(Ox,Oy,corrMaskLinear,dx,dy,3)).*greenMapMask; % Nonlinear transformation
         % corrMasked(:,:,iDate,iRun)  = imresize(corrMaskNonLinear,1/scaleRatio,'OutputSize',[size(processedDat,1) size(processedDat,2)]); % Resize the image

         
    end
    if noSeedSig; continue; end
end

% Reshape the correlation map
corrMapR = reshape(corrMap,[size(corrMap,1) size(corrMap,2) size(corrMap,3)*size(corrMap,4)]);

% Check if any entity in the matrix is empty
checkMat = squeeze(sum(corrMapR,[1 2])== 0);
corrMapR(:,:,checkMat) = [];

% Average out the correlation maps
avgCorrMap = mean(corrMapR,3,'omitnan');

% Compare the run-wise correlation maps with the average correlation map
% Get the pearson correlation between them
imSize = size(corrMapR);
corrMapMat    = reshape(corrMapR,[imSize(1)*imSize(2) imSize(3)]);
corrMapMatAvg = reshape(avgCorrMap,[imSize(1)*imSize(2) 1]);

for iL = 1:size(corrMapMat,2)
    rho(iL,1) = corr(corrMapMatAvg,corrMapMat(:,iL)); % Get the correlation coefficient
end

% Compute the mean and std dev and reject the runs whose values are greater
% or lesser than 1.5 std dev
rhoM   = mean(rho,'omitnan');
rhoStd = std(rho,'omitnan');
lowThresh = rhoM - (1.5*rhoStd); % rhoM - (1.5*rhoStd);
highThresh = rhoM + (1.5*rhoStd); % rhoM + (1.5*rhoStd);
goodRuns = ((rho >= lowThresh) & (rho <= highThresh));
corrMapMat(:,~goodRuns) = []; % Remove bad runs 

% Get the average of good runs 
corrMapGood = reshape(corrMapMat,[imSize(1) imSize(2) size(corrMapMat,2)]);
rsConnMatrix = mean(corrMapGood,3,'omitnan'); avgCorrMapFinal = rsConnMatrix;
rsConnMatrix = imresize(rsConnMatrix,scaleRatio,'OutputSize',[greenMapSize(1) greenMapSize(2)]);

% Only include pixels which are significantly greater than 0
[~,p]     = ttest(corrMapMat');
reshapedP = reshape(p,[imSize(1) imSize(2)]);
sigPixels = (reshapedP<0.001) & (avgCorrMapFinal>0);
avgCorrMapFinal(~sigPixels) = 0; 

avgCorrMapFinalR = imresize(avgCorrMapFinal,scaleRatio,'OutputSize',[greenMapSize(1) greenMapSize(2)]);

end
