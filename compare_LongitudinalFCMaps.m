%% This code compares how FC maps are changed over time 

clc; clear;
commonDir = 'C:\Users\KEM294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\Imaging']));
addpath(genpath([commonDir '\Codes\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\fly_track\videoIO']));

%% Store/Retrieve imaging data for all experiments and all sessions...
monkeyName = 'Whiskey'; hemisphere = 'Left';
[allDates,allRuns, refDate, refDir,lensCombo, roiSize, ephysFileNameAll, serverPath,probeLabel,...
    chInCortexNotes, greenMapRef] = getMonkeyParams_Imaging_Ephys(monkeyName,commonDir, hemisphere);
dateVals = [1 3 4];
checkDates = ['04_25_2022';'08_08_2022'; allDates(dateVals,:)];

runNames00 = ['run00'; 'run00';'run00'; 'run00'; 'run00'];
spatialBin = 3;

% Get average RS map for 5 seeds
if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\avgCorrMapControls.mat'],'file')
    figure('units','normalized','outerposition',[0 0 1 1]);
    imagesc(greenMapRef);axis image off; colormap gray; hold on;
    for iSeed = 1:20

        title(strrep(['Select the location where the probe is located for ' monkeyName ' - seed :' num2str(iSeed)],'_','\_'));
        refSeed(iSeed,:,:) = ginput(1);
    end
    close gcf;

    for iSeed = 1:20
        [~,corrMapFinal(:,:,iSeed)] = getRSConnectivityMaps(squeeze(refSeed(iSeed,:,:))',monkeyName);
    end

    save(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\avgCorrMapControls.mat'],'corrMapFinal','refSeed');

else
    corrMapFinal = load(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\avgCorrMapControls.mat']);
    refSeed       = corrMapFinal.refSeed;
    corrMapFinal = corrMapFinal.corrMapFinal;
end

 %% Load the recording-wise data and compute FC maps for different seeds 
for iDate = 1:size(checkDates,1)
    clear expDate dataDir templateDir 
    expDate    = checkDates(iDate,:);

    % Load all the information
    runName    = runNames00(iDate,:);
    fileNum    = 0;

    % Get the directory and filenames

    dataDir  = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
    templateDir = ['X:\Data\' monkeyName '_SqM\Left Hemisphere\'];
    cDriveFlag  = 0; 

    if isempty(dir(dataDir))        
        if (exist(['C:\Users\KEM294\Documents\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ],'dir'))
            dataDir     = ['C:\Users\KEM294\Documents\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
            templateDir = ['C:\Users\KEM294\Documents\Data\' monkeyName '_SqM\Left Hemisphere\'];
            cDriveFlag  = 1; 
        else
            [~,~] = mkdir(dataDir);
        end
    end


    numFiles    = length(dir([dataDir '\Spatial Downsample SS3' '/*.mat'])); % loading the downsampled data only
    datName     = 'Data_RS_10Hz_SS3_';
    
    fileInfo    = dir(dataDir);

    % Get greens - Check if greens exist in the data folder
    if isempty(find(strcmp({fileInfo.name},['green' runName(end-1:end) '.png']), 1)) && isempty(find(strcmp({fileInfo.name},['green' runName(end-1:end) '.bmp']), 1))
        try
            copyfile([serverPath '\' expDate '\' runName '\green' runName(end-1:end) '.png'], dataDir);
        catch
            copyfile([serverPath '\' expDate '\' runName '\green' runName(end-1:end) '.bmp'], dataDir);
        end
    end

    % Check if greens are edited
    if exist([templateDir  expDate '\' runNames00(iDate,:) '\green0' runNames00(iDate,5) '_Edited.png'],'file')
        greenTemp = imread([dataDir '\green0' runNames00(iDate,5) '_Edited.png']);

    elseif exist([templateDir  expDate '\' runNames00(iDate,:) '\green0' runNames00(iDate,5) '_Edited.bmp'],'file')
        greenTemp = imread([dataDir '\green0' runNames00(iDate,5) '_Edited.bmp']);

    else
        error('Greens are not edited...');
    end

    greenIm{iDate,1}  = greenTemp(:,:,1);

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
        processedDat{iDate,1} = matfile([dataDir '\processedFrames.mat']); clear tempBandPass;

    else
        % Retrieve processed data
        disp(['Loading imaging data for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)]);
        processedDat{iDate,1} = matfile([dataDir '\processedFrames.mat']);
    end

 cd(commonDir);
 clear pDatTemp fcMap seedLoc seedSig 
 pDatTemp = processedDat{iDate,1}.tempBandPass;

 % Get 20 points that can be compared between recordings
 if ~exist([dataDir '\controlSeedLocs.mat'],'file') 
     figure('units','normalized','outerposition',[0 0 1 1]);
     imagesc(greenIm{iDate,1}); colormap gray; axis image off;
     for iSeed = 1:20
         title(strrep(['Pick a seed to get a FC map for ' monkeyName ' ' expDate ' Seed: ' num2str(iSeed)],'_','\_')); hold on;
         seedLoc(iSeed,:) = ginput(1);
     end
     seedLoc = fliplr(round(seedLoc));
     close gcf;

     % Get recording wise map for 5 seeds
     for iSeed = 1: size(seedLoc,1)
         clear seedSig
         seedSig = calculateSeedSignal(imresize(greenIm{iDate,1},1/spatialBin), clipMask,fliplr(round(seedLoc(iSeed,:)./spatialBin)),12,pDatTemp); % Get Gaussian weighted seed signal
         fcMap(:,:,iSeed) = plotCorrMap(seedSig,pDatTemp,0);
     end
     save([dataDir '\controlSeedLocs.mat'],'fcMap','seedLoc');

 else
     seedLoc = load([dataDir '\controlSeedLocs.mat']);
     fcMap    = seedLoc.fcMap;
     seedLoc  = seedLoc.seedLoc;
 end

% Nonlinear registration to master green
clear regisPoints
if exist([dataDir '\imRegisMaster.mat'],'file') || exist(['C:\Users\KEM294\Documents\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate  '\imageRegisPoints_RS.mat'],'file')

    if ~cDriveFlag
        regisPoints = load([dataDir '\imRegisMaster.mat'],'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp');
    else
        regisPoints = load(['C:\Users\KEM294\Documents\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate  '\imageRegisPoints_RS.mat'],...
            'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp');
    end

    OTemp = regisPoints.OTemp;
    spacingTemp = regisPoints.spacingTemp;

else

    clear movPointsTemp fixedPointsTemp OTemp spacingTemp

    [movPointsTemp,fixedPointsTemp] = cpselect(greenTemp,greenMapRef,'Wait',true);
    [OTemp,spacingTemp,~] = point_registration(size(greenMapRef),fliplr(fixedPointsTemp),fliplr(movPointsTemp));
   
    if ~cDriveFlag
        save([dataDir '\imRegisMaster.mat'],'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp');
    else
        save(['C:\Users\KEM294\Documents\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate  '\imageRegisPoints_RS.mat'],...
            'movPointsTemp','fixedPointsTemp','OTemp','spacingTemp');

    end

end

greenSize = size(greenMapRef);
tempMask  = zeros(size(greenMapRef));
tempMask(1:size(greenTemp,1),1:size(greenTemp,2)) = imresize(corrMask,[1082 1312]);
[imMaskTemp,tMaskTemp]   = bspline_transform(OTemp,tempMask,spacingTemp,3); % bicubic spline interpolation

for iSeed = 1:size(seedLoc,1)
    tempFCMap = NaN(size(greenMapRef));
    tempFCMap(1:size(greenTemp,1),1:size(greenTemp,2)) = imresize(squeeze(fcMap(:,:,iSeed)),[1082 1312]);

    [imFCMapTemp,~] = bspline_transform(OTemp,tempFCMap,spacingTemp,3);

    imMaskTempR = logical(reshape(imMaskTemp,[greenSize(1)*greenSize(2) 1]));
    imFCMapTempR = reshape(imFCMapTemp,[greenSize(1)*greenSize(2) 1]);
    imFCMapTempR(~imMaskTempR) = NaN;

    corrMapFinalR = reshape(squeeze(corrMapFinal(:,:,iSeed)),[greenSize(1)*greenSize(2) 1]);
    corrMapFinalR(~imMaskTempR) = NaN;

    corrFC_Avg_SessionWise(iSeed,iDate) = corr(imFCMapTempR,corrMapFinalR,'rows','complete');
end
end

seedClass = [1 1 2 1 1 2 2 2 2 2 1 1 2 1 2 1 2 1 2 1 ];
recSites = corrFC_Avg_SessionWise((seedClass==2),:);
nonRecSites = corrFC_Avg_SessionWise((seedClass==1),:);
figure; 
scatter(1:5,recSites,[],[0.9922 0.2392 0.7098],'filled'); hold on;
scatter(1:5,nonRecSites,[],[0.0431 0.8549 0.3176],'filled');
ylim([0 1]); box off;
x = reshape(repmat(1:5,10,1),[50 1]);

[f,gof] = fit(x,reshape(recSites,[50 1]),'poly1'); % Line fitting
c = coeffvalues(f); r2 = gof.rsquare;
xFit = linspace(1, 5, 1000);
yFit = c(1)*xFit + c(2);
plot(xFit,yFit,'-','Color',[0.9922 0.2392 0.7098],'LineWidth',1);
mdl = fitlm(x,reshape(recSites,[50 1]));
text(4, 0.95,['R^2 : ' num2str(gof.rsquare*100) '%'],'Color',[0.9922 0.2392 0.7098]);
text(4,0.9,['p-val: ' num2str(mdl.Coefficients.pValue(2))],'Color',[0.9922 0.2392 0.7098]);

[f,gof] = fit(x,reshape(nonRecSites,[50 1]),'poly1'); % Line fitting
c = coeffvalues(f); r2 = gof.rsquare;
xFit = linspace(1, 5, 1000);
yFit = c(1)*xFit + c(2);
 plot(xFit,yFit,'-','Color',[0.0431 0.8549 0.3176],'LineWidth',1);
 mdl = fitlm(x,reshape(nonRecSites,[50 1]));
text(4, 0.85,['R^2 : ' num2str(gof.rsquare*100) '%'],'Color',[0.0431 0.8549 0.3176]);
text(4,0.8,['p-val: ' num2str(mdl.Coefficients.pValue(2))],'Color',[0.0431 0.8549 0.3176]);

xticks(1:5); xlim([0.8 5.2])
xticklabels({'04\_25\_2022'; '08\_08\_2022';'08\_14\_2023';'12\_04\_2023';'02\_20\_2024'});
xlabel('Procedure dates'); ylabel('Correlation of session-wise FC to average FC');

