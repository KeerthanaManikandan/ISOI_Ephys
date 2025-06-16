% ephysImaging_compilation_vFinal - MASTER SCRIPT
% This script compiles data processed from ephys_Imaging_v3 for ONE MONKEY
% at a time and saves the processed data. 
% See ephysImaging_v3 for reference
% see ephysImaging_compilation.m,ephysImaging_compilation_v4 and
% ephysImaging_v5 for previous versions of this script

% Final MATLAB script for ISOI_Ephys paper

% June 12,2025 - KM
% Set paths
clc; clear;
commonDir = 'C:\Users\kem294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\ISOI_Ephys\neuroshare']));
addpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12\fly_track\videoIO']));

%% Initialize variables and get monkey data
hemisphere = 'Left'; spatialBin = 3;
iM = 1; % 1 - Charlie Sheen, 2 - Whiskey

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
[processedDat,greenIm,probe,badCh,badTimesLFP,badTimeThresh,estChInCortex] = ...
    getAllData_Imaging_Ephys(monkeyName,hemisphere,commonDir,serverPath,allDates,allRuns,...
    ephysFileNameAll,greenMapRef,chInCortexNotes,probeLabel,spatialBin);

clc; disp(['All physiology and imaging data for ' monkeyName ' loaded']);


%% Temporal correlations: Correlate rs-ISOI with LFP at ROI
% This section calculates cross-correlations between rs-ISOI and LFP
% (different frequency bands) within an ROI 

clear tempProfileNoRef tempProfileSuperNoRef tempProfileDeepNoRef tempProfileMidNoRef...
    tempProfile_6Ch_NoRef tempProfileSuper_6Ch_NoRef tempProfileDeep_6Ch_NoRef tempProfileMid_6Ch_NoRef

% Obtain cross-correlations for all recordings
for iDate = 1:size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:); % Get experiment dates 

    for iRun = 1:size(allRuns{iDate,1}) % Get individual runs within an experiment 
        clear runName dataDir clipMask elecMask clipMaskCortex corrMask x negIdx lowIdx
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

        clc; disp(['Analyzing data for ' monkeyName ' '  expDate ' run: ' runName]);

        if ~exist([dataDir '\ROIAllVars.mat'],'file')

            % IMAGING: Load the appropriate masks for the imaging data
            [clipMaskCortex, corrMask] = getMasks(dataDir,runName);

            % Cross-correlate rs-ISOI and LFP within the ROI
            % Pick the ROI (1mm x 1mm) - 1mm diameter/ 500um radius;
            disp('Getting temporal profile...');
            seedRad  = round(roiSize{iDate}(iRun)*2./spatialBin); % seed radius for FC map 
            imSize   = size(imresize(greenIm{iDate,iRun},1/spatialBin)); % image size
            greenFig = imresize(greenIm{iDate,iRun},1/spatialBin,'OutputSize',[imSize(1) imSize(2)]); 

            % Check if ROI center and electrode location are saved
            clear varInfo;
            if exist(fullfile(dataDir,'roiCenterLoc.mat'),'file')
                varInfo = who('-file', fullfile(dataDir,'roiCenterLoc.mat'));
            else
                varInfo =[];
            end
            if find(ismember(varInfo,'seedLocProbe')); probeLocFlag = 1; else; probeLocFlag = 0; end % check if electrode location are saved

            % Save ROI center and probe location if not saved already
            if ~exist(fullfile(dataDir,'roiCenterLoc.mat'),'file') || ~probeLocFlag

                % Check if the ROI center and electrode location are saved
                if ~exist(fullfile(dataDir,'roiCenterLoc.mat'),'file') 
                    disp('Picking ROI near the electrode...');
                    figure; imagesc(greenFig); colormap gray; axis image off;
                    title('Pick a seed to get the ROI');
                    seedLocIn = ginput(1); seedLocIn = (round(seedLocIn)); close gcf;
                else
                    seedLocIn = load([dataDir '\roiCenterLoc.mat'],'seedLocIn');
                    seedLocIn = seedLocIn.seedLocIn;
                end

                if ~probeLocFlag % Check if electrode location has been saved
                    figure; imagesc(greenFig); colormap gray; axis image off; hold on;
                    title ('Pick the location of the probe...')
                    seedLocProbe = ginput(1); seedLocProbe = (round(seedLocProbe)); close gcf;
                end

                save([dataDir '\roiCenterLoc.mat'],'seedLocProbe','seedLocIn'); % Save electrode and ROI center location

            else
                disp('Seed for ROI already picked...');

                seedLocIn = load([dataDir '\roiCenterLoc.mat'],'seedLocIn');
                seedLocIn = seedLocIn.seedLocIn;

                seedLocProbe = load([dataDir '\roiCenterLoc.mat'],'seedLocProbe');
                seedLocProbe = seedLocProbe.seedLocProbe;
            end

            % Show the ROI on the blood vessel map
            if ~exist(fullfile(dataDir,'ROI.png'),'file')
                figure; imagesc(greenFig); axis image off; colormap gray;
                rectangle('Position',[seedLocIn(1)-round(seedRad/2),seedLocIn(2)-round(seedRad/2),...
                    seedRad,seedRad],'EdgeColor','r','LineWidth',2);
                title('ROI near the electrode');
                f = gcf; exportgraphics(f,[dataDir '\ROI.png'],'Resolution',300); close gcf;
            end
            
            % Generate a recording-wise FC map
            if ~exist(fullfile(dataDir,'FCMap_ROI.png'),'file')
                clear pDatTemp;
                circleRad = round(roiSize{1}(1)/(spatialBin)); % seed radius
                pDatTemp  = processedDat{iDate,iRun}.tempBandPass; % rs-ISOI data

                seedSigT  = calculateSeedSignal(greenFig,corrMask,...
                    seedLocIn,circleRad,pDatTemp); % Get Gaussian weighted seed signal

                corrMapT   = plotCorrMap(seedSigT,pDatTemp,0); % Generate FC map
                
                % Visualize FC map 
                greenImRGB = ind2rgb(greenFig,gray(256));
                figure; imagesc(greenImRGB);axis image; colormap jet; axis image off;
                hold on; imagesc(corrMapT,'AlphaData',corrMapT.*0.8);clim([0 1]);colorbar;
                f = gcf; exportgraphics(f,[dataDir '\FCMap_ROI.png'],'Resolution',300); close gcf;
            end

            % Get the ROI for gamma cross correlations
            clipMaskROI = clipMaskCortex(seedLocIn(2)-round(seedRad/2):seedLocIn(2)+round(seedRad/2),...
                seedLocIn(1)-round(seedRad/2):seedLocIn(1)+round(seedRad/2));

            clear tempProfileNoRefRun tempProfileSuperNoRefRun tempProfileDeepNoRefRun tempProfileMidNoRefRun...
                tempProfile_6Ch_NoRefRun tempProfileSuper_6Ch_NoRefRun tempProfileDeep_6Ch_NoRefRun tempProfileMid_6Ch_NoRefRun

            % Perform cross-correlations for different re-referencing
            % schemes

            [tempProfileNoRefRun,tempProfileSuperNoRefRun,... % 10-channel split No-Ref
                tempProfileDeepNoRefRun,tempProfileMidNoRefRun,~] = ...
                getCrossCorrROI(dataDir,monkeyName,expDate,runName,...
                processedDat{iDate,iRun}.tempBandPass,probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,...
                badTimesLFP{iDate,iRun},badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),seedLocIn,...
                seedRad,clipMaskROI,10,'NoRef');

            [tempProfile_6Ch_NoRefRun,tempProfileSuper_6Ch_NoRefRun,... % 6-channel split No-Ref
                tempProfileDeep_6Ch_NoRefRun,tempProfileMid_6Ch_NoRefRun,~] = ...
                getCrossCorrROI(dataDir,monkeyName,expDate,runName,...
                processedDat{iDate,iRun}.tempBandPass,probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,...
                badTimesLFP{iDate,iRun},badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),seedLocIn,...
                seedRad,clipMaskROI,6,'NoRef');

            % Save the variables 
            save([dataDir '\ROIAllVars.mat'],'tempProfileNoRefRun','tempProfileSuperNoRefRun',...
                'tempProfileDeepNoRefRun','tempProfile_6Ch_NoRefRun','tempProfileDeep_6Ch_NoRefRun',...
                'tempProfileSuper_6Ch_NoRefRun','tempProfileMidNoRefRun','tempProfileMid_6Ch_NoRefRun');
            
            % No re-referencing scheme. splitting channels into two
            % compartments. Middle compartment will be a matrix of NaNs
            tempProfileNoRef(iDate,iRun)           = tempProfileNoRefRun; % All channels 
            tempProfileSuperNoRef(iDate,iRun)      = tempProfileSuperNoRefRun; % Superficial compartment
            tempProfileMidNoRef(iDate,iRun)        = tempProfileMidNoRefRun; % Middle compartment
            tempProfileDeepNoRef(iDate,iRun)       = tempProfileDeepNoRefRun; % Deep compartment

            % No re-referencing scheme,splitting channels into three
            % compartments
            tempProfile_6Ch_NoRef(iDate,iRun)      = tempProfile_6Ch_NoRefRun; % All channels
            tempProfileSuper_6Ch_NoRef(iDate,iRun) = tempProfileSuper_6Ch_NoRefRun; % Superficial
            tempProfileMid_6Ch_NoRef(iDate,iRun)   = tempProfileMid_6Ch_NoRefRun; % Middle compartment
            tempProfileDeep_6Ch_NoRef(iDate,iRun)  = tempProfileDeep_6Ch_NoRefRun; % Deep compartment

        else % Retrieve the temporal correlations for different re-referencing schemes
            clear allVars
            allVars                                 = load([dataDir '\ROIAllVars.mat']);
            tempProfileNoRef(iDate,iRun)            = allVars.tempProfileNoRefRun;
            tempProfileSuperNoRef(iDate,iRun)       = allVars.tempProfileSuperNoRefRun;
            tempProfileMidNoRef(iDate,iRun)         = allVars.tempProfileMidNoRefRun;
            tempProfileDeepNoRef(iDate,iRun)        = allVars.tempProfileDeepNoRefRun;
            tempProfile_6Ch_NoRef(iDate,iRun)       = allVars.tempProfile_6Ch_NoRefRun;
            tempProfileSuper_6Ch_NoRef(iDate,iRun)  = allVars.tempProfileSuper_6Ch_NoRefRun;
            tempProfileMid_6Ch_NoRef(iDate,iRun)    = allVars.tempProfileMid_6Ch_NoRefRun;
            tempProfileDeep_6Ch_NoRef(iDate,iRun)   = allVars.tempProfileDeep_6Ch_NoRefRun;          
        end
    end
end

% Compile temporal correlations, remove bad runs, and save them
clear bandLabels tempProfilesAll medTempProfileAll allChCorr allChLagVal...
    allChCorrSuper allChCorrMid allChCorrDeep
bandLabels = {'Theta'; 'Alpha'; 'Beta'; 'Gamma'; 'Spiking'};

% Temporal profiles for all frequency bands
tempProfilesAll(:,:,1)  = [tempProfileNoRef.profileTheta]; 
tempProfilesAll(:,:,2)  = [tempProfileNoRef.profileAlpha]; 
tempProfilesAll(:,:,3)  = [tempProfileNoRef.profileBeta]; 
tempProfilesAll(:,:,4)  = [tempProfileNoRef.profile];      
tempProfilesAll(:,:,5)  = [tempProfileNoRef.profileRaw];   

tempProfilesAll(:,~goodRuns,:) = [];

% Get the distribution of lags and peak negative correlations 
allChCorr(:,1) = [tempProfileNoRef.magLowTheta]'; 
allChCorr(:,2) = [tempProfileNoRef.magLowAlpha]'; 
allChCorr(:,3) = [tempProfileNoRef.magLowBeta]'; 
allChCorr(:,4) = [tempProfileNoRef.magLow]'; 
allChCorr(:,5) = [tempProfileNoRef.magLowRaw]'; 
allChCorr(~goodRuns,:) = []; 

allChLagVal(:,1)= [tempProfileNoRef.lagLowTheta]'; 
allChLagVal(:,2)= [tempProfileNoRef.lagLowAlpha]'; 
allChLagVal(:,3)= [tempProfileNoRef.lagLowBeta]'; 
allChLagVal(:,4) = [tempProfileNoRef.lagLow]'; 
allChLagVal(:,5) = [tempProfileNoRef.lagLowRaw]'; 
allChLagVal(~goodRuns,:) = []; 

% Get the distribution of lags and peak negative correlations - Superficial
allChCorrSuper(:,1) = [tempProfileSuper_6Ch_NoRef.magLowTheta]'; 
allChCorrSuper(:,2) = [tempProfileSuper_6Ch_NoRef.magLowAlpha]'; 
allChCorrSuper(:,3) = [tempProfileSuper_6Ch_NoRef.magLowBeta]'; 
allChCorrSuper(:,4) = [tempProfileSuper_6Ch_NoRef.magLow]'; 
allChCorrSuper(:,5) = [tempProfileSuper_6Ch_NoRef.magLowRaw]'; 
allChCorrSuper(~goodRuns,:) = []; 

% Get the distribution of lags and peak negative correlations - Middle
allChCorrMid(:,1) = [tempProfileMid_6Ch_NoRef.magLowTheta]'; 
allChCorrMid(:,2) = [tempProfileMid_6Ch_NoRef.magLowAlpha]'; 
allChCorrMid(:,3) = [tempProfileMid_6Ch_NoRef.magLowBeta]'; 
allChCorrMid(:,4) = [tempProfileMid_6Ch_NoRef.magLow]'; 
allChCorrMid(:,5) = [tempProfileMid_6Ch_NoRef.magLowRaw]'; 
allChCorrMid(~goodRuns,:) = []; 

% Get the distribution of lags and peak negative correlations - Deep
allChCorrDeep(:,1) = [tempProfileDeep_6Ch_NoRef.magLowTheta]'; 
allChCorrDeep(:,2) = [tempProfileDeep_6Ch_NoRef.magLowAlpha]'; 
allChCorrDeep(:,3) = [tempProfileDeep_6Ch_NoRef.magLowBeta]'; 
allChCorrDeep(:,4) = [tempProfileDeep_6Ch_NoRef.magLow]'; 
allChCorrDeep(:,5) = [tempProfileDeep_6Ch_NoRef.magLowRaw]'; 
allChCorrDeep(~goodRuns,:) = []; 


% Get the distributions of lags and correlations for 10/10 and 6/6 split
% for superficial and deep channels
superCorr_10Ch_Gamma = [tempProfileSuperNoRef.magLow]'; 
deepCorr_10Ch_Gamma  = [tempProfileDeepNoRef.magLow]';
superCorr_6Ch_Gamma  = [tempProfileSuper_6Ch_NoRef.magLow]'; 
deepCorr_6Ch_Gamma   = [tempProfileDeep_6Ch_NoRef.magLow]'; 

superCorr_10Ch_Gamma(~(goodRuns & ~singleChFlag)) = []; 
deepCorr_10Ch_Gamma(~(goodRuns & ~singleChFlag))  = []; 
superCorr_6Ch_Gamma(~(goodRuns & ~singleChFlag))  = []; 
deepCorr_6Ch_Gamma(~(goodRuns & ~singleChFlag))  = [];  

superLag_10Ch_Gamma = [tempProfileSuperNoRef.lagLow]'./10; 
deepLag_10Ch_Gamma  = [tempProfileDeepNoRef.lagLow]'./10;
superLag_6Ch_Gamma  = [tempProfileSuper_6Ch_NoRef.lagLow]'./10; 
deepLag_6Ch_Gamma   = [tempProfileDeep_6Ch_NoRef.lagLow]'./10;

superLag_10Ch_Gamma(~(goodRuns & ~singleChFlag)) = []; 
deepLag_10Ch_Gamma(~(goodRuns & ~singleChFlag))  = []; 
superLag_6Ch_Gamma(~(goodRuns & ~singleChFlag))  = []; 
deepLag_6Ch_Gamma(~(goodRuns & ~singleChFlag))   = [];  

% Get the median +/- mad lag 
x = -200:200;
allChLag = [tempProfileNoRef.lagLow]'; allChLag(~(goodRuns & ~singleChFlag)) = []; 
lagFrameRange = find(x == round(median(allChLag)- mad(allChLag))) : find(x == round(median(allChLag)+ mad(allChLag)));
%169:189; % Median +/- MAD across two monkeys

smFlagFOV = smFlag; smFlagFOV(~(goodRunsSpatial & ~singleChFlag)) = []; % good run flag for FOV and ROI 
smFlagROI = smFlag; smFlagROI(~(goodRuns)) = [];

% Save all compiled variables
if ~exist(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat'],'file') 
    save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat'],...
        'smFlagFOV','tempProfilesAll','allChCorr','allChLagVal','superCorr_10Ch_Gamma','deepCorr_10Ch_Gamma',...
        'superCorr_6Ch_Gamma','deepCorr_6Ch_Gamma','superLag_10Ch_Gamma','deepLag_10Ch_Gamma',...
        'superLag_6Ch_Gamma','deepLag_6Ch_Gamma','lagFrameRange','allChLag','allChCorrSuper',...
        'allChCorrMid','allChCorrDeep','smFlagROI','-append');
end

%% Spatial correlations: Correlate cross-modal and FC map
% Pre-allocate variables variables
videoFlag      = 0; % To save videos of hybrids - change if you need to
corrFCHybrid   = NaN(size(processedDat,1),size(processedDat,2),5,5,401);
corrFCSuper    = NaN(size(processedDat,1),size(processedDat,2),5,5,401);
corrFCDeep     = NaN(size(processedDat,1),size(processedDat,2),5,5,401);
super_DeepCorr = NaN(size(processedDat,1),size(processedDat,2),5,5,401);
corrFCMid      = NaN(size(processedDat,1),size(processedDat,2),5,5,401);
super_MidCorr  = NaN(size(processedDat,1),size(processedDat,2),5,5,401);
mid_DeepCorr   = NaN(size(processedDat,1),size(processedDat,2),5,5,401);

super_DeepAvgFrames = NaN(size(processedDat,1),size(processedDat,2),5,5);
super_MidAvgFrames  = NaN(size(processedDat,1),size(processedDat,2),5,5);
deep_MidAvgFrames   = NaN(size(processedDat,1),size(processedDat,2),5,5);
peakNegValsAll      = NaN(size(processedDat,1),size(processedDat,2),5,5);
peakNegTimesAll     = NaN(size(processedDat,1),size(processedDat,2),5,5);

superHybridAllBands = NaN(size(processedDat,1),size(processedDat,2),5,5,5);
midHybridAllBands   = NaN(size(processedDat,1),size(processedDat,2),5,5,5);
deepHybridAllBands  = NaN(size(processedDat,1),size(processedDat,2),5,5,5);
bandNames           = {'Theta'; 'Alpha';'Beta';'Gamma';'Spiking'};

% Obtain spatial correlations for all recordings
tic;
for iDate = 1:size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);

    for iRun = 1: size(allRuns{iDate,1},1)       
        clear runName dataDir clipMask elecMask clipMaskCortex corrMask...
            x negIdx lowIdx serverDir
        runName   = allRuns{iDate,1}(iRun,:);
        dataDir   = ['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName];
        serverDir = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\kem294\Data\' ...
            monkeyName '_SqM\' hemisphere ' Hemisphere\'  expDate '\' runName ];
        
        if ~exist(serverDir,'dir') 
            [~,~] = mkdir(serverDir); 
        end 

        clc; disp(['Analyzing data for ' monkeyName ' '  expDate ' run: ' runName]);
        
        % Cross-correlate cross-modal map to FC map for different
        % re-referencing schemes

        % No re-referencing, 10 superficial, 10 deep channels
        if ~exist([serverDir '\crossCorrFOV_10_NoRef.mat'],'file')
            disp('No reference - channel split 10(superficial)/10(Deep)');
            getCrossCorrFOV(monkeyName,expDate,runName,serverDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,10,'NoRef');
        end

        % No re-referencing, 6 superficial, 6 deep, 6 middle channels
        if ~exist([serverDir '\crossCorrFOV_6_NoRef.mat'],'file') 
            disp('No reference - channel split top 6(superficial)/bottom 6(Deep)');
            getCrossCorrFOV(monkeyName,expDate,runName,serverDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,6,'NoRef');
        end

         % Bipolar re-referencing, 6 superficial, 6 deep, 6 middle channels
        if ~exist([serverDir '\crossCorrFOV_6_BipolarRef.mat'],'file')
            disp('Bipolar reference -  channel split top 6(superficial)/bottom 6(Deep)');
            getCrossCorrFOV(monkeyName,expDate,runName,serverDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,6,'BipolarRef');
        end

        % Average re-referencing, 6 superficial, 6 deep, 6 middle channels
        if ~exist([serverDir '\crossCorrFOV_6_AvgRef.mat'],'file')
            disp('Avg reference - channel split top 6(superficial)/bottom 6(Deep)');
            getCrossCorrFOV(monkeyName,expDate,runName,serverDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,6,'AvgRef');
        end

        % CSD re-referencing, 6 superficial, 6 deep, 6 middle channels
        if ~exist([serverDir '\crossCorrFOV_6_CSDRef.mat'],'file') 
            disp('CSD reference - channel split top 6(superficial)/bottom 6(Deep)');
            getCrossCorrFOV(monkeyName,expDate,runName,serverDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,6,'CSDRef');
        end
      
        % Check if spatial correlation values have been compiled and saved
        clear varInfo;
        varFlag = 0; 
        try matfile([dataDir '\processedHybridMapVars.mat']);
            varInfo = who('-file',[dataDir '\processedHybridMapVars.mat']);
            if sum(ismember(varInfo,'crossFreqCrossLayerHybridT'))==0
                varFlag = 1;
            end
        catch
            varFlag = 1;
        end 

        if varFlag
            % IMAGING: Load the appropriate masks for the imaging data
            % Load clipmask
            [clipMaskCortex, corrMask] = getMasks(dataDir,runName);          
            
            imSize    = size(clipMaskCortex);       
            corrMaskT = reshape(corrMask,[imSize(1)*imSize(2) 1]);

            % Get ROI location and FC map
            seedLocIn    = load([dataDir '\roiCenterLoc.mat']);
            seedLocProbe = seedLocIn.seedLocProbe;
            seedLocIn    = seedLocIn.seedLocIn;
            circleRad    = round(roiSize{iDate}(iRun)/(spatialBin)); % 500um radius
            greenFig     = imresize(greenIm{iDate,iRun},1/spatialBin,'OutputSize',[imSize(1) imSize(2)]);

            pDatTemp = processedDat{iDate,iRun}.tempBandPass; % Load rs-ISOI data
            seedSigT = calculateSeedSignal(greenFig,corrMask,...
                seedLocIn,circleRad,pDatTemp); % Get Gaussian weighted seed signal

            fcMap             = plotCorrMap(seedSigT,pDatTemp,0); % Get FC map 
            fcMap             = reshape(fcMap,[361*438 1]);
            fcMap(~corrMaskT) = NaN;

            chLen = estChInCortex{1,iDate}(iRun,2) - estChInCortex{1,iDate}(iRun,1); % Estimate # of channels in cortex

            % Initialize run-related variables
            % Correlations within the same frequencies            
            corrFCHybridT = NaN(5,4,401); corrFCSuperT    = NaN(5,4,401);
            corrFCDeepT   = NaN(5,4,401); super_DeepCorrT = NaN(5,4,401);

            corrFCMidT          = NaN(5,5,401); super_MidCorrT       = NaN(5,5,401);
            mid_DeepCorrT       = NaN(5,5,401); peakNegValsAllT      = NaN(5,5);
            peakNegTimesAllT    = NaN(5,5);     super_DeepAvgFramesT = NaN(5,5);
            super_MidAvgFramesT = NaN(5,5);     deep_MidAvgFramesT   = NaN(5,5);
           
            % Correlations between frequencies    
            superHybridAllBandsT = NaN(5,5,5);  deepHybridAllBandsT = NaN(5,5,5);
            midHybridAllBandsT   = NaN(5,5,5); crossFreqCrossLayerHybridT = NaN(5,5,5,3,3);

            for iType = 1:5 % 10/10 split or 6/6 split of superficial/deep channels for different re-referencing schemes
                clear crossCorrFOV allXCorr superXCorr deepXCorr allLags fileName
                switch iType
                    case 1 % No re-ref, 10/10 split between compartments
                        crossCorrFOV = matfile([serverDir '\crossCorrFOV_10_NoRef.mat']);
                        fileName     = '10_NoRef';

                    case 2 % No re-ref, 6/6/6 split between compartments
                        crossCorrFOV = matfile([serverDir '\crossCorrFOV_6_NoRef.mat']);
                        fileName     = '6_NoRef';

                    case 3 % Bipolar re-referencing, 6/6/6 split between compartments
                        crossCorrFOV = matfile([serverDir '\crossCorrFOV_6_BipolarRef.mat']);
                        fileName     = '6_BipolarRef';

                    case 4 % Average re-referencing, 6/6/6 split between compartments
                        crossCorrFOV = matfile([serverDir '\crossCorrFOV_6_AvgRef.mat']);
                        fileName     = '6_AvgRef';

                    case 5 % CSD re-referencing, 6/6/6 split between compartments
                        crossCorrFOV = matfile([serverDir '\crossCorrFOV_6_CSDRef.mat']);
                        fileName     = '6_CSDRef';
                end
                
                % Load the correlations
                allXcorr     = crossCorrFOV.spatialProfile;
                superXcorr   = crossCorrFOV.spatialProfileSuper;
                deepXcorr    = crossCorrFOV.spatialProfileDeep;

                if chLen~=0 && iType~=1 % Check if you have data for the middle compartments
                    midXCorr     = crossCorrFOV.spatialProfileMid;
                end
                allLags      = crossCorrFOV.lagFull;

                x = allLags; % Get lag information
                negIdx = x<0 & x>=-150; negVals = x(negIdx);
                lowIdx = x<0 & x>= -80; xLow = x(lowIdx);

                % Cross-modal maps for different frequencies

                % Theta
                crossCorrTheta      = reshape(allXcorr.ccFullTheta,[401 imSize(1)*imSize(2)]);   crossCorrTheta(:,~corrMaskT)      = NaN;
                crossCorrSuperTheta = reshape(superXcorr.ccFullTheta,[401 imSize(1)*imSize(2)]); crossCorrSuperTheta(:,~corrMaskT) = NaN;
                crossCorrDeepTheta  = reshape(deepXcorr.ccFullTheta,[401 imSize(1)*imSize(2)]);  crossCorrDeepTheta(:,~corrMaskT)  = NaN;
                lagLowTheta         = tempProfileNoRef(iDate,iRun).lagLowTheta;
            
                % Alpha
                crossCorrAlpha      = reshape(allXcorr.ccFullAlpha,[401 imSize(1)*imSize(2)]);   crossCorrAlpha(:,~corrMaskT)      = NaN;
                crossCorrSuperAlpha = reshape(superXcorr.ccFullAlpha,[401 imSize(1)*imSize(2)]); crossCorrSuperAlpha(:,~corrMaskT) = NaN;
                crossCorrDeepAlpha  = reshape(deepXcorr.ccFullAlpha,[401 imSize(1)*imSize(2)]);  crossCorrDeepAlpha(:,~corrMaskT)  = NaN;
                lagLowAlpha         = tempProfileNoRef(iDate,iRun).lagLowAlpha;
              
                % Beta
                crossCorrBeta      = reshape(allXcorr.ccFullBeta,[401 imSize(1)*imSize(2)]);   crossCorrBeta(:,~corrMaskT)      = NaN;
                crossCorrSuperBeta = reshape(superXcorr.ccFullBeta,[401 imSize(1)*imSize(2)]); crossCorrSuperBeta(:,~corrMaskT) = NaN;
                crossCorrDeepBeta  = reshape(deepXcorr.ccFullBeta,[401 imSize(1)*imSize(2)]);  crossCorrDeepBeta(:,~corrMaskT)  = NaN;
                lagLowBeta         = tempProfileNoRef(iDate,iRun).lagLowBeta;
               
                % Gamma
                crossCorrGamma      = reshape(allXcorr.ccFull,[401 imSize(1)*imSize(2)]);   crossCorrGamma(:,~corrMaskT)      = NaN;
                crossCorrSuperGamma = reshape(superXcorr.ccFull,[401 imSize(1)*imSize(2)]); crossCorrSuperGamma(:,~corrMaskT) = NaN;
                crossCorrDeepGamma  = reshape(deepXcorr.ccFull,[401 imSize(1)*imSize(2)]);  crossCorrDeepGamma(:,~corrMaskT)  = NaN;
                lagLowGamma         = tempProfileNoRef(iDate,iRun).lagLow;
               
                % Spiking
                crossCorrSpiking      = reshape(allXcorr.ccFullRaw,[401 imSize(1)*imSize(2)]);   crossCorrSpiking(:,~corrMaskT)      = NaN;
                crossCorrSuperSpiking = reshape(superXcorr.ccFullRaw,[401 imSize(1)*imSize(2)]); crossCorrSuperSpiking(:,~corrMaskT) = NaN;
                crossCorrDeepSpiking  = reshape(deepXcorr.ccFullRaw,[401 imSize(1)*imSize(2)]);  crossCorrDeepSpiking(:,~corrMaskT)  = NaN;
                lagLowSpiking         = tempProfileNoRef(iDate,iRun).lagLowRaw;               

                if iType >= 2 && chLen~=0
                    crossCorrMidTheta = reshape(midXCorr.ccFullTheta,[401 imSize(1)*imSize(2)]); crossCorrMidTheta(:,~corrMaskT)   = NaN;
                    crossCorrMidAlpha = reshape(midXCorr.ccFullAlpha,[401 imSize(1)*imSize(2)]); crossCorrMidAlpha(:,~corrMaskT)   = NaN;
                    crossCorrMidBeta = reshape(midXCorr.ccFullBeta,[401 imSize(1)*imSize(2)]);   crossCorrMidBeta(:,~corrMaskT)    = NaN;
                    crossCorrMidGamma = reshape(midXCorr.ccFull,[401 imSize(1)*imSize(2)]);      crossCorrMidGamma(:,~corrMaskT)   = NaN;
                    crossCorrMidSpiking = reshape(midXCorr.ccFullRaw,[401 imSize(1)*imSize(2)]); crossCorrMidSpiking(:,~corrMaskT) = NaN;
                end               
                   
                % Correlate all cross-modal maps to FC map
                for iBand = 1:5
                    clear bandName mapsAll crossCorrSuperR crossCorrDeepR crossCorrMidR
                    switch iBand
                        case 1 % Theta
                            bandName        = 'Theta';
                            mapsAll         = crossCorrTheta;
                            crossCorrSuperR = crossCorrSuperTheta;
                            crossCorrDeepR  = crossCorrDeepTheta;
                            lagLow          = lagLowTheta;

                            if iType >= 2 && chLen~=0
                                crossCorrMidR   = crossCorrMidTheta;                                
                            end

                        case 2 % Alpha
                            bandName        = 'Alpha';
                            mapsAll         = crossCorrAlpha;
                            crossCorrSuperR = crossCorrSuperAlpha;
                            crossCorrDeepR  = crossCorrDeepAlpha;
                            lagLow          = lagLowAlpha;

                            if iType >= 2 && chLen~=0
                                crossCorrMidR   = crossCorrMidAlpha;
                            end

                        case 3 % Beta
                            bandName        = 'Beta';
                            mapsAll         = crossCorrBeta;
                            crossCorrSuperR = crossCorrSuperBeta;
                            crossCorrDeepR  = crossCorrDeepBeta;
                            lagLow          = lagLowBeta;

                            if iType >= 2 && chLen~=0
                                crossCorrMidR   = crossCorrMidBeta;
                            end

                        case 4 % Gamma
                            bandName        = 'Gamma';
                            mapsAll         = crossCorrGamma;
                            crossCorrSuperR = crossCorrSuperGamma;
                            crossCorrDeepR  = crossCorrDeepGamma;
                            lagLow          = lagLowGamma;

                            if iType >= 2 && chLen~=0
                                crossCorrMidR   = crossCorrMidGamma;
                            end

                        case 5 % Spiking
                            bandName        = 'Spiking';
                            mapsAll         = crossCorrSpiking;
                            crossCorrSuperR = crossCorrSuperSpiking;
                            crossCorrDeepR  = crossCorrDeepSpiking;
                            lagLow          = lagLowSpiking;

                            if iType >= 2 && chLen~=0
                                crossCorrMidR   = crossCorrMidSpiking;
                            end
                    end
                    
                    % Correlate cross-modal maps to FC map
                    for iMap = 1:size(crossCorrGamma,1)
                        corrFCHybridT(iBand,iType,iMap)   = corr(fcMap,mapsAll(iMap,:)','rows','complete');
                        corrFCSuperT(iBand,iType,iMap)    = corr(fcMap,crossCorrSuperR(iMap,:)','rows','complete');
                        corrFCDeepT(iBand,iType,iMap)     = corr(fcMap,crossCorrDeepR(iMap,:)','rows','complete');
                        super_DeepCorrT(iBand,iType,iMap) = corr(crossCorrSuperR(iMap,:)',crossCorrDeepR(iMap,:)','rows','complete');

                        if iType >= 2 && chLen~=0
                            corrFCMidT(iBand,iType,iMap)     = corr(fcMap,crossCorrMidR(iMap,:)','rows','complete');
                            super_MidCorrT(iBand,iType,iMap) = corr(crossCorrSuperR(iMap,:)',crossCorrMidR(iMap,:)','rows','complete');
                            mid_DeepCorrT(iBand,iType,iMap)  = corr(crossCorrMidR(iMap,:)',crossCorrDeepR(iMap,:)','rows','complete');
                        end
                    end

                    % Get the peak negative correlation value and the lags
                    % for all channel hybrids
                    [peakNegValsAllT(iBand,iType),peakNegTimesAllT(iBand,iType)] = min(squeeze(corrFCHybridT(iBand,iType,negIdx)));
                    peakNegTimesAllT(iBand,iType) = negVals(peakNegTimesAllT(iBand,iType))./10;
                    
                    % Average the super/deep/mid cross-modal maps for a lag
                    % range and correlate them for 10/10 and 6/6 split
                    clear superHybridAvg deepHybridAvg
                    superHybridAvg = mean(crossCorrSuperR(lagFrameRange,:),1,'omitnan');
                    deepHybridAvg  = mean(crossCorrDeepR(lagFrameRange,:),1,'omitnan');
                    super_DeepAvgFramesT(iBand,iType) = corr(superHybridAvg',deepHybridAvg','rows','complete');

                    if iType >= 2 && chLen~=0
                        midHybridAvg = mean(crossCorrMidR(lagFrameRange,:),1,'omitnan');
                        super_MidAvgFramesT(iBand,iType) = corr(superHybridAvg',midHybridAvg','rows','complete');
                        deep_MidAvgFramesT(iBand,iType)  = corr(deepHybridAvg',midHybridAvg','rows','complete');
                    end
                    
                    % Show and save cross-modal maps
                    % All channels
                    if ~exist([dataDir '\HybridMapFOV_' fileName '_' bandName '.png'],'file') 
                        clear frameLow grayImFull crossCorr
                        frameLow       = (x == lagLow);
                        crossCorr      = reshape(mapsAll,[401 imSize(1) imSize(2)]);
                        figure; imagesc(squeeze(crossCorr(frameLow,:,:))); hold on; axis image off;
                        colormap(flipud(jet)); clim([-0.5 0.5]);
                        grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));
                        h = imagesc(grayImFull); hold off; set(h,'AlphaData',~corrMask);
                        title(['Peak negative at ' num2str(lagLow/10) 's']);
                        f = gcf; exportgraphics(f,[dataDir '\HybridMapFOV_' fileName '_' bandName '.png'],'Resolution',300); close gcf;
                    end

                   % Superficial compartment 
                    if ~exist([dataDir '\HybridMapFOV_Super_' fileName '_' bandName '.png'],'file')  
                        clear frameLow grayImFull crossCorrSuper
                        frameLow       = (x == lagLow);
                        crossCorrSuper = reshape(crossCorrSuperR,[401 imSize(1) imSize(2)]);
                        
                        figure; imagesc(squeeze(crossCorrSuper(frameLow,:,:))); hold on; axis image off;
                        colormap(flipud(jet)); clim([-0.5 0.5]);
                        grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));
                        h = imagesc(grayImFull); hold off; set(h,'AlphaData',~corrMask);
                        
                        title(['Superficial channels-Peak negative at ' num2str(lagLow/10) 's']);
                        f = gcf; exportgraphics(f,[dataDir '\HybridMapFOV_Super_' fileName '_' bandName '.png'],'Resolution',300); close gcf;
                    end

                    % Deep compartment
                    if ~exist([dataDir '\HybridMapFOV_Deep_' fileName '_' bandName '.png'],'file')  
                        clear frameLow grayImFull crossCorrDeep
                        frameLow      = (x == lagLow);
                        crossCorrDeep = reshape(crossCorrDeepR,[401 imSize(1) imSize(2)]);
                        
                        figure; imagesc(squeeze(crossCorrDeep(frameLow,:,:))); hold on; axis image off;
                        colormap(flipud(jet)); clim([-0.5 0.5]);
                        grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));
                        h = imagesc(grayImFull); hold off; set(h,'AlphaData',~corrMask);
                        
                        title(['Deep channels-Peak negative at ' num2str(lagLow/10) 's']);
                        f = gcf; exportgraphics(f,[dataDir '\HybridMapFOV_Deep_' fileName '_' bandName '.png'],'Resolution',300); close gcf;
                    end

                    % Middle compartment
                    if (~exist([dataDir '\HybridMapFOV_Mid_' fileName '_' bandName '.png'],'file')) && iType >= 2 && chLen~=0
                        clear frameLow grayImFull crossCorrMid 
                        frameLow     = (x == lagLow);
                        crossCorrMid = reshape(crossCorrMidR,[401 imSize(1) imSize(2)]);

                        figure; imagesc(squeeze(crossCorrMid(frameLow,:,:))); hold on; axis image off;                       
                        colormap(flipud(jet)); clim([-0.5 0.5]);
                        grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));                        
                        h = imagesc(grayImFull); hold off; set(h,'AlphaData',~corrMask);

                        title(['Middle channels-Peak negative at ' num2str(lagLow/10) 's']);
                        f = gcf; exportgraphics(f,[dataDir '\HybridMapFOV_Mid_' fileName '_' bandName '.png'],'Resolution',300); close gcf;
                    end
                end

                % Correlate compartment maps between frequencies 
                for iBand1 = 1:5
                    clear super1 deep1 mid1
                    switch iBand1
                        case 1
                            super1 = mean(crossCorrSuperTheta(lagFrameRange,:),1,'omitnan');
                            deep1  = mean(crossCorrDeepTheta(lagFrameRange,:),1,'omitnan');
                            if iType >= 2 && chLen~=0
                                mid1 = mean(crossCorrMidTheta(lagFrameRange,:),1,'omitnan');
                            end
                        
                        case 2
                            super1 = mean(crossCorrSuperAlpha(lagFrameRange,:),1,'omitnan');
                            deep1  = mean(crossCorrDeepAlpha(lagFrameRange,:),1,'omitnan');
                            if iType >= 2 && chLen~=0
                                mid1 = mean(crossCorrMidAlpha(lagFrameRange,:),1,'omitnan');
                            end
                        
                        case 3
                            super1 = mean(crossCorrSuperBeta(lagFrameRange,:),1,'omitnan');
                            deep1  = mean(crossCorrDeepBeta(lagFrameRange,:),1,'omitnan');
                            if iType >= 2 && chLen~=0
                                mid1 = mean(crossCorrMidBeta(lagFrameRange,:),1,'omitnan');
                            end
                       
                        case 4
                            super1 = mean(crossCorrSuperGamma(lagFrameRange,:),1,'omitnan');
                            deep1  = mean(crossCorrDeepGamma(lagFrameRange,:),1,'omitnan');
                            if iType >= 2 && chLen~=0
                                mid1 = mean(crossCorrMidGamma(lagFrameRange,:),1,'omitnan');
                            end
                        
                        case 5
                            super1 = mean(crossCorrSuperSpiking(lagFrameRange,:),1,'omitnan');
                            deep1  = mean(crossCorrDeepSpiking(lagFrameRange,:),1,'omitnan');
                            if iType >= 2 && chLen~=0
                                mid1 = mean(crossCorrMidSpiking(lagFrameRange,:),1,'omitnan');
                            end
                    end 
                        
                    for iBand2 = 1:5
                        clear super2 deep2 mid2
                        switch iBand2
                            case 1
                                super2 = mean(crossCorrSuperTheta(lagFrameRange,:),1,'omitnan');
                                deep2  = mean(crossCorrDeepTheta(lagFrameRange,:),1,'omitnan');
                                if iType >= 2 && chLen~=0
                                    mid2 = mean(crossCorrMidTheta(lagFrameRange,:),1,'omitnan');
                                end

                            case 2
                                super2 = mean(crossCorrSuperAlpha(lagFrameRange,:),1,'omitnan');
                                deep2  = mean(crossCorrDeepAlpha(lagFrameRange,:),1,'omitnan');
                                if iType >= 2 && chLen~=0
                                    mid2 = mean(crossCorrMidAlpha(lagFrameRange,:),1,'omitnan');
                                end

                            case 3
                                super2 = mean(crossCorrSuperBeta(lagFrameRange,:),1,'omitnan');
                                deep2  = mean(crossCorrDeepBeta(lagFrameRange,:),1,'omitnan');
                                if iType >= 2 && chLen~=0
                                    mid2 = mean(crossCorrMidBeta(lagFrameRange,:),1,'omitnan');
                                end

                            case 4
                                super2 = mean(crossCorrSuperGamma(lagFrameRange,:),1,'omitnan');
                                deep2  = mean(crossCorrDeepGamma(lagFrameRange,:),1,'omitnan');
                                if iType >= 2 && chLen~=0
                                    mid2 = mean(crossCorrMidGamma(lagFrameRange,:),1,'omitnan');
                                end

                            case 5
                                super2 = mean(crossCorrSuperSpiking(lagFrameRange,:),1,'omitnan');
                                deep2  = mean(crossCorrDeepSpiking(lagFrameRange,:),1,'omitnan');
                                if iType >= 2 && chLen~=0
                                    mid2 = mean(crossCorrMidSpiking(lagFrameRange,:),1,'omitnan');
                                end
                        end
                        
                        % Correlate compartment maps between frequencies
                        superHybridAllBandsT(iType,iBand1,iBand2) = corr(super1',super2','rows','complete');
                        deepHybridAllBandsT(iType,iBand1,iBand2)  = corr(deep1',deep2','rows','complete');
                       
                        if iType>=2 && chLen~=0
                            midHybridAllBandsT(iType,iBand1,iBand2) = corr(mid1',mid2','rows','complete');
                        end

                        if iBand1~=iBand2
                            crossFreqCrossLayerHybridT(iType,iBand1,iBand2,1,1) = corr(super1',super2','rows','complete');
                            crossFreqCrossLayerHybridT(iType,iBand1,iBand2,1,3) = corr(super1',deep2','rows','complete');
                            crossFreqCrossLayerHybridT(iType,iBand1,iBand2,3,1) = corr(deep1',super2','rows','complete');
                            crossFreqCrossLayerHybridT(iType,iBand1,iBand2,3,3) = corr(deep1',deep2','rows','complete');

                            if  iType>=2 && chLen~=0
                                crossFreqCrossLayerHybridT(iType,iBand1,iBand2,1,2) = corr(super1',mid2','rows','complete');
                                crossFreqCrossLayerHybridT(iType,iBand1,iBand2,2,1) = corr(mid1',super2','rows','complete');
                                crossFreqCrossLayerHybridT(iType,iBand1,iBand2,2,2) = corr(mid1',mid2','rows','complete');
                                crossFreqCrossLayerHybridT(iType,iBand1,iBand2,2,3) = corr(mid1',deep2','rows','complete');
                                crossFreqCrossLayerHybridT(iType,iBand1,iBand2,3,2) = corr(deep1',mid2','rows','complete');
                            end
                        end                      
                    end
                end
            end

            % Save all correlations for a single recording
            save([dataDir '\processedHybridMapVars.mat'],'corrFCHybridT','corrFCSuperT',...
                'corrFCDeepT','super_DeepCorrT','corrFCMidT','super_MidCorrT','mid_DeepCorrT',...
                'peakNegValsAllT','peakNegTimesAllT','super_DeepAvgFramesT','super_MidAvgFramesT',...
                'deep_MidAvgFramesT','superHybridAllBandsT','deepHybridAllBandsT','midHybridAllBandsT',...
                'crossFreqCrossLayerHybridT','fcMap');
            
            % Compile all variables for post-processing
            corrFCHybrid(iDate,iRun,:,:,:)   = corrFCHybridT; % Correlations between cross-modal and FC maps 
            corrFCSuper(iDate,iRun,:,:,:)    = corrFCSuperT;
            corrFCDeep(iDate,iRun,:,:,:)     = corrFCDeepT;
            corrFCMid(iDate,iRun,:,:,:)      = corrFCMidT;
           
            super_DeepCorr(iDate,iRun,:,:,:) = super_DeepCorrT; % Correlations between compartments           
            super_MidCorr(iDate,iRun,:,:,:)  = super_MidCorrT;
            mid_DeepCorr(iDate,iRun,:,:,:)   = mid_DeepCorrT ;

            peakNegValsAll(iDate,iRun,:,:)  = peakNegValsAllT; % Peak correlations and lags
            peakNegTimesAll(iDate,iRun,:,:) = peakNegTimesAllT;

            super_DeepAvgFrames(iDate,iRun,:,:) = super_DeepAvgFramesT; % Avg frames for correlations between compartments
            super_MidAvgFrames(iDate,iRun,:,:)  = super_MidAvgFramesT;
            deep_MidAvgFrames(iDate,iRun,:,:)   = deep_MidAvgFramesT;

            superHybridAllBands(iDate,iRun,:,:,:) = superHybridAllBandsT; % Correlations between frequencies 
            midHybridAllBands(iDate,iRun,:,:,:)   = midHybridAllBandsT;
            deepHybridAllBands(iDate,iRun,:,:,:)  = deepHybridAllBandsT;

            crossFreqCrossLayerHybrid(iDate,iRun,:,:,:,:,:) = crossFreqCrossLayerHybridT; 
        
        else % Load saved variables 
            clear allHybridVars
            allHybridVars                    = matfile([dataDir '\processedHybridMapVars.mat']);
            corrFCHybrid(iDate,iRun,:,:,:)   = allHybridVars.corrFCHybridT;% Correlations between cross-modal and FC maps 
            corrFCSuper(iDate,iRun,:,:,:)    = allHybridVars.corrFCSuperT;
            corrFCDeep(iDate,iRun,:,:,:)     = allHybridVars.corrFCDeepT;
            super_DeepCorr(iDate,iRun,:,:,:) = allHybridVars.super_DeepCorrT;

            corrFCMid(iDate,iRun,:,:,:)     = allHybridVars.corrFCMidT;
            super_MidCorr(iDate,iRun,:,:,:) = allHybridVars.super_MidCorrT;
            mid_DeepCorr(iDate,iRun,:,:,:)  = allHybridVars.mid_DeepCorrT ;

            peakNegValsAll(iDate,iRun,:,:)  = allHybridVars.peakNegValsAllT;% Peak correlations and lags
            peakNegTimesAll(iDate,iRun,:,:) = allHybridVars.peakNegTimesAllT;

            super_DeepAvgFrames(iDate,iRun,:,:) = allHybridVars.super_DeepAvgFramesT; % Avg frames for correlations between compartments
            super_MidAvgFrames(iDate,iRun,:,:)  = allHybridVars.super_MidAvgFramesT;
            deep_MidAvgFrames(iDate,iRun,:,:)   = allHybridVars.deep_MidAvgFramesT;

            superHybridAllBands(iDate,iRun,:,:,:) = allHybridVars.superHybridAllBandsT;% Correlations between frequencies
            midHybridAllBands(iDate,iRun,:,:,:)   = allHybridVars.midHybridAllBandsT;
            deepHybridAllBands(iDate,iRun,:,:,:)  = allHybridVars.deepHybridAllBandsT;
           
            crossFreqCrossLayerHybrid(iDate,iRun,:,:,:,:,:) = allHybridVars.crossFreqCrossLayerHybridT;

        end
    end
end
toc;

% Compile all spatial correlation variables

% Correlations between cross-modal and FC maps 
corrFCHybridT = reshape(corrFCHybrid,[size(corrFCHybrid,1)*size(corrFCHybrid,2) size(corrFCHybrid,3) size(corrFCHybrid,4) size(corrFCHybrid,5)]);
nanRow        = isnan(corrFCHybridT(:,1,1,1));
corrFCHybridT(nanRow,:,:,:) = [];
corrFCHybridT(~goodRunsSpatial,:,:,:) = []; 

% Correlations between compartments and FC map
corrFCSuperT = mean(corrFCSuper(:,:,:,:,lagFrameRange),5,'omitnan'); % Superficial
corrFCSuperT = reshape(corrFCSuperT,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]); 
corrFCSuperT(nanRow,:,:) = []; 
corrFCSuperT(~(goodRunsSpatial & ~singleChFlag),:,:) = []; 

corrFCMidT = mean(corrFCMid(:,:,:,:,lagFrameRange),5,'omitnan'); % Middle
corrFCMidT = reshape(corrFCMidT,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]); 
corrFCMidT(nanRow,:,:) = []; 
corrFCMidT(~(goodRunsSpatial & ~singleChFlag),:,:) = []; 

corrFCDeepT = mean(corrFCDeep(:,:,:,:,lagFrameRange),5,'omitnan'); % Deep
corrFCDeepT = reshape(corrFCDeepT,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]); 
corrFCDeepT(nanRow,:,:) = []; 
corrFCDeepT(~(goodRunsSpatial & ~singleChFlag),:,:) = []; 

% Peak negative correlations and lags for good runs
peakNegValsAllT = reshape(peakNegValsAll,[size(peakNegValsAll,1)*size(peakNegValsAll,2) size(peakNegValsAll,3) size(peakNegValsAll,4)]);
peakNegValsAllT(nanRow,:,:) = [];
peakNegValsAllT(~goodRunsSpatial,:,:) = [];

% Correlations between compartments ie; S/D, S/M and M/D 
super_DeepAvgFramesT = reshape(super_DeepAvgFrames,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]);
super_DeepAvgFramesT(nanRow,:,:)           = [];
super_DeepAvgFramesT(~(goodRunsSpatial & ~singleChFlag),:,:) = [];

super_MidAvgFramesT = reshape(super_MidAvgFrames,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]);
super_MidAvgFramesT(nanRow,:,:)           = [];
super_MidAvgFramesT(~(goodRunsSpatial & ~singleChFlag),:,:) = [];

deep_MidAvgFramesT = reshape(deep_MidAvgFrames,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]);
deep_MidAvgFramesT(nanRow,:,:)           = [];
deep_MidAvgFramesT(~(goodRunsSpatial & ~singleChFlag),:,:) = [];

% Cross-frequency correlations between hybrid maps for each layer 
superHybridAllBandsT = reshape(superHybridAllBands,[size(superHybridAllBands,1)*size(superHybridAllBands,2) size(superHybridAllBands,3) size(superHybridAllBands,4) size(superHybridAllBands,5)]);
nanRow = isnan(superHybridAllBandsT(:,2,1,1));
superHybridAllBandsT(nanRow,:,:,:) = [];
superHybridAllBandsT(~(goodRunsSpatial & ~singleChFlag),:,:,:)   = [];

midHybridAllBandsT = reshape(midHybridAllBands,[size(superHybridAllBands,1)*size(superHybridAllBands,2) size(superHybridAllBands,3) size(superHybridAllBands,4) size(superHybridAllBands,5)]);
midHybridAllBandsT(nanRow,:,:,:) = [];
midHybridAllBandsT(~(goodRunsSpatial & ~singleChFlag),:,:,:)   = [];

deepHybridAllBandsT = reshape(deepHybridAllBands,[size(superHybridAllBands,1)*size(superHybridAllBands,2) size(superHybridAllBands,3) size(superHybridAllBands,4) size(superHybridAllBands,5)]);
deepHybridAllBandsT(nanRow,:,:,:) = [];
deepHybridAllBandsT(~(goodRunsSpatial & ~singleChFlag),:,:,:)   = [];

% Check if the complied data has been saved
clear varInfo;
varFlag = 0;
try matfile(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat']);
    varInfo = who('-file',['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat']);
    if sum(ismember(varInfo,'corrFCHybridT'))==0 || sum(ismember(varInfo,'superHybridAllBandsT'))==0
        varFlag = 1;
    end
catch
    varFlag = 1;
end
        
if varFlag % Save variables
    save(['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat'],...
       'corrFCHybridT','corrFCSuperT','corrFCMidT','corrFCDeepT','peakNegValsAllT',...
       'super_DeepAvgFramesT','super_MidAvgFramesT','deep_MidAvgFramesT',...
       'superHybridAllBandsT','midHybridAllBandsT','deepHybridAllBandsT','-append');
end

