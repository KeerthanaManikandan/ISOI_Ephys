% ephysImaging_compilation_v4
% This script obtains the data procesed from ephysImaging_v3 and compiles
% the following:
% 1. Cross correlations between physiology and imaging for ROI and FOV
% 2. Temporal and spatial controls.
% See ephysImaging_v3 for reference
% see ephysImaging_compilation.m and ephysImaging_compilation_v3 for
% previous versions of this script
% September 19, 2024 - KM
% Set paths
clc; clear;
commonDir = 'C:\Users\KEM294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\Imaging']));
addpath(genpath([commonDir '\Codes\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\fly_track\videoIO']));

%% Initialize variables and get monkey data 
hemisphere = 'Left'; spatialBin = 3;
iM = 2;
switch iM
    case 1
        monkeyName = 'CharlieSheen';
        goodRuns  = ([1 1 NaN NaN NaN;... % Date x run
            1 1 NaN NaN NaN; ...
            1 1 1 0 1; ...
            1 1 1 1 NaN]);

        singleChFlag  = ([0 0 NaN NaN NaN;... % Date x run
            0 0 NaN NaN NaN; ...
            0 0 0 0 0; ...
            0 0 0 0 NaN]);

        goodRunsSpatial  = ([1 1 NaN NaN NaN;... % Date x run
            1 1 NaN NaN NaN; ...
            1 1 1 0 1; ...
            1 1 1 1 NaN]);

        smFlag  = (['S' 'M' '#' '#' '#';... % Date x run
            'S' 'M' '#' '#' '#'; ...
            'S' 'S' 'M' 'M' 'M'; ...
            'S' 'M' 'M' 'M' '#']); 

    case 2
        monkeyName = 'Whiskey';
        goodRuns   = ([1 1 1 NaN NaN NaN NaN; ... % Date x run
            1 1 1 1 1 1 1 ; ...
            1 0 1 1 NaN NaN NaN; ...
            1 1 0 1 1 1 0; ...
            1 0 0 1 1 1 1]);

        singleChFlag = ([1 1 1 NaN NaN NaN NaN; ... % Date x run
            0 0 0 0 0 0 0 ; ...
            0 0 0 0 NaN NaN NaN; ...
            0 0 0 0 0 0 0; ...
            0 0 0 0 0 0 0]);

        goodRunsSpatial = ([1 1 1 NaN NaN NaN NaN; ... % Date x run
            1 1 1 1 1 1 1 ; ...
            1 0 1 1 NaN NaN NaN; ...
            1 1 0 1 1 1 1; ...
            1 0 0 1 1 1 0]);

         smFlag = (['M' 'M' 'M' '#' '#' '#' '#'; ... % Date x run; # - in place of NaNs
            'M' 'M' 'M' 'M' 'S' 'S' 'S' ; ...
            'S' 'M' 'M' 'S' '#' '#' '#'; ...
            'S' 'M' 'M' 'S' 'M' 'M' 'M'; ...
            'S' 'S' 'M' 'M' 'M' 'S' 'M']);
end
goodRuns = reshape(goodRuns,[size(goodRuns,1)*size(goodRuns,2) 1]);
goodRuns(isnan(goodRuns)) = []; goodRuns = logical(goodRuns);

singleChFlag = reshape(singleChFlag,[size(singleChFlag,1)*size(singleChFlag,2) 1]);
singleChFlag(isnan(singleChFlag)) = []; singleChFlag = logical(singleChFlag); 

goodRunsSpatial = reshape(goodRunsSpatial,[size(goodRunsSpatial,1)*size(goodRunsSpatial,2) 1]);
goodRunsSpatial(isnan(goodRunsSpatial)) = []; goodRunsSpatial = logical(goodRunsSpatial);

smFlag = reshape(smFlag,[size(smFlag,1)*size(smFlag,2) 1]);
smFlag((smFlag == '#')) = [];
sensoryGoodRuns = (smFlag =='S') & goodRuns; 
sensoryGoodSpatialRuns = (smFlag == 'S') & goodRunsSpatial;
motorGoodRuns = (smFlag == 'M') & goodRuns; 
motorGoodSpatialRuns = (smFlag == 'M') & goodRuns; 

% Get monkey data and parameters
[allDates,allRuns, refDate, refDir,lensCombo, roiSize, ephysFileNameAll, serverPath,probeLabel,...
    chInCortexNotes, greenMapRef] = getMonkeyParams_Imaging_Ephys(monkeyName,commonDir, hemisphere);

% Get monkey data....
[processedDat,greenIm,probe,badCh,badTimesLFP,badTimeThresh,estChInCortex] = ...
    getAllData_Imaging_Ephys(monkeyName,hemisphere,commonDir,serverPath,allDates,allRuns,...
    ephysFileNameAll,greenMapRef,chInCortexNotes,probeLabel,spatialBin);

clc; disp(['All physiology and imaging data for ' monkeyName ' loaded']);
% 
% %% Move data to the server
% for iDate = 1:size(allDates,1)
%     clear expDate;
%     expDate = allDates(iDate,:);
%     for iRun = 1:size(allRuns{iDate,1})
%         clear runName dataDir clipMask elecMask clipMaskCortex corrMask x negIdx lowIdx
%         runName = allRuns{iDate,1}(iRun,:);
%         dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
%     end
% end


%% Cross correlating physiology and imaging at ROI 
% tempProfileNoRef           = NaN(iM,size(probe,2), size(probe,1));
% tempProfileSuperNoRef      = NaN(iM,size(probe,2), size(probe,1));
% tempProfileDeepNoRef       = NaN(iM,size(probe,2), size(probe,1));
% tempProfile_6Ch_NoRef      = NaN(iM,size(probe,2), size(probe,1));
% tempProfileSuper_6Ch_NoRef = NaN(iM,size(probe,2), size(probe,1));
% tempProfileDeep_6Ch_NoRef  = NaN(iM,size(probe,2), size(probe,1));

for iDate = 1:size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);
    for iRun = 1:size(allRuns{iDate,1})
        clear runName dataDir clipMask elecMask clipMaskCortex corrMask x negIdx lowIdx
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

        clc; disp(['Analyzing data for ' monkeyName ' '  expDate ' run: ' runName]);

        if ~exist([dataDir '\ROIAllVars.mat'],'file')
            % IMAGING: Load the appropriate masks for the imaging data
            % Load clipmask
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

            if corrMaskFlag
                corrMask = imresize(corrMask,1/3); % Resize the mask
                corrMask = corrMask(:,:,1)>0;
                corrMask = corrMask & ~elecMask;
            else
                corrMask = clipMaskCortex;
            end

            % Get the cross correlations for the ROI and FOV
            % 1. Get the temporal profile for all recordings...
            % Pick the ROI (1mm x 1mm) - 1mm diameter/ 500um radius;
            disp('Getting temporal profile...');
            seedRad = round(roiSize{iDate}(iRun)*2./spatialBin);
            imSize = size(imresize(greenIm{iDate,iRun},1/spatialBin));
            greenFig = imresize(greenIm{iDate,iRun},1/spatialBin,'OutputSize',[imSize(1) imSize(2)]);

            % Check if ROI center and location of probe are saved
            clear varInfo;
            if exist(fullfile(dataDir,'roiCenterLoc.mat'),'file')
                varInfo = who('-file', fullfile(dataDir,'roiCenterLoc.mat'));
            else
                varInfo =[];
            end
            if find(ismember(varInfo,'seedLocProbe')); probeLocFlag = 1; else; probeLocFlag = 0; end

            % Save ROI center and probe location if not saved already
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

                if ~probeLocFlag % Check if probe location has been saved
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

            % Show the ROI on the green map
            if ~exist(fullfile(dataDir,'ROI.png'),'file')
                figure; imagesc(greenFig); axis image off; colormap gray;
                rectangle('Position',[seedLocIn(1)-round(seedRad/2),seedLocIn(2)-round(seedRad/2),...
                    seedRad,seedRad],'EdgeColor','r','LineWidth',2);
                title('ROI near the electrode');
                f = gcf; exportgraphics(f,[dataDir '\ROI.png'],'Resolution',300); close gcf;
            end

            if ~exist(fullfile(dataDir,'FCMap_ROI.png'),'file')
                clear pDatTemp;
                circleRad = round(roiSize{1}(1)/(spatialBin));
                pDatTemp  = processedDat{iDate,iRun}.tempBandPass;
                seedSigT  = calculateSeedSignal(greenFig,corrMask,...
                    seedLocIn,circleRad,pDatTemp); % Get Gaussian weighted seed signal
                corrMapT   = plotCorrMap(seedSigT,pDatTemp,0);
                greenImRGB = ind2rgb(greenFig,gray(256));
                figure; imagesc(greenImRGB);axis image; colormap jet; axis image off;
                hold on; imagesc(corrMapT,'AlphaData',corrMapT.*0.8);caxis([0 1]);colorbar;
                f = gcf; exportgraphics(f,[dataDir '\FCMap_ROI.png'],'Resolution',300); close gcf;
            end

            % Get the ROI for gamma cross correlations
            clipMaskROI = clipMaskCortex(seedLocIn(2)-round(seedRad/2):seedLocIn(2)+round(seedRad/2),...
                seedLocIn(1)-round(seedRad/2):seedLocIn(1)+round(seedRad/2));

            clear tempProfileNoRefRun tempProfileSuperNoRefRun tempProfileDeepNoRefRun

            [tempProfileNoRefRun,tempProfileSuperNoRefRun,... % 10-channel split No-Ref
                tempProfileDeepNoRefRun,~] = ...
                getCrossCorrROI(dataDir,monkeyName,expDate,runName,...
                processedDat{iDate,iRun}.tempBandPass,probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,...
                badTimesLFP{iDate,iRun},badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),seedLocIn,...
                seedRad,clipMaskROI,10,'NoRef');

            [tempProfile_6Ch_NoRefRun,tempProfileSuper_6Ch_NoRefRun,... % 6-channel split No-Ref
                tempProfileDeep_6Ch_NoRefRun,~] = ...
                getCrossCorrROI(dataDir,monkeyName,expDate,runName,...
                processedDat{iDate,iRun}.tempBandPass,probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,...
                badTimesLFP{iDate,iRun},badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),seedLocIn,...
                seedRad,clipMaskROI,6,'NoRef');

            save([dataDir '\ROIAllVars.mat'],'tempProfileNoRefRun','tempProfileSuperNoRefRun',...
                'tempProfileDeepNoRefRun','tempProfile_6Ch_NoRefRun','tempProfileDeep_6Ch_NoRefRun',...
                'tempProfileSuper_6Ch_NoRefRun');

            tempProfileNoRef(iM,iDate,iRun)           = tempProfileNoRefRun; %#ok<*SAGROW> 
            tempProfileSuperNoRef(iM,iDate,iRun)      = tempProfileSuperNoRefRun;
            tempProfileDeepNoRef(iM,iDate,iRun)       = tempProfileDeepNoRefRun;
            tempProfile_6Ch_NoRef(iM,iDate,iRun)      = tempProfile_6Ch_NoRefRun;
            tempProfileSuper_6Ch_NoRef(iM,iDate,iRun) = tempProfileSuper_6Ch_NoRefRun;
            tempProfileDeep_6Ch_NoRef(iM,iDate,iRun)  = tempProfileDeep_6Ch_NoRefRun;

        else
            clear allVars
            allVars                                    = load([dataDir '\ROIAllVars.mat']);
            tempProfileNoRef(iM,iDate,iRun)            = allVars.tempProfileNoRefRun;
            tempProfileSuperNoRef(iM,iDate,iRun)       = allVars.tempProfileSuperNoRefRun;
            tempProfileDeepNoRef(iM,iDate,iRun)        = allVars.tempProfileDeepNoRefRun;
            tempProfile_6Ch_NoRef(iM,iDate,iRun)       = allVars.tempProfile_6Ch_NoRefRun;
            tempProfileSuper_6Ch_NoRef(iM,iDate,iRun)  = allVars.tempProfileSuper_6Ch_NoRefRun;
            tempProfileDeep_6Ch_NoRef(iM,iDate,iRun)   = allVars.tempProfileDeep_6Ch_NoRefRun;
        end


    end
end

%% Compile the ROI level variables and plot
bandLabels = {'Theta'; 'Alpha'; 'Beta'; 'Gamma'; 'Spiking'};

% Temporal profiles for all frequency bands
profileAlpha    = [tempProfileNoRef(iM,:,:).profileAlpha]; 
profileBeta     = [tempProfileNoRef(iM,:,:).profileBeta]; 
profileGamma    = [tempProfileNoRef(iM,:,:).profile];      
profileTheta    = [tempProfileNoRef(iM,:,:).profileTheta]; 
profileSpiking  = [tempProfileNoRef(iM,:,:).profileRaw];   

profileAlpha(:,~goodRuns) = [];
profileBeta(:,~goodRuns) = [];
profileGamma(:,~goodRuns) = [];
profileTheta(:,~goodRuns) = [];
profileSpiking(:,~goodRuns) = [];

% Taking the medians of the temporal profiles across all runs
medTempProfileAll(:,1) = median(profileTheta,2,'omitnan');
medTempProfileAll(:,2) = median(profileAlpha,2,'omitnan');
medTempProfileAll(:,3) = median(profileBeta,2,'omitnan');
medTempProfileAll(:,4) = median(profileGamma,2,'omitnan');
medTempProfileAll(:,5) = median(profileSpiking,2,'omitnan');

% Plot median temporal profiles for all runs
figure;
for iPlot = 1:5
    clear semProfile
    switch iPlot
        case 1
            semProfile = std(profileTheta,0,2)./sqrt(size(profileTheta,2));
        case 2
            semProfile = std(profileAlpha,0,2)./sqrt(size(profileAlpha,2));
        case 3
            semProfile = std(profileBeta,0,2)./sqrt(size(profileBeta,2));
        case 4
            semProfile = std(profileGamma,0,2)./sqrt(size(profileGamma,2));
        case 5
            semProfile = std(profileSpiking,0,2)./sqrt(size(profileSpiking,2));
    end

    subplot(2,3,iPlot);
    plot(-200:200,smooth(medTempProfileAll(:,iPlot),3),'k','LineWidth',1); hold on;
    patch([-200:200 fliplr(-200:200)], [(medTempProfileAll(:,iPlot) - 2.*semProfile);...
        flipud((medTempProfileAll(:,iPlot) + 2.*semProfile))],'blue','FaceAlpha',0.3,'EdgeColor','none')
    title(bandLabels{iPlot}); box off; xline(0);
    xticks(-200:50:200);xticklabels(-20:5:20);ylim([-0.35 0.2]); yticks(-0.5:0.1:0.5);
end
xlabel('Lags (s)'); ylabel('Cross correlation between ISOI and LFP Powers');

% Show all the median temporal profiles of all bands on the same plot
figure; plot(-200:200,movmedian(medTempProfileAll,3),'LineWidth',1.5);
xticks(-200:50:200);xticklabels(-20:5:20);ylim([-0.5 0.5]);
xlabel('Lags (s)'); ylabel('Cross correlation between ISOI and LFP Powers');
legend(bandLabels,'AutoUpdate','off','Location','northeast','box', 'off');
xline(0); box off;


% Show the distribution of lags and peak negative correlations 
allChCorr(:,1)= [tempProfileNoRef(iM,:,:).magLowTheta]'; 
allChCorr(:,2)= [tempProfileNoRef(iM,:,:).magLowAlpha]'; 
allChCorr(:,3)= [tempProfileNoRef(iM,:,:).magLowBeta]'; 
allChCorr(:,4) = [tempProfileNoRef(iM,:,:).magLow]'; 
allChCorr(:,5) = [tempProfileNoRef(iM,:,:).magLowRaw]'; 
allChCorr(~goodRuns,:) = []; 

allChLagVal(:,1)= [tempProfileNoRef(iM,:,:).lagLowTheta]'; 
allChLagVal(:,2)= [tempProfileNoRef(iM,:,:).lagLowAlpha]'; 
allChLagVal(:,3)= [tempProfileNoRef(iM,:,:).lagLowBeta]'; 
allChLagVal(:,4) = [tempProfileNoRef(iM,:,:).lagLow]'; 
allChLagVal(:,5) = [tempProfileNoRef(iM,:,:).lagLowRaw]'; 
allChLagVal(~goodRuns,:) = []; 

figure; 
subplot(121); boxplot(allChCorr,{'Theta';'Alpha';'Beta';'Gamma';'Spiking'});
ylim([-0.5 0.15]); box off;
subplot(122); boxplot(allChLagVal,{'Theta';'Alpha';'Beta';'Gamma';'Spiking'}); 
ylim([-100 20]); yticks(-100:10:20); yticklabels(-10:1:2); box off;

% Check if there is an effect of frequency on peak negative correlations
[pCorrROI,tblCorrROI,statsCorrROI] = anova1(allChCorr,{'Theta';'Alpha';'Beta';'Gamma';'Spiking'},'off');
[rCorrROI,~,~,gnamesCorrROI] = multcompare(statsCorrROI,"CriticalValueType","bonferroni","Display","off");

tblCorrMROI = array2table(rCorrROI,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblCorrMROI.("Group") = gnamesCorrROI(tblCorrMROI.("Group"));
tblCorrMROI.("Control Group") = gnamesCorrROI(tblCorrMROI.("Control Group"));

% Check if there is an effect of frequency on lags 
[pLagROI,tblLagROI,statsLagROI] = anova1(allChCorr,{'Theta';'Alpha';'Beta';'Gamma';'Spiking'},'off');
[rLagROI,~,~,gnamesLagROI] = multcompare(statsLagROI,"CriticalValueType","bonferroni","Display","off");

tblLagMROI = array2table(rLagROI,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblLagMROI.("Group") = gnamesLagROI(tblLagMROI.("Group"));
tblLagMROI.("Control Group") = gnamesLagROI(tblLagMROI.("Control Group"));

% Show the distributions of lags and correlations for 10/10 and 6/6 split
% for superficial and deep channels
superCorr_10Ch_Gamma = [tempProfileSuperNoRef(iM,:,:).magLow]'; 
deepCorr_10Ch_Gamma  = [tempProfileDeepNoRef(iM,:,:).magLow]';
superCorr_6Ch_Gamma  = [tempProfileSuper_6Ch_NoRef(iM,:,:).magLow]'; 
deepCorr_6Ch_Gamma   = [tempProfileDeep_6Ch_NoRef(iM,:,:).magLow]'; 

superCorr_10Ch_Gamma(~goodRuns) = []; superCorr_10Ch_Gamma(singleChFlag) = []; 
deepCorr_10Ch_Gamma(~goodRuns)  = []; deepCorr_10Ch_Gamma(singleChFlag)  = []; 
superCorr_6Ch_Gamma(~goodRuns)  = []; superCorr_6Ch_Gamma(singleChFlag)  = [];
deepCorr_6Ch_Gamma(~goodRuns)  = [];  deepCorr_6Ch_Gamma(singleChFlag)  = []; 

superLag_10Ch_Gamma = [tempProfileSuperNoRef(iM,:,:).lagLow]'./10; 
deepLag_10Ch_Gamma  = [tempProfileDeepNoRef(iM,:,:).lagLow]'./10;
superLag_6Ch_Gamma  = [tempProfileSuper_6Ch_NoRef(iM,:,:).lagLow]'./10; 
deepLag_6Ch_Gamma   = [tempProfileDeep_6Ch_NoRef(iM,:,:).lagLow]'./10;

superLag_10Ch_Gamma(~goodRuns) = []; superLag_10Ch_Gamma(singleChFlag) = []; 
deepLag_10Ch_Gamma(~goodRuns)  = []; deepLag_10Ch_Gamma(singleChFlag)  = [];
superLag_6Ch_Gamma(~goodRuns)  = []; superLag_6Ch_Gamma(singleChFlag)  = [];
deepLag_6Ch_Gamma(~goodRuns)   = []; deepLag_6Ch_Gamma(singleChFlag)   = []; 

% Show the distributions of the lags and correlations for super/deep ch
figure; 
% Show the correlation distributions
subplot(121); boxplot([superCorr_10Ch_Gamma, deepCorr_10Ch_Gamma, superCorr_6Ch_Gamma, ...
    deepCorr_6Ch_Gamma],{'Super_10Ch';'Deep_10Ch'; 'Super_6Ch'; 'Deep_6Ch'}); 
ylim([-0.8 0]);yticks(-0.8:0.1:0); ylabel('ROI - Peak negative Correlation');

% Show the lag distributions
subplot(122); boxplot([superLag_10Ch_Gamma deepLag_10Ch_Gamma superLag_6Ch_Gamma ...
    deepLag_6Ch_Gamma],{'Super_10Ch';'Deep_10Ch'; 'Super_6Ch'; 'Deep_6Ch'}); 
ylim([-10 2]);yticks(-15:1:2); ylabel('ROI - Peak negative Lag');

% Check if channel split has an effect on correlations 
[pCorr,~,statsCorr] = anova1([superCorr_10Ch_Gamma deepCorr_10Ch_Gamma superCorr_6Ch_Gamma ...
    deepCorr_6Ch_Gamma],{'Super_10Ch';'Deep_10Ch'; 'Super_6Ch'; 'Deep_6Ch'},'off');
[rCorr,~,~,gnamesCorr] = multcompare(statsCorr,"CriticalValueType","bonferroni","Display","off");

tblCorrM = array2table(rCorr,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblCorrM.("Group") = gnamesCorr(tblCorrM.("Group"));
tblCorrM.("Control Group") = gnamesCorr(tblCorrM.("Control Group"));

% Check if channel split has an effect on lags  
[pLag,tblLag,statsLag] = anova1([superLag_10Ch_Gamma deepLag_10Ch_Gamma superLag_6Ch_Gamma ...
    deepLag_6Ch_Gamma],{'Super_10Ch';'Deep_10Ch'; 'Super_6Ch'; 'Deep_6Ch'},'off');
[rLag,~,~,gnamesLag] = multcompare(statsLag,"CriticalValueType","bonferroni","Display","off");

tblLagM = array2table(rLag,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblLagM.("Group") = gnamesLag(tblLagM.("Group"));
tblLagM.("Control Group") = gnamesLag(tblLagM.("Control Group"));

% Get the median +/- mad lag 
x = -200:200;
allChLag = [tempProfileNoRef(iM,:,:).lagLow]'; allChLag(~goodRuns) = [];  allChLag(singleChFlag) = []; 
lagFrameRange = find(x == round(median(allChLag)- mad(allChLag))) : find(x == round(median(allChLag)+ mad(allChLag)));

%% Move files to server
tic;
for iDate = 1:size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);

    for iRun = 1: size(allRuns{iDate,1})
        clear runName dataDir clipMask elecMask clipMaskCortex corrMask...
            x negIdx lowIdx
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName];
        serverDir = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\kem294\Data\' ...
            monkeyName '_SqM\' hemisphere ' Hemisphere\'  expDate '\' runName ];
        
        disp(['Moving  data for ' monkeyName ' '  expDate ' run: ' runName]);
        if ~exist(serverDir,'dir')
            [~,~] = mkdir(serverDir);
        end
        if ~exist([serverDir '\crossCorrFOV_10_NoRef.mat'],'file')
            movefile([dataDir '\crossCorrFOV_10_NoRef.mat'],serverDir);
        end

        if ~exist([serverDir '\crossCorrFOV_6_NoRef.mat'],'file')
            movefile([dataDir '\crossCorrFOV_6_NoRef.mat'],serverDir);
        end

        if ~exist([serverDir '\crossCorrFOV_6_BipolarRef.mat'],'file')
            movefile([dataDir '\crossCorrFOV_6_BipolarRef.mat'],serverDir);
        end

        if ~exist([serverDir '\crossCorrFOV_6_AvgRef.mat'],'file')
            movefile([dataDir '\crossCorrFOV_6_AvgRef.mat'],serverDir);
        end
    end
end
toc;

%%  Cross correlating physiology and imaging for FOV
videoFlag      = 0; % To save videos of hybrids - change if you need to
corrFCHybrid   = NaN(size(processedDat,1),size(processedDat,2),5,4,401);
corrFCSuper    = NaN(size(processedDat,1),size(processedDat,2),5,4,401);
corrFCDeep     = NaN(size(processedDat,1),size(processedDat,2),5,4,401);
super_DeepCorr = NaN(size(processedDat,1),size(processedDat,2),5,4,401);
corrFCMid      = NaN(size(processedDat,1),size(processedDat,2),5,4,401);
super_MidCorr  = NaN(size(processedDat,1),size(processedDat,2),5,4,401);
mid_DeepCorr   = NaN(size(processedDat,1),size(processedDat,2),5,4,401);

super_DeepAvgFrames = NaN(size(processedDat,1),size(processedDat,2),5,4);
super_MidAvgFrames  = NaN(size(processedDat,1),size(processedDat,2),5,4);
deep_MidAvgFrames   = NaN(size(processedDat,1),size(processedDat,2),5,4);
peakNegValsAll      = NaN(size(processedDat,1),size(processedDat,2),5,4);
peakNegTimesAll     = NaN(size(processedDat,1),size(processedDat,2),5,4);

superHybridAllBands = NaN(size(processedDat,1),size(processedDat,2),4,5,5);
midHybridAllBands   = NaN(size(processedDat,1),size(processedDat,2),4,5,5);
deepHybridAllBands  = NaN(size(processedDat,1),size(processedDat,2),4,5,5);
bandNames           = {'Theta'; 'Alpha';'Beta';'Gamma';'Spiking'};

tic;
for iDate = 1:size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);
    for iRun = 1: size(allRuns{iDate,1})
        clear runName dataDir clipMask elecMask clipMaskCortex corrMask...
            x negIdx lowIdx serverDir
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName];
        serverDir = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\kem294\Data\' ...
            monkeyName '_SqM\' hemisphere ' Hemisphere\'  expDate '\' runName ];

        clc; disp(['Analyzing data for ' monkeyName ' '  expDate ' run: ' runName]);

        if ~exist([serverDir '\crossCorrFOV_10_NoRef.mat'],'file')
            disp('No reference - channel split 10(superficial)/10(Deep)');
            getCrossCorrFOV(monkeyName,expDate,runName,serverDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,10,'NoRef');
        end

        if ~exist([serverDir '\crossCorrFOV_6_NoRef.mat'],'file')
            disp('No reference - channel split top 6(superficial)/bottom 6(Deep)');
            getCrossCorrFOV(monkeyName,expDate,runName,serverDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,6,'NoRef');
        end

        if ~exist([serverDir '\crossCorrFOV_6_BipolarRef.mat'],'file')
            disp('Bipolar reference -  channel split top 6(superficial)/bottom 6(Deep)');
            getCrossCorrFOV(monkeyName,expDate,runName,serverDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,6,'BipolarRef');
        end

        if ~exist([serverDir '\crossCorrFOV_6_AvgRef.mat'],'file')
            disp('Avg reference - channel split top 6(superficial)/bottom 6(Deep)');
            getCrossCorrFOV(monkeyName,expDate,runName,serverDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,6,'AvgRef');
        end
        
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

            if corrMaskFlag
                corrMask = imresize(corrMask,1/3); % Resize the mask
                corrMask = corrMask(:,:,1)>0;
                corrMask = corrMask & ~elecMask;
            else
                corrMask = clipMaskCortex;
            end

            corrMaskT = reshape(corrMask,[imSize(1)*imSize(2) 1]);

            % Get ROI location and FC map
            seedLocIn    = load([dataDir '\roiCenterLoc.mat']);
            seedLocProbe = seedLocIn.seedLocProbe;
            seedLocIn    = seedLocIn.seedLocIn;
            circleRad    = round(roiSize{iDate}(iRun)/(spatialBin)); % 500um radius
            greenFig     = imresize(greenIm{iDate,iRun},1/spatialBin,'OutputSize',[imSize(1) imSize(2)]);

            pDatTemp = processedDat{iDate,iRun}.tempBandPass;
            seedSigT = calculateSeedSignal(greenFig,corrMask,...
                seedLocIn,circleRad,pDatTemp); % Get Gaussian weighted seed signal

            fcMap             = plotCorrMap(seedSigT,pDatTemp,0);
            fcMap             = reshape(fcMap,[361*438 1]);
            fcMap(~corrMaskT) = NaN;

            chLen = estChInCortex{1,iDate}(iRun,2) - estChInCortex{1,iDate}(iRun,1);

            % Initialize run-related variables
            % Correlations within the same frequencies            
            corrFCHybridT = NaN(5,4,401); corrFCSuperT    = NaN(5,4,401);
            corrFCDeepT   = NaN(5,4,401); super_DeepCorrT = NaN(5,4,401);

            corrFCMidT          = NaN(5,4,401); super_MidCorrT       = NaN(5,4,401);
            mid_DeepCorrT       = NaN(5,4,401); peakNegValsAllT      = NaN(5,4);
            peakNegTimesAllT    = NaN(5,4);     super_DeepAvgFramesT = NaN(5,4);
            super_MidAvgFramesT = NaN(5,4);     deep_MidAvgFramesT   = NaN(5,4);
           
            % Correlations between frequencies    
            superHybridAllBandsT = NaN(4,5,5);  deepHybridAllBandsT = NaN(4,5,5);
            midHybridAllBandsT   = NaN(4,5,5); crossFreqCrossLayerHybridT = NaN(4,5,5,3,3);

            for iType = 1:4 % 10/10 split or 6/6 split of superficial/deep channels
                clear crossCorrFOV allXCorr superXCorr deepXCorr allLags fileName
                switch iType
                    case 1
                        crossCorrFOV = matfile([serverDir '\crossCorrFOV_10_NoRef.mat']);
                        allXcorr     = crossCorrFOV.spatialProfile;
                        superXcorr   = crossCorrFOV.spatialProfileSuper;
                        deepXcorr    = crossCorrFOV.spatialProfileDeep;
                        allLags      = crossCorrFOV.lagFull;
                        fileName     = '10_NoRef';

                    case 2
                        crossCorrFOV = matfile([serverDir '\crossCorrFOV_6_NoRef.mat']);
                        allXcorr     = crossCorrFOV.spatialProfile;
                        superXcorr   = crossCorrFOV.spatialProfileSuper;
                        deepXcorr    = crossCorrFOV.spatialProfileDeep;
                        if chLen~=0
                            midXCorr     = crossCorrFOV.spatialProfileMid;
                        end
                        allLags      = crossCorrFOV.lagFull;
                        fileName     = '6_NoRef';

                    case 3
                        crossCorrFOV = matfile([serverDir '\crossCorrFOV_6_BipolarRef.mat']);
                        allXcorr     = crossCorrFOV.spatialProfile;
                        superXcorr   = crossCorrFOV.spatialProfileSuper;
                        deepXcorr    = crossCorrFOV.spatialProfileDeep;
                        if chLen~=0
                            midXCorr     = crossCorrFOV.spatialProfileMid;
                        end
                        allLags      = crossCorrFOV.lagFull;
                        fileName     = '6_BipolarRef';

                    case 4
                        crossCorrFOV = matfile([serverDir '\crossCorrFOV_6_AvgRef.mat']);
                        allXcorr     = crossCorrFOV.spatialProfile;
                        superXcorr   = crossCorrFOV.spatialProfileSuper;
                        deepXcorr    = crossCorrFOV.spatialProfileDeep;
                        if chLen~=0
                            midXCorr     = crossCorrFOV.spatialProfileMid;
                        end
                        allLags      = crossCorrFOV.lagFull;
                        fileName     = '6_AvgRef';

                end

                x = allLags;
                negIdx = x<0 & x>=-150; negVals = x(negIdx);
                lowIdx = x<0 & x>= -80; xLow = x(lowIdx);

                % Show and save the full field of view for the hybrid map
                % Within frequency comparison

                % Theta
                crossCorrTheta      = reshape(allXcorr.ccFullTheta,[401 imSize(1)*imSize(2)]);   crossCorrTheta(:,~corrMaskT)      = NaN;
                crossCorrSuperTheta = reshape(superXcorr.ccFullTheta,[401 imSize(1)*imSize(2)]); crossCorrSuperTheta(:,~corrMaskT) = NaN;
                crossCorrDeepTheta  = reshape(deepXcorr.ccFullTheta,[401 imSize(1)*imSize(2)]);  crossCorrDeepTheta(:,~corrMaskT)  = NaN;
                lagLowTheta         = tempProfileNoRef(iM,iDate,iRun).lagLowTheta;
            
                % Alpha
                crossCorrAlpha      = reshape(allXcorr.ccFullAlpha,[401 imSize(1)*imSize(2)]);   crossCorrAlpha(:,~corrMaskT)      = NaN;
                crossCorrSuperAlpha = reshape(superXcorr.ccFullAlpha,[401 imSize(1)*imSize(2)]); crossCorrSuperAlpha(:,~corrMaskT) = NaN;
                crossCorrDeepAlpha  = reshape(deepXcorr.ccFullAlpha,[401 imSize(1)*imSize(2)]);  crossCorrDeepAlpha(:,~corrMaskT)  = NaN;
                lagLowAlpha         = tempProfileNoRef(iM,iDate,iRun).lagLowAlpha;
              
                % Beta
                crossCorrBeta      = reshape(allXcorr.ccFullBeta,[401 imSize(1)*imSize(2)]);   crossCorrBeta(:,~corrMaskT)      = NaN;
                crossCorrSuperBeta = reshape(superXcorr.ccFullBeta,[401 imSize(1)*imSize(2)]); crossCorrSuperBeta(:,~corrMaskT) = NaN;
                crossCorrDeepBeta  = reshape(deepXcorr.ccFullBeta,[401 imSize(1)*imSize(2)]);  crossCorrDeepBeta(:,~corrMaskT)  = NaN;
                lagLowBeta         = tempProfileNoRef(iM,iDate,iRun).lagLowBeta;
               
                % Gamma
                crossCorrGamma      = reshape(allXcorr.ccFull,[401 imSize(1)*imSize(2)]);   crossCorrGamma(:,~corrMaskT)      = NaN;
                crossCorrSuperGamma = reshape(superXcorr.ccFull,[401 imSize(1)*imSize(2)]); crossCorrSuperGamma(:,~corrMaskT) = NaN;
                crossCorrDeepGamma  = reshape(deepXcorr.ccFull,[401 imSize(1)*imSize(2)]);  crossCorrDeepGamma(:,~corrMaskT)  = NaN;
                lagLowGamma         = tempProfileNoRef(iM,iDate,iRun).lagLow;
               
                % Spiking
                crossCorrSpiking      = reshape(allXcorr.ccFullRaw,[401 imSize(1)*imSize(2)]);   crossCorrSpiking(:,~corrMaskT)      = NaN;
                crossCorrSuperSpiking = reshape(superXcorr.ccFullRaw,[401 imSize(1)*imSize(2)]); crossCorrSuperSpiking(:,~corrMaskT) = NaN;
                crossCorrDeepSpiking  = reshape(deepXcorr.ccFullRaw,[401 imSize(1)*imSize(2)]);  crossCorrDeepSpiking(:,~corrMaskT)  = NaN;
                lagLowSpiking         = tempProfileNoRef(iM,iDate,iRun).lagLowRaw;
               

                if iType >= 2 && chLen~=0
                    crossCorrMidTheta = reshape(midXCorr.ccFullTheta,[401 imSize(1)*imSize(2)]); crossCorrMidTheta(:,~corrMaskT)   = NaN;
                    crossCorrMidAlpha = reshape(midXCorr.ccFullAlpha,[401 imSize(1)*imSize(2)]); crossCorrMidAlpha(:,~corrMaskT)   = NaN;
                    crossCorrMidBeta = reshape(midXCorr.ccFullBeta,[401 imSize(1)*imSize(2)]);   crossCorrMidBeta(:,~corrMaskT)    = NaN;
                    crossCorrMidGamma = reshape(midXCorr.ccFull,[401 imSize(1)*imSize(2)]);      crossCorrMidGamma(:,~corrMaskT)   = NaN;
                    crossCorrMidSpiking = reshape(midXCorr.ccFullRaw,[401 imSize(1)*imSize(2)]); crossCorrMidSpiking(:,~corrMaskT) = NaN;
                end               
                   
                % Correlate each hybrid map(all,super,deep,mid maps)
                % with FC map and correlate super/deep/mid maps to one
                % another
                for iBand = 1:5
                    clear bandName mapsAll crossCorrSuperR crossCorrDeepR crossCorrMidR
                    switch iBand
                        case 1
                            bandName        = 'Theta';
                            mapsAll         = crossCorrTheta;
                            crossCorrSuperR = crossCorrSuperTheta;
                            crossCorrDeepR  = crossCorrDeepTheta;
                            lagLow          = lagLowTheta;

                            if iType >= 2 && chLen~=0
                                crossCorrMidR   = crossCorrMidTheta;                                
                            end

                        case 2
                            bandName        = 'Alpha';
                            mapsAll         = crossCorrAlpha;
                            crossCorrSuperR = crossCorrSuperAlpha;
                            crossCorrDeepR  = crossCorrDeepAlpha;
                            lagLow          = lagLowAlpha;

                            if iType >= 2 && chLen~=0
                                crossCorrMidR   = crossCorrMidAlpha;
                            end

                        case 3
                            bandName        = 'Beta';
                            mapsAll         = crossCorrBeta;
                            crossCorrSuperR = crossCorrSuperBeta;
                            crossCorrDeepR  = crossCorrDeepBeta;
                            lagLow          = lagLowBeta;

                            if iType >= 2 && chLen~=0
                                crossCorrMidR   = crossCorrMidBeta;
                            end

                        case 4
                            bandName        = 'Gamma';
                            mapsAll         = crossCorrGamma;
                            crossCorrSuperR = crossCorrSuperGamma;
                            crossCorrDeepR  = crossCorrDeepGamma;
                            lagLow          = lagLowGamma;

                            if iType >= 2 && chLen~=0
                                crossCorrMidR   = crossCorrMidGamma;
                            end

                        case 5
                            bandName        = 'Spiking';
                            mapsAll         = crossCorrSpiking;
                            crossCorrSuperR = crossCorrSuperSpiking;
                            crossCorrDeepR  = crossCorrDeepSpiking;
                            lagLow          = lagLowSpiking;

                            if iType >= 2 && chLen~=0
                                crossCorrMidR   = crossCorrMidSpiking;
                            end
                    end

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

                    % Average the super/deep hybrid maps for a certain range and correlate them
                    % with one another for 10/10 and 6/6
                    clear superHybridAvg deepHybridAvg
                    superHybridAvg = mean(crossCorrSuperR(lagFrameRange,:),1,'omitnan');
                    deepHybridAvg  = mean(crossCorrDeepR(lagFrameRange,:),1,'omitnan');
                    super_DeepAvgFramesT(iBand,iType) = corr(superHybridAvg',deepHybridAvg','rows','complete');

                    if iType >= 2 && chLen~=0
                        midHybridAvg = mean(crossCorrMidR(lagFrameRange,:),1,'omitnan');
                        super_MidAvgFramesT(iBand,iType) = corr(superHybridAvg',midHybridAvg','rows','complete');
                        deep_MidAvgFramesT(iBand,iType)  = corr(deepHybridAvg',midHybridAvg','rows','complete');
                    end

                    if ~exist([dataDir '\HybridMapFOV_' fileName '_' bandName '.png'],'file')
                        clear frameLow grayImFull crossCorr
                        frameLow       = (x == lagLow);
                        crossCorr      = reshape(mapsAll,[401 imSize(1) imSize(2)]);
                        figure; imagesc(squeeze(crossCorr(frameLow,:,:))); hold on; axis image off;
                        colormap(flipud(jet)); caxis([-0.5 0.5]);
                        grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));
                        h = imagesc(grayImFull); hold off; set(h,'AlphaData',~corrMask);
                        title(['Peak negative at ' num2str(lagLow/10) 's']);
                        f = gcf; exportgraphics(f,[dataDir '\HybridMapFOV_' fileName '_' bandName '.png'],'Resolution',300); close gcf;
                    end

                    if ~exist([dataDir '\HybridMapFOV_Super_' fileName '_' bandName '.png'],'file')
                        clear frameLow grayImFull crossCorrSuper
                        frameLow       = (x == lagLow);
                        crossCorrSuper = reshape(crossCorrSuperR,[401 imSize(1) imSize(2)]);
                        figure; imagesc(squeeze(crossCorrSuper(frameLow,:,:))); hold on; axis image off;
                        colormap(flipud(jet)); caxis([-0.5 0.5]);
                        grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));
                        h = imagesc(grayImFull); hold off; set(h,'AlphaData',~corrMask);
                        title(['Superficial channels-Peak negative at ' num2str(lagLow/10) 's']);
                        f = gcf; exportgraphics(f,[dataDir '\HybridMapFOV_Super_' fileName '_' bandName '.png'],'Resolution',300); close gcf;
                    end

                    if ~exist([dataDir '\HybridMapFOV_Deep_' fileName '_' bandName '.png'],'file')
                        clear frameLow grayImFull crossCorrDeep
                        frameLow      = (x == lagLow);
                        crossCorrDeep = reshape(crossCorrDeepR,[401 imSize(1) imSize(2)]);
                        figure; imagesc(squeeze(crossCorrDeep(frameLow,:,:))); hold on; axis image off;
                        colormap(flipud(jet)); caxis([-0.5 0.5]);
                        grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));
                        h = imagesc(grayImFull); hold off; set(h,'AlphaData',~corrMask);
                        title(['Deep channels-Peak negative at ' num2str(lagLow/10) 's']);
                        f = gcf; exportgraphics(f,[dataDir '\HybridMapFOV_Deep_' fileName '_' bandName '.png'],'Resolution',300); close gcf;
                    end

                    if ~exist([dataDir '\HybridMapFOV_Mid_' fileName '_' bandName '.png'],'file') && iType >= 2 && chLen~=0
                        clear frameLow grayImFull crossCorrMid
                        frameLow     = (x == lagLow);
                        crossCorrMid = reshape(crossCorrMidR,[401 imSize(1) imSize(2)]);
                        figure; imagesc(squeeze(crossCorrMid(frameLow,:,:))); hold on; axis image off;
                        colormap(flipud(jet)); caxis([-0.5 0.5]);
                        grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));
                        h = imagesc(grayImFull); hold off; set(h,'AlphaData',~corrMask);
                        title(['Middle channels-Peak negative at ' num2str(lagLow/10) 's']);
                        f = gcf; exportgraphics(f,[dataDir '\HybridMapFOV_Mid_' fileName '_' bandName '.png'],'Resolution',300); close gcf;
                    end
                end

                % Correlate maps across frequencies for each layer
                % compartment
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

            save([dataDir '\processedHybridMapVars.mat'],'corrFCHybridT','corrFCSuperT',...
                'corrFCDeepT','super_DeepCorrT','corrFCMidT','super_MidCorrT','mid_DeepCorrT',...
                'peakNegValsAllT','peakNegTimesAllT','super_DeepAvgFramesT','super_MidAvgFramesT',...
                'deep_MidAvgFramesT','superHybridAllBandsT','deepHybridAllBandsT','midHybridAllBandsT',...
                'crossFreqCrossLayerHybridT','fcMap');

            corrFCHybrid(iDate,iRun,:,:,:)   = corrFCHybridT;
            corrFCSuper(iDate,iRun,:,:,:)    = corrFCSuperT;
            corrFCDeep(iDate,iRun,:,:,:)     = corrFCDeepT;
            corrFCMid(iDate,iRun,:,:,:)      = corrFCMidT;
           
            super_DeepCorr(iDate,iRun,:,:,:) = super_DeepCorrT;           
            super_MidCorr(iDate,iRun,:,:,:)  = super_MidCorrT;
            mid_DeepCorr(iDate,iRun,:,:,:)   = mid_DeepCorrT ;

            peakNegValsAll(iDate,iRun,:,:)  = peakNegValsAllT;
            peakNegTimesAll(iDate,iRun,:,:) = peakNegTimesAllT;

            super_DeepAvgFrames(iDate,iRun,:,:) = super_DeepAvgFramesT;
            super_MidAvgFrames(iDate,iRun,:,:)  = super_MidAvgFramesT;
            deep_MidAvgFrames(iDate,iRun,:,:)   = deep_MidAvgFramesT;

            superHybridAllBands(iDate,iRun,:,:,:) = superHybridAllBandsT;
            midHybridAllBands(iDate,iRun,:,:,:)   = midHybridAllBandsT;
            deepHybridAllBands(iDate,iRun,:,:,:)  = deepHybridAllBandsT;

            crossFreqCrossLayerHybrid(iDate,iRun,:,:,:,:,:) = crossFreqCrossLayerHybridT; 

        else
            clear allHybridVars
            allHybridVars                    = matfile([dataDir '\processedHybridMapVars.mat']);
            corrFCHybrid(iDate,iRun,:,:,:)   = allHybridVars.corrFCHybridT;
            corrFCSuper(iDate,iRun,:,:,:)    = allHybridVars.corrFCSuperT;
            corrFCDeep(iDate,iRun,:,:,:)     = allHybridVars.corrFCDeepT;
            super_DeepCorr(iDate,iRun,:,:,:) = allHybridVars.super_DeepCorrT;

            corrFCMid(iDate,iRun,:,:,:)     = allHybridVars.corrFCMidT;
            super_MidCorr(iDate,iRun,:,:,:) = allHybridVars.super_MidCorrT;
            mid_DeepCorr(iDate,iRun,:,:,:)  = allHybridVars.mid_DeepCorrT ;

            peakNegValsAll(iDate,iRun,:,:)  = allHybridVars.peakNegValsAllT;
            peakNegTimesAll(iDate,iRun,:,:) = allHybridVars.peakNegTimesAllT;

            super_DeepAvgFrames(iDate,iRun,:,:) = allHybridVars.super_DeepAvgFramesT;
            super_MidAvgFrames(iDate,iRun,:,:)  = allHybridVars.super_MidAvgFramesT;
            deep_MidAvgFrames(iDate,iRun,:,:)   = allHybridVars.deep_MidAvgFramesT;

            superHybridAllBands(iDate,iRun,:,:,:) = allHybridVars.superHybridAllBandsT;
            midHybridAllBands(iDate,iRun,:,:,:)   = allHybridVars.midHybridAllBandsT;
            deepHybridAllBands(iDate,iRun,:,:,:)  = allHybridVars.deepHybridAllBandsT;
           
            crossFreqCrossLayerHybrid(iDate,iRun,:,:,:,:,:) = allHybridVars.crossFreqCrossLayerHybridT;

        end
    end
end
toc;

%% Compile all variables and plot
corrFCHybridT = reshape(corrFCHybrid,[size(corrFCHybrid,1)*size(corrFCHybrid,2) size(corrFCHybrid,3) size(corrFCHybrid,4) size(corrFCHybrid,5)]);
nanRow        = isnan(corrFCHybridT(:,1,1,1));
corrFCHybridT(nanRow,:,:,:) = [];
corrFCHybridT(~goodRunsSpatial,:,:,:) = [];

% Get the median profiles for all frequency bands
medTempProfileHybrid =  squeeze(median(squeeze(corrFCHybridT(:,:,2,:)),1,'omitnan')); 

% Plot the median +/- 2*sem
figure;
for iPlot = 1:5
    subplot(2,3,iPlot);
    semProfile = std(squeeze(corrFCHybridT(:,iPlot,2,:)),0,1)./sqrt(size(squeeze(corrFCHybridT(:,iPlot,2,:)),1));
    plot(-200:200,smooth(medTempProfileHybrid(iPlot,:),7),'k','LineWidth',1); hold on;
    patch([-200:200 fliplr(-200:200)], [(medTempProfileHybrid(iPlot,:) - 2.*semProfile)...
        fliplr((medTempProfileHybrid(iPlot,:) + 2.*semProfile))],'blue','FaceAlpha',0.3,'EdgeColor','none')
    title(bandLabels{iPlot}); box off; xline(0);
    xticks(-200:50:200);xticklabels(-20:5:20);ylim([-0.9 0.8]); yticks(-1:0.2:0.8);
end 

% Compile the peak negative correlations and lags for good runs
peakNegValsAllT = reshape(peakNegValsAll,[size(peakNegValsAll,1)*size(peakNegValsAll,2) size(peakNegValsAll,3) size(peakNegValsAll,4)]);
peakNegValsAllT(nanRow,:,:) = [];
peakNegValsAllT(~goodRunsSpatial,:,:) = [];
figure; boxplot(squeeze(peakNegValsAllT(:,:,2)),{'Theta';'Alpha';'Beta';'Gamma';'Spiking'});
ylim([-1 0.6]); box off; ylabel('Peak negative correlations with FC map');

% Check if there is an effect of frequency on peak negative correlations
[pBand,~,statsBand] = anova1(peakNegValsAllT(:,:,2),...
    {'Theta';'Alpha';'Beta';'Gamma';'Spiking'},'off');
[rBand,~,~,gnamesBand] = multcompare(statsBand,"Display","off");

tblBand = array2table(rBand,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblBand.("Group") = gnamesBand(tblBand.("Group"));
tblBand.("Control Group") = gnamesBand(tblBand.("Control Group"));

% Compile the correlations between superficial and deep hybrid maps 
super_DeepCorrT = reshape(super_DeepCorr,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4) size(super_DeepCorr,5)]);
nanRow         = find(isnan(super_DeepCorrT(:,1,1,1)));
super_DeepCorrT(nanRow,:,:,:) = [];
super_DeepCorrT(~goodRunsSpatial,:,:,:) = [];
super_DeepCorrT(singleChFlag,:,:,:) = [];

% Plot cross correlations as a function of lag for all runs and for all
% frequency bands and check if 10/10 vs 6/6 correlation splits are similar
for iBand = 2:4
    figure;
    for iType = 1:2
        subplot(1,2,iType); plot(-200:200,squeeze(super_DeepCorrT(:,iBand,iType,:)),'Color',[0.65 0.65 0.65]);
        xline(0,'LineWidth',1); xlabel('Lags (s)'); ylabel('Correlations between superficial and deep maps');
        xticks(-200:50:200); xticklabels(-20:5:20); ylim([-0.5 1]);
        sgtitle(bandNames{iBand});
    end

    corrBetween6_10 = diag(corr(squeeze(super_DeepCorrT(:,iBand,1,:))',squeeze(super_DeepCorrT(:,iBand,2,:))','rows','complete'));
    figure; boxplot(corrBetween6_10,{'Correlation between 10/10 vs 6/6'});
    title(bandNames{iBand}); ylim([0.9 1]);yticks(0.85:0.01:1);
end 

% Compile the Superficial/Deep, Superficial/Middle and Middle/Deep
% correlations by averaging frames/correlations based on lag obtained from ROI
avgCorrSuperDeepHybrids = mean(super_DeepCorr(:,:,:,:,lagFrameRange),5,'omitnan'); 
avgCorrSuperDeepHybrids = reshape(avgCorrSuperDeepHybrids,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]); 
avgCorrSuperDeepHybrids(nanRow,:,:)           = []; 
avgCorrSuperDeepHybrids(~goodRunsSpatial,:,:) = []; 
avgCorrSuperDeepHybrids(singleChFlag,:,:)     = []; 

super_DeepAvgFramesT = reshape(super_DeepAvgFrames,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]);
super_DeepAvgFramesT(nanRow,:,:)           = [];
super_DeepAvgFramesT(~goodRunsSpatial,:,:) = [];
super_DeepAvgFramesT(singleChFlag,:,:)     = [];

super_MidAvgFramesT = reshape(super_MidAvgFrames,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]);
super_MidAvgFramesT(nanRow,:,:)           = [];
super_MidAvgFramesT(~goodRunsSpatial,:,:) = [];
super_MidAvgFramesT(singleChFlag,:,:)     = [];

deep_MidAvgFramesT = reshape(deep_MidAvgFrames,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]);
deep_MidAvgFramesT(nanRow,:,:)           = [];
deep_MidAvgFramesT(~goodRunsSpatial,:,:) = [];
deep_MidAvgFramesT(singleChFlag,:,:)     = [];

% Show the S/D, S/M, D/M for the different channel splits and averaging
% methods - Conclusion - Everything gives the same results
for iBand = 1:5
    figure; subplot(121);
    boxplot([squeeze(super_DeepAvgFramesT(:,iBand,1:2)) squeeze(avgCorrSuperDeepHybrids(:,iBand,1:2))],...
        {'10/10 Frame Avg';'6/6 Frame Avg';'10/10 Corr Avg'; ...
        '6/6 Corr Average'}); ylim([0 1]); box off;
    ylabel('Correlation between hybrids'); ylim([0 1.1]); yticks(0:0.1:1);

    % Show the distributions of super/deep; super/mid and mid/deep correlations
   subplot(122);
    boxplot([super_DeepAvgFramesT(:,iBand,2) super_MidAvgFramesT(:,iBand,2) deep_MidAvgFramesT(:,iBand,2)],...
        {'Super_Deep'; 'Super_Mid'; 'Deep_Mid'}); box off;
    ylabel('Correlation between hybrids'); ylim([-1.1 1.1]); yticks(-1:0.1:1);
    sgtitle(bandNames{iBand});
end

% Plot the S/D, S/M and M/D correlations for different referencing schemes
% - Change the band-value (second column of variables) if you want to check
% for other frequency bands
figure;
for  iPlot = 1:4
    subplot(2,2,iPlot);
    clear tempVar;
    tempVar = [super_DeepAvgFramesT(:,4,iPlot) super_MidAvgFramesT(:,4,iPlot) ...
        deep_MidAvgFramesT(:,4,iPlot)];

    boxplot(tempVar,{'Super_Deep'; 'Super_Mid';'Mid_Deep'});
    ylabel('Correlation between hybrids'); ylim([-0.3 1.1]); yticks(-0.5:0.1:1);
    switch iPlot
        case 1
            title('No Ref - 10/10 channel split');
        case 2
            title('No Ref - 6/6 channel split');
        case 3
            title('Bipolar Ref - 6/6 channel split');
        case 4
            title('Avg Ref - 6/6 channel split');
    end
    axis square; box off;
end
% Check if there is an effect of layer compartments on correlations for
% gamma band 
[~,~,statsType] = anova1([super_DeepAvgFramesT(:,4,2) super_MidAvgFramesT(:,4,2) deep_MidAvgFramesT(:,4,2)],...
    {'Super_Deep'; 'Super_Mid'; 'Deep_Mid'},'off');
[rType,~,~,gnamesType] = multcompare(statsType,"Display","off");

tblType = array2table(rType,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblType.("Group") = gnamesType(tblType.("Group"));
tblType.("Control Group") = gnamesType(tblType.("Control Group"));


% Compile correlations between layer compartments and FC map
corrFCSuperT = mean(corrFCSuper(:,:,:,:,lagFrameRange),5,'omitnan'); 
corrFCSuperT = reshape(corrFCSuperT,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]); 
corrFCSuperT(nanRow,:,:) = []; 
corrFCSuperT(~goodRunsSpatial,:,:) = []; 
corrFCSuperT(singleChFlag,:,:) = []; 

corrFCDeepT = mean(corrFCDeep(:,:,:,:,lagFrameRange),5,'omitnan'); 
corrFCDeepT = reshape(corrFCDeepT,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]); 
corrFCDeepT(nanRow,:,:) = []; 
corrFCDeepT(~goodRunsSpatial,:,:) = []; 
corrFCDeepT(singleChFlag,:,:) = []; 

corrFCMidT = mean(corrFCMid(:,:,:,:,lagFrameRange),5,'omitnan'); 
corrFCMidT = reshape(corrFCMidT,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]); 
corrFCMidT(nanRow,:,:) = []; 
corrFCMidT(~goodRunsSpatial,:,:) = []; 
corrFCMidT(singleChFlag,:,:) = [];

% Plot the correlations between hybrids and FC maps for different reference schemes
for iBand = 2: 4
    figure;
    for  iPlot = 1:4
        subplot(2,2,iPlot);
        boxplot([corrFCSuperT(:,iBand,iPlot) corrFCDeepT(:,iBand,iPlot) corrFCMidT(:,iBand,iPlot)],...
            {'Super_FC'; 'Deep_FC'; 'Middle_FC'});
        ylabel('Correlation between hybrids'); ylim([-1 1]); yticks(-1:0.1:1);
        switch iPlot
            case 1
                title('No Ref - 10/10 channel split');
            case 2
                title('No Ref - 6/6 channel split');
            case 3
                title('Bipolar Ref - 6/6 channel split');
            case 4
                title('Avg Ref - 6/6 channel split');
        end
        axis square; box off;
    end
    sgtitle(bandNames{iBand});

    figure;
    subplot(121);
    boxplot([corrFCSuperT(:,iBand,2) corrFCMidT(:,iBand,2) corrFCDeepT(:,iBand,2)],...
        {'Super_FC'; 'Middle_FC'; 'Deep_FC'});box off;
    ylabel('Correlation between hybrids and FC map'); ylim([-1 1]); yticks(-1:0.1:1);

    subplot(122);
    boxplot([squeeze(corrFCSuperT(:,iBand,1:2)) squeeze(corrFCDeepT(:,iBand,1:2))],...
        {'10/10 Super'; '6/6 Super';'10/10 Deep'; '6/6 Deep'}); box off; 
    ylabel('Correlation between hybrids and FC map'); ylim([-1 0.5]); yticks(-1:0.1:1);

    sgtitle(bandNames{iBand}); 
end

[~,~,statsType] = anova1([corrFCSuperT(:,4,2) corrFCMidT(:,4,2) corrFCDeepT(:,4,2)],...
    {'Super_FC'; 'Mid_FC'; 'Deep_FC'},'off');
[rType,~,~,gnamesType] = multcompare(statsType,"Display","off");
tblType = array2table(rType,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblType.("Group") = gnamesType(tblType.("Group"));
tblType.("Control Group") = gnamesType(tblType.("Control Group"));

%% Plot the distributions of hybrid map correlations with FC map for
% different frequencies for S/M/D layer compartments separately for sensory
% and motor/pre-motor areas and everything together

for iArea = 3% 1:3
    clear corrFCSuperArea corrFCDeepArea corrFCMidArea
    switch iArea
        case 1
            runFlag   = (smFlag == 'S') & goodRunsSpatial & ~singleChFlag;
            areaLabel = 'Sensory areas';
        case 2
            runFlag   = (smFlag == 'M') & goodRunsSpatial & ~singleChFlag;
            areaLabel = 'Motor/Pre-motor areas';
        case 3
            runFlag   =  goodRunsSpatial & ~singleChFlag;
            areaLabel = 'All sites';

    end

    corrFCSuperArea = mean(corrFCSuper(:,:,:,:,lagFrameRange),5,'omitnan');
    corrFCSuperArea = reshape(corrFCSuperArea,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]);
    corrFCSuperArea(nanRow,:,:) = [];
    corrFCSuperArea(~runFlag,:,:) = [];

    corrFCDeepArea  = mean(corrFCDeep(:,:,:,:,lagFrameRange),5,'omitnan');
    corrFCDeepArea  = reshape(corrFCDeepArea,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]);
    corrFCDeepArea(nanRow,:,:) = [];
    corrFCDeepArea(~runFlag,:,:) = [];

    corrFCMidArea   = mean(corrFCMid(:,:,:,:,lagFrameRange),5,'omitnan');
    corrFCMidArea   = reshape(corrFCMidArea,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]);
    corrFCMidArea(nanRow,:,:) = [];
    corrFCMidArea(~runFlag,:,:) = [];

    figure;
    for iLayer = 1:3
        clear fcCorr fcTitle fcCorrT
        switch iLayer
            case 1            
                fcCorr  = corrFCSuperArea(:,:,2);
                fcTitle = 'Superficial-FC';
            case 2             
                fcCorr  = corrFCMidArea(:,:,2);
               fcTitle = 'Middle-FC';
            case 3               
                fcCorr  = corrFCDeepArea(:,:,2);
                fcTitle = 'Deep-FC';
        end

        subplot(1,3,iLayer); boxplot(fcCorr,{'Theta';'Alpha';'Beta';'Gamma';'Spiking'});
        title([areaLabel ' ' fcTitle]); box off; ylim([-1 1]); yticks(-1:0.2:1);

        [pFC,~,statsFC] = anova1(fcCorr,{'Theta';'Alpha';'Beta';'Gamma';'Spiking'});
        [rFC,~,~,gnamesFC] = multcompare(statsFC,"Display","off");
        tblFC = array2table(rFC,"VariableNames",["Group","Control Group","Lower Limit",...
            "Difference","Upper limit","p-val"]);
        tblFC.("Group") = gnamesFC(tblFC.("Group"));
        tblFC.("Control Group") = gnamesFC(tblFC.("Control Group"));

    end
end


%% Compile cross-frequency correlations and plot the correlations
clear superHybridAllBandsT midHybridAllBandsT deepHybridAllBandsT corrFCSuperArea corrFCDeepArea corrFCMidArea
smFlagNew = smFlag; smFlagNew(~(goodRunsSpatial & ~singleChFlag)) = [];

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

corrFCSuperArea = mean(corrFCSuper(:,:,:,:,lagFrameRange),5,'omitnan');
corrFCSuperArea = reshape(corrFCSuperArea,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]);
corrFCSuperArea(nanRow,:,:) = [];
corrFCSuperArea(~(goodRunsSpatial & ~singleChFlag),:,:) = [];

corrFCMidArea   = mean(corrFCMid(:,:,:,:,lagFrameRange),5,'omitnan');
corrFCMidArea   = reshape(corrFCMidArea,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]);
corrFCMidArea(nanRow,:,:) = [];
corrFCMidArea(~(goodRunsSpatial & ~singleChFlag),:,:) = [];

corrFCDeepArea  = mean(corrFCDeep(:,:,:,:,lagFrameRange),5,'omitnan');
corrFCDeepArea  = reshape(corrFCDeepArea,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]);
corrFCDeepArea(nanRow,:,:) = [];
corrFCDeepArea(~(goodRunsSpatial & ~singleChFlag),:,:) = [];

sH = reshape(squeeze(superHybridAllBandsT(:,2,:,:)),[size(superHybridAllBandsT,1) 25]);
mH = reshape(squeeze(midHybridAllBandsT(:,2,:,:)),[size(superHybridAllBandsT,1) 25]); 
dH = reshape(squeeze(deepHybridAllBandsT(:,2,:,:)),[size(superHybridAllBandsT,1) 25]); 


fcS = squeeze(corrFCSuperArea(:,:,2)); 
fcM = squeeze(corrFCMidArea(:,:,2));
fcD = squeeze(corrFCDeepArea(:,:,2));

for iArea = 1:3
    clear super mid deep fcSuper fcMid fcDeep
    switch iArea
        case 1
            super   = superHybridAllBandsT((smFlagNew == 'S'),:,:,:);
            mid     = midHybridAllBandsT((smFlagNew == 'S'),:,:,:);
            deep    = deepHybridAllBandsT((smFlagNew == 'S'),:,:,:);
            fcSuper = corrFCSuperArea((smFlagNew == 'S'),:,:);
            fcMid   = corrFCMidArea((smFlagNew == 'S'),:,:);
            fcDeep  = corrFCDeepArea((smFlagNew == 'S'),:,:);
            areaLabel = 'Sensory areas';
        case 2
            super     = superHybridAllBandsT((smFlagNew == 'M'),:,:,:);
            mid       = midHybridAllBandsT((smFlagNew == 'M'),:,:,:);
            deep      = deepHybridAllBandsT((smFlagNew == 'M'),:,:,:);  
            fcSuper   = corrFCSuperArea((smFlagNew == 'M'),:,:);
            fcMid     = corrFCMidArea((smFlagNew == 'M'),:,:);
            fcDeep    = corrFCDeepArea((smFlagNew == 'M'),:,:);
            areaLabel = 'Motor/Pre-motor areas';
        case 3
            super     = superHybridAllBandsT;
            mid       = midHybridAllBandsT;
            deep      = deepHybridAllBandsT;
            fcSuper   = corrFCSuperArea;
            fcMid     = corrFCMidArea;
            fcDeep    = corrFCDeepArea;
            areaLabel = 'All areas';            
    end

    figure; subplot(131); % All frequencies for superficial compartment
    imagesc(squeeze(median(squeeze(super(:,2,:,:)),1,'omitnan')));
    xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels');
    caxis([0 1]);colormap jet;colorbar; title('Superficial'); axis square;

    subplot(132); % All frequencies for middle compartment
    imagesc(squeeze(median(squeeze(mid(:,2,:,:)),1,'omitnan')));
    xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels');
    caxis([0 1]);colormap jet;colorbar;title('Middle');axis square;

    subplot(133); % All frequencies for deep compartment
    imagesc(squeeze(median(squeeze(deep(:,2,:,:)),1,'omitnan')));
    xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels');
    caxis([0 1]);colormap jet;colorbar; title('Deep');axis square;

    sgtitle(areaLabel);
   
    % Show the FC map
    figure; subplot(131); imagesc(median(fcSuper(:,:,2),1,'omitnan')'); axis square image; xticks();
    caxis([-0.5 0]);colormap(flipud(jet));colorbar; yticks(1:5); yticklabels(bandLabels');title('Superficial');
    
    subplot(132); imagesc(median(fcMid(:,:,2),1,'omitnan')'); axis square image; xticks();
    caxis([-0.5 0]);colormap(flipud(jet));colorbar;yticks(1:5);  yticklabels(bandLabels');title('Middle')
    
    subplot(133); imagesc(median(fcDeep(:,:,2),1,'omitnan')'); axis square image; xticks();
    caxis([-0.5 0]);colormap(flipud(jet));colorbar; yticks(1:5); yticklabels(bandLabels');title('Deep');
    
    sgtitle(areaLabel);   

    % Show FC map v2
    figure; clear fcAll
    fcAll = [median(fcSuper(:,:,2),1,'omitnan') ; median(fcMid(:,:,2),1,'omitnan'); median(fcDeep(:,:,2),1,'omitnan')] ;
    imagesc(fcAll); axis image square; xticks(1:5); xticklabels(bandLabels);
    yticks(1:3); yticklabels({'Superficial'; 'Middle' ; 'Deep'}); 
    caxis([-0.5 0]);colormap(flipud(jet));colorbar; title([areaLabel ' - Hybrids vs FC']);
end

%% Controls - Disambiguating channel position vs channel count to solidify the results
posControlVsRunCorr = NaN(size(processedDat,1),size(processedDat,2),3,5);
superDeepPosCheck   = NaN(size(processedDat,1),size(processedDat,2),3, 401,5);
fs                  = 1e3;
for iDate = 1: size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);

    for iRun = 1:size(allRuns{iDate,1})
        clear runName dataDir
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

        clc; disp([monkeyName ' '  expDate ' run: ' runName]);

        if ~exist([dataDir '\superDeepControls.mat'],'file')
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

            if corrMaskFlag
                corrMask = imresize(corrMask,1/3); % Resize the mask
                corrMask = corrMask(:,:,1)>0;
                corrMask = corrMask & ~elecMask;
            else
                corrMask = clipMaskCortex;
            end

            corrMaskT = reshape(corrMask,[imSize(1)*imSize(2) 1]);

            % Obtaining ephys data
            clear probeCh timeStamp badTimes badTimesThIm chInCortex szLFP...
                szIm szMin badTimeIm badTimes10Hz
            probeCh      = probe{iRun,iDate}.probeCh;
            timeStamp    = probe{iRun,iDate}.timeStamp;
            badTimes     = badTimesLFP{iDate,iRun};
            badTimesThIm = badTimeThresh{iDate,iRun};
            chInCortex   = estChInCortex{1,iDate}(iRun,:);

            % Removing bad channels
            probeCh(:,badCh{iDate,iRun}) = [];
            szLFP   = size(probeCh,1);

            if chInCortex(1)-chInCortex(2) ~= 0

                % Get imaging data
                clear pDatTemp imSize processedDat10
                pDatTemp = processedDat{iDate,iRun}.tempBandPass;
                imSize   = size(pDatTemp);
                pDatTemp = reshape(pDatTemp,[imSize(1)*imSize(2) imSize(3)]);

                % Upsampling imaging data to 10 Hz
                disp('Upsampling imaging data to 10 Hz... ');
                parfor iP = 1:size(pDatTemp,1)
                    processedDat10(iP,:) = interp(pDatTemp(iP,:),5);
                end

                szIm = size(processedDat10,2)*100;

                if ~(szLFP == szIm)
                    szMin          = min([szLFP, szIm]);
                    probeCh        = probeCh(1:szMin,:);
                    processedDat10 = processedDat10(:,1:floor(szMin/100));

                    badTimes(badTimes>szMin)         = [];
                    badTimesThIm(badTimesThIm>szMin) = [];
                else
                    szMin = szLFP;
                end

                timeStampSorted = timeStamp- timeStamp(1);
                badTimes10Hz    = unique(badTimesThIm./1000);
                badTimeIm       = [];

                for iT = 1: length(badTimes10Hz)
                    badTimeIm(iT) = find((floor(abs(timeStampSorted - badTimes10Hz(iT))*100)./100)<=0.05,1,'first');
                end

                badTimeIm = unique(badTimeIm);
                badTimeIm(badTimeIm>size(processedDat10,2)) = [];
                processedDat10(:,badTimeIm) = [];
                probeCh(badTimes,:) = []; % Remove bad times from LFP

                % Remove bad times determined visually from spectrogram
                [probeCh,~,processedDat10] = removeBadTimesFromSpec(monkeyName,expDate,runName,probeCh,...
                    [],processedDat10);

                clear ch
                ch = chInCortex(1):chInCortex(2);

                % 1. Pick 6 random channels 5 times and get the mean hybrid maps
                % Pick 5 iterations of groups of 6 channels - Check if position of
                % channels affect the super/deep hybrids
                clear superDeepPosCheckTemp posControlVsRunCorrTemp
                for iBand = 1:3
                    for iIter = 1:5
                        clear superCh deepCh
                        rng("shuffle"); superCh = ch(randperm(6,6));
                        rng("shuffle"); deepCh  = ch(randperm(6,6));

                        % Get infraslow ephys of gamma bands
                        clear bandVals ephysSuper ephysDeep
                        switch iBand
                            case 1
                                bandVals = [8 12];
                            case 2
                                bandVals = [13 30];
                            case 3
                                bandVals = [30 90]; 
                        end
                        [bT,aT]   = butter(3,bandVals./(1e3/2),'bandpass');  % Gamma band filtering parameters

                        ephysSuper = abs(single(filtfilt(bT,aT,double(probeCh(:,superCh)))));
                        ephysDeep  = abs(single(filtfilt(bT,aT,double(probeCh(:,deepCh)))));

                        envelopeDatSuper = envelope(ephysSuper,5);
                        envelopeDatDeep  = envelope(ephysDeep,5);

                        % Get filter parameters...
                        clear z p k sos g enSize                      

                        % Bandpass - 0.01 Hz - 0.1 Hz
                        [z,p,k] = butter(3,[0.01 0.1]./(fs/2),'bandpass');
                        [sos,g] = zp2sos(z,p,k);
                        enSize  = size(envelopeDatSuper);

                        envelopeFiltSuper = filtfilt(sos,g,double([envelopeDatSuper; envelopeDatSuper; envelopeDatSuper]));
                        envelopeFiltSuper = envelopeFiltSuper(enSize(1)+1:(end-enSize(1)),:);

                        envelopeFiltDeep = filtfilt(sos,g,double([envelopeDatDeep; envelopeDatDeep; envelopeDatDeep]));
                        envelopeFiltDeep = envelopeFiltDeep(enSize(1)+1:(end-enSize(1)),:);

                        infraEphysSuper = mean(single(downsample(envelopeFiltSuper,100)),2,'omitnan');
                        infraEphysDeep  = mean(single(downsample(envelopeFiltDeep,100)),2,'omitnan');

                        % Check size of timeseries of both modalities
                        clear szIm szLFP processedDat10Temp
                        szIm  = size(processedDat10,2);
                        szLFP = size(infraEphysSuper,1);

                        if ~(szLFP == szIm)
                            szMin               = min([szLFP, szIm]);
                            infraEphysSuper     = infraEphysSuper(1:szMin,:);
                            infraEphysDeep      = infraEphysDeep(1:szMin,:);
                            processedDat10Temp  = processedDat10(:,1:szMin);
                        else
                            processedDat10Temp  = processedDat10;
                        end

                        tic;
                        disp('Performing cross correlations...');

                        parfor iP = 1:size(processedDat10Temp,1)
                            [ccFullSuper(:,iP),~] = xcorr(infraEphysSuper',processedDat10Temp(iP,:),200,'normalized');
                            [ccFullDeep(:,iP),~]  = xcorr(infraEphysDeep',processedDat10Temp(iP,:),200,'normalized');
                        end

                        disp('Correlating superficial and deep hybrids for all lags...');

                        parfor iMap = 1:size(ccFullSuper,1)
                            superDeepPosCheckTemp(iBand,iMap,iIter) = corr(ccFullSuper(iMap,:)',ccFullDeep(iMap,:)','rows','complete');
                        end
                        toc;
                    end

                    runSuperDeepCorr(:,iBand)        = squeeze(super_DeepCorr(iDate,iRun,iBand,2,:));
                    posControlVsRunCorrTemp(iBand,:) = corr(runSuperDeepCorr(:,iBand) ,squeeze(superDeepPosCheckTemp(iBand,:,:)),'rows','complete');
                end
               
                save([dataDir '\superDeepControls.mat'],'superDeepPosCheckTemp','posControlVsRunCorrTemp');

                posControlVsRunCorr(iDate,iRun,:,:) = posControlVsRunCorrTemp;
                superDeepPosCheck(iDate,iRun,:,:,:) = superDeepPosCheckTemp;
            end
            
        else
            clear controlVars
            controlVars                       = load([dataDir '\superDeepControls.mat']);
            posControlVsRunCorr(iDate,iRun,:,:) = controlVars.posControlVsRunCorrTemp;
            superDeepPosCheck(iDate,iRun,:,:,:) = controlVars.superDeepPosCheckTemp;
            runSuperDeepCorr                    = squeeze(super_DeepCorr(iDate,iRun,:,2,:))';
        end

        if ~exist([dataDir '\Super_DeepPosControls.png'],'file') && (chInCortex(1)-chInCortex(2) ~= 0)
            figure; plot(-200:200,runSuperDeepCorr(:,3),'r','LineWidth',1); hold on;
            plot(-200:200,squeeze(superDeepPosCheck(iDate,iRun,3,:,:)),'k'); ylim([0 1.2]);
            xlabel('Lags (s)'); ylabel('Correlation between superficial and deep hybrids');
            f = gcf; exportgraphics(f,[dataDir '\Super_DeepPosControls.png'],'Resolution',300); close gcf;
        end
    end
end

%% Compiling the data and plotting S/D, S/M, M/D with Ideal condition
x = -200:200; % Lag distribution
lagROI = NaN(size(super_DeepCorr,[1 2 5]));
lagFOV = NaN(size(super_DeepCorr,[1 2 5]));
for iDate = 1: size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);
    for iRun = 1:size(allRuns{iDate,1})
        clear runName
        runName = allRuns{iDate,1}(iRun,:);
        clc; disp([monkeyName ' '  expDate ' run: ' runName]);
        lagROI(iDate,iRun,:) = tempProfileNoRef(iM,iDate,iRun).lagLow;
        lagROI(iDate,iRun,:) = (x == lagROI(iDate,iRun)); % lag from ROI - all channels
        lagFOV(iDate,iRun,:) = peakNegTimesAll(iDate,iRun,1)*10;
        lagFOV(iDate,iRun,:) = (x == lagFOV(iDate,iRun)); % lag from FOV - all channels
    end
end
lagROI = reshape(lagROI,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,5)]);
lagFOV = reshape(lagFOV,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,5)]);
lagROI(isnan(lagROI(:,1)),:) = []; lagROI(~goodRunsSpatial,:) = []; lagROI(singleChFlag,:) = [];
lagFOV(isnan(lagFOV(:,1)),:) = []; lagFOV(~goodRunsSpatial,:) = []; lagFOV(singleChFlag,:) = [];

superDeepPosCheckT = reshape(superDeepPosCheck,[size(super_DeepCorr,1)*size(super_DeepCorr,2) 3 401 5]); 
superDeepPosCheckT(nanRow,:,:,:)           = []; 
superDeepPosCheckT(~goodRunsSpatial,:,:,:) = []; 
superDeepPosCheckT(singleChFlag,:,:,:)     = []; 
superDeepPosCheckT                         = squeeze(median(superDeepPosCheckT,4,'omitnan'));
superDeepPosCheckT                         = permute(superDeepPosCheckT,[2 1 3]);
superDeepPosCheckT                         = superDeepPosCheckT(:,logical(lagFOV))'; 

for iBand = 1:3
    figure;
    for  iPlot = 1:4
        subplot(2,2,iPlot);
        clear tempVar;
        tempVar = [super_DeepAvgFramesT(:,iBand+1,iPlot) super_MidAvgFramesT(:,iBand+1,iPlot) ...
            deep_MidAvgFramesT(:,iBand+1,iPlot) superDeepPosCheckT(:,iBand) ];

        boxplot(tempVar,{'Super_Deep'; 'Super_Mid';'Mid_Deep'; 'Randomized channels'});
        ylabel('Correlation between hybrids'); ylim([-0.7 1.1]); yticks(-0.7:0.1:1);

        switch iPlot
            case 1
                title('No Ref - 10/10 channel split');
            case 2
                title('No Ref - 6/6 channel split');
            case 3
                title('Bipolar Ref - 6/6 channel split');
            case 4
                title('Avg Ref - 6/6 channel split');
        end
        axis square; box off;
    end
    sgtitle(bandNames{iBand+1});
end


% [pType,~,statsType] = anova1([allSuperDeep allSuperMid allDeepMid allCtrl],{'Super_Deep'; 'Super_Mid';'Mid_Deep'; 'Randomized channels'});
[pType,~,statsType] = anova1([super_DeepAvgFramesT(:,4,2) super_MidAvgFramesT(:,4,2) deep_MidAvgFramesT(:,4,2) ...
    superDeepPosCheckT(:,3)],...
    {'Super_Deep'; 'Super_Mid'; 'Deep_Mid'; 'Randomized channels'});
[rType,~,~,gnamesType] = multcompare(statsType,'Alpha',0.008);


tblType = array2table(rType,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblType.("Group") = gnamesType(tblType.("Group"));
tblType.("Control Group") = gnamesType(tblType.("Control Group"));

%% Correlate between superficial/middle/deep layer compartments
% Initialize variables
lfpLayerCorr = NaN(size(probe,2),size(probe,1),3,5);
infraLayerCorr = NaN(size(probe,2),size(probe,1),3,5);
powerLayerCorr = NaN(size(probe,2),size(probe,1),3,5);

superInfraAllBands = NaN(size(probe,2),size(probe,1),5,5);
deepInfraAllBands  = NaN(size(probe,2),size(probe,1),5,5);
midInfraAllBands   = NaN(size(probe,2),size(probe,1),5,5);
superPowerAllBands = NaN(size(probe,2),size(probe,1),5,5);
deepPowerAllBands  = NaN(size(probe,2),size(probe,1),5,5);
midPowerAllBands   = NaN(size(probe,2),size(probe,1),5,5);

crossFreqCrossLayerInfra = NaN(size(probe,2),size(probe,1),5,5,3,3);
crossFreqCrossLayerPower = NaN(size(probe,2),size(probe,1),5,5,3,3);

% Get filter parameters...
fs = 1e3; 
gammaBand   = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); % Gamma band filtering parameters
alphaBand   = [8 12];  [bA,aA] = butter(3,alphaBand./(fs/2),'bandpass'); % Alpha band filtering parameters
betaBand   = [13 30]; [bB,aB] = butter(3,betaBand./(fs/2),'bandpass');  % Beta band filtering parameters
thetaBand  = [6 8];   [bT,aT] = butter(3,thetaBand./(fs/2),'bandpass'); % Theta band filtering parameters

tic;
for iDate = 1:size(allDates,1)
    clear expDate
    expDate = allDates(iDate,:);

    for iRun = 1:size(allRuns{iDate,1})
        clear runName dataDir probeCh rawCh lfpBadTimes lfpBadCh chInCortex 
      
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

        clc; disp(['Analyzing data for ' monkeyName ' '  expDate ' run: ' runName]);
        
        % Get run related variables
        probeCh     = probe{iRun,iDate}.probeCh;
        rawCh       = probe{iRun,iDate}.rawCh; 
        lfpBadTimes = badTimesLFP{iDate,iRun};
        lfpBadCh    = badCh{iDate,iRun};
        chInCortex  = estChInCortex{1,iDate}(iRun,:);
        
        probeCh(:,lfpBadCh)    = []; % Remove bad channels from LFP
        probeCh(lfpBadTimes,:) = []; % Remove bad times from LFP  

        rawCh(:,lfpBadCh)    = []; % Remove bad channels from raw
        rawCh(lfpBadTimes,:) = []; % Remove bad times from raw

        % Remove bad times determined visually from spectrogram
        [probeCh,rawCh,~] = removeBadTimesFromSpec(monkeyName,expDate,runName,probeCh,rawCh,[]);
      
        if size(probeCh,2)~=1
            % Get frequency band based information
            clear thetaLFP alphaLFP betaLFP gammaLFP infraEphysTheta infraEphysAlpha...
                infraEphysBeta infraEphysGamma infraEphysSpiking powerTheta powerAlpha...
                powerBeta powerGamma powerSpiking

            thetaLFP = single(filtfilt(bT,aT,double(probeCh(:,chInCortex(1):chInCortex(2)))));
            alphaLFP = single(filtfilt(bA,aA,double(probeCh(:,chInCortex(1):chInCortex(2)))));
            betaLFP  = single(filtfilt(bB,aB,double(probeCh(:,chInCortex(1):chInCortex(2)))));
            gammaLFP = single(filtfilt(bG,aG,double(probeCh(:,chInCortex(1):chInCortex(2)))));

            powerTheta   = envelope(abs(thetaLFP),5); 
            powerAlpha   = envelope(abs(alphaLFP),5); 
            powerBeta    = envelope(abs(betaLFP),5);
            powerGamma   = envelope(abs(gammaLFP),5); 
            powerSpiking = envelope(abs(rawCh(:,chInCortex(1):chInCortex(2))),5);

            infraEphysTheta   = getInfraSlowPowerLFP(probeCh,bT,aT,chInCortex);
            infraEphysAlpha   = getInfraSlowPowerLFP(probeCh,bA,aA,chInCortex);
            infraEphysBeta    = getInfraSlowPowerLFP(probeCh,bB,aB,chInCortex);
            infraEphysGamma   = getInfraSlowPowerLFP(probeCh,bG,aG,chInCortex);
            infraEphysSpiking = getInfraSlowPowerLFP(rawCh,[],[],chInCortex);

            probeCh  = single(probeCh(:,chInCortex(1):chInCortex(2)));
            rawCh    = single(rawCh(:,chInCortex(1):chInCortex(2)));

            % Separate into superficial/middle/deep layer compartments for the
            % different frequency bands
            clear lfpSuper lfpMid lfpMid thetaSuper thetaMid thetaDeep alphaSuper...
                alphaMid alphaDeep betaSuper betaMid betaDeep gammaSuper gammaMid...
                gammaDeep spikingSuper spikingMid spikingDeep infraThetaSuper...
                infraThetaMid infraThetaDeep infraAlphaSuper infraAlphaMid...
                infraAlphaDeep infraBetaSuper infraBetaMid infraBetaDeep...
                infraGammaSuper infraGammaMid infraGammaDeep infraSpikingSuper...
                infraSpikingMid infraSpikingDeep thetaPowerSuper thetaPowerMid...
                thetaPowerDeep alphaPowerSuper alphaPowerMid alphaPowerDeep...
                betaPowerSuper betaPowerMid betaPowerDeep gammaPowerSuper...
                gammaPowerMid gammaPowerDeep spikingPowerSuper spikingPowerMid...
                spikingPowerDeep

            thetaSuper = mean(thetaLFP(:,1:6),2,'omitnan');
            thetaMid   = mean(thetaLFP(:,7:(end-6)),2,'omitnan');
            thetaDeep  = mean(thetaLFP(:,end-6+1:end),2,'omitnan');

            alphaSuper = mean(alphaLFP(:,1:6),2,'omitnan');
            alphaMid   = mean(alphaLFP(:,7:(end-6)),2,'omitnan');
            alphaDeep  = mean(alphaLFP(:,end-6+1:end),2,'omitnan');

            betaSuper = mean(betaLFP(:,1:6),2,'omitnan');
            betaMid   = mean(betaLFP(:,7:(end-6)),2,'omitnan');
            betaDeep  = mean(betaLFP(:,end-6+1:end),2,'omitnan');

            gammaSuper = mean(gammaLFP(:,1:6),2,'omitnan');
            gammaMid   = mean(gammaLFP(:,7:(end-6)),2,'omitnan');
            gammaDeep  = mean(gammaLFP(:,end-6+1:end),2,'omitnan');

            spikingSuper = mean(rawCh(:,1:6),2,'omitnan');
            spikingMid   = mean(rawCh(:,7:(end-6)),2,'omitnan');
            spikingDeep  = mean(rawCh(:,end-6+1:end),2,'omitnan');

            thetaPowerSuper = mean(powerTheta(:,1:6),2,'omitnan');
            thetaPowerMid   = mean(powerTheta(:,7:(end-6)),2,'omitnan');
            thetaPowerDeep  = mean(powerTheta(:,end-6+1:end),2,'omitnan');

            alphaPowerSuper = mean(powerAlpha(:,1:6),2,'omitnan');
            alphaPowerMid   = mean(powerAlpha(:,7:(end-6)),2,'omitnan');
            alphaPowerDeep  = mean(powerAlpha(:,end-6+1:end),2,'omitnan');

            betaPowerSuper = mean(powerBeta(:,1:6),2,'omitnan');
            betaPowerMid   = mean(powerBeta(:,7:(end-6)),2,'omitnan');
            betaPowerDeep  = mean(powerBeta(:,end-6+1:end),2,'omitnan');

            gammaPowerSuper = mean(powerGamma(:,1:6),2,'omitnan');
            gammaPowerMid   = mean(powerGamma(:,7:(end-6)),2,'omitnan');
            gammaPowerDeep  = mean(powerGamma(:,end-6+1:end),2,'omitnan');

            spikingPowerSuper = mean(powerSpiking(:,1:6),2,'omitnan');
            spikingPowerMid   = mean(powerSpiking(:,7:(end-6)),2,'omitnan');
            spikingPowerDeep  = mean(powerSpiking(:,end-6+1:end),2,'omitnan');

            infraThetaSuper = mean(infraEphysTheta(:,1:6),2,'omitnan');
            infraThetaMid   = mean(infraEphysTheta(:,7:(end-6)),2,'omitnan');
            infraThetaDeep  = mean(infraEphysTheta(:,end-6+1:end),2,'omitnan');

            infraAlphaSuper = mean(infraEphysAlpha(:,1:6),2,'omitnan');
            infraAlphaMid   = mean(infraEphysAlpha(:,7:(end-6)),2,'omitnan');
            infraAlphaDeep  = mean(infraEphysAlpha(:,end-6+1:end),2,'omitnan');

            infraBetaSuper = mean(infraEphysBeta(:,1:6),2,'omitnan');
            infraBetaMid   = mean(infraEphysBeta(:,7:(end-6)),2,'omitnan');
            infraBetaDeep  = mean(infraEphysBeta(:,end-6+1:end),2,'omitnan');

            infraGammaSuper = mean(infraEphysGamma(:,1:6),2,'omitnan');
            infraGammaMid   = mean(infraEphysGamma(:,7:(end-6)),2,'omitnan');
            infraGammaDeep  = mean(infraEphysGamma(:,end-6+1:end),2,'omitnan');

            infraSpikingSuper = mean(infraEphysSpiking(:,1:6),2,'omitnan');
            infraSpikingMid   = mean(infraEphysSpiking(:,7:(end-6)),2,'omitnan');
            infraSpikingDeep  = mean(infraEphysSpiking(:,end-6+1:end),2,'omitnan');

            infraThetaSuperAllCh = (infraEphysTheta(:,1:6));
            infraThetaMidAllCh   = (infraEphysTheta(:,7:(end-6)));
            infraThetaDeepAllCh  = (infraEphysTheta(:,end-6+1:end));

            infraAlphaSuperAllCh = (infraEphysAlpha(:,1:6));
            infraAlphaMidAllCh   = (infraEphysAlpha(:,7:(end-6)));
            infraAlphaDeepAllCh  = (infraEphysAlpha(:,end-6+1:end));

            infraBetaSuperAllCh = (infraEphysBeta(:,1:6));
            infraBetaMidAllCh   = (infraEphysBeta(:,7:(end-6)));
            infraBetaDeepAllCh  = (infraEphysBeta(:,end-6+1:end));

            infraGammaSuperAllCh = (infraEphysGamma(:,1:6));
            infraGammaMidAllCh   = (infraEphysGamma(:,7:(end-6)));
            infraGammaDeepAllCh  = (infraEphysGamma(:,end-6+1:end));

            infraSpikingSuperAllCh = (infraEphysSpiking(:,1:6));
            infraSpikingMidAllCh   = (infraEphysSpiking(:,7:(end-6)));
            infraSpikingDeepAllCh  = (infraEphysSpiking(:,end-6+1:end));


            % Correlate between superficial/middle/deep layers
            % Time-series
            lfpLayerCorr(iDate,iRun,1,1) = corr(thetaSuper,thetaDeep,'rows','complete');
            lfpLayerCorr(iDate,iRun,2,1) = corr(thetaSuper,thetaMid,'rows','complete');
            lfpLayerCorr(iDate,iRun,3,1) = corr(thetaMid,thetaDeep,'rows','complete');

            lfpLayerCorr(iDate,iRun,1,2) = corr(alphaSuper,alphaDeep,'rows','complete');
            lfpLayerCorr(iDate,iRun,2,2) = corr(alphaSuper,alphaMid,'rows','complete');
            lfpLayerCorr(iDate,iRun,3,2) = corr(alphaMid,alphaDeep,'rows','complete');

            lfpLayerCorr(iDate,iRun,1,3) = corr(betaSuper,betaDeep,'rows','complete');
            lfpLayerCorr(iDate,iRun,2,3) = corr(betaSuper,betaMid,'rows','complete');
            lfpLayerCorr(iDate,iRun,3,3) = corr(betaMid,betaDeep,'rows','complete');

            lfpLayerCorr(iDate,iRun,1,4) = corr(gammaSuper,gammaDeep,'rows','complete');
            lfpLayerCorr(iDate,iRun,2,4) = corr(gammaSuper,gammaMid,'rows','complete');
            lfpLayerCorr(iDate,iRun,3,4) = corr(gammaMid,gammaDeep,'rows','complete');

            lfpLayerCorr(iDate,iRun,1,5) = corr(spikingSuper,spikingDeep,'rows','complete');
            lfpLayerCorr(iDate,iRun,2,5) = corr(spikingSuper,spikingMid,'rows','complete');
            lfpLayerCorr(iDate,iRun,3,5) = corr(spikingMid,spikingDeep,'rows','complete');

            [gammaLayerXCorr(iDate,iRun,:,1),lags] = xcorr(detrend(gammaSuper,0),detrend(gammaDeep,0),1000,'coeff');
            [gammaLayerXCorr(iDate,iRun,:,2),~]    = xcorr(detrend(gammaSuper,0),detrend(gammaMid,0),1000,'coeff');
            [gammaLayerXCorr(iDate,iRun,:,3),~]    = xcorr(detrend(gammaMid,0),detrend(gammaDeep,0),1000,'coeff');
%              
%             clear superH midH deepH
%             superH = hilbert(gammaSuper); midH = hilbert(gammaMid); deepH = hilbert(gammaDeep); 
%             phaseRad{iDate,iRun}(1) = angle(superH./deepH);
%             phaseRad{iDate,iRun}(2) = angle(superH./midH);
%             phaseRad{iDate,iRun}(3)= angle(midH./deepH); 

            % Power
            powerLayerCorr(iDate,iRun,1,1) = corr(thetaPowerSuper,thetaPowerDeep,'rows','complete');
            powerLayerCorr(iDate,iRun,2,1) = corr(thetaPowerSuper,thetaPowerMid,'rows','complete');
            powerLayerCorr(iDate,iRun,3,1) = corr(thetaPowerMid,thetaPowerDeep,'rows','complete');

            powerLayerCorr(iDate,iRun,1,2) = corr(alphaPowerSuper,alphaPowerDeep,'rows','complete');
            powerLayerCorr(iDate,iRun,2,2) = corr(alphaPowerSuper,alphaPowerMid,'rows','complete');
            powerLayerCorr(iDate,iRun,3,2) = corr(alphaPowerMid,alphaPowerDeep,'rows','complete');

            powerLayerCorr(iDate,iRun,1,3) = corr(betaPowerSuper,betaPowerDeep,'rows','complete');
            powerLayerCorr(iDate,iRun,2,3) = corr(betaPowerSuper,betaPowerMid,'rows','complete');
            powerLayerCorr(iDate,iRun,3,3) = corr(betaPowerMid,betaPowerDeep,'rows','complete');

            powerLayerCorr(iDate,iRun,1,4) = corr(gammaPowerSuper,gammaPowerDeep,'rows','complete');
            powerLayerCorr(iDate,iRun,2,4) = corr(gammaPowerSuper,gammaPowerMid,'rows','complete');
            powerLayerCorr(iDate,iRun,3,4) = corr(gammaPowerMid,gammaPowerDeep,'rows','complete');

            powerLayerCorr(iDate,iRun,1,5) = corr(spikingPowerSuper,spikingPowerDeep,'rows','complete');
            powerLayerCorr(iDate,iRun,2,5) = corr(spikingPowerSuper,spikingPowerMid,'rows','complete');
            powerLayerCorr(iDate,iRun,3,5) = corr(spikingPowerMid,spikingPowerDeep,'rows','complete');

            [gammaPowerLayerXCorr(iDate,iRun,:,1),~] = xcorr(detrend(gammaPowerSuper,0),detrend(gammaPowerDeep,0),1000,'coeff');
            [gammaPowerLayerXCorr(iDate,iRun,:,2),~] = xcorr(detrend(gammaPowerSuper,0),detrend(gammaPowerMid,0),1000,'coeff');
            [gammaPowerLayerXCorr(iDate,iRun,:,3),~] = xcorr(detrend(gammaPowerMid,0),detrend(gammaPowerDeep,0),1000,'coeff');
            
%             clear superH midH deepH
%             superH = hilbert(gammaPowerSuper); midH = hilbert(gammaPowerMid); deepH = hilbert(gammaPowerDeep);
%             phaseRadPower(iDate,iRun,:,1) = angle(superH./deepH);
%             phaseRadPower(iDate,iRun,:,2) = angle(superH./midH);
%             phaseRadPower(iDate,iRun,:,3) = angle(midH./deepH);

            % Infraslow power
            infraLayerCorr(iDate,iRun,1,1) = corr(infraThetaSuper,infraThetaDeep,'rows','complete');
            infraLayerCorr(iDate,iRun,2,1) = corr(infraThetaSuper,infraThetaMid,'rows','complete');
            infraLayerCorr(iDate,iRun,3,1) = corr(infraThetaMid,infraThetaDeep,'rows','complete');

            infraLayerCorr(iDate,iRun,1,2) = corr(infraAlphaSuper,infraAlphaDeep,'rows','complete');
            infraLayerCorr(iDate,iRun,2,2) = corr(infraAlphaSuper,infraAlphaMid,'rows','complete');
            infraLayerCorr(iDate,iRun,3,2) = corr(infraAlphaMid,infraAlphaDeep,'rows','complete');

            infraLayerCorr(iDate,iRun,1,3) = corr(infraBetaSuper,infraBetaDeep,'rows','complete');
            infraLayerCorr(iDate,iRun,2,3) = corr(infraBetaSuper,infraBetaMid,'rows','complete');
            infraLayerCorr(iDate,iRun,3,3) = corr(infraBetaMid,infraBetaDeep,'rows','complete');

            infraLayerCorr(iDate,iRun,1,4) = corr(infraGammaSuper,infraGammaDeep,'rows','complete');
            infraLayerCorr(iDate,iRun,2,4) = corr(infraGammaSuper,infraGammaMid,'rows','complete');
            infraLayerCorr(iDate,iRun,3,4) = corr(infraGammaMid,infraGammaDeep,'rows','complete');

            infraLayerCorr(iDate,iRun,1,5) = corr(infraSpikingSuper,infraSpikingDeep,'rows','complete');
            infraLayerCorr(iDate,iRun,2,5) = corr(infraSpikingSuper,infraSpikingMid,'rows','complete');
            infraLayerCorr(iDate,iRun,3,5) = corr(infraSpikingMid,infraSpikingDeep,'rows','complete');

            [infraGammaLayerXCorr(iDate,iRun,:,1),~] = xcorr(detrend(infraGammaSuper,0),detrend(infraGammaDeep,0),1000,'coeff');
            [infraGammaLayerXCorr(iDate,iRun,:,2),~] = xcorr(detrend(infraGammaSuper,0),detrend(infraGammaMid,0),1000,'coeff');
            [infraGammaLayerXCorr(iDate,iRun,:,3),~] = xcorr(detrend(infraGammaMid,0),detrend(infraGammaDeep,0),1000,'coeff');

%             clear superH midH deepH
%             superH = hilbert(infraGammaSuper); midH = hilbert(infraGammaMid); deepH = hilbert(infraGammaDeep);
%             phaseRadInfra(iDate,iRun,:,1) = angle(superH./deepH);
%             phaseRadInfra(iDate,iRun,:,2) = angle(superH./midH);
%             phaseRadInfra(iDate,iRun,:,3) = angle(midH./deepH);

            % Correlate between frequencies given a layer compartment
            clear chLen
            chLen = estChInCortex{1,iDate}(iRun,2) - estChInCortex{1,iDate}(iRun,1);
            for iBand1 = 1:5
                clear super1 deep1 mid1 superPower1 midPower1 deepPower1
                switch iBand1
                    case 1
                        super1 = infraThetaSuperAllCh;
                        deep1  = infraThetaDeepAllCh;

                        superPower1 = thetaPowerSuper;
                        deepPower1  = thetaPowerDeep;

                        if  chLen~=0
                            mid1      = infraThetaMidAllCh;
                            midPower1 = thetaPowerMid;
                        end

                    case 2
                        super1 = infraAlphaSuperAllCh;
                        deep1  = infraAlphaDeepAllCh;

                        superPower1 = alphaPowerSuper;
                        deepPower1 = alphaPowerDeep;

                        if chLen~=0
                            mid1 = infraAlphaMidAllCh;
                            midPower1 = alphaPowerMid;
                        end

                    case 3
                        super1 = infraBetaSuperAllCh;
                        deep1  = infraBetaDeepAllCh;

                        superPower1 = betaPowerSuper;
                        deepPower1  = betaPowerDeep;

                        if chLen~=0
                            mid1      = infraBetaMidAllCh;
                            midPower1 = betaPowerMid;
                        end

                    case 4
                        super1 = infraGammaSuperAllCh;
                        deep1  = infraGammaDeepAllCh;

                        superPower1 = gammaPowerSuper;
                        deepPower1  = gammaPowerDeep;

                        if  chLen~=0
                            mid1      = infraGammaMidAllCh;
                            midPower1 = gammaPowerMid;
                        end

                    case 5
                        super1 = infraSpikingSuperAllCh;
                        deep1  = infraSpikingDeepAllCh;

                        superPower1 = spikingPowerSuper;
                        deepPower1  = spikingPowerDeep;

                        if  chLen~=0
                            mid1      = infraSpikingMidAllCh;
                            midPower1 = spikingPowerMid;
                        end
                end

                for iBand2 = 1:5
                    clear super2 deep2 mid2  superPower2 midPower2 deepPower2
                    switch iBand2
                        case 1
                            super2 = infraThetaSuperAllCh;
                            deep2  = infraThetaDeepAllCh;

                            superPower2 = thetaPowerSuper;
                            deepPower2  = thetaPowerDeep;

                            if  chLen~=0
                                mid2      = infraThetaMidAllCh;
                                midPower2 = thetaPowerMid;
                            end

                        case 2
                            super2 = infraAlphaSuperAllCh;
                            deep2  = infraAlphaDeepAllCh;

                            superPower2 = alphaPowerSuper;
                            deepPower2  = alphaPowerDeep;

                            if chLen~=0
                                mid2      = infraAlphaMidAllCh;
                                midPower2 = alphaPowerMid;
                            end

                        case 3
                            super2 = infraBetaSuperAllCh;
                            deep2  = infraBetaDeepAllCh;

                            superPower2 = betaPowerSuper;
                            deepPower2  = betaPowerDeep;

                            if chLen~=0
                                mid2      = infraBetaMidAllCh;
                                midPower2 = betaPowerMid;
                            end

                        case 4
                            super2 = infraGammaSuperAllCh;
                            deep2  = infraGammaDeepAllCh;

                            superPower2 = gammaPowerSuper;
                            deepPower2  = gammaPowerDeep;

                            if  chLen~=0
                                mid2      = infraGammaMidAllCh;
                                midPower2 = gammaPowerMid;
                            end

                        case 5
                            super2 = infraSpikingSuperAllCh;
                            deep2  = infraSpikingDeepAllCh;

                            superPower2 = spikingPowerSuper;
                            deepPower2  = spikingPowerDeep;

                            if chLen~=0
                                mid2      = infraSpikingMidAllCh;
                                midPower2 = spikingPowerMid;
                            end
                    end

                    superInfraAllBands(iDate,iRun,iBand1,iBand2,:,:) = ((corr(super1,super2,'rows','complete')));
                    deepInfraAllBands(iDate,iRun,iBand1,iBand2,:,:)  = ((corr(deep1,deep2,'rows','complete')));

                    superPowerAllBands(iDate,iRun,iBand1,iBand2) = corr(superPower1,superPower2,'rows','complete');
                    deepPowerAllBands(iDate,iRun,iBand1,iBand2)  = corr(deepPower1,deepPower2,'rows','complete');

                    if  chLen~=0
                        midInfraAllBands(iDate,iRun,iBand1,iBand2,:,:) = ((corr(mid1,mid2,'rows','complete')));
                        midPowerAllBands(iDate,iRun,iBand1,iBand2) = corr(midPower1,midPower2,'rows','complete');
                    end

%                     if iBand1~=iBand2
%                         crossFreqCrossLayerInfra(iDate,iRun,iBand1,iBand2,1,1) = corr(super1,super2,'rows','complete');                         
%                         crossFreqCrossLayerInfra(iDate,iRun,iBand1,iBand2,1,3) = corr(super1,deep2,'rows','complete');                        
%                         crossFreqCrossLayerInfra(iDate,iRun,iBand1,iBand2,3,1) = corr(deep1,super2,'rows','complete');                      
%                         crossFreqCrossLayerInfra(iDate,iRun,iBand1,iBand2,3,3) = corr(deep1,deep2,'rows','complete'); 
% 
%                         crossFreqCrossLayerPower(iDate,iRun,iBand1,iBand2,1,1) = corr(superPower1,superPower2,'rows','complete');
%                         crossFreqCrossLayerPower(iDate,iRun,iBand1,iBand2,1,3) = corr(superPower1,deepPower2,'rows','complete');
%                         crossFreqCrossLayerPower(iDate,iRun,iBand1,iBand2,3,1) = corr(deepPower1,superPower2,'rows','complete');
%                         crossFreqCrossLayerPower(iDate,iRun,iBand1,iBand2,3,3) = corr(deepPower1,deepPower2,'rows','complete');
% 
%                          if  chLen~=0
%                              crossFreqCrossLayerInfra(iDate,iRun,iBand1,iBand2,1,2) = corr(super1,mid2,'rows','complete');
%                              crossFreqCrossLayerInfra(iDate,iRun,iBand1,iBand2,2,1) = corr(mid1,super2,'rows','complete');
%                              crossFreqCrossLayerInfra(iDate,iRun,iBand1,iBand2,2,2) = corr(mid1,mid2,'rows','complete');
%                              crossFreqCrossLayerInfra(iDate,iRun,iBand1,iBand2,2,3) = corr(mid1,deep2,'rows','complete');
%                              crossFreqCrossLayerInfra(iDate,iRun,iBand1,iBand2,3,2) = corr(deep1,mid2,'rows','complete');
% 
%                              crossFreqCrossLayerPower(iDate,iRun,iBand1,iBand2,1,2) = corr(superPower1,midPower2,'rows','complete');
%                              crossFreqCrossLayerPower(iDate,iRun,iBand1,iBand2,2,1) = corr(midPower1,superPower2,'rows','complete');
%                              crossFreqCrossLayerPower(iDate,iRun,iBand1,iBand2,2,2) = corr(midPower1,midPower2,'rows','complete');
%                              crossFreqCrossLayerPower(iDate,iRun,iBand1,iBand2,2,3) = corr(midPower1,deepPower2,'rows','complete');
%                              crossFreqCrossLayerPower(iDate,iRun,iBand1,iBand2,3,2) = corr(deepPower1,midPower2,'rows','complete');
%                          end
%                     end 
                end
            end
          
          
        else
            lfpLayerCorr(iDate,iRun,:,:)   = 0;
            infraLayerCorr(iDate,iRun,:,:) = 0;
            powerLayerCorr(iDate,iRun,:,:) = 0;

            superInfraAllBands(iDate,iRun,:,:) = 0;
            deepInfraAllBands(iDate,iRun,:,:)  = 0;
            midInfraAllBands(iDate,iRun,:,:)   = 0;
            superPowerAllBands(iDate,iRun,:,:) = 0;
            deepPowerAllBands(iDate,iRun,:,:)  = 0;
            midPowerAllBands(iDate,iRun,:,:)   = 0;
% 
%             crossFreqCrossLayerInfra(iDate,iRun,:,:,:,:) = 0;
%             crossFreqCrossLayerPower(iDate,iRun,:,:,:,:) = 0;
            
        end
    end
end
toc;

%% Compile and plot 
lfpLayerCorrT   = reshape(lfpLayerCorr,[size(probe,2)*size(probe,1) 3 5]);
infraLayerCorrT = reshape(infraLayerCorr,[size(probe,2)*size(probe,1) 3 5]);
powerLayerCorrT = reshape(powerLayerCorr,[size(probe,2)*size(probe,1) 3 5]);

superInfraAllBandsT = reshape(superInfraAllBands,[size(probe,2)*size(probe,1) 5 5]);
deepInfraAllBandsT  = reshape(deepInfraAllBands,[size(probe,2)*size(probe,1) 5 5]);
midInfraAllBandsT   = reshape(midInfraAllBands,[size(probe,2)*size(probe,1) 5 5]);
superPowerAllBandsT = reshape(superPowerAllBands,[size(probe,2)*size(probe,1) 5 5]);
deepPowerAllBandsT  = reshape(deepPowerAllBands,[size(probe,2)*size(probe,1) 5 5]);
midPowerAllBandsT   = reshape(midPowerAllBands,[size(probe,2)*size(probe,1) 5 5]);

gammaLayerXCorrT      = reshape(gammaLayerXCorr,[size(probe,2)*size(probe,1)  2001 3]);
gammaPowerLayerXCorrT = reshape(gammaPowerLayerXCorr,[size(probe,2)*size(probe,1) 2001 3]);
infraGammaLayerXCorrT = reshape(infraGammaLayerXCorr,[size(probe,2)*size(probe,1)  2001 3]);
% 
% crossFreqCrossLayerInfraT = reshape(crossFreqCrossLayerInfra,[size(probe,2)*size(probe,1) 5 5 3 3]);
% crossFreqCrossLayerPowerT = reshape(crossFreqCrossLayerPower,[size(probe,2)*size(probe,1) 5 5 3 3]);

nanLocs = (isnan(lfpLayerCorrT(:,1,1)));
lfpLayerCorrT(nanLocs,:,:)         = [];
infraLayerCorrT(nanLocs,:,:)       = []; 
powerLayerCorrT(nanLocs,:,:)       = [];
gammaLayerXCorrT(nanLocs,:,:)      = [];
gammaPowerLayerXCorrT(nanLocs,:,:) = []; 
infraGammaLayerXCorrT(nanLocs,:,:) = [];
superInfraAllBandsT(nanLocs,:,:)   = [];
deepInfraAllBandsT(nanLocs,:,:)    = [];
midInfraAllBandsT(nanLocs,:,:)     = [];
superPowerAllBandsT(nanLocs,:,:)   = [];
deepPowerAllBandsT(nanLocs,:,:)    = [];
midPowerAllBandsT(nanLocs,:,:)     = [];
% crossFreqCrossLayerInfraT(nanLocs,:,:,:,:) = [];
% crossFreqCrossLayerPowerT(nanLocs,:,:,:,:) = [];

lfpLayerCorrT(~(goodRunsSpatial & ~singleChFlag),:,:)         = [];
infraLayerCorrT(~(goodRunsSpatial & ~singleChFlag),:,:)       = []; 
powerLayerCorrT(~(goodRunsSpatial & ~singleChFlag),:,:)       = [];
gammaLayerXCorrT(~(goodRunsSpatial & ~singleChFlag),:,:)      = [];
gammaPowerLayerXCorrT(~(goodRunsSpatial & ~singleChFlag),:,:) = []; 
infraGammaLayerXCorrT(~(goodRunsSpatial & ~singleChFlag),:,:) = [];
% crossFreqCrossLayerInfraT(~(goodRunsSpatial & ~singleChFlag),:,:,:,:) = [];
% crossFreqCrossLayerPowerT(~(goodRunsSpatial & ~singleChFlag),:,:,:,:) = [];

superInfraAllBandsT(~(goodRunsSpatial & ~singleChFlag),:,:)   = [];
deepInfraAllBandsT(~(goodRunsSpatial & ~singleChFlag),:,:)    = [];
midInfraAllBandsT(~(goodRunsSpatial & ~singleChFlag),:,:)     = [];
superPowerAllBandsT(~(goodRunsSpatial & ~singleChFlag),:,:)   = [];
deepPowerAllBandsT(~(goodRunsSpatial & ~singleChFlag),:,:)    = [];
midPowerAllBandsT(~(goodRunsSpatial & ~singleChFlag),:,:)     = [];

smFlagNew = smFlag; smFlagNew(~(goodRunsSpatial & ~singleChFlag)) = [];

% Check phase differences between the layer compartments
for iType = 1:3
    clear typeVar typeLabel
    switch iType
        case 1
            typeLabel = 'LFP';
            typeVar   = gammaLayerXCorrT;
        case 2
            typeLabel = 'Power';
            typeVar   = gammaPowerLayerXCorrT; 
        case 3
            typeLabel = 'Infraslow';
            typeVar   = infraGammaLayerXCorrT;  
    end

    figure;
    subplot(1,3,1); plot(lags,squeeze(typeVar(:,:,1))); title('S/D'); ylim([-0.75 1]); xticks(-1000:200:1000); box off;
    subplot(1,3,2); plot(lags,squeeze(typeVar(:,:,2))); title('S/M'); ylim([-0.75 1]); xticks(-1000:200:1000); box off;
    subplot(1,3,3); plot(lags,squeeze(typeVar(:,:,3))); title('M/D'); ylim([-0.75 1]); xticks(-1000:200:1000); box off;
    sgtitle(typeLabel); ylabel('Correlation'); xlabel('Lag (samples)');
end

% Plot the corelations between the layer compartments 
for iType = 1:3
    clear typeLabel typeVar
    switch iType
        case 1
            typeLabel = 'LFP';
            typeVar   = lfpLayerCorrT;
        case 2
            typeLabel = 'Power';
            typeVar   = powerLayerCorrT; 
        case 3
            typeLabel = 'Infraslow';
            typeVar   = infraLayerCorrT; 
    end
    
    figure;
    for iBand = 1:5
        switch iBand
            case 1
                bandLabel = 'Theta';
            case 2
                bandLabel = 'Alpha';
            case 3
                bandLabel = 'Beta';
            case 4
                bandLabel = 'Gamma';
            case 5
                bandLabel = 'Spiking';
        end

        subplot(2,3,iBand);
        boxplot(squeeze(typeVar(:,:,iBand)),'Labels',{'S/D'; 'S/M'; 'M/D'}); 
        title(bandLabel); ylim([-0.5 1]); yticks(-0.5:0.1:1); box off; clear s rTemp;
        [pAll(iType,iBand),~,s] = anova1(squeeze(typeVar(:,:,iBand)),{'S/D'; 'S/M'; 'M/D'},'off');
        rTemp = multcompare(s,"Display","off");
        rAll(:,iType,iBand) = squeeze(rTemp(:,6));
   
    end
    sgtitle([typeLabel ' Correlations']);
end

corrSD = [lfpLayerCorrT(:,1,4) powerLayerCorrT(:,1,4) infraLayerCorrT(:,1,4)];
[pCorrSD,~,statsCorrSD] = anova1(corrSD,{'Time-series';'Power';'InfraSlow power'},'off');
[rCorrSD,~,~,gnamesCorrSD] = multcompare(statsCorrSD,"CriticalValueType","bonferroni","Display","off");

tblCorrSD = array2table(rCorrSD,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblCorrSD.("Group") = gnamesCorrSD(tblCorrSD.("Group"));
tblCorrSD.("Control Group") = gnamesCorrSD(tblCorrSD.("Control Group"));


%% Plot cross-frequency power/infraslow power correlation
for iArea = 1:3
    clear infraSuper infraMid infraDeep powerSuper powerMid powerDeep areaLabel
    switch iArea
        case 1
            infraSuper = superInfraAllBandsT(smFlagNew == 'S',:,:);
            infraMid   = midInfraAllBandsT(smFlagNew == 'S',:,:);
            infraDeep  = deepInfraAllBandsT(smFlagNew == 'S',:,:);
            powerSuper = superPowerAllBandsT(smFlagNew == 'S',:,:);
            powerMid   = midPowerAllBandsT(smFlagNew == 'S',:,:);
            powerDeep  = deepPowerAllBandsT(smFlagNew == 'S',:,:);
            areaLabel  = 'Sensory';
        case 2
            infraSuper = superInfraAllBandsT(smFlagNew == 'M',:,:);
            infraMid   = midInfraAllBandsT(smFlagNew == 'M',:,:);
            infraDeep  = deepInfraAllBandsT(smFlagNew == 'M',:,:);
            powerSuper = superPowerAllBandsT(smFlagNew == 'M',:,:);
            powerMid   = midPowerAllBandsT(smFlagNew == 'M',:,:);
            powerDeep  = deepPowerAllBandsT(smFlagNew == 'M',:,:);
            areaLabel  = 'Motor/Pre-motor';
        case 3
            infraSuper = superInfraAllBandsT;
            infraMid   = midInfraAllBandsT;
            infraDeep  = deepInfraAllBandsT;
            powerSuper = superPowerAllBandsT;
            powerMid   = midPowerAllBandsT;
            powerDeep  = deepPowerAllBandsT;
            areaLabel  = 'All areas'; 
    end 
    for iSig = 1:2
        clear super deep mid sigLabel       
        switch iSig
            case 1
                super    = infraSuper;
                mid      = infraMid;
                deep     = infraDeep;
                sigLabel = 'Infraslow power';
                cLim     = [0 1];

            case 2
                super    = powerSuper;
                mid      = powerMid;
                deep     = powerDeep;
                sigLabel = 'Power';
                cLim     = [0 0.5];
        end
        figure;
        subplot(131); imagesc(squeeze(median(super,1,'omitnan'))); colormap jet;
        caxis(cLim); title('Superficial'); axis image; colorbar;
        xticks(1:5); yticks(1:5); xticklabels(bandLabels); yticklabels(bandLabels);

        subplot(132); imagesc(squeeze(median(mid,1,'omitnan'))); colormap jet;
        caxis(cLim); title('Middle');axis image; colorbar;
        xticks(1:5); yticks(1:5); xticklabels(bandLabels); yticklabels(bandLabels);

        subplot(133); imagesc(squeeze(median(deep,1,'omitnan'))); colormap jet;
        caxis(cLim); title('Deep');axis image; colorbar;
        xticks(1:5); yticks(1:5); xticklabels(bandLabels); yticklabels(bandLabels);
        sgtitle([areaLabel ' ' sigLabel]);

    end
end

%% Show the distribution of correlations
for  iArea = 3%1:3
    clear infraSuper infraMid infraDeep powerSuper powerMid powerDeep areaLabel
    switch iArea
        case 1
            infraSuper = superInfraAllBandsT(smFlagNew == 'S',:,:);
            infraMid   = midInfraAllBandsT(smFlagNew == 'S',:,:);
            infraDeep  = deepInfraAllBandsT(smFlagNew == 'S',:,:);
            powerSuper = superPowerAllBandsT(smFlagNew == 'S',:,:);
            powerMid   = midPowerAllBandsT(smFlagNew == 'S',:,:);
            powerDeep  = deepPowerAllBandsT(smFlagNew == 'S',:,:);
            areaLabel  = 'Sensory';
        case 2
            infraSuper = superInfraAllBandsT(smFlagNew == 'M',:,:);
            infraMid   = midInfraAllBandsT(smFlagNew == 'M',:,:);
            infraDeep  = deepInfraAllBandsT(smFlagNew == 'M',:,:);
            powerSuper = superPowerAllBandsT(smFlagNew == 'M',:,:);
            powerMid   = midPowerAllBandsT(smFlagNew == 'M',:,:);
            powerDeep  = deepPowerAllBandsT(smFlagNew == 'M',:,:);
            areaLabel  = 'Motor/Pre-motor';
        case 3
            infraSuper = superInfraAllBandsT;
            infraMid   = midInfraAllBandsT;
            infraDeep  = deepInfraAllBandsT;
            powerSuper = superPowerAllBandsT;
            powerMid   = midPowerAllBandsT;
            powerDeep  = deepPowerAllBandsT;
            areaLabel  = 'All areas';
    end

    for iBand = 1:5
        figure;
        subplot(231); boxplot(squeeze(infraSuper(:,iBand,:)),bandLabels); 
        ylim([-0.5 1.2]); title('Superficial - infraslow power'); box off; 
        subplot(232); boxplot(squeeze(infraMid(:,iBand,:)),bandLabels);  
        ylim([-0.5 1.2]); title('Middle - infraslow power'); box off; 
        subplot(233); boxplot(squeeze(infraDeep(:,iBand,:)),bandLabels);  
        ylim([-0.5 1.2]);title('Deep - infraslow power'); box off; 

        subplot(234); boxplot(squeeze(powerSuper(:,iBand,:)),bandLabels); 
        ylim([-0.5 1.2]); title('Superficial - power'); box off; 
        subplot(235); boxplot(squeeze(powerMid(:,iBand,:)),bandLabels);  
        ylim([-0.5 1.2]);title('Middle - power'); box off; 
        subplot(236); boxplot(squeeze(powerDeep(:,iBand,:)),bandLabels);  
        ylim([-0.5 1.2]);title('Deep - power'); box off; 
       
        sgtitle(['Correlation with ' bandLabels{iBand} ' for ' areaLabel]); 
    end
end


%% Show the cross-frequency and layer correlations 

crossFreqCrossLayerInfraT = reshape(crossFreqCrossLayerInfraT,[size(crossFreqCrossLayerPowerT,1) 25 3 3]); 
crossFreqCrossLayerPowerT = reshape(crossFreqCrossLayerPowerT,[size(crossFreqCrossLayerPowerT,1) 25 3 3]); 

goodLabels = [0 1 1 1 1; 0 0 1 1 1; 0 0 0 1 1; 0 0 0 0 1;0 0 0 0 0];
goodLabels = reshape(goodLabels,[25 1]); % To avoid repeats in correlations

crossFreqCrossLayerInfraT(:,~goodLabels,:,:) = []; 
crossFreqCrossLayerPowerT(:,~goodLabels,:,:) = []; 

freqCombs = {'Theta-alpha'; 'Theta-beta'; 'Alpha-beta'; 'Theta-gamma'; 'Alpha-gamma';...
    'Beta-gamma';'Theta-spiking';'Alpha-spiking'; 'Beta-spiking'; 'Gamma-spiking'}; 
layerLabels  = {'Superficial', 'Middle', 'Deep'};

%%
figure;
for iComb = 1:10
    subplot(2,5,iComb);
    imagesc(squeeze(median(crossFreqCrossLayerInfraT(:,iComb,:,:),1,'omitnan'))); axis square; 
    xticks(1:3); yticks(1:3); xticklabels(layerLabels); yticklabels(layerLabels);
    caxis([0 0.75]); colorbar; colormap jet; 
    title(freqCombs{iComb}); 
end 

%% Hybrids 
crossFreqCrossLayerHybridT = reshape(crossFreqCrossLayerHybrid,[size(corrFCHybrid,1)*size(corrFCHybrid,2) ...
    4 5 5 3 3]);
nanRow        = ~isnan(crossFreqCrossLayerHybridT(:,1,1,1,1,1));
crossFreqCrossLayerHybridT(nanRow,:,:,:,:,:) = [];
crossFreqCrossLayerHybridT(~(goodRunsSpatial & ~singleChFlag),:,:,:,:,:) = [];


crossFreqCrossLayerHybridT = reshape(crossFreqCrossLayerHybridT(:,2,:,:,:,:),[size(crossFreqCrossLayerHybridT,1) 25 3 3]);
crossFreqCrossLayerHybridT(:,~goodLabels,:,:) = []; 

figure;
for iComb = 1:10
    subplot(2,5,iComb);
    imagesc(squeeze(median(crossFreqCrossLayerHybridT(:,iComb,:,:),1,'omitnan'))); axis square; 
    xticks(1:3); yticks(1:3); xticklabels(layerLabels); yticklabels(layerLabels);
    caxis([0 0.75]); colorbar; colormap jet; 
    title(freqCombs{iComb}); 
end 

% Take the 10 cross-frequency comparisons
% crossFreqCrossLayerInfraT(:,)

% %
% infraSuper = reshape(superInfraAllBandsT,[size(superInfraAllBandsT,1) 25]);
% infraMid   = reshape(midInfraAllBandsT,[size(superInfraAllBandsT,1) 25]);
% infraDeep  = reshape(deepInfraAllBandsT,[size(superInfraAllBandsT,1) 25]);
% powerSuper = reshape(superPowerAllBandsT,[size(superInfraAllBandsT,1) 25]);
% powerMid   = reshape(midPowerAllBandsT,[size(superInfraAllBandsT,1) 25]);
% powerDeep  = reshape(deepPowerAllBandsT,[size(superInfraAllBandsT,1) 25]);
% 
% figure;
% subplot(131);imagesc(reshape(median(infraSuper,1,'omitnan'),[5 5])); caxis([0 1]); colormap jet;   xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels'); axis square; colorbar;
% subplot(132);imagesc(reshape(median(infraMid,1,'omitnan'),[5 5]));caxis([0 1]); colormap jet;  xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels');axis square; colorbar;
% subplot(133);imagesc(reshape(median(infraDeep,1,'omitnan'),[5 5]));caxis([0 1]); colormap jet;  xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels');axis square;colorbar;
% 
% figure;
% subplot(131);imagesc(reshape(median(infraSuper(smFlagNew == 'S',:),1,'omitnan'),[5 5])); caxis([0 1]); colormap jet;   xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels'); axis square; colorbar;
% subplot(132);imagesc(reshape(median(infraMid(smFlagNew == 'S',:),1,'omitnan'),[5 5]));caxis([0 1]); colormap jet;  xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels');axis square; colorbar;
% subplot(133);imagesc(reshape(median(infraDeep(smFlagNew == 'S',:),1,'omitnan'),[5 5]));caxis([0 1]); colormap jet;  xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels');axis square;colorbar;
% 
% figure;
% subplot(131);imagesc(reshape(median(infraSuper(smFlagNew == 'M',:),1,'omitnan'),[5 5])); caxis([0 1]); colormap jet;   xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels'); axis square; colorbar;
% subplot(132);imagesc(reshape(median(infraMid(smFlagNew == 'M',:),1,'omitnan'),[5 5]));caxis([0 1]); colormap jet;  xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels');axis square; colorbar;
% subplot(133);imagesc(reshape(median(infraDeep(smFlagNew == 'M',:),1,'omitnan'),[5 5]));caxis([0 1]); colormap jet;  xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels');axis square;colorbar;
% 
% 
% figure;
% subplot(131);imagesc(reshape(median(powerSuper,1,'omitnan'),[5 5])); caxis([0 0.5]); colormap jet;   xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels'); axis square; colorbar;
% subplot(132);imagesc(reshape(median(powerMid,1,'omitnan'),[5 5]));caxis([0 0.5]); colormap jet;  xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels');axis square; colorbar;
% subplot(133);imagesc(reshape(median(powerDeep,1,'omitnan'),[5 5]));caxis([0 0.5]); colormap jet;  xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels');axis square;colorbar;
% 
% figure;
% subplot(131);imagesc(reshape(median(powerSuper(smFlagNew == 'S',:),1,'omitnan'),[5 5])); caxis([0 0.5]); colormap jet;   xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels'); axis square; colorbar;
% subplot(132);imagesc(reshape(median(powerMid(smFlagNew == 'S',:),1,'omitnan'),[5 5]));caxis([0 0.5]); colormap jet;  xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels');axis square; colorbar;
% subplot(133);imagesc(reshape(median(powerDeep(smFlagNew == 'S',:),1,'omitnan'),[5 5]));caxis([0 0.5]); colormap jet;  xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels');axis square;colorbar;
% 
% figure;
% subplot(131);imagesc(reshape(median(powerSuper(smFlagNew == 'M',:),1,'omitnan'),[5 5])); caxis([0 0.5]); colormap jet;   xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels'); axis square; colorbar;
% subplot(132);imagesc(reshape(median(powerMid(smFlagNew == 'M',:),1,'omitnan'),[5 5]));caxis([0 0.5]); colormap jet;  xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels');axis square; colorbar;
% subplot(133);imagesc(reshape(median(powerDeep(smFlagNew == 'M',:),1,'omitnan'),[5 5]));caxis([0 0.5]); colormap jet;  xticks(1:5); yticks(1:5); xticklabels(bandLabels');yticklabels(bandLabels');axis square;colorbar;
% 
