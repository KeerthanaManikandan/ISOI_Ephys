% ephysImaging_compilation_v2
% This script compiles the data procesed from ephysImaging_v3 and compiles
% the following:
% 1. Cross correlations between physiology and imaging for ROI and FOV
% 2. Temporal and spatial controls.
% 3. Distributions of cross correlations and lags as a function of distance
% See ephysImaging_v3 for reference
% see ephysImaging_compilation.m for previous version of this script
% July 3, 2024 - KM
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
rmpath(genpath([commonDir '\Codes\chronux_2_12\spectral_analysis\continuous\dupes']));

%% Retrieve the hybrid map for both animals
hemisphere = 'Left'; spatialBin = 3;
iM = 2;
switch iM
    case 1
        monkeyName = 'CharlieSheen';
        goodRuns  = ([1 1 NaN NaN NaN;... % Date x run
            1 1 NaN NaN NaN; ...
            1 1 1 0 1; ...
            1 1 1 1 NaN]);

        goodRunsSpatial  = ([1 1 NaN NaN NaN;... % Date x run
            1 1 NaN NaN NaN; ...
            1 1 1 0 1; ...
            1 1 1 1 NaN]);

    case 2
        monkeyName = 'Whiskey';
        goodRuns   = ([1 1 1 NaN NaN NaN NaN; ... % Date x run
            1 1 1 1 1 1 1 ; ...
            1 0 1 1 NaN NaN NaN; ...
            1 1 0 1 1 1 0; ...
            1 0 0 1 1 1 1]);

        goodRunsSpatial = ([1 1 1 NaN NaN NaN NaN; ... % Date x run
            1 1 1 1 1 1 1 ; ...
            1 0 1 1 NaN NaN NaN; ...
            1 1 0 1 1 1 1; ...
            1 0 0 1 1 1 0]);
end
goodRuns = reshape(goodRuns,[size(goodRuns,1)*size(goodRuns,2) 1]);
goodRuns(isnan(goodRuns)) = []; goodRuns = logical(goodRuns);

goodRunsSpatial = reshape(goodRunsSpatial,[size(goodRunsSpatial,1)*size(goodRunsSpatial,2) 1]);
goodRunsSpatial(isnan(goodRunsSpatial)) = []; goodRunsSpatial = logical(goodRunsSpatial);

%% Get monkey parameters
[allDates,allRuns, refDate, refDir,lensCombo, roiSize, ephysFileNameAll, serverPath,probeLabel,...
    chInCortexNotes, greenMapRef] = getMonkeyParams_Imaging_Ephys(monkeyName,commonDir, hemisphere);

% Get monkey data....
[processedDat,greenIm,probe,badCh,badTimesLFP,badTimeThresh,estChInCortex] = ...
    getAllData_Imaging_Ephys(monkeyName,hemisphere,commonDir,serverPath,allDates,allRuns,...
    ephysFileNameAll,greenMapRef,chInCortexNotes,probeLabel,spatialBin);

clc; disp(['All physiology and imaging data for ' monkeyName ' loaded']);

%% Get cross correlation for ROI 
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
                circleRad =round(roiSize{1}(1)/(spatialBin));
                pDatTemp = processedDat{iDate,iRun}.tempBandPass;
                seedSigT = calculateSeedSignal(greenFig,corrMask,...
                    seedLocIn,circleRad,pDatTemp); % Get Gaussian weighted seed signal
                corrMapT           =plotCorrMap(seedSigT,pDatTemp,0);
                greenImRGB = ind2rgb(greenFig,gray(256));
                figure; imagesc(greenImRGB);axis image; colormap jet; axis image off;
                hold on; imagesc(corrMapT,'AlphaData',corrMapT.*0.8);caxis([0 1]);colorbar;
                f = gcf; exportgraphics(f,[dataDir '\FCMap_ROI.png'],'Resolution',300); close gcf;
            end

            % Get the ROI for gamma cross correlations
            clipMaskROI = clipMaskCortex(seedLocIn(2)-round(seedRad/2):seedLocIn(2)+round(seedRad/2),...
                seedLocIn(1)-round(seedRad/2):seedLocIn(1)+round(seedRad/2));

            clear tempProfileNoRefRun tempProfileSuperNoRefRun tempProfileDeepNoRefRun


            [tempProfileNoRefRun,tempProfileSuperNoRefRun,...
                tempProfileDeepNoRefRun,~] = ...
                getCrossCorrROI(dataDir,monkeyName,expDate,runName,...
                processedDat{iDate,iRun}.tempBandPass,probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,...
                badTimesLFP{iDate,iRun},badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),seedLocIn,...
                seedRad,clipMaskROI,10,'NoRef');
            save([dataDir '\ROIAllVars.mat'],'tempProfileNoRefRun','tempProfileSuperNoRefRun','tempProfileDeepNoRefRun');

        else
            clear allVars
            allVars = load([dataDir '\ROIAllVars.mat']);
            tempProfileNoRefRun      = allVars.tempProfileNoRefRun;
            tempProfileSuperNoRefRun = allVars.tempProfileSuperNoRefRun;
            tempProfileDeepNoRefRun  = allVars.tempProfileDeepNoRefRun;
        end

        tempProfileNoRef(iM,iDate,iRun)      = tempProfileNoRefRun;
        tempProfileSuperNoRef(iM,iDate,iRun) = tempProfileSuperNoRefRun;
        tempProfileDeepNoRef(iM,iDate,iRun)  = tempProfileDeepNoRefRun;

        %
        %         if ~exist([dataDir '\AvgRef'],'dir')
        %             [tempProfileAvg(iM,iDate,iRun),tempProfileSuperAvg(iM,iDate,iRun),...
        %                 tempProfileDeepAvg(iM,iDate,iRun),...
        %                 lags] = getCrossCorrROI(dataDir,monkeyName,expDate,runName,...
        %                 processedDat{iDate,iRun}.tempBandPass,probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,...
        %                 badTimesLFP{iDate,iRun},badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),seedLocIn,...
        %                 seedRad,clipMaskROI,10,'AvgRef');
        %         end
        %
        %         if ~exist([dataDir '\BipolarRef'],'dir')
        %             [tempProfileBipolar(iM,iDate,iRun),tempProfileSuperBipolar(iM,iDate,iRun),...
        %                 tempProfileDeepBipolar(iM,iDate,iRun),~] ...
        %                 = getCrossCorrROI(dataDir,monkeyName,expDate,runName,...
        %                 processedDat{iDate,iRun}.tempBandPass,probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,...
        %                 badTimesLFP{iDate,iRun},badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),seedLocIn,...
        %                 seedRad,clipMaskROI,10,'BipolarRef');
        %         end
        %
        %         if ~exist([dataDir '\AvgRefTop5_Bottom5'],'dir')
        %             [tempProfileAvg_5(iM,iDate,iRun),tempProfileSuperAvg_5(iM,iDate,iRun),...
        %                 tempProfileDeepAvg_5(iM,iDate,iRun),...
        %                 lags] = getCrossCorrROI(dataDir,monkeyName,expDate,runName,...
        %                 processedDat{iDate,iRun}.tempBandPass,probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,...
        %                 badTimesLFP{iDate,iRun},badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),seedLocIn,...
        %                 seedRad,clipMaskROI,5,'AvgRefTop5_Bottom5');
        %         end
        %
        %         if ~exist([dataDir '\BipolarRef_Top5_Bottom5'],'dir')
        %             [tempProfileBipolar_5(iM,iDate,iRun),tempProfileSuperBipolar_5(iM,iDate,iRun),...
        %                 tempProfileDeepBipolar_5(iM,iDate,iRun),~] = ...
        %                 getCrossCorrROI(dataDir,monkeyName,expDate,runName,...
        %                 processedDat{iDate,iRun}.tempBandPass,probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,...
        %                 badTimesLFP{iDate,iRun},badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),seedLocIn,...
        %                 seedRad,clipMaskROI,5,'BipolarRef_Top5_Bottom5');
        %         end
        %
        %         if ~exist([dataDir '\NoRefTop5_Bottom5'],'dir')
        %             [tempProfileNoRef_5(iM,iDate,iRun),tempProfileSuperNoRef_5(iM,iDate,iRun),...
        %                 tempProfileDeepNoRef_5(iM,iDate,iRun),~]...
        %                 = getCrossCorrROI(dataDir,monkeyName,expDate,runName,...
        %                 processedDat{iDate,iRun}.tempBandPass,probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,...
        %                 badTimesLFP{iDate,iRun},badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),seedLocIn,...
        %                 seedRad,clipMaskROI,5,'NoRefTop5_Bottom5');
        %         end
    end
end


%%
% Plot all the temporal profiles... (median+/- sem) Refer ephysImaging_v2
clear medTempProfile bandWiseLagsAll bandWiseMagsAll bandWiseLagsSuper...
    bandWiseMagsSuper bandWiseLagsDeep bandWiseMagsDeep
x = -200:200;%lags;
negIdx = x<0 & x>=-150; xNew = x(negIdx);
 bandLabels = {'Gamma band'; 'Alpha band' ;'Beta band'; 'Theta band';'Spike band'};

for iRef = 3%1:3
    clear allCh superCh deepCh avgLabel
    switch iRef
        case 1
            allCh    = tempProfileAvg;
            superCh  = tempProfileSuperAvg;
            deepCh   = tempProfileDeepAvg;
            refLabel = 'Avg';
        case 2
            allCh    = tempProfileBipolar;
            superCh  = tempProfileSuperBipolar;
            deepCh   = tempProfileDeepBipolar;
            refLabel = 'Bipolar';
        case 3
            allCh    = tempProfileNoRef;
            superCh  = tempProfileSuperNoRef;
            deepCh   = tempProfileDeepNoRef;
            refLabel = 'NoRef';
        case 4
            allCh    = tempProfileAvg_5;
            superCh  = tempProfileSuperAvg_5;
            deepCh   = tempProfileDeepAvg_5;
            refLabel = 'Avg_5';
        case 5
            allCh    = tempProfileBipolar_5;
            superCh  = tempProfileSuperBipolar_5;
            deepCh   = tempProfileDeepBipolar_5;
            refLabel = 'Bipolar_5';
        case 6
            allCh    = tempProfileNoRef_5;
            superCh  = tempProfileSuperNoRef_5;
            deepCh   = tempProfileDeepNoRef_5;
            refLabel = 'NoRef_5';
    end

    if iRef<4
        avgLabel = '10_Ch';
    else
        avgLabel ='Top5_Bottom5';
    end

    for iChType = 1%:3
        clear profileAll chLabel allBandLags allBandMags
        switch iChType
            case 1
                profileAll = allCh;
                chLabel    = 'All_Ch';
            case 2
                profileAll = superCh;
                chLabel    = 'Super_Ch';
            case 3
                profileAll = deepCh;
                chLabel    = 'Deep_Ch';
        end
        figure;
        for iBand = 1:5
            switch iBand
                case 1
                    profileT  = [profileAll(iM,:,:).profile];
                    lagVals   = [profileAll(iM,:,:).lagLow];
                    magVals   = [profileAll(iM,:,:).magLow];
                    bandLabel = 'Gamma band';
                case 2
                    profileT  = [profileAll(iM,:,:).profileAlpha];
                    lagVals   = [profileAll(iM,:,:).lagLowAlpha];
                    magVals   = [profileAll(iM,:,:).magLowAlpha];
                    bandLabel = 'Alpha band';
                case 3
                    profileT  = [profileAll(iM,:,:).profileBeta];
                    lagVals   = [profileAll(iM,:,:).lagLowBeta];
                    magVals   = [profileAll(iM,:,:).magLowBeta];
                    bandLabel = 'Beta band';
                case 4
                    profileT  = [profileAll(iM,:,:).profileTheta];
                    lagVals   = [profileAll(iM,:,:).lagLowTheta];
                    magVals   = [profileAll(iM,:,:).magLowTheta];
                    bandLabel = 'Theta band';
                case 5
                    profileT  = [profileAll(iM,:,:).profileRaw];
                    lagVals   = [profileAll(iM,:,:).lagLowRaw];
                    magVals   = [profileAll(iM,:,:).magLowRaw];
                    bandLabel = 'Spike band';
            end

            % Remove rows of zeros
            tempProfileR = profileT';
            zeroR = find(all(tempProfileR == 0,2));
            if ~isempty(zeroR);tempProfileR(zeroR,:) = [];end

            medTempProfileA = squeeze(median(tempProfileR(goodRuns,:),1,'omitnan'));
            semA  = std(tempProfileR(goodRuns,:),0,1)./sqrt(size(tempProfileR(goodRuns,:),1));
% 
%             medTempProfileR = squeeze(median(tempProfileR(~goodRuns,:),1,'omitnan'));
%             semR  = std(tempProfileR(~goodRuns,:),0,1)./sqrt(size(tempProfileR(~goodRuns,:),1));
% 
            subplot(2,3,iBand); %plot(x,tempProfileR(goodRuns,:),'Color',[0.65 0.65 0.65]);
            plot(x,smooth(medTempProfileA,5),'k','LineWidth',1);xlim([-200 200]); hold on;
            patch([x fliplr(x)],[(medTempProfileA-2.*semA) fliplr((medTempProfileA+2.*semA))],...
                'blue','FaceAlpha',0.3,'EdgeColor','none'); box off;
            xline(0);xticks(-200:50:200);xticklabels(-20:5:20);
            ylim([-0.35 0.2]);  yticks(-0.3:0.1:0.2); title(bandLabels{iBand});
% 
%             subplot(122); %plot(x,tempProfileR(~goodRuns,:),'Color',[0.65 0.65 0.65]); hold on;
%             plot(x,movmean(medTempProfileR,1),'k','LineWidth',1);xlim([-200 200]);hold on;
%             patch([x' fliplr(x')],[(medTempProfileR-2.*semR) fliplr((medTempProfileR+2.*semR))],...
%                 'blue','FaceAlpha',0.3,'EdgeColor','none'); grid on;
%             xline(0);xticklabels(-20:10:20);ylim([-0.5 0.5]); title('Rejected runs');
%             sgtitle(strrep([bandLabel ' ' chLabel ' Ref Type: ' refLabel],'_','\_'));
% 
%             f = gcf; exportgraphics(f,['X:\Data\' monkeyName '_SqM\Left Hemisphere\Results\ISOI_Ephys\tempProfile_'...
%                 bandLabel(1:end-4) '_' chLabel '_' refLabel '_' avgLabel '.png'],'Resolution',300); close gcf;

            medTempProfile(iRef,iChType,iBand,:) = medTempProfileA;

            % Lag magnitudes
            lagR = lagVals';
            lagR(lagR==0) = []; allBandLags(:,iBand) = lagR./10;

            magR = magVals';
            magR(magR==0) = []; allBandMags(:,iBand) = magR;
        end
%{
        boxLabels = {'Gamma band'; 'Alpha band' ;'Beta band'; 'Theta band';'Spike band'};
        figure; subplot(141); boxplot(allBandLags(goodRuns,:),'Labels',boxLabels); ylim([-15 2]);
        ylabel('Peak negative Lags (s)'); title('Accepted runs');
        subplot(142); boxplot(allBandMag(goodRuns,:),'Labels',boxLabels); ylim([-0.5 0.5]);
        ylabel('Peak negative correlation magnitude'); title('Accepted runs');

        subplot(143); boxplot(allBandLags(~goodRuns,:),'Labels',boxLabels); ylim([-15 2]);
        ylabel('Peak negative Lags (s)'); title('Rejected runs');
        subplot(144); boxplot(allBandMag(~goodRuns,:),'Labels',boxLabels); ylim([-0.5 0.5]);
        ylabel('Peak negative correlation magnitude'); title('Rejected runs');

        sgtitle(strrep([chLabel ' Ref Type: ' refLabel],'_','\_'));
        f = gcf; exportgraphics(f,['X:\Data\' monkeyName '_SqM\Left Hemisphere\Results\ISOI_Ephys\tempMag_' ...
            chLabel '_' avgLabel '_' refLabel '.png'],'Resolution',300); close gcf;

%}
        bandWiseLagsAll(iRef,iChType,:,:) = allBandLags(goodRuns,:);
        bandWiseMagsAll(iRef,iChType,:,:) = allBandMags(goodRuns,:);

    end
end

% Show temporal profiles for all frequency bands
profileAlpha    = [tempProfileNoRef(iM,:,:).profileAlpha]; profileAlpha(:,~goodRuns) = []; 
profileBeta     =  [tempProfileNoRef(iM,:,:).profileBeta]; profileBeta(:,~goodRuns) = []; 
profileGamma    =  [tempProfileNoRef(iM,:,:).profile]; profileGamma(:,~goodRuns) = []; 
profileTheta    =  [tempProfileNoRef(iM,:,:).profileTheta]; profileTheta(:,~goodRuns) = []; 
profileSpiking  =  [tempProfileNoRef(iM,:,:).profileRaw]; profileSpiking(:,~goodRuns) = []; 

medTempProfileAll(:,1) = median(profileGamma,2,'omitnan');
medTempProfileAll(:,2) = median(profileAlpha,2,'omitnan');
medTempProfileAll(:,3) = median(profileBeta,2,'omitnan');
medTempProfileAll(:,4) = median(profileTheta,2,'omitnan');
medTempProfileAll(:,5) = median(profileSpiking,2,'omitnan');

figure;
for iBand = 1:5
    plot(x,smooth(medTempProfileAll(:,iBand),10),'LineWidth',1.5);hold on;
end

xlim([-200 200]); xticks(-200:50:200); xticklabels(-20:5:20); ylim([-0.5 0.5]); yticks(-0.5:0.1:0.5);
xlabel('Lag (s)'); ylabel('Cross correlations between FC map and hybrids'); xline(0); box off;

figure; 
for iPlot = 1:5
    subplot(2,3,iPlot);
    plot(x,smooth(medTempProfileAll(:,iPlot),3),'LineWidth',1.5);
    title(bandLabels{iPlot});
    xticks(-200:50:200);xticklabels(-20:5:20);ylim([-0.5 0.5]);
end 
xlabel('Lags (s)'); ylabel('Cross correlation between ISOI and LFP Powers');
figure; plot(x,movmedian(medTempProfileAll,3),'LineWidth',1.5); 
xticks(-200:50:200);xticklabels(-20:5:20);ylim([-0.5 0.5]);
xlabel('Lags (s)'); ylabel('Cross correlation between ISOI and LFP Powers');
legend(bandLabels,'AutoUpdate','off','Location','northeast','box', 'off');xline(0); box off; 

% Show box plot of correlations for good runs
allBandMags(~goodRuns,:) = []; 

figure;
allBandSort1 = allBandMags(:,[4 2 3 1 5]);
 subplot(121); boxplot(allBandSort1,{'Theta';'Alpha';'Beta';'Gamma';'Spiking'}); 
ylim([-0.5 0.5]); ylabel('Cross Correlation between ISOI and LFP Powers');box off; 

allBandSort2 = allBandMags(:,[1 3 2 5 4]);
subplot(122);  boxplot(allBandSort2,{'Gamma,';'Beta';'Alpha'; 'Spiking';'Theta'});
ylim([-0.5 0.5]); ylabel('Cross Correlation between ISOI and LFP Powers'); box off;


[pCorr,tblCorr,statsCorr] = anova1(allBandMags,bandLabels);
[rCorr,mCorr,~,gnamesCorr] = multcompare(statsCorr,"CriticalValueType","bonferroni");
tblCorrM = array2table(rCorr,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblCorrM.("Group") = gnamesCorr(tblCorrM.("Group"));
tblCorrM.("Control Group") = gnamesCorr(tblCorrM.("Control Group"));


%% Show the lags and correlation magnitudes for all bands,references and channel type
bandLabels = {'Gamma band'; 'Alpha band'; 'Beta band'; 'Theta band'; 'Spiking'};
layerLabels = {'All_Ch'; 'Superficial_Ch'; 'Deep_Ch'};
refLabels = {'Avg';'Bipolar';'No Ref'};
for iRef = 3%1:3
    for iChType = 1:3
        bandLags = squeeze(bandWiseLagsAll(iRef,iChType,:,:));
        bandMags = squeeze(bandWiseMagsAll(iRef,iChType,:,:));

        figure; subplot(121); boxplot(bandLags,'Labels',bandLabels);ylim([-15 2]);
        ylabel('Peak negative Lags (s)');

        subplot(122); boxplot(bandMags,'Labels',bandLabels); ylim([-0.5 0.5]);
        ylabel('Correlation at peak negative');
        sgtitle(strrep([layerLabels{iChType} ' ' refLabels{iRef}],'_','\_'));
    end
end

%%  Show the temporal profile for the three types of references
for iBand = 1:5
    allRefProfile = squeeze(medTempProfile(:,:,iBand,:));

    figure;
    for iLayer = 1:3
        subplot(1,3,iLayer); plot(x,squeeze(allRefProfile(1:3,iLayer,:)),'LineWidth',2);
        title(strrep(layerLabels{iLayer},'_','\_')); xlim([-200 200]); grid on;
        xline(0);xticks(-200:100:200);xticklabels(-20:10:20);ylim([-0.3 0.3]);
    end
    legend(refLabels,'Location','northeast');
    sgtitle(strrep(['All Channels ' bandLabels{iBand}],'_','\_'));

    figure;
    for iLayer = 1:3
        subplot(1,3,iLayer); plot(x,squeeze(allRefProfile(4:6,iLayer,:)),'LineWidth',2);
        title(strrep(layerLabels{iLayer},'_','\_'));xlim([-200 200]);
        xline(0);xticks(-200:100:200);xticklabels(-20:10:20);ylim([-0.3 0.3]); grid on;
    end
    legend(refLabels,'Location','northeast');
    sgtitle(strrep([ 'Top5_Bottom5 ' bandLabels{iBand}],'_','\_'));
end

%% Show the distribution of lag for different frequencies and references
for iBand = 1:5
    figure;
    for iRef = 1:3
        lagDistAll   = squeeze(bandWiseLagsAll(iRef,:,goodRuns,iBand));
        subplot(1,3,iRef); boxplot(lagDistAll','Labels',layerLabels); ylim([-10 0]);
        ylabel('Peak negative Lags (s)'); title(refLabels{iRef});
    end
    sgtitle(bandLabels{iBand});

    figure;
    for iRef = 1:3
        magDistAll   = squeeze(bandWiseMagsAll(iRef,:,goodRuns,iBand));
        subplot(1,3,iRef); boxplot(magDistAll','Labels',layerLabels); ylim([-1 0]);
        ylabel('Magnitude at peak negative'); title(refLabels{iRef});
    end
    sgtitle(bandLabels{iBand});
end

%{
%% Show temporal profile for different references for one example run
iDate = 2; iRun = 7;
for iChType = 1:3
    switch iChType
        case 1
            avgProfile     = tempProfileAvg;
            bipProfile     = tempProfileBipolar;
            noRefProfile   = tempProfileNoRef;
            avgProfile_5   = tempProfileAvg_5;
            bipProfile_5   = tempProfileBipolar_5;
            noRefProfile_5 = tempProfileNoRef_5;
        case 2
            avgProfile     = tempProfileSuperAvg;
            bipProfile     = tempProfileSuperBipolar;
            noRefProfile   = tempProfileSuperNoRef;
            avgProfile_5   = tempProfileSuperAvg_5;
            bipProfile_5   = tempProfileSuperBipolar_5;
            noRefProfile_5 = tempProfileSuperNoRef_5;
        case 3
            avgProfile     = tempProfileDeepAvg;
            bipProfile     = tempProfileBipolar;
            noRefProfile   = tempProfileNoRef;
            avgProfile_5   = tempProfileAvg_5;
            bipProfile_5   = tempProfileBipolar_5;
            noRefProfile_5 = tempProfileNoRef_5;
    end

    figure;
    for iType = 1:2
        switch iType
            case 1
                profileDat(:,1) = avgProfile(iM,iDate,iRun).profile;
                profileDat(:,2) = bipProfile(iM,iDate,iRun).profile;
                profileDat(:,3) = noRefProfile(iM,iDate,iRun).profile;
                typeLabel       = '10 channels (super/deep)';
            case 2
                profileDat(:,1) = avgProfile_5(iM,iDate,iRun).profile;
                profileDat(:,2) = bipProfile_5(iM,iDate,iRun).profile;
                profileDat(:,3) = noRefProfile_5(iM,iDate,iRun).profile;
                typeLabel       = 'Top 5 - Bottom 5 channels';
        end

        subplot(1,2,iType);plot(x,profileDat,'LineWidth',2); title(typeLabel);
        ylim([-0.5 0.5]); grid on; xlim([-200 200]);
        xline(0);xticks(-200:100:200);xticklabels(-20:10:20);
        legend(refLabels);
    end
    sgtitle(strrep(layerLabels{iChType},'_','\_'));
end
%}
%% Full FOV
videoFlag = 1; % To save videos
for iDate = 3:5%size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);
    for iRun = 1:size(allRuns{iDate,1})
        clear runName dataDir clipMask elecMask clipMaskCortex corrMask x negIdx lowIdx
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

        clc; disp(['Analyzing data for ' monkeyName ' '  expDate ' run: ' runName]);

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

        corrMaskT = reshape(corrMask,[361*438 1]);

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
        
%         % Get the cross correlations for the entire FOV
%         if ~exist(fullfile(dataDir,'crossCorrFOV.mat'),'file') 
%             fovFlag = 1;
%         else
%             crossCorrFOV{iRun,iDate} = matfile([dataDir '\crossCorrFOV.mat']);
%             try allXcorr = crossCorrFOV{iRun,iDate}.spatialProfile;
%                 fovFlag = 0;
%             catch
%                 fovFlag = 1;
%             end
%         end

        if ~exist([dataDir '\crossCorrFOV_10_NoRef.mat'],'file')
            disp('No reference - channel split 10(superficial)/10(Deep)');
            getCrossCorrFOV(monkeyName,expDate,runName,dataDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,10,'NoRef');

            disp('No reference - channel split top 6(superficial)/bottom 6(Deep)');
            getCrossCorrFOV(monkeyName,expDate,runName,dataDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,6,'NoRef');

            disp('Adjacent channel avg - channel split 10(superficial)/10(Deep)');
            getCrossCorrFOV(monkeyName,expDate,runName,dataDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,10,'AdjacentAvg');
          
            disp('Adjacent channel avg - channel split top 6(superficial)/bottom 6(Deep)');
              getCrossCorrFOV(monkeyName,expDate,runName,dataDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,6,'AdjacentAvg');


        end
        for iType = 1:4
            clear crossCorrFOV allXCorr superXCorr deepXCorr allLags fileName
            switch iType
                case 1
                    crossCorrFOV = matfile([dataDir '\crossCorrFOV_10_NoRef.mat']);
                    allXcorr     = crossCorrFOV.spatialProfile;
                    superXcorr   = crossCorrFOV.spatialProfileSuper;
                    deepXcorr    = crossCorrFOV.spatialProfileDeep;
                    allLags      = crossCorrFOV.lagFull;
                    fileName     = '10_NoRef';
                case 2
                    crossCorrFOV = matfile([dataDir '\crossCorrFOV_6_NoRef.mat']);
                    allXcorr     = crossCorrFOV.spatialProfile;
                    superXcorr   = crossCorrFOV.spatialProfileSuper;
                    deepXcorr    = crossCorrFOV.spatialProfileDeep;
                    allLags      = crossCorrFOV.lagFull;
                    fileName     = '6_NoRef';
                case 3
                    crossCorrFOV = matfile([dataDir '\crossCorrFOV_10_AdjacentAvg.mat']);
                    allXcorr     = crossCorrFOV.spatialProfile;
                    superXcorr   = crossCorrFOV.spatialProfileSuper;
                    deepXcorr    = crossCorrFOV.spatialProfileDeep;
                    allLags      = crossCorrFOV.lagFull;
                    fileName     = '10_AdjacentAvg';
                case 4
                    crossCorrFOV = matfile([dataDir '\crossCorrFOV_6_AdjacentAvg.mat']);
                    allXcorr     = crossCorrFOV.spatialProfile;
                    superXcorr   = crossCorrFOV.spatialProfileSuper;
                    deepXcorr    = crossCorrFOV.spatialProfileDeep;
                    allLags      = crossCorrFOV.lagFull;
                    fileName     = '6_AdjacentAvg';
            end

            x = allLags;
            negIdx = x<0 & x>=-150; negVals = x(negIdx);
            lowIdx = x<0 & x>= -80; xLow = x(lowIdx);

            %         crossCorrFOV{iRun,iDate} = matfile([dataDir '\crossCorrFOV.mat']);
            %         allXcorr                 = crossCorrFOV{iRun,iDate}.spatialProfile;
            %         superXcorr               = crossCorrFOV{iRun,iDate}.spatialProfileSuper;
            %         deepXcorr                = crossCorrFOV{iRun,iDate}.spatialProfileDeep;
            %         allLags{iRun,iDate}      = crossCorrFOV{iRun,iDate}.lagFull;
            %
            %
            %         x = allLags{iRun,iDate};
            %         negIdx = x<0 & x>=-150; xNew = x(negIdx);
            %         lowIdx = x<0 & x>= -80; xLow = x(lowIdx);
            %
            %         % Show and save the full field of view for the hybrid map
            for iBand = 4%1:5
                clear bandName crossCorr lagLow lagLowSuper lagLowDeep
                switch iBand
                    case 1
                        bandName       = 'Theta';
                        crossCorr      = allXcorr.ccFullTheta;
                        crossCorrSuper = superXcorr.ccFullTheta;
                        crossCorrDeep  = deepXcorr.ccFullTheta;
                        lagLow         = tempProfileNoRef(iM,iDate,iRun).lagLowTheta;
                        lagLowSuper    = tempProfileSuperNoRef(iM,iDate,iRun).lagLowTheta;
                        lagLowDeep     = tempProfileDeepNoRef(iM,iDate,iRun).lagLowTheta;

                    case 2
                        bandName       = 'Alpha';
                        crossCorr      = allXcorr.ccFullAlpha;
                        crossCorrSuper = superXcorr.ccFullAlpha;
                        crossCorrDeep  = deepXcorr.ccFullAlpha;
                        lagLow         = tempProfileNoRef(iM,iDate,iRun).lagLowAlpha;
                        lagLowSuper    = tempProfileSuperNoRef(iM,iDate,iRun).lagLowAlpha;
                        lagLowDeep     = tempProfileDeepNoRef(iM,iDate,iRun).lagLowAlpha;
                   
                    case 3
                        bandName       = 'Beta';
                        crossCorr      = allXcorr.ccFullBeta;
                        crossCorrSuper = superXcorr.ccFullBeta;
                        crossCorrDeep  = deepXcorr.ccFullBeta;
                        lagLow         = tempProfileNoRef(iM,iDate,iRun).lagLowBeta;
                        lagLowSuper    = tempProfileSuperNoRef(iM,iDate,iRun).lagLowBeta;
                        lagLowDeep     = tempProfileDeepNoRef(iM,iDate,iRun).lagLowBeta;
                    
                    case 4
                        bandName       = 'Gamma';
                        crossCorr      = allXcorr.ccFull;
                        crossCorrSuper = superXcorr.ccFull;
                        crossCorrDeep  = deepXcorr.ccFull;
                        lagLow         = tempProfileNoRef(iM,iDate,iRun).lagLow;
                        lagLowSuper    = tempProfileSuperNoRef(iM,iDate,iRun).lagLow;
                        lagLowDeep     = tempProfileDeepNoRef(iM,iDate,iRun).lagLow;
                    
                    case 5
                        bandName       = 'Spike';
                        crossCorr      = allXcorr.ccFullRaw;
                        crossCorrSuper = superXcorr.ccFullRaw;
                        crossCorrDeep  = deepXcorr.ccFullRaw;
                        lagLow         = tempProfileNoRef(iM,iDate,iRun).lagLowRaw;
                        lagLowSuper    = tempProfileSuperNoRef(iM,iDate,iRun).lagLowRaw;
                        lagLowDeep     = tempProfileDeepNoRef(iM,iDate,iRun).lagLowRaw;
                end

                if ~exist(fullfile(dataDir,['HybridMapFOV_' fileName '_' bandName '.png']),'file')
                    clear frameLow grayImFull
                    frameLow       = (x == lagLow);
                    figure; imagesc(squeeze(crossCorr(frameLow,:,:))); hold on; axis image off;
                    colormap(flipud(jet)); caxis([-0.5 0.5]);
                    grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));
                    h = imagesc(grayImFull); hold off; set(h,'AlphaData',~corrMask);
                    title(['Peak negative at ' num2str(lagLow/10) 's']);
                    f = gcf; exportgraphics(f,[dataDir '\HybridMapFOV_' fileName '_' bandName '.png'],'Resolution',300); close gcf;
                end

                if ~exist(fullfile(dataDir,['HybridMapFOV_Super_' fileName '_' bandName '.png']),'file')
                    clear frameLow grayImFull
                    frameLow       = (x == lagLowSuper);
                    figure; imagesc(squeeze(crossCorrSuper(frameLow,:,:))); hold on; axis image off;
                    colormap(flipud(jet)); caxis([-0.5 0.5]);
                    grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));
                    h = imagesc(grayImFull); hold off; set(h,'AlphaData',~corrMask);
                    title(['Superficial channels-Peak negative at ' num2str(lagLowSuper/10) 's']);
                    f = gcf; exportgraphics(f,[dataDir '\HybridMapFOV_Super_' fileName '_' bandName '.png'],'Resolution',300); close gcf;
                end

                if ~exist(fullfile(dataDir,['HybridMapFOV_Deep_' fileName '_' bandName '.png']),'file')
                    clear frameLow grayImFull
                    frameLow       = (x == lagLowDeep);
                    figure; imagesc(squeeze(crossCorrDeep(frameLow,:,:))); hold on; axis image off;
                    colormap(flipud(jet)); caxis([-0.5 0.5]);
                    grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));
                    h = imagesc(grayImFull); hold off; set(h,'AlphaData',~corrMask);
                    title(['Deep channels-Peak negative at ' num2str(lagLowDeep/10) 's']);
                    f = gcf; exportgraphics(f,[dataDir '\HybridMapFOV_Deep_' fileName '_' bandName '.png'],'Resolution',300); close gcf;
                end

                % Correlate FC map with hybrids
                clear mapsAll
                if iBand == 4
                    mapsAll   = reshape(crossCorr,[401 imSize(1)*imSize(2)]);
                    mapsAll(:,~corrMaskT) = NaN; 
                    for iMap = 1:size(mapsAll,1)
                        corrFCHybrid(iDate,iRun,iType,iMap,:) = corr(fcMap,mapsAll(iMap,:)','rows','complete');
                    end

                    [peakNegValsAll(iDate,iRun,iType),peakNegTimesAll(iDate,iRun,iType)] = min(squeeze(corrFCHybrid(iDate,iRun,iType,negIdx)));
                    peakNegTimesAll(iDate,iRun,iType) = negVals(peakNegTimesAll(iDate,iRun,iType))./10;

                    superHybridMap = reshape(squeeze(crossCorrSuper((x == lagLowSuper),:,:)),[361*438 1]);
                    deepHybridMap  = reshape(squeeze(crossCorrDeep((x == lagLowDeep),:,:)),[361*438 1]);
                    superHybridMap(~corrMaskT) = NaN;
                    deepHybridMap(~corrMaskT) = NaN;
                    super_DeepCorr(iDate,iRun,iType) = corr(superHybridMap,deepHybridMap,'rows','complete');
                end

                %
                %
                %             % Save videos of hybrid maps
                % %             if videoFlag
                % %                 saveVideoFOV(dataDir,bandName,monkeyName,expDate,runName,crossCorr,clipMaskCortex,x);
                % %             end
                %
                %
            end
        end
    end
end
% 1. Is there an effect of reference scheme on correlation or lag
gR= ([1 1 1 NaN NaN NaN NaN; ... % Date x run
            1 1 1 1 1 1 1 ; ...
            1 0 1 1 NaN NaN NaN; ...
            1 1 0 1 1 1 1; ...
            1 0 0 1 1 1 0]);
gR(1,:) = []; gR = reshape(gR,[28 1]);
nanValsGR = isnan(gR); gR(isnan(gR)) = [];
peakNegValsAll(1,:,:) = [];
peakNegTimesAll(1,:,:) = []; 
peakNegValsT = reshape(peakNegValsAll,[28 4]); peakNegTimesT = reshape(peakNegTimesAll,[28 4]); 
peakNegValsT(nanValsGR,:) = []; peakNegTimesT(nanValsGR,:) = []; 
peakNegValsT(~gR,:) = []; peakNegTimesT(~gR,:) = [];
figure;subplot(121); boxplot(peakNegValsT(:,[1 3]),{'No Ref'; 'Adj channel ref'});
ylabel('Max correlation with FC map'); ylim([-1 -0.2]);
subplot(122); boxplot(peakNegTimesT(:,[1 3]),{'No Ref'; 'Adj channel ref'});
ylabel('Lag at peak negative'); ylim([-16 1]);


% 1. Is there an effect of channel picking on correlation between
% superficial/deep hybrid maps
super_DeepCorr(1,:,:) = [];
super_DeepCorrT = reshape(super_DeepCorr,[28 4]); super_DeepCorrT(nanValsGR,:) = []; 
super_DeepCorrT(~gR,:) = [];
figure; boxplot(super_DeepCorrT(:,[1 2]),{'10/10'; 'Top 3rd/bottom 3rd'}); ylim([-0.2 1.1]);
ylabel('Correlation between peak negative superficial and deep hybrid maps');


%% Correlate superficial vs deep hybrids 
x = -200:200;
tic;
for iDate = 1:size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);
    
    for iRun = 1:size(allRuns{iDate,1})
        clear runName dataDir clipMask elecMask clipMaskCortex corrMask negIdx lowIdx
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
        crossCorrFOV{iRun,iDate} = matfile([dataDir '\crossCorrFOV.mat']);
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
        corrMaskT = reshape(corrMask,[361*438 1]);
       
        superXcorr     = crossCorrFOV{iRun,iDate}.spatialProfileSuper;
        deepXcorr      = crossCorrFOV{iRun,iDate}.spatialProfileDeep;

        crossCorrSuper = superXcorr.ccFull;
        crossCorrDeep  = deepXcorr.ccFull;

        lagLowSuper    = tempProfileSuperNoRef(iM,iDate,iRun).lagLow;
        lagLowDeep     = tempProfileDeepNoRef(iM,iDate,iRun).lagLow;

        superHybridMap = reshape(squeeze(crossCorrSuper((x == lagLowSuper),:,:)),[361*438 1]);
        deepHybridMap  = reshape(squeeze(crossCorrDeep((x == lagLowDeep),:,:)),[361*438 1]);
        superHybridMap(~corrMaskT) = NaN;
        deepHybridMap(~corrMaskT) = NaN;

        super_DeepCorrNew(iDate,iRun) = corr(superHybridMap,deepHybridMap,'rows','complete');
        
    end
end 
super_DeepCorrNew(1,:)= [];
super_DeepCorrNewT = reshape(super_DeepCorrNew,[28 1]);
super_DeepCorrNewT(super_DeepCorrNewT==0)=[];
super_DeepCorrNewT(~gR) = [];
% super_DeepCorrNewT(super_DeepCorrNewT ==1)=[];
edges = -0.2:0.1:1;
toc;
figure;histogram(super_DeepCorrNewT,edges);
xlabel('Correlation'); ylabel('Count'); grid on;
xlim([-0.3 1.1]); xticks(-0.2:0.1:1);yticks(0:14);


%% Find the distance between probe and ROI center
clear eucDist
for iDate = 1:size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);
    for iRun = 1:size(allRuns{iDate,1})
        clear runName dataDir
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

        imSize = size(imresize(greenIm{iDate,iRun},1/spatialBin));
        greenFig = imresize(greenIm{iDate,iRun},1/spatialBin,'OutputSize',[imSize(1) imSize(2)]);

        % Load probe and ROI center location
        seedLocIn    = load([dataDir '\roiCenterLoc.mat']);
        seedLocProbe = seedLocIn.seedLocProbe;
        seedLocIn    = seedLocIn.seedLocIn;

        %         figure; imagesc(greenFig); colormap gray; axis image off; hold on;
        %         title ('Pick the location of the probe...')
        %         seedLocProbe = ginput(1); seedLocProbe = (round(seedLocProbe)); close gcf;
        %         save([dataDir '\roiCenterLoc.mat'],'seedLocProbe','seedLocIn');

        % Calculate euclidean distance between the two locations
        eucDist(iDate,iRun) = pdist2(seedLocProbe,seedLocIn,'euclidean');
        eucDist(iDate,iRun) = eucDist(iDate,iRun)*(500/round(roiSize{iDate}(iRun)/spatialBin)); 
    end
end

eucDist(eucDist == 0) = NaN;
eucDist = reshape(eucDist,[size(eucDist,1)*size(eucDist,2) 1]);
eucDist(isnan(eucDist)) = [];
eucDist(~goodRunsSpatial) = [];
edges = 0:100:2000;
figure; histogram(eucDist,edges);
xlabel('Distance of center of ROI from probe (\mum)');
ylabel('# of recordings'); box off;
xlim([0 2000]); xticks(0:200:2000);
ylim([0 5]);yticks(0:5);

%% Show probe location and ROI center on green
for iDate = 4%:size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);
    for iRun =7%:size(allRuns{iDate,1})
        clear runName dataDir
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

        imSize = size(imresize(greenIm{iDate,iRun},1/spatialBin));
        greenFig = imresize(greenIm{iDate,iRun},1/spatialBin,'OutputSize',[imSize(1) imSize(2)]);

        % Load probe and ROI center location
        seedLocIn    = load([dataDir '\roiCenterLoc.mat']);
        seedLocProbe = seedLocIn.seedLocProbe;
        seedLocIn    = seedLocIn.seedLocIn;

        figure; imagesc(greenFig); colormap gray; hold on;
        plot(seedLocIn(1),seedLocIn(2),'.w','MarkerSize',15);
        plot(seedLocProbe(1),seedLocProbe(2),'.r','MarkerSize',15);
        axis image off;
        hold off;
    end
end
%% Correlate the FC map with the hybrid map for all lags
%  figure; imagesc(greenFig(seedLocInROI(2)-round(seedRad/2):seedLocInROI(2)+round(seedRad/2),...
%     seedLocInROI(1)-round(seedRad/2):seedLocInROI(1)+round(seedRad/2)));
x = -200:200;
negIdx = (-100<=x)&(x<=0); negVals = x(negIdx);
videoFlag = 0; 
clear peakNegVals peakNegTimes allSpatialCorrMap
% iDate = 2; iRun = 6;
tic;
for iDate = 1:size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);
    for iRun = 1:size(allRuns{iDate,1})
        clear runName dataDir allHybridMaps pDatTemp processedDat10 ...
            allHybridMaps 
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
        clc;disp(['Processing - ' monkeyName ' ' expDate ': ' runName]);
        pDatTemp = processedDat{iDate,iRun}.tempBandPass;

        imSize = size(imresize(greenIm{iDate,iRun},1/spatialBin));
        greenFig = imresize(greenIm{iDate,iRun},1/spatialBin,'OutputSize',[imSize(1) imSize(2)]);

        allHybridMaps = matfile([dataDir '\crossCorrFOV.mat']);

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

        % Get ROI location and FC map
        seedLocIn    = load([dataDir '\roiCenterLoc.mat']);
        seedLocProbe = seedLocIn.seedLocProbe;
        seedLocIn    = seedLocIn.seedLocIn;
        circleRad    = round(roiSize{iDate}(iRun)/(spatialBin)); % 500um radius

        seedSigT     = calculateSeedSignal(greenFig,corrMask,...
            seedLocIn,circleRad,pDatTemp); % Get Gaussian weighted seed signal

        fcMap            = plotCorrMap(seedSigT,pDatTemp,0);
        corrMaskT        = reshape(corrMask,[imSize(1)*imSize(2) 1]);
       
        grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));

        hybridMapsAll = allHybridMaps.spatialProfile;
        mapsAll       = hybridMapsAll.ccFull;
        mapsAlpha     = hybridMapsAll.ccFullAlpha;
        mapsBeta      = hybridMapsAll.ccFullBeta;
        mapsTheta     = hybridMapsAll.ccFullTheta;
        mapsRaw       = hybridMapsAll.ccFullRaw;

%         hybridMapsSuper = allHybridMaps.spatialProfileSuper;
%         mapsAll       = hybridMapsAll.ccFull;
%         mapsAlpha     = hybridMapsAll.ccFullAlpha;
%         mapsBeta      = hybridMapsAll.ccFullBeta;
%         mapsTheta     = hybridMapsAll.ccFullTheta;
%         mapsRaw       = hybridMapsAll.ccFullRaw;
% 
% 
%         hybridMapsDeep = allHybridMaps.spatialProfileDeep;
%         mapsAll       = hybridMapsAll.ccFull;
%         mapsAlpha     = hybridMapsAll.ccFullAlpha;
%         mapsBeta      = hybridMapsAll.ccFullBeta;
%         mapsTheta     = hybridMapsAll.ccFullTheta;
%         mapsRaw       = hybridMapsAll.ccFullRaw;

        if videoFlag
            saveVideoFOV(dataDir,'Gamma',monkeyName,expDate,runName,mapsAll,clipMaskCortex,x); %#ok<UNRCH> 
        end

        % Get the seed location
        %         clear seedLocIn
        %         figure; imagesc(greenIm{iDate,iRun}); colormap gray; axis image off;
        %         title('Pick seed to get FC map');
        %         [seedLocIn(1),seedLocIn(2)] = ginput(1); close gcf;
        %         seedLocIn = floor(seedLocIn./spatialBin);
        %
        %         save([dataDir '\roiCenterLoc.mat'],'seedLocIn','seedLocProbe');
        % Get the FC map for the seed
        % hybridMapsAllNew = reshape(hybridMapsAll.ccFull,[imSize(1)*imSize(2) length(x)]);
        %         figure('units','normalized','outerposition',[0 0 1 1]);
        %         subplot(121); imagesc(fcMap);hold on;
        %         h = imagesc(grayImFull); hold off;
        %         set(h,'AlphaData',~corrMask); hold on;
        %         plot(seedLocIn(1),seedLocIn(2),'.w','MarkerSize',10);
        %         colormap jet;caxis([0 1]); axis image off; colorbar;
        %         title('FC map for seed');

        fcMap     = reshape(fcMap,[361*438 1]);
        mapsAll   = reshape(mapsAll,[401 361*438]);
        mapsAlpha = reshape(mapsAlpha,[401 361*438]);
        mapsBeta  = reshape(mapsBeta ,[401 361*438]);
        mapsTheta = reshape(mapsTheta,[401 361*438]);
        mapsRaw   = reshape(mapsRaw,[401 361*438]);

        mapsAll(:,~corrMaskT) = NaN;
        mapsAlpha(:,~corrMaskT) = NaN;
        mapsBeta(:,~corrMaskT) = NaN;
        mapsTheta(:,~corrMaskT) = NaN;
        mapsRaw(:,~corrMaskT) = NaN;
        fcMap(~corrMaskT)     = NaN;

        for iMap = 1:401
            mapVals(iMap)      = corr(fcMap,mapsAll(iMap,:)','rows','complete');
%             mapValsAlpha(iMap) = corr(fcMap,mapsAlpha(iMap,:)','rows','complete');
%             mapValsBeta(iMap)  = corr(fcMap,mapsBeta(iMap,:)','rows','complete');
%             mapValsTheta(iMap) = corr(fcMap,mapsTheta(iMap,:)','rows','complete');
%             mapValsRaw(iMap)   = corr(fcMap,mapsRaw(iMap,:)','rows','complete');
        end

        if ~exist([dataDir '\FCMap_vs_HybridMaps.png'],'file') 
            figure; subplot(121); imagesc(reshape(fcMap,[361 438])); axis image off; hold on;
            grayIm = cat(3, 0.25.*ones((imSize)),0.25.*ones((imSize)), 0.25.*ones((imSize)));
            h = imagesc(grayIm); hold off; set(h,'AlphaData',~(reshape(corrMaskT,[361 438])));
            colormap(jet); colorbar; caxis([0 1]);
            title('FC map for seed');
            subplot(122);  plot(-200:200,mapVals,'.');ylim([-1 1]); xticks(-200:50:200);xticklabels(-20:5:20);
            xlabel('Lags (s)'); ylabel('Correlation between FC map and hybrid map'); grid on; xline(0);
            sgtitle(strrep([ monkeyName ' ' expDate ' ' runName],'_','\_'));
            f = gcf; exportgraphics(f,[dataDir '\FCMap_vs_HybridMaps.png'],'Resolution',300); close gcf;
        end

        allSpatialCorrMap(iDate,iRun,1,:) = mapVals;
        allSpatialCorrMap(iDate,iRun,2,:) = mapValsAlpha;
        allSpatialCorrMap(iDate,iRun,3,:) = mapValsBeta;
        allSpatialCorrMap(iDate,iRun,4,:) = mapValsTheta;
        allSpatialCorrMap(iDate,iRun,5,:) = mapValsRaw;

        [peakNegVals(iDate,iRun,1),peakNegTimes(iDate,iRun,1)] = min(mapVals(negIdx));
        peakNegTimes(iDate,iRun,1) = negVals(peakNegTimes(iDate,iRun,1))./10;

        [peakNegVals(iDate,iRun,2),peakNegTimes(iDate,iRun,2)] = min(mapValsAlpha(negIdx));
        peakNegTimes(iDate,iRun,2) = negVals(peakNegTimes(iDate,iRun,2))./10;

        [peakNegVals(iDate,iRun,3),peakNegTimes(iDate,iRun,3)] = min(mapValsBeta(negIdx));
        peakNegTimes(iDate,iRun,3) = negVals(peakNegTimes(iDate,iRun,3))./10;

        [peakNegVals(iDate,iRun,4),peakNegTimes(iDate,iRun,4)] = min(mapValsTheta(negIdx));
        peakNegTimes(iDate,iRun,4) = negVals(peakNegTimes(iDate,iRun,4))./10;

        [peakNegVals(iDate,iRun,5),peakNegTimes(iDate,iRun,5)] = min(mapValsRaw(negIdx));
        peakNegTimes(iDate,iRun,5) = negVals(peakNegTimes(iDate,iRun,5))./10;
    end
end
toc;

allSpatialCorrMapR = reshape(allSpatialCorrMap,[size(allSpatialCorrMap,1)*size(allSpatialCorrMap,2) 5 401]);
zeroRows = ~all(squeeze(allSpatialCorrMapR(:,1,:)),2);
allSpatialCorrMapR(zeroRows,:,:) = []; % Recheck this piece of code
allSpatialCorrMapR(~goodRunsSpatial,:,:) = [];

peakNegValsR = reshape (peakNegVals,[size(allSpatialCorrMap,1)*size(allSpatialCorrMap,2) 5]);
peakNegValsR(zeroRows,:) = [];
peakNegValsR(~goodRunsSpatial,:) = [];

peakNegTimesR = reshape(peakNegTimes,[size(allSpatialCorrMap,1)*size(allSpatialCorrMap,2) 5]);
peakNegTimesR(zeroRows,:) = [];
peakNegTimesR(~goodRunsSpatial,:) = [];

% Show the distribution of correlations
figure; boxplot(peakNegValsR(:,1),'Label','Distance of FC map seed from probe (um)');
ylim([-1 0.15]); yticks(-1:0.1:0.1); ylabel('Peak negative correlation');

% Check if correlations vary with distance between ROI and probe
figure; scatter(eucDist,peakNegValsR(:,1),'filled'); hold on;
xlabel('Distance of FC map seed from probe (\mum)');
ylabel('Peak negative correlation');grid on;
xlim([0 2000]); ylim([-1 1]); yticks(-1:0.1:1); xticks(0:200:32000);
[f,gof] = fit(eucDist,double(peakNegValsR(:,1)),'poly1');  % Line fitting
c = coeffvalues(f); r2 = gof.rsquare;
xFit  = linspace(min(eucDist), max(eucDist), 1000);
yFit =  c(1)*xFit + c(2);
plot(xFit,yFit,'-','Color',[0 0.4470 0.7410],'LineWidth',1);
mdl = fitlm(eucDist,double(peakNegValsR(:,1)));
text(1500, 0.8,['R^2 : ' num2str(gof.rsquare*100) '%']);
text(1500,0.7,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);

%% Plot all the cross correlation vs lag for the full FOV
bandLabels = {'Gamma band'; 'Alpha band' ; 'Beta band'; 'Theta band'; 'Spiking'};
figure;
for iBand = 1:5
    subplot(2,3,iBand);
    medTempProfile = squeeze(median(squeeze(allSpatialCorrMapR(:,iBand,:)),1,'omitnan'));
    semAll  = mad(squeeze(allSpatialCorrMapR(:,iBand,:)),0,1)./sqrt(size(squeeze(allSpatialCorrMapR(:,iBand,:)),1));
    plot(x,smooth(medTempProfile,5),'LineWidth',1.5);xlim([-200 200]); hold on;
    patch([x fliplr(x)],[(medTempProfile-2.*semAll) fliplr((medTempProfile+2.*semAll))],...
        'blue','FaceAlpha',0.3,'EdgeColor','none');xline(0); box off; 
    xticks(-200:50:200); xticklabels(-20:5:20); ylim([-1 1]); yticks(-1:0.2:1);
    xlabel('Lag (s)'); ylabel('Cross correlations between FC map and hybrids');
   title(bandLabels{iBand}); 
end

figure;
for iBand = 1:5
    medTempProfile = squeeze(median(squeeze(allSpatialCorrMapR(:,iBand,:)),1,'omitnan'));
    plot(x,smooth(medTempProfile,10),'LineWidth',1.5);hold on;
end

xlim([-200 200]); xticks(-200:50:200); xticklabels(-20:5:20); ylim([-1 1]); yticks(-1:0.2:1);
xlabel('Lag (s)'); ylabel('Cross correlations between FC map and hybrids'); xline(0); box off; 

figure;
for iBand = 1:5
    subplot(2,4,iBand);
    plot(x,squeeze(allSpatialCorrMapR(:,iBand,:)),'Color',[0.65 0.65 0.65],'LineWidth',1);
    xlim([-200 200]); hold on;
    medTempProfile = squeeze(median(squeeze(allSpatialCorrMapR(:,iBand,:)),1,'omitnan'));
    plot(x, smooth(medTempProfile,10),'k','LineWidth',2); % 10 point moving avg filter
    xline(0);xticks(-200: 50: 200); xticklabels(-20:5:20);
    xlabel('Lag (s)'); ylabel('Cross correlation for full FOV');
    ylim([-0.5 0.5]); title(bandLabels{iBand});
    ylim([-1 1]);grid on;
end

figure;
allBandFOVSort1 = peakNegValsR(:,[4 2 3 1 5]);
 subplot(121); boxplot(allBandFOVSort1,{'Theta';'Alpha';'Beta';'Gamma';'Spiking'}); 
ylim([-1 1]); ylabel('Cross correlations between FC map and hybrids');box off; 

allBandFOVSort2 = peakNegValsR(:,[1 3 2 5 4]);
subplot(122);  boxplot(allBandFOVSort2,{'Gamma,';'Beta';'Alpha'; 'Spiking';'Theta'});
ylim([-1 1]); ylabel('Cross correlations between FC map and hybrids'); box off;


% Show distributions of negative correlations
figure;subplot(121); boxplot(peakNegValsR,'Labels',bandLabels);
ylim([-1 1]); ylabel('Minimum Cross correlations between FC map and hybrids');

% Show distributions of lags for different frequency bands
subplot(122); boxplot(peakNegTimesR,'Labels',bandLabels);
ylim([-15 10]); yticks(-15:3:10); ylabel('Lag at max cross correlation');

[pCorr,tblCorr,statsCorr] = anova1(peakNegValsR,bandLabels);
[rCorr,mCorr,~,gnamesCorr] = multcompare(statsCorr,"CriticalValueType","bonferroni");
tblCorrM = array2table(rCorr,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblCorrM.("Group") = gnamesCorr(tblCorrM.("Group"));
tblCorrM.("Control Group") = gnamesCorr(tblCorrM.("Control Group"));

[pLags,tblLags,statsLags] = anova1(peakNegTimesR,bandLabels,'off');
[rLags,mLags,~,gnamesLags] = multcompare(statsLags,"CriticalValueType","bonferroni");
tblLagsM = array2table(rLags,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblLagsM.("Group") = gnamesLags(tblLagsM.("Group"));
tblLagsM.("Control Group") = gnamesLags(tblLagsM.("Control Group"));

%% Calculate FC maps for seeds placed <1mm near the probe - Full FOV controls
%%% and correlate with hybrid map
x = -200:200;
negIdx = (-100<=x)&(x<=0); negVals = x(negIdx);

for iDate = 1:size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);
    for iRun = 1: size(allRuns{iDate,1},1)
        clear runName dataDir clipMask elecMask clipMaskCortex corrMask lowIdx ...
            pDatTemp greenFig seedLocIn crossCorrTemp allHybridMaps mapsAll...
            peakNegHybridMap mapsAllTemp probeCh badTimes szLFP skullMask ...
            inDatSize infraEphys allCortexMask

        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

        % IMAGING: Load the appropriate masks for the imaging data
        % Load clipmask
        if exist([dataDir '\clipMask0' runName(end) '.BMP'],'file') == 0
            clipMask = imread([dataDir '\clipMask0' runName(end) '.png']);
        else
            clipMask = imread([dataDir '\clipMask0' runName(end) '.bmp']);
        end

        % Load the cortex mask
        if exist([dataDir '\skullMask.bmp'],'file') == 0
            skullMask = imread([dataDir '\skullMask.png']);
        else
            skullMask = imread([dataDir '\skullMask.bmp']);
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
        skullMask      = imresize(skullMask,1/3);
        skullMask      = skullMask(:,:,1)>0; 
        allCortexMask  = skullMask & ~elecMask; % Mask of cortex with vessels
        clipMaskCortex = clipMask & ~elecMask; % Mask of cortex without vessels
        
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

        % Get imaging data for the run
        pDatTemp = processedDat{iDate,iRun}.tempBandPass;
        imSize   = size(pDatTemp);
        greenFig = imresize(greenIm{iDate,iRun},1/spatialBin,'OutputSize',[imSize(1) imSize(2)]);
        
        % Get ROI location and FC map
        seedLocIn    = load([dataDir '\roiCenterLoc.mat']);
        seedLocProbe = seedLocIn.seedLocProbe;
        seedLocIn    = seedLocIn.seedLocIn;
        circleRad    = 6;

        seedSigT     = calculateSeedSignal(greenFig,corrMask,...
            seedLocIn,circleRad,pDatTemp); % Get Gaussian weighted seed signal

        fcMap            = plotCorrMap(seedSigT,pDatTemp,0);
        corrMaskT        = reshape(corrMask,[imSize(1)*imSize(2) 1]);
        fcMap            = reshape(fcMap,[361*438 1]);
        fcMap(~corrMaskT) = NaN;

        % Get the hybrid maps and determine lag  for the run
        try
            allHybridMaps  = matfile([dataDir '\crossCorrFOV.mat']);
        catch
            getCrossCorrFOV(monkeyName,expDate,runName,dataDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp);

            allHybridMaps = matfile([dataDir '\crossCorrFOV.mat']);

        end

        mapsAllTemp           = allHybridMaps.spatialProfile;
        mapsAll               = mapsAllTemp.ccFull;
        mapsAll               = reshape(mapsAll,[401 361*438]);
        mapsAll(:,~corrMaskT) = NaN;
        peakNegHybridMap      = squeeze(mapsAll((x == peakNegTimes(iDate,iRun,1).*10),:));

        clc; disp(['Analyzing data for ' monkeyName ' '  expDate ' run: ' runName]);

%         %%% 3. TEMPORAL CONTROL: Shuffle Ephys using variable windows and
%         %%% determine the correlation at peak negative
% 
        if ~exist([dataDir '\tempControlVarsFOV.mat'],'file') 
            % Get infraslow ephys powers and imaging data
            clear  runWiseCorrShuffled runWiseLagShuffled gammaEphys probeCh...
                processedDat10 ch badChannels badTimes szLFP timeStamp...
                badTimeThreshTemp badTimes

            disp('Obtaining temporal controls...');

            % Upsampling imaging data to 10 Hz           
            pDatTemp   = reshape(pDatTemp,[imSize(1)*imSize(2) imSize(3)]);
            
            parfor iP = 1:size(pDatTemp,1)
                processedDat10(iP,:) = interp(pDatTemp(iP,:),5);
            end

            szIm = size(processedDat10,2)*100;


            % Get ephys data
            probeCh               = probe{iRun,iDate}.probeCh;
            ch                     = estChInCortex{iDate}(iRun,:);
            badChannels            = badCh{iDate,iRun};
            badTimes               = badTimesLFP{iDate,iRun};
            probeCh(:,badChannels) = [];
            szLFP                  = size(probeCh,1);

            % Make both matrices equal...
            badTimeThreshTemp = badTimeThresh{iDate,iRun};

            if ~(szLFP == szIm)
                szMin          = min([szLFP, szIm]);
                probeCh        = probeCh(1:szMin,:);
                processedDat10 = processedDat10(:,1:floor(szMin/100));

                badTimes(badTimes>szMin)           = [];
                badTimeThreshTemp(badTimeThreshTemp>szMin) = [];
            else
                szMin = szLFP;
            end

            timeStamp       = probe{iRun,iDate}.timeStamp;
            timeStampSorted = timeStamp- timeStamp(1);
            badTimes10Hz    = unique(badTimeThreshTemp./1000);
            badTimeIm       = [];

            % Identifying frames to be removed from RS-ISOI
            for iT = 1: length(badTimes10Hz)
                badTimeIm(iT) = find((floor(abs(timeStampSorted - badTimes10Hz(iT))*100)./100)<=0.05,1,'first');
            end

            % Remove bad times determined from LFP
            badTimeIm = unique(badTimeIm);
            badTimeIm(badTimeIm>size(processedDat10,2)) = [];
            processedDat10(:,badTimeIm) = [];
            probeCh(badTimes,:) = [];

            % Remove bad times determined visually from spectrogram
            [probeCh,~,processedDat10] = removeBadTimesFromSpec(monkeyName,expDate,runName,probeCh,[],processedDat10);

            % Get infraslow powers
            gammaBand  = [30 90]; [bG,aG] = butter(3,gammaBand./(1e3/2),'bandpass');% Gamma band filtering parameters
            gammaEphys = single(filtfilt(bG,aG,double(probeCh(:,ch(1):ch(2)))));

            % Set the time windows to perform cross correlations
            % between ISOI and Ephys
            tic;
            timeLen = min([size(gammaEphys,1) size(processedDat10,2)*1e2]);
            winLen  = [0.001 0.1 1 5 10 50 100 500 timeLen./1e3].*1e3;

            for iWin = 1:length(winLen) % shuffle time series in windows of 10s
                clear comb1 infraEphysS roiS ccT ccFull z p k sos...
                    enSize envelopeFiltered gammaEphysS
                for iRep = 1: 5
                    clear comb1 infraEphysS roiS ccT ccFull z p k sos...
                        enSize envelopeFiltered gammaEphysS envelopeDat
                    rng('shuffle');
                    comb1 = randperm(round(timeLen/winLen(iWin)));

                    if iWin == 1
                        gammaEphysS  = gammaEphys(comb1,:);
                    elseif iWin == 9
                        gammaEphysS  = gammaEphys;
                    else
                        gammaEphysS = [];
                        for iL = 1:length(comb1)
                            clear win1
                            win1 = (comb1(iL)-1)*winLen(iWin)+1 : (comb1(iL)-1)*winLen(iWin)+winLen(iWin);
                            win1(win1>timeLen) = [];
                            gammaEphysS = [gammaEphysS; gammaEphys(win1,:)];
                        end
                    end

                    % Get infraslow ephys
                    envelopeDat = envelope(abs(gammaEphysS),5);

                    % Bandpass - 0.01 Hz - 0.1 Hz
                    [z,p,k] = butter(3,[0.01 0.1]./(1e3/2),'bandpass');
                    [sos,g] = zp2sos(z,p,k);
                    enSize  = size(envelopeDat);

                    envelopeFiltered = filtfilt(sos,g,double([envelopeDat; envelopeDat; envelopeDat ]));
                    envelopeFiltered = envelopeFiltered(enSize(1)+1:(end-enSize(1)),:);
                    infraEphysS      = mean(single(downsample(envelopeFiltered,100)),2);

                    % Check size of imaging and ephys after this operation
                    clear processedDat10R szMin
                    szIm  = size(processedDat10,2);
                    szLFP = size(infraEphysS,1);
                    if ~(szLFP == szIm)
                        szMin        = min([  szLFP, szIm]);
                        infraEphysS   = infraEphysS(1:szMin,:);
                        processedDat10R = processedDat10(:,1:szMin);
                    else
                        processedDat10R = processedDat10;
                    end

                    parfor iP = 1:size(processedDat10R,1)
                        if iP == 1
                            [ccFull(:,iP),lagFull(:,iP)]  = xcorr(infraEphysS',processedDat10R(iP,:),200,'normalized');
                        else
                            [ccFull(:,iP),~]      = xcorr(infraEphysS',processedDat10R(iP,:),200,'normalized');
                        end
                    end

                    ccFull(:,~corrMaskT) = NaN;

                    for iMap = 1:401; mapVals(iMap) = corr(fcMap,ccFull(iMap,:)','rows','complete'); end

                    [runWiseCorrShuffled(iWin,iRep),runWiseLagShuffled(iWin,iRep)] = min(mapVals(negIdx));
                    runWiseLagShuffled(iWin,iRep) = negVals(runWiseLagShuffled(iWin,iRep))./10;
                end
            end
            toc;
            runWiseCorrAllShuffled{iDate,iRun} = runWiseCorrShuffled;
            runWiseLagAllShuffled{iDate,iRun}  = runWiseLagShuffled;
            corrNegShuffle(iDate,iRun,:)       = median(runWiseCorrShuffled,2,'omitnan');%runWiseCorrShuffled;% 
            corrNegTimes(iDate,iRun,:)         = median(runWiseLagShuffled,2,'omitnan');%runWiseLagShuffled;%

            save([dataDir '\tempControlVarsFOV.mat'],'runWiseCorrShuffled','runWiseLagShuffled');

        else
            clear allVars runWiseCorrShuffled
            allVars = load([dataDir '\tempControlVarsFOV.mat']);
            runWiseCorrAllShuffled{iDate,iRun} = allVars.runWiseCorrShuffled;
            runWiseLagAllShuffled{iDate,iRun}  = allVars.runWiseLagShuffled;
            runWiseCorrShuffled                = allVars.runWiseCorrShuffled;
            corrNegShuffle(iDate,iRun,:)       = median(runWiseCorrShuffled,2,'omitnan');
            corrNegTimes(iDate,iRun,:)         = median(allVars.runWiseLagShuffled,2,'omitnan');
        end

        if ~exist([dataDir '\temporalControlFOV_v2.png'],'file')
           winLen  = {0.001 0.1 1 5 10 50 100 500 'Unshuffled'}; 
            figure; plot(movmean(runWiseCorrShuffled,[0 1]),'Color',[0.65 0.65 0.65],'LineWidth',1); hold on;
            plot(movmean(squeeze(corrNegShuffle(iDate,iRun,:)),[0 1]),'k','LineWidth',2);
%             semAll  = std(runWiseCorrShuffled,0,2)./sqrt(size(runWiseCorrShuffled,2));
%             xVar = [(1:length(winLen)) fliplr((1:length(winLen)))];
%             yVar = [(squeeze(corrNegShuffle(iDate,iRun,:))-2.*semAll)' ...
%                 flipud((squeeze(corrNegShuffle(iDate,iRun,:))+2.*semAll))'];
%             patch(xVar,yVar,'blue','FaceAlpha',0.3,'EdgeColor','none');

            xticklabels(winLen);  hold on;  ylim([-1 0.5]); yticks(-1:0.2:0.5);
            xlabel('Length of window for shuffling (s)'); ylabel('Cross correlation'); grid on;
            title(strrep(['Full FOV temporal control for ' monkeyName ' ' expDate ' ' runName],'_','\_'));
            f = gcf; exportgraphics(f,[dataDir '\temporalControlFOV_v2.png'],'Resolution',300); close gcf;
        end

        %%% 4.SPATIAL CONTROL (FC MAP BASED): Correlating FC maps
        %%% obtained from seeds placed at varying distances from the
        %%% probe location with the hybrid map
        disp('Obtaining spatial controls...');
        pDatTemp = reshape(pDatTemp,[imSize(1) imSize(2) imSize(3)]);
        tic;
        if ~exist([dataDir '\spatialControlVarsFOV.mat'],'file') 
            clear  lagLow  hybridLowMap runWiseSpatialCorr runWiseSpatialTimes locAll corrHybridMap
            minShift  = round(roiSize{iDate}(iRun)/spatialBin);
            theta     = linspace(0,2*pi, round(pi*minShift)); % number of angles
            maxPoints = length(theta);
            distShift = {'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'};

            for iShift = 1:5
                clear loc locShift row col numPoints pixelLoc
                if iShift == 1
                    locShift = round(roiSize{iDate}(iRun)/spatialBin);
                else
                    locShift = round(roiSize{iDate}(iRun)*2*(iShift-1)/spatialBin);
                end

                % Get coordinates for a circle
                loc(:,1) = round(locShift * cos(theta) + seedLocProbe(1));
                loc(:,2) = round(locShift * sin(theta) + seedLocProbe(2));

                % Get the FC map for the locations and correlate with hybrid map
                for iPoint = 1:size(loc,1)
                    clear seedSigT corrMapT mapVals seedRad seed 
                    seedRad =  round(roiSize{iDate}(iRun)/(spatialBin)); % seed radius for FC map- 500um

                    if any((fliplr(loc(iPoint,:))+seedRad)>size(greenFig)) || any(loc(iPoint,:)-seedRad<= 0)
                        runWiseSpatialCorr(iShift,iPoint)  = NaN;
                        runWiseSpatialTimes(iShift,iPoint) = NaN;
                        loc(iPoint,:) = NaN;
                    
                    else
                        seed = loc(iPoint,:);
                        clipMask_seed = ~corrMask(seed(2)-seedRad:seed(2)+seedRad,seed(1)-seedRad:seed(1)+seedRad);
                        allCortex_seed = ~allCortexMask(seed(2)-seedRad:seed(2)+seedRad,seed(1)-seedRad:seed(1)+seedRad);
                        maskSize = size(clipMask_seed,1) *size(clipMask_seed,2);

                        % Check if the masks occupy more than 50% of seed
                        % or if the seeds are close to the edge (the edges
                        % of the image should at max occupy 25% of seed)
                        if ((sum(clipMask_seed,'all')/maskSize)>0.5) || ((sum(allCortex_seed,'all')/maskSize)>0.25)
                            runWiseSpatialCorr(iShift,iPoint)  = NaN;
                            runWiseSpatialTimes(iShift,iPoint) = NaN;
                            loc(iPoint,:) = NaN;
                            continue;

                        else
                            % Get Gaussian weighted seed signal
                            seedSigT = calculateSeedSignal(greenFig,clipMaskCortex,loc(iPoint,:),seedRad,pDatTemp);
                            corrMapT = reshape(plotCorrMap(seedSigT,pDatTemp,0),[imSize(1)*imSize(2) 1]);

                            % Correlate FC map with all hybrid maps. This
                            % correlation makes sure that we remain agnostic
                            % about the peak negative lag and this gives a
                            % distribution of lags and/or correlations. The
                            % variables that change are lags, distance between
                            % probe and seed locations.
                            for iMap = 1:401; mapVals(iMap) = corr(corrMapT,mapsAll(iMap,:)','rows','complete'); end

                            [runWiseSpatialCorr(iShift,iPoint),runWiseSpatialTimes(iShift,iPoint)] = min(mapVals(negIdx));
                            runWiseSpatialTimes(iShift,iPoint) = negVals(runWiseSpatialTimes(iShift,iPoint))./10;

                            % Correlate FC map with peak negative hybrid map -
                            % the lag is fixed here, the distance between the
                            % probe and the seed is the variable that is
                            % changing here.
                            corrHybridMap(iShift,iPoint) = corr(corrMapT,peakNegHybridMap','rows','complete');
                        end
                    end
                    locAll{iShift} = loc;
                end
            end

            save([dataDir '\spatialControlVarsFOV.mat'],'runWiseSpatialCorr','runWiseSpatialTimes','locAll','corrHybridMap');
            spCorrControl{iRun,iDate}      = runWiseSpatialCorr;
            spCorrControlTimes{iRun,iDate} = runWiseSpatialTimes;
            spCorrHybridMap{iRun,iDate}    = corrHybridMap;

            toc;
        else
            clear allSpatialVars
            allSpatialVars = load([dataDir '\spatialControlVarsFOV.mat']);
            spCorrControl{iRun,iDate}      = allSpatialVars.runWiseSpatialCorr;
            spCorrControlTimes{iRun,iDate} = allSpatialVars.runWiseSpatialTimes;
            spCorrHybridMap{iRun,iDate}    = allSpatialVars.corrHybridMap;
            locAll                         = allSpatialVars.locAll;
        end

        if ~exist([dataDir '\spatialControlFOV_v2.png'],'file')|| 1
            cVals = {'w','k','b','g','m'};
            figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(131); imagesc(greenFig); hold on; colormap gray; axis image off;
            plot(seedLocProbe(1),seedLocProbe(2),'Marker','pentagram','MarkerSize',15,...
                'MarkerFaceColor','r','MarkerEdgeColor','none');

            for iShift = 1:5
                if ~isempty(locAll{iShift})
                    plot(locAll{iShift}(:,1),locAll{iShift}(:,2),'.','Color',cVals{iShift},'MarkerSize',10);
                end
            end

            subplot(132);boxplot(spCorrControl{iRun,iDate}',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});
            xlabel('Distance from probe (mm)'); ylabel('Correlation between FC map and peak negative map');
            ylim([-1 0.3]);

            subplot(133);boxplot(spCorrControlTimes{iRun,iDate}',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});
            xlabel('Distance from probe (mm)'); ylabel('Lag at peak negative correlation'); ylim([-11 3]);

            sgtitle(strrep(['FC maps vs Hybrid maps (varying lag) for ' monkeyName ' ' expDate ' ' runName],'_','\_'));
            f = gcf; exportgraphics(f,[dataDir '\spatialControlFOV_v3.png'],'Resolution',300); close gcf;
        end

        if ~exist([dataDir '\spatialControl_HybridMap_v2.png'],'file')
            cVals = {'w','k','b','g','m'};
            figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(121); imagesc(greenFig); hold on; colormap gray; axis image off;
            plot(seedLocProbe(1),seedLocProbe(2),'Marker','pentagram','MarkerSize',15,...
                'MarkerFaceColor','r','MarkerEdgeColor','none');

            for iShift = 1:5
                if ~isempty(locAll{iShift})
                    plot(locAll{iShift}(:,1),locAll{iShift}(:,2),'.','Color',cVals{iShift},'MarkerSize',10);
                end
            end

            subplot(122);boxplot(spCorrHybridMap{iRun,iDate}',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});
            xlabel('Distance from probe (mm)'); ylabel('Correlation between FC map and peak negative map');
            sgtitle(strrep(['FC maps vs Hybrid map (fixed lag) for ' monkeyName ' ' expDate ' ' runName],'_','\_')); ylim([-1 0.3]);

            f = gcf; exportgraphics(f,[dataDir '\spatialControl_HybridMap_v3.png'],'Resolution',300); close gcf;
        end
    end
end

%% Group population data from controls
% Get medians (for each distance group) for each run and combine across all
% runs
clear spCorrNorm medSPCorrControl medSpCorrAll spCorrMinT
spCorrMin = cellfun(@(x) x./abs(min(x,[],'all','omitnan')),spCorrControl,'un',0);  
spCorrMin = reshape(spCorrMin,[size(spCorrControl,1)*size(spCorrControl,2) 1]);
zeroInd   = cell2mat(cellfun(@(x) isempty(x),spCorrMin,'un',0));
spCorrMin(zeroInd) = [];
spCorrMin(~goodRunsSpatial) = [];
spCorrMin = cellfun(@(x) median(x,2,'omitnan'),spCorrMin,'un',0);
spCorrMin = (cat(2,spCorrMin{:}))'; 

figure; boxplot(spCorrMin,{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'}); hold on; 
 ylim([-1 0.5]); yticks(-1:0.2:1);
xlabel('Distance from probe (mm)'); ylabel('Correlation between FC map and peak negative map');

[pSpCorr,tblSpCorr,statsSpCorr] = anova1(spCorrMin,{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'},'off');
[rSpCorr,mSpCorr,~,gnamesSpCorr] = multcompare(statsSpCorr,"CriticalValueType","bonferroni");

tblSpCorrM = array2table(rSpCorr,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblSpCorrM.("Group") = gnamesSpCorr(tblSpCorrM.("Group"));
tblSpCorrM.("Control Group") = gnamesSpCorr(tblSpCorrM.("Control Group"));

spCorrMinT = reshape(cellfun(@(x) x./abs(min(x,[],'all','omitnan')),spCorrControl,'un',0),[size(spCorrControl,1)*size(spCorrControl,2) 1]);
spCorrMinT(zeroInd) = [];
spCorrMinT(~goodRunsSpatial) = [];
spCorrMinT = cat(2,spCorrMinT{:});

figure; boxplot(spCorrMinT',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'}); hold on; 
ylim([-1 0.5]);yticks(-1:0.2:0.5);
xlabel('Distance from probe (mm)'); ylabel('Correlation between FC map and peak negative map');

medSpCorrControl = cellfun(@(x) median(x,2,'omitnan'),spCorrControl,'un',0);
medSpCorrControl = reshape(medSpCorrControl,[size(spCorrControl,1)*size(spCorrControl,2) 1]);
zeroInd = cell2mat(cellfun(@(x) (isempty(x)), medSpCorrControl,'un',0)); 
medSpCorrControl(zeroInd) = [];
medSpCorrControl(~goodRunsSpatial) = [];

for iShift = 1:5
    medSpCorrAll(:,iShift) = cell2mat(cellfun(@(x) x(iShift),medSpCorrControl,'un',0));
end

figure; boxplot(medSpCorrAll,{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'}); hold on; 
med = median(medSpCorrAll); plot(med,'ro-'); ylim([-1 0.5]); yticks(-1:0.2:0.5);
xlabel('Distance from probe (mm)'); ylabel('Correlation between FC map and peak negative map');

[pSpCorr,tblSpCorr,statsSpCorr] = anova1(medSpCorrAll,{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'},'off');
[rSpCorr,mSpCorr,~,gnamesSpCorr] = multcompare(statsSpCorr,"CriticalValueType","bonferroni");

tblSpCorrM = array2table(rSpCorr,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblSpCorrM.("Group") = gnamesSpCorr(tblSpCorrM.("Group"));
tblSpCorrM.("Control Group") = gnamesSpCorr(tblSpCorrM.("Control Group"));

% Another method
medSpControlT = reshape(spCorrControl,[35 1]);
medSpControlT(zeroInd) = [];
medSpControlT(~goodRunsSpatial) = [];
medSpControlT = cat(2,medSpControlT{:});

figure; boxplot(medSpControlT',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'}); hold on; 
med = median(medSpControlT,2,'omitnan'); plot(med,'ro-');ylim([-1 0.5]);yticks(-1:0.2:0.5);
xlabel('Distance from probe (mm)'); ylabel('Correlation between FC map and peak negative map');

%% Grouping temporal controls

clear runWiseCorrAllShuffledT 
runWiseCorrAllShuffledT = reshape(runWiseCorrAllShuffled,[size(runWiseCorrAllShuffled,1)*size(runWiseCorrAllShuffled,2) 1]);
zeroInd   = cell2mat(cellfun(@(x) isempty(x),runWiseCorrAllShuffledT,'un',0));
runWiseCorrAllShuffledT(zeroInd) = [];
runWiseCorrAllShuffledT(~goodRunsSpatial) = [];
runWiseCorrAllShuffledT = cellfun(@(x) x./abs(min(x(9,:))),runWiseCorrAllShuffledT,'un',0)';
runWiseCorrAllShuffledT = cell2mat(cellfun(@(x) median(x,2,'omitnan'),runWiseCorrAllShuffledT,'un',0));

% figure; plot(smoothdata(runWiseCorrAllShuffledT,2,'movmean',10),'Color',[0.65 0.65 0.65],'LineWidth',1); hold on;
figure; plot(smoothdata(runWiseCorrAllShuffledT,2,'movmean',10),'Color',[0.65 0.65 0.65],'LineWidth',1); hold on;
plot(median(smoothdata(runWiseCorrAllShuffledT,2,'movmean',2),2),'k','LineWidth',2);
xticklabels(winLen);  hold on;  ylim([-1 0.5]); yticks(-1:0.2:0.5);box off; ylim([-0.8 -0.1]);
xlabel('Length of window for shuffling (s)'); ylabel('Cross correlation'); grid off;
box off; ylim([-1 -0.3])


clear corrNegShuffleT

corrNegShuffleT = reshape(corrNegShuffle,[size(corrNegShuffle,1)*size(corrNegShuffle,2) 9]);
zeroInd =  find((corrNegShuffleT(:,1)) == 0); 
corrNegShuffleT(zeroInd,:) = [];
corrNegShuffleT(~goodRunsSpatial,:) = [];

winLen  = {0.001 0.1 1 5 10 50 100 500 'Unshuffled'};
figure; plot(smoothdata(corrNegShuffleT',2,'movmean',8),'Color',[0.65 0.65 0.65],'LineWidth',1); hold on;
plot(median(smoothdata(corrNegShuffleT,2,'movmean',3),1),'k','LineWidth',2);
xticklabels(winLen);  hold on;  ylim([-1 0.5]); yticks(-1:0.2:0.5);
xlabel('Length of window for shuffling (s)'); ylabel('Cross correlation'); grid off;
box off; ylim([-0.8 -0.1]);






%% Calculate FC maps for seeds placed <1mm near the probe - ROI controls
%%% and correlate with hybrid map
for iDate = 1:2%size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);
    for iRun = 1: size(allRuns{iDate,1})
        clear runName dataDir clipMask elecMask clipMaskCortex corrMask x negIdx lowIdx
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
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
        clc; disp(['Analyzing data for ' monkeyName ' '  expDate ' run: ' runName]);
        clear rad col row circlePixels pixelLoc spCorrAllT
        rad          = round(roiSize{iDate}(iRun)*2/spatialBin); % 1mm
        [col, row]   = meshgrid(1:imSize(2), 1:imSize(1));
        circlePixels = (row - seedLocProbe(2)).^2 + (col - seedLocProbe(1)).^2 <= rad.^2;
        [pixelLoc(:,2), pixelLoc(:,1)] = find(circlePixels);

        % Remove pixels that are within seed (seed diameter  = 500 um)
        numPoints = length(pixelLoc); iPixel = 1;
        circleRad =  round(roiSize{iDate}(iRun)/(spatialBin*2)); % seed radius for FC map

        while iPixel <= numPoints
            clear locs circlePixelsTemp
            circlePixelsTemp       = (row - pixelLoc(iPixel,2)).^2 + (col - pixelLoc(iPixel,1)).^2 <= circleRad.^2;
            [locs(:,2), locs(:,1)] = find(circlePixelsTemp);
            repPoints              = ismember(pixelLoc,locs,'rows');
            repPoints(iPixel)      = 0; % Reference pixel should not be removed
            pixelLoc(repPoints,:)  = [];
            numPoints              = length(pixelLoc);
            iPixel                 = iPixel+1;
        end

        clear sz1 sz2 lagLow frameNumMedLow hybridLowMap
        sz1 = imSize(1); sz2 = imSize(2);
        x = allLags{iRun,iDate};
        lagLow         = tempProfile(iM,iDate,iRun).lagLow;
        frameNumMedLow = (x == lagLow);
        hybridLowMap    = reshape(squeeze(allXcorr.ccFull(frameNumMedLow,:,:)),[imSize(1)*imSize(2) 1]);

        clear pDatTemp
        pDatTemp = processedDat{iDate,iRun}.tempBandPass;

        disp('Correlating hybrid maps with FC maps...');
        tic;
        parfor iMap = 1:size(pixelLoc,1) % Pixels inside 1mm radius
            seedSigT           = calculateSeedSignal(greenFig,corrMask,...
                pixelLoc(iMap,:),circleRad,pDatTemp); % Get Gaussian weighted seed signal
            corrMapT           = reshape(plotCorrMap(seedSigT,pDatTemp,0),[sz1*sz2 1]);
            spCorrAllT(iMap,1) = corr(corrMapT,hybridLowMap,'rows','complete');
        end
        toc;
        spCorr(iRun,iDate,iM)    = median(spCorrAllT,'all','omitnan'); % Take median correlations
        spCorrSort               = sort(spCorrAllT,'ascend');
        spCorrTen(iRun,iDate,iM) = median(spCorrSort(ceil(length(spCorrSort)*0.1)),'all','omitnan'); % Top 10%

        %%% 3. TEMPORAL CONTROL: Shuffle Ephys/ISOI or both and
        %%% determine the correlation at peak negative
        if ~exist([dataDir '\temporalControl.png'],'file')
            % Get infraslow ephys powers and imaging data
            clear probeCh badTimes szLFP pDatTemp imSize processedDat10 inDatSize

            probeCh = probe{iRun,iDate}.probeCh;
            probeCh(:,badCh{iDate,iRun}) = [];
            badTimes = badTimesLFP{iDate,iRun};

            szLFP    = size(probeCh,1);
            pDatTemp = processedDat{iDate,iRun}.tempBandPass;
            imSize   = size(pDatTemp);
            pDatTemp = reshape(pDatTemp,[imSize(1)*imSize(2) imSize(3)]);

            % Upsampling imaging data to 10 Hz
            disp('Upsampling imaging data to 10 Hz... ');
            parfor iP = 1:size(pDatTemp,1)
                processedDat10(iP,:) = interp(pDatTemp(iP,:),5);
            end

            % Get the ROI
            processedDat10R = reshape(processedDat10,[imSize(1) imSize(2) size(processedDat10,2)]);

            inDat = processedDat10R(seedLocIn(2)-round(seedRad/2):seedLocIn(2)+round(seedRad/2),...
                seedLocIn(1)-round(seedRad/2):seedLocIn(1)+round(seedRad/2),:);

            inDatSize = size(inDat);
            inDat     = reshape(inDat,[inDatSize(1)*inDatSize(2) inDatSize(3)]);

            % Make both matrices equal...
            szIm = size(processedDat10,2)*100;
            if ~(szLFP == szIm)
                szMin     = min([szLFP, szIm]);
                probeCh     = probeCh(1:szMin,:);
                badTimes(badTimes>szMin) = [];
            else
                szMin = szLFP;
            end

            probeCh(badTimes,:) = [];

            % IMAGING: Upsample ROI (x100 times) and remove concurrent bad
            % frames
            clear inDatUp
            tic;
            parfor iP = 1:size(inDat,1)
                inDatTemp = interp(inDat(iP,:),100);
                if ~(szLFP == szIm); inDatTemp = inDatTemp(:,1:szMin); end
                inDatTemp(:,badTimes) = []; % Remove bad frames
                inDatUp(iP,:) = downsample(inDatTemp,100); % Downsample imaging data
            end
            toc;

            % Remove bad times determined visually from spectrogram
            [probeCh,~,inDatUp] = removeBadTimesFromSpec(monkeyName,expDate,runName,probeCh,[],inDatUp);
            gammaBand  = [30 90]; [bG,aG] = butter(3,gammaBand./(1e3/2),'bandpass');% Gamma band filtering parameters
            infraEphys = getInfraSlowPowerLFP(probeCh,bG,aG,ch); % Gamma band
            infraEphys = mean(infraEphys,2);

            % Set the time windows to perform cross correlations
            % between ISOI and Ephys
            timeLen = length(infraEphys);
            winLen  = [1 5 10 50 100 500 1000 5000 timeLen];

            lagLow         = tempProfile(iM,iDate,iRun).lagLow; % Get peak
            x = allLags{iRun,iDate};
            frameNumMedLow = (x == lagLow);

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

        %%% SPATIAL CONTROL (FC MAP BASED): Correlating FC maps
        %%% obtained from seeds placed at varying distances from the
        %%% probe location with the hybrid map

        if ~exist([dataDir '\spatialControl.png'],'file')
            clear pDatTemp lagLow frameNumMedLow hybridLowMap
            pDatTemp = processedDat{iDate,iRun}.tempBandPass;

            lagLow         = tempProfile(iM,iDate,iRun).lagLow;
            frameNumMedLow = (x == lagLow);
            hybridLowMap    = reshape(squeeze(allXcorr.ccFull(frameNumMedLow,:,:)),[imSize(1)*imSize(2) 1]);

            cVals = {'w','k','b','g','m'};
            figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(131); imagesc(greenFig); hold on; colormap gray; axis image off;
            plot(seedLocProbe(1),seedLocProbe(2),'Marker','pentagram','MarkerSize',15,...
                'MarkerFaceColor','r','MarkerEdgeColor','none');

            minShift  = round(roiSize{iDate}(iRun)/spatialBin);
            theta     = linspace(0,2*pi, round(pi*minShift)); % number of angles
            maxPoints = length(theta);

            for iShift = 1:5
                clear locShift row col numPoints pixelLoc
                if iShift == 1
                    locShift = round(roiSize{iDate}(iRun)/spatialBin);
                else
                    locShift = round(roiSize{iDate}(iRun)*2*(iShift-1)/spatialBin);
                end

                % Get coordinates for a circle
                loc(:,1) = round(locShift * cos(theta) + seedLocProbe(1));
                loc(:,2) = round(locShift * sin(theta) + seedLocProbe(2));

                % Get the FC map for the locations and correlate with hybrid map
                for iPoint = 1:size(loc,1)
                    clear seedSigT corrMapT
                    if sum(fliplr(loc(iPoint,:))+12>size(greenFig)) || (sum(loc(iPoint,:)-12<= 0))
                        spCorrControl{iRun,iDate}(iShift,iPoint) = NaN;
                        loc(iPoint,:) = NaN;
                    else
                        % Get Gaussian weighted seed signal
                        seedSigT = calculateSeedSignal(greenFig,clipMaskCortex,loc(iPoint,:),12,pDatTemp);
                        corrMapT = reshape(plotCorrMap(seedSigT,pDatTemp,0),[imSize(1)*imSize(2) 1]);

                        spCorrControl{iRun,iDate}(iShift,iPoint) =  corr(corrMapT,hybridLowMap,'rows','complete');
                        spCorrDiff{iRun,iDate}(iShift,iPoint)    = spCorr(iRun,iDate) - corr(corrMapT,hybridLowMap,...
                            'rows','complete');
                    end
                end
                % Plot the points sampled on the green blood vessel map
                plot(loc(:,1),loc(:,2),'.','Color',cVals{iShift},'MarkerSize',10);
            end

            subplot(132);boxplot(spCorrControl{iRun,iDate}',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});
            xlabel('Distance from probe (mm)'); ylabel('Correlation between FC map and peak negative map');

            subplot(133);boxplot(spCorrDiff{iRun,iDate}',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});
            xlabel('Distance from probe (mm)'); ylabel(['Deviation from the actual correlation between FC map' ...
                ' and peak negative map']);

            f = gcf; exportgraphics(f,[dataDir '\spatialControl.png'],'Resolution',300); close gcf;
        end

        %%% 5.SPATIAL CONTROL (ROI BASED): Mismatch ROI across sessions
        %%% and compute the temporal profile

    end
end

%%  Show the distribution of correlation between FC map and hybrid map
spCorrR = reshape(squeeze(spCorr(:,:,iM)),[size(spCorr,1)*size(spCorr,2) 1]);
spCorrR(spCorrR == 0) = [];

spCorrTenR = reshape(squeeze(spCorrTen(:,:,iM)),[size(spCorrTen,1)*size(spCorrTen,2) 1]);
spCorrTenR(spCorrTenR == 0) = [];

figure; subplot(121);boxplot([spCorrR(goodRuns) spCorrTenR(goodRuns)],'Labels',{'Median correlation','Top 10% correlation'});
ylim([-1 1]); title('Accepted runs'); ylabel('Correlation between FC map and hybrid maps');
subplot(122); boxplot([spCorrR(~goodRuns) spCorrTenR(~goodRuns)],'Labels',{'Median correlation','Top 10% correlation'});
ylim([-1 1]); title('Rejected runs');

% end

%% Functions used in this script %%
%% Function 1: Save hybrid maps as a video
