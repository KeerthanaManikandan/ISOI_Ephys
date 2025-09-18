% Plots all compiled data for spatial, temporal, laminar analysis and other
% figures in the paper. Modified from ephysImaging_plotting.m
% June 13,2025 - KM
% ISOI vs LFP paper figures

clc; clear;
commonDir = 'C:\Users\kem294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\ISOI_Ephys\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12\fly_track\videoIO']));

%% Load all compiled variables
clear allMonkeyVars hemisphere spatialBin
% Initialize variables 
hemisphere = 'Left'; spatialBin = 3;
monkeys = {'CharlieSheen';'Whiskey'};

% Get experiment dates
dateVals1 = char(caldiff(datetime({'11/29/2021' ; '01/11/2022'; '08/07/2023'; '11/20/2023'},'InputFormat','MM/dd/yyyy'),'Days')); 
dateVals2 = char(caldiff(datetime({'08/14/2023';'10/16/2023';'12/04/2023';'02/20/2024';'04/29/2024'},'InputFormat','MM/dd/yyyy'),'Days'));

dateVals1 = cumsum([1; str2num(dateVals1(:,1:end-1))]); %#ok<*ST2NM> 
dateVals2 = cumsum([1; str2num(dateVals2(:,1:end-1))]);

% Load monkey variables
for iM = 1:2
    allMonkeyVars(iM) =  load(['D:\Data\' monkeys{iM} '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat']); 
end
bandLabels = {'Theta';'Alpha';'Beta';'Gamma';'Spiking'}; 

%% Temporal correlations 
% Show the median cross-correlation vs lag for ROI 
clear tempProfiles

for iM = 1:3
    if iM~=3
        tempProfiles = allMonkeyVars(iM).tempProfilesAll; % Single animal
        monkeyName   = monkeys{iM};
    else
        tempProfiles = [allMonkeyVars(:).tempProfilesAll]; % Both animals
        monkeyName   = 'Combined data'; 
    end

    figure; 
    for iBand = 1:5 % Frequency bands 
        clear semProfile medTempProfileAll
        subplot(2,3,iBand);   
        semProfile = std(tempProfiles(:,:,iBand),0,2)./sqrt(size(tempProfiles,2)); 
        medTempProfileAll = median(tempProfiles(:,:,iBand),2,'omitnan');
       
        plot(-200:200,smooth(medTempProfileAll,3),'k','LineWidth',1); hold on;
        patch([-200:200 fliplr(-200:200)], [(medTempProfileAll- 2.*semProfile);...
            flipud((medTempProfileAll+ 2.*semProfile))],'blue','FaceAlpha',0.3,'EdgeColor','none')
        
        title(bandLabels{iBand}); box off; xline(0);
        xticks(-200:50:200);xticklabels(-20:5:20);ylim([-0.35 0.2]); yticks(-0.5:0.1:0.5);
    end
     sgtitle(monkeyName);
end 

%% Show distributions of temporal correlations 
clear roiCorr
for iM = 1:3
clear roiCorr
    if iM~=3
        roiCorr = allMonkeyVars(iM).allChCorr; % Both monkeys
        monkeyName   = monkeys{iM};
    else
        roiCorr = [allMonkeyVars(1).allChCorr ; allMonkeyVars(2).allChCorr]; % Single monkey data
        monkeyName   = 'Combined data'; 
    end

    figure; 
    boxplot(roiCorr,bandLabels); ylim([-0.6 0.15]); box off;   
    title(monkeyName);
end

%% Show distribution of temporal correlations separated by cortical area
clear allROICorr

allROICorr   = [allMonkeyVars(1).allChCorr ; allMonkeyVars(2).allChCorr];
smFlagT      = [allMonkeyVars(1).smFlagROI; allMonkeyVars(2).smFlagROI];
ssROICorr    = allROICorr(smFlagT == 'S',:);
motorROICorr = allROICorr(smFlagT == 'M',:);

% Bar plot version 1
figure; bar(-[median(ssROICorr,1,'omitnan')' median(motorROICorr,1)']);
box off; xticklabels(bandLabels); legend('Somatosensory','Motor');
ylabel('Temporal correlations'); ylim([0 0.4]);

% Boxplot version 1
figure; subplot(121); boxplot(-ssROICorr,bandLabels);  ylim([-0.15 0.6]); ylabel('Temporal correlations'); box off; title('Somatosensory');
subplot(122); boxplot(-motorROICorr,bandLabels); ylim([-0.15 0.6]); ylabel('Temporal correlations'); box off; title('Motor');

% Boxplot version 2
ssROICorr(15:21,:) = NaN;
fullMat = reshape([ssROICorr; motorROICorr],size(ssROICorr,1),[]); 
figure; boxplot(-fullMat); box off;  ylim([-0.15 0.6]); ylabel('Temporal correlations'); 
% xticks(1:5); xticklabels(bandLabels);



%% Show temporal correlation distributions for ROI for each compartment 
for iM = 1:3
clear super mid deep
    if iM~=3 % Single animal
        super      = allMonkeyVars(iM).allChCorrSuper;
        mid        = allMonkeyVars(iM).allChCorrMid;
        deep       = allMonkeyVars(iM).allChCorrDeep;
        smFlagTemp = allMonkeyVars(iM).smFlagROI;
        monkeyName = monkeys{iM};

    else % Both animals
        super      = [allMonkeyVars(1).allChCorrSuper; allMonkeyVars(2).allChCorrSuper];
        mid        = [allMonkeyVars(1).allChCorrMid;   allMonkeyVars(2).allChCorrMid];
        deep       = [allMonkeyVars(1).allChCorrDeep;  allMonkeyVars(2).allChCorrDeep];
        smFlagTemp = [allMonkeyVars(1).smFlagROI;      allMonkeyVars(2).smFlagROI];
        monkeyName = 'Combined data'; 
    end

    nanRows             = find(isnan(mid(:,1)));
    super(nanRows,:)    = [];
    mid(nanRows,:)      = [];
    deep(nanRows,:)     = [];
    smFlagTemp(nanRows) = [];

 figure;
    for iType = 1:3 % Separate by cortical area 
        switch iType 
            case 1 % Somatosensory 
                superTemp = super(smFlagTemp =='S',:);
                midTemp   = mid(smFlagTemp =='S',:);
                deepTemp  = deep(smFlagTemp=='S',:);
                typeLabel = 'Sensory';

            case 2 % Motor
                superTemp = super(smFlagTemp =='M',:);
                midTemp   = mid(smFlagTemp =='M',:);
                deepTemp  = deep(smFlagTemp=='M',:);
                typeLabel = 'Motor';

            case 3 % Both
                superTemp = super;
                midTemp   = mid;
                deepTemp  = deep;
                typeLabel = 'All sites';
        end
        
        % Plot the heatmap of median correlations
        figPlot = [median(superTemp,1,'omitnan');median(midTemp,1,'omitnan' ); median(deepTemp,1,'omitnan')];
        subplot(1,3,iType); imagesc(figPlot); xticks(1:5); xticklabels(bandLabels); yticks(1:3);
        yticklabels({'Superficial';'Middle';'Deep'}); colorbar;colormap(flipud(jet)); 
        clim([-0.35 0]); title(typeLabel); sgtitle(monkeyName); axis square;
    end
end



%% Spatial correlations
% Show the median cross-correlation vs lag for FOV 
clear fovProfile
iType = 2; % Re-referencing type

for iM = 1:3
clear roiCorr
    if iM~=3 % Single animal
        fovProfile = -squeeze(allMonkeyVars(iM).corrFCHybridT(:,:,iType,:)); 
        monkeyName   = monkeys{iM};

    else %  Both animals
        fovProfile = -[squeeze(allMonkeyVars(1).corrFCHybridT(:,:,iType,:)) ; squeeze(allMonkeyVars(2).corrFCHybridT(:,:,iType,:))];
        monkeyName   = 'Combined data'; 
    end

    figure;
    for iBand = 1:5
        clear semProfile medTempProfileAll
        subplot(2,3,iBand);   
        semProfile = (std(squeeze(fovProfile(:,iBand,:)),0,1)./sqrt(size(fovProfile,1)))';
        medTempProfileAll = median(squeeze(fovProfile(:,iBand,:)),1,'omitnan')';
       
        plot(-200:200,smooth(medTempProfileAll,7),'k','LineWidth',1); hold on;
        patch([-200:200 fliplr(-200:200)], [(medTempProfileAll- 2.*semProfile);...
            flipud((medTempProfileAll+ 2.*semProfile))],'blue','FaceAlpha',0.3,'EdgeColor','none')
        
        title(bandLabels{iBand}); box off; xline(0);
        xticks(-200:50:200);xticklabels(-20:5:20);ylim([-1 1]); yticks(-1:0.1:1);
    end
     sgtitle(monkeyName);
end

%% Show distribution of spatial correlations
clear fovCorr iType
iType = 2; % Re-referencing type
for iM = 1:3
clear fovCorr
    if iM~=3 % Both animals
        fovCorr = allMonkeyVars(iM).peakNegValsAllT(:,:,iType);
        monkeyName   = monkeys{iM};
    else % Single animal
        fovCorr = [allMonkeyVars(1).peakNegValsAllT(:,:,iType) ; allMonkeyVars(2).peakNegValsAllT(:,:,iType)];
        monkeyName   = 'Combined data'; 
    end

    figure; 
    boxplot(fovCorr,bandLabels); ylim([-1 1]); box off;    
    title(monkeyName);
end

%% Show distribution of spatial correlations separated by cortical area
clear allFOVCorr
iType = 2;
allFOVCorr   = [allMonkeyVars(1).peakNegValsAllT(:,:,iType) ; allMonkeyVars(2).peakNegValsAllT(:,:,iType)];
smFlagT      = [allMonkeyVars(1).smFlagROI; allMonkeyVars(2).smFlagROI];
ssFOVCorr    = allFOVCorr(smFlagT == 'S',:);
motorFOVCorr = allFOVCorr(smFlagT == 'M',:);

% Bar plot version
figure; bar(-[median(ssFOVCorr,1)' median(motorFOVCorr,1)']);
box off; xticklabels(bandLabels); legend('Somatosensory','Motor');
ylabel('Spatial correlations'); ylim([0 0.8]);

% Boxplot version
figure; subplot(121); boxplot(-ssFOVCorr,bandLabels); ylim([-0.3 1]); ylabel('Spatial correlations'); box off; title('Somatosensory');
subplot(122); boxplot(-motorFOVCorr,bandLabels); ylim([-0.3 1]); ylabel('Spatial correlations'); box off; title('Motor');

% Boxplot version 2
ssFOVCorr(15:21,:) = NaN;
fullMat = reshape([ssFOVCorr; motorFOVCorr],size(ssFOVCorr,1),[]); 
figure; boxplot(-fullMat); box off;  ylim([-0.3 1]); ylabel('Spatial correlations'); 
% xticks(1:5); xticklabels(bandLabels);


%% Show spatial correlation distributions for each compartment
% Boxplots of frequencies for each layer
clear super mid deep
iType =2; % Re-referencing type
for iM = 1:3
    clear super mid deep
    if iM~=3 % Single monkey
        super      = allMonkeyVars(iM).corrFCSuperT(:,:,iType);
        mid        = allMonkeyVars(iM).corrFCMidT(:,:,iType);
        deep       = allMonkeyVars(iM).corrFCDeepT(:,:,iType);
        smFlagTemp = allMonkeyVars(iM).smFlagFOV;
        monkeyName = monkeys{iM};

    else % Both monkeys
        super      = [allMonkeyVars(1).corrFCSuperT(:,:,iType); allMonkeyVars(2).corrFCSuperT(:,:,iType)];
        mid        = [allMonkeyVars(1).corrFCMidT(:,:,iType); allMonkeyVars(2).corrFCMidT(:,:,iType)];
        deep       = [allMonkeyVars(1).corrFCDeepT(:,:,iType); allMonkeyVars(2).corrFCDeepT(:,:,iType)];
        smFlagTemp = [allMonkeyVars(1).smFlagFOV; allMonkeyVars(2).smFlagFOV];
        monkeyName = 'Combined data';
    end

    for iArea = 1:3 % Cortical area
        switch iArea
            case 1 % Somatosensory
                superTemp = super(smFlagTemp =='S',:);
                midTemp   = mid(smFlagTemp =='S',:);
                deepTemp  = deep(smFlagTemp=='S',:);
                typeLabel = 'Sensory';

            case 2 % Motor 
                superTemp = super(smFlagTemp =='M',:);
                midTemp   = mid(smFlagTemp =='M',:);
                deepTemp  = deep(smFlagTemp=='M',:);
                typeLabel = 'Motor';

            case 3 % Both
                superTemp = super;
                midTemp   = mid;
                deepTemp  = deep;
                typeLabel = 'All sites'; 
        end

        figure;
        subplot(131); boxplot(superTemp,bandLabels); ylim([-1 1]); box off; title('Superficial');
        subplot(132); boxplot(midTemp,bandLabels); ylim([-1 1]); box off; title('Middle');
        subplot(133); boxplot(deepTemp,bandLabels); ylim([-1 1]); box off; title('Deep');
        sgtitle([monkeyName ' - ' typeLabel]);
    end
end

%% Show correlation with FC map for superficial/middle/deep layers for Sensory/Motor areas
% Matrix of Layers vs frequencies for S/M sites
clear super mid deep
iType = 2;
for iM = 1:3
clear super mid deep
    if iM~=3
        super      = allMonkeyVars(iM).corrFCSuperT(:,:,iType);
        mid        = allMonkeyVars(iM).corrFCMidT(:,:,iType);
        deep       = allMonkeyVars(iM).corrFCDeepT(:,:,iType);
        smFlagTemp = allMonkeyVars(iM).smFlagFOV;
        monkeyName = monkeys{iM};
    else
        super      = [allMonkeyVars(1).corrFCSuperT(:,:,iType); allMonkeyVars(2).corrFCSuperT(:,:,iType)];
        mid        = [allMonkeyVars(1).corrFCMidT(:,:,iType); allMonkeyVars(2).corrFCMidT(:,:,iType)];
        deep       = [allMonkeyVars(1).corrFCDeepT(:,:,iType); allMonkeyVars(2).corrFCDeepT(:,:,iType)];
        smFlagTemp = [allMonkeyVars(1).smFlagFOV; allMonkeyVars(2).smFlagFOV];
        monkeyName = 'Combined data';
    end

    medAll     = [median(super,1,'omitnan'); median(mid,1,'omitnan'); median(deep,1,'omitnan')];
    sensoryAll = [median(super(smFlagTemp=='S',:),1,'omitnan'); median(mid(smFlagTemp=='S',:),1,'omitnan'); median(deep(smFlagTemp=='S',:),1,'omitnan')];
    motorAll   = [median(super(smFlagTemp=='M',:),1,'omitnan'); median(mid(smFlagTemp=='M',:),1,'omitnan'); median(deep(smFlagTemp=='M',:),1,'omitnan')];

    figure; 
    subplot(131); imagesc(-sensoryAll); xticks(1:5); yticks(1:3);
    xticklabels(bandLabels); yticklabels({'Superficial';'Middle';'Deep'});
    colormap((jet));clim([0 0.65 ]); colorbar;title('Sensory'); axis image square; 

    subplot(132); imagesc(-motorAll);xticks(1:5); yticks(1:3);
    xticklabels(bandLabels); yticklabels({'Superficial';'Middle';'Deep'});
    colormap((jet));clim([0 0.65]); colorbar;title('Motor'); axis image square;

    subplot(133); imagesc(-medAll); xticks(1:5); yticks(1:3);
    xticklabels(bandLabels); yticklabels({'Superficial';'Middle';'Deep'});
    colormap((jet));clim([0 0.65]); colorbar;title('All sites');axis image square;
    sgtitle(monkeyName);

    % Show the distributions for all frequencies in the same plot
    clear sensoryOrg motorOrg allOrg
    sensoryOrg = [super(smFlagTemp=='S',1) mid(smFlagTemp=='S',1) deep(smFlagTemp=='S',1) super(smFlagTemp=='S',2) mid(smFlagTemp=='S',2) deep(smFlagTemp=='S',2) ...
        super(smFlagTemp=='S',3) mid(smFlagTemp=='S',3) deep(smFlagTemp=='S',3) super(smFlagTemp=='S',4) mid(smFlagTemp=='S',4) deep(smFlagTemp=='S',4) ...
        super(smFlagTemp=='S',5) mid(smFlagTemp=='S',5) deep(smFlagTemp=='S',5)];

    motorOrg = [super(smFlagTemp=='M',1) mid(smFlagTemp=='M',1) deep(smFlagTemp=='M',1) super(smFlagTemp=='M',2) mid(smFlagTemp=='M',2) deep(smFlagTemp=='M',2) ...
        super(smFlagTemp=='M',3) mid(smFlagTemp=='M',3) deep(smFlagTemp=='M',3) super(smFlagTemp=='M',4) mid(smFlagTemp=='M',4) deep(smFlagTemp=='M',4) ...
        super(smFlagTemp=='M',5) mid(smFlagTemp=='M',5) deep(smFlagTemp=='M',5)];

    allOrg = [super(:,1) mid(:,1) deep(:,1) super(:,2) mid(:,2) deep(:,2) ...
        super(:,3) mid(:,3) deep(:,3) super(:,4) mid(:,4) deep(:,4) ...
        super(:,5) mid(:,5) deep(:,5)];

    labels = {'Theta-Super'; 'Theta-Middle';'Theta-Deep'; 'Alpha-Super'; 'Alpha-Mid'; 'Alpha-Deep'; 'Beta-Super'; 'Beta-Mid'; 'Beta-deep';...
        'Gamma-Super'; 'Gamma-Mid';'Gamma-Deep'; 'Spiking-Super'; 'Spiking-Mid';'Spiking-Deep'};
    
    figure; 
    subplot(131); boxplot(sensoryOrg,labels); box off; ylim([-1 0.8]); title('Sensory');
    subplot(132); boxplot(motorOrg,labels);   box off; ylim([-1 0.8]); title('Motor');
    subplot(133); boxplot(allOrg,labels);     box off; ylim([-1 0.8]); title('All areas');

end

%% Show distributions of correlations between compartment maps 
clear superDeep superMid midDeep 
iType = 2;
for iM = 1:3
    if iM~=3 % Single monkey
        superDeep  = allMonkeyVars(iM).super_DeepAvgFramesT(:,:,iType);
        superMid   = allMonkeyVars(iM).super_MidAvgFramesT(:,:,iType);
        midDeep    = allMonkeyVars(iM).deep_MidAvgFramesT(:,:,iType);
        controls   = allMonkeyVars(iM). superDeepPosCheckT;
        monkeyName = monkeys{iM};

    else % Both monkeys
        superDeep  = [allMonkeyVars(1).super_DeepAvgFramesT(:,:,iType); allMonkeyVars(2).super_DeepAvgFramesT(:,:,iType)];
        superMid   = [allMonkeyVars(1).super_MidAvgFramesT(:,:,iType); allMonkeyVars(2).super_MidAvgFramesT(:,:,iType)];
        midDeep    = [allMonkeyVars(1).deep_MidAvgFramesT(:,:,iType); allMonkeyVars(2).deep_MidAvgFramesT(:,:,iType)];
        controls   = [allMonkeyVars(1). superDeepPosCheckT; allMonkeyVars(2). superDeepPosCheckT];
        monkeyName = 'Combined data';
    end 
    
   figure;
   for iBand = 1:5
       if iBand>= 2 && iBand<= 4
           subplot(2,3,iBand); boxplot([superMid(:,iBand) midDeep(:,iBand) superDeep(:,iBand) controls(:,iBand-1)],{'S/M';'M/D';'S/D';'Controls'});
       else
           subplot(2,3,iBand); boxplot([superMid(:,iBand) midDeep(:,iBand) superDeep(:,iBand)],{'S/M';'M/D';'S/D'});
       end
       title(bandLabels{iBand}); box off; ylim([-0.8 1.2]);
   end
   sgtitle(monkeyName); 
end 


%% Show correlations between cross-modal and FC for each compartment as a heatmap
clear super mid deep
iType = 2;
for iM = 1:3
clear super mid deep
    if iM~=3 % Single animal
        super      = squeeze(allMonkeyVars(iM).superHybridAllBandsT(:,iType,:,:));
        mid        = squeeze(allMonkeyVars(iM).midHybridAllBandsT(:,iType,:,:));
        deep       = squeeze(allMonkeyVars(iM).deepHybridAllBandsT(:,iType,:,:));
        smFlagTemp = allMonkeyVars(iM).smFlagFOV;
        monkeyName = monkeys{iM};

    else % Both animals
        super      = [squeeze(allMonkeyVars(1).superHybridAllBandsT(:,iType,:,:)); squeeze(allMonkeyVars(2).superHybridAllBandsT(:,iType,:,:))];
        mid        = [squeeze(allMonkeyVars(1).midHybridAllBandsT(:,iType,:,:)); squeeze(allMonkeyVars(2).midHybridAllBandsT(:,iType,:,:))];
        deep       = [squeeze(allMonkeyVars(1).deepHybridAllBandsT(:,iType,:,:)); squeeze(allMonkeyVars(2).deepHybridAllBandsT(:,iType,:,:))];
        smFlagTemp = [allMonkeyVars(1).smFlagFOV; allMonkeyVars(2).smFlagFOV];
        monkeyName = 'Combined data';
    end

    figure;
    for iArea = 1:3 % Cortical area
        switch iArea 
            case 1 % Somatosensory
                superTemp = squeeze(median(super(smFlagTemp =='S',:,:),1,'omitnan'));
                midTemp   = squeeze(median(mid(smFlagTemp =='S',:,:),1,'omitnan'));
                deepTemp  = squeeze(median(deep(smFlagTemp=='S',:,:),1,'omitnan'));
                typeLabel = 'Sensory';

            case 2 % Motor 
                superTemp = squeeze(median(super(smFlagTemp =='M',:,:),1,'omitnan'));
                midTemp   = squeeze(median(mid(smFlagTemp =='M',:,:),1,'omitnan'));
                deepTemp  = squeeze(median(deep(smFlagTemp=='M',:,:),1,'omitnan'));
                typeLabel = 'Motor';

            case 3 % Both
                superTemp = squeeze(median(super,1,'omitnan'));
                midTemp   = squeeze(median(mid,1,'omitnan'));
                deepTemp  = squeeze(median(deep,1,'omitnan'));
                typeLabel = 'All sites';
        end

        % Superficial
        subplot(3,3,(3*(iArea-1)+1)); imagesc(superTemp); xticks(1:5); yticks(1:5); xticklabels(bandLabels);
        yticklabels(bandLabels); colormap jet; axis image; colorbar; clim([0 1]); box off; title([typeLabel ' - Superficial']);
        
        % Middle
        subplot(3,3,(3*(iArea-1)+2)); imagesc(midTemp);   xticks(1:5); yticks(1:5); xticklabels(bandLabels);
        yticklabels(bandLabels); colormap jet;axis image; colorbar;  clim([0 1]); box off; title([typeLabel ' - Middle']);
        
        % Deep
        subplot(3,3,(3*(iArea-1)+3)); imagesc(deepTemp);  xticks(1:5); yticks(1:5); xticklabels(bandLabels);
        yticklabels(bandLabels); colormap jet; axis image; colorbar; clim([0 1]); box off; title([typeLabel ' - Deep']);
        
        sgtitle(monkeyName);
    end

end

%% Show distributions of pairwise correlations of infraslow ephys powers
clear super mid deep
for iM = 1:3
    if iM~=3 % Single animal
        super      = squeeze(allMonkeyVars(iM).superInfraAllBandsT);
        mid        = squeeze(allMonkeyVars(iM).midInfraAllBandsT);
        deep       = squeeze(allMonkeyVars(iM).deepInfraAllBandsT);
        smFlagTemp = allMonkeyVars(iM).smFlagFOV;
        monkeyName = monkeys{iM};

    else % Both animals 
        super      = [squeeze(allMonkeyVars(1).superInfraAllBandsT); squeeze(allMonkeyVars(2).superInfraAllBandsT)];
        mid        = [squeeze(allMonkeyVars(1).midInfraAllBandsT); squeeze(allMonkeyVars(2).midInfraAllBandsT)];
        deep       = [squeeze(allMonkeyVars(1).deepInfraAllBandsT); squeeze(allMonkeyVars(2).deepInfraAllBandsT)];
        smFlagTemp = [allMonkeyVars(1).smFlagFOV; allMonkeyVars(2).smFlagFOV];
        monkeyName = 'Combined data';
    end

    smMat = logical(diag([1 1 1 1 1 1])); smMat = reshape(smMat,[36 1]);
    dMat  = logical(diag([1 1 1 1 1 1 1 1 1])); dMat = reshape(dMat,[81 1]);



    figure;
    for iType = 1:3 % Cortical area
        clear superTemp midTemp deepTemp
        switch iType
            case 1 % Somatosensory
                superTemp = squeeze(median(super(smFlagTemp =='S',:,:,:,:),[1 4 5],'omitnan'));
                midTemp   = squeeze(median(mid(smFlagTemp =='S',:,:,:,:),[1 4 5],'omitnan'));
                deepTemp  = squeeze(median(deep(smFlagTemp=='S',:,:,:,:),[1 4 5],'omitnan'));
                typeLabel = 'Sensory';
                 
            case 2 % Motor
                superTemp = squeeze(median(super(smFlagTemp =='M',:,:,:,:),[1 4 5],'omitnan'));
                midTemp   = squeeze(median(mid(smFlagTemp =='M',:,:,:,:),[1 4 5],'omitnan'));
                deepTemp  = squeeze(median(deep(smFlagTemp=='M',:,:,:,:),[1 4 5],'omitnan'));
                typeLabel = 'Motor';

            case 3 % Both
                superTemp = squeeze(median(super,[1 4 5],'omitnan'));
                midTemp   = squeeze(median(mid,[1 4 5],'omitnan'));
                deepTemp  = squeeze(median(deep,[1 4 5],'omitnan')); 
                typeLabel = 'All sites';
        end

        
        subplot(3,3,(3*(iType-1)+1)); imagesc(superTemp); xticks(1:5); yticks(1:5); xticklabels(bandLabels);
        yticklabels(bandLabels); colormap jet; axis image; colorbar; clim([0 1]); box off; title([typeLabel ' - Superficial']);
        
        subplot(3,3,(3*(iType-1)+2)); imagesc(midTemp);   xticks(1:5); yticks(1:5); xticklabels(bandLabels);
        yticklabels(bandLabels); colormap jet;axis image; colorbar;  clim([0 1]); box off; title([typeLabel ' - Middle']);
        
        subplot(3,3,(3*(iType-1)+3)); imagesc(deepTemp);  xticks(1:5); yticks(1:5); xticklabels(bandLabels);
        yticklabels(bandLabels); colormap jet; axis image; colorbar; clim([0 1]); box off; title([typeLabel ' - Deep']);
        
        sgtitle(monkeyName);
    end
end


%% Sort the correlations based on date of recording

% Charlie Sheen
dateValsCh = ([dateVals1(1) dateVals1(1) NaN NaN NaN;... % Date x run
    dateVals1(2) dateVals1(2) NaN NaN NaN; ...
    dateVals1(3) dateVals1(3) dateVals1(3) 0 dateVals1(3); ...
    dateVals1(4) dateVals1(4) dateVals1(4) dateVals1(4) NaN]);

dateValsCh = reshape(dateValsCh,[size(dateValsCh,1)*size(dateValsCh,2) 1]);
dateValsCh(isnan(dateValsCh)) = [];
dateValsCh(dateValsCh==0) = [];

% Whiskey
dateValsW = ([dateVals2(1) dateVals2(1) dateVals2(1) NaN NaN NaN NaN; ... % Date x run
    dateVals2(2) dateVals2(2) dateVals2(2) dateVals2(2) dateVals2(2) dateVals2(2) dateVals2(2) ; ...
    dateVals2(3) 0 dateVals2(3) dateVals2(3) NaN NaN NaN; ...
    dateVals2(4) dateVals2(4) 0 dateVals2(4) dateVals2(4) dateVals2(4) dateVals2(4); ...
    dateVals2(5) 0 0 dateVals2(5) dateVals2(5) dateVals2(5) 0]);

dateValsW = reshape(dateValsW,[size(dateValsW,1)*size(dateValsW,2) 1]);
dateValsW(isnan(dateValsW)) = [];
dateValsW(dateValsW==0) = [];

% Both animals
dateValsAll = [dateValsCh; dateValsW];
fovCorr = [allMonkeyVars(1).peakNegValsAllT(:,:,2) ; allMonkeyVars(2).peakNegValsAllT(:,:,2)];     
roiCorr = [allMonkeyVars(1).allChCorr ; allMonkeyVars(2).allChCorr];
smFlag  = ['M';'M';'S';'S';'S';'M';'M';'M';'M';'M';'M';'M';'S';'S';'M';'S';'M';'M';'S';'M';'S';'S';'M'];

smFlagFOV = [allMonkeyVars(1).smFlagFOV; smFlag]; % Cortical area flag
smFlagROI = [allMonkeyVars(1).smFlagROI; allMonkeyVars(2).smFlagROI];

[tAll,loc] = sort(dateValsAll,'ascend'); % Sort in chronological order
fovCorrS = fovCorr(loc,:);
roiCorrS = roiCorr(loc,:);

smFlagFOVS = smFlagFOV(loc); smFlagFOVT = NaN(size(smFlagFOVS));
smFlagROIS = smFlagROI(loc); smFlagROIT = NaN(size(smFlagROIS));

smFlagFOVT(smFlagFOVS=='S') = 1; smFlagFOVT(smFlagFOVS=='M') = 2;
smFlagROIT(smFlagROIS=='S') = 1; smFlagROIT(smFlagROIS=='M') = 2;

comb = cumsum([0; groupcounts(tAll)]);

% Initializing variables
fovGroups   = NaN(length(comb)-1,7,5);
roiGroups   = NaN(length(comb)-1,7,5);
smFOVGroups = NaN(7,length(comb)-1) ;
smROIGroups = NaN(7,length(comb)-1) ;

% Organize the data in chronological order
for iL = 1:length(comb)-1
    fovGroups(iL,1:length(comb(iL)+1:comb(iL+1)),:) = fovCorrS(comb(iL)+1:comb(iL+1),:);
    roiGroups(iL,1:length(comb(iL)+1:comb(iL+1)),:) = roiCorrS(comb(iL)+1:comb(iL+1),:);
    smFOVGroups(1:length(comb(iL)+1:comb(iL+1)),iL) = smFlagFOVT(comb(iL)+1:comb(iL+1));
    smROIGroups(1:length(comb(iL)+1:comb(iL+1)),iL) = smFlagROIT(comb(iL)+1:comb(iL+1));
    timeline(iL,1) = tAll(comb(iL)+1);
end

sensoryCountROI = sum(smROIGroups== 1,1); 
motorCountROI   = sum(smROIGroups== 2,1); 
sensoryCountFOV = sum(smFOVGroups== 1,1); 
motorCountFOV   = sum(smFOVGroups== 2,1); 

% Histogram of data recorded from cortical areas
figure; bar(1:8,[sensoryCountROI;  motorCountROI]'); xticks(1:8); xticklabels(unique(tAll));
ylim([0 5]); xlim([-0 9]); xlabel('Days'); ylabel('Count'); title('ROI');
legend({'Sensory';'Motor'},'Location','northeast','AutoUpdate','off'); box off;

figure; bar(1:8,[sensoryCountFOV;  motorCountFOV]'); xticks(1:8); xticklabels(unique(tAll));
ylim([0 5]); xlim([-0 9]); xlabel('Days'); ylabel('Count'); title('FOV');
legend({'Sensory';'Motor'},'Location','northeast','AutoUpdate','off'); box off;

timelineNew  = [1 30 44 64 90 113 191 240 260 300 400 500 617 722]; 
fovGroupsNew = NaN(length(timelineNew),7,5);
roiGroupsNew = NaN(length(timelineNew),7,5);

for iT = 1: length(timelineNew)
    if any(timeline == timelineNew(iT))
        newLoc = find(timeline == timelineNew(iT));
        fovGroupsNew(iT,:,:) = fovGroups(newLoc,:,:);
        roiGroupsNew(iT,:,:) = roiGroups(newLoc,:,:);
    end
end

timeLineRep = log2(reshape(repmat(timelineNew,[7 1]),[98 1]));

% Show the temporal and spatial correlations over time
% Spatial correlations
figure;
for iBand = 1:5
    clear coeff xFit yFit mdl fVal idx
    subplot(2,3,iBand);
    fVal = (reshape(-fovGroupsNew(:,:,iBand)',[98 1])); 
    plot(timeLineRep,fVal,'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on;

    idx   = isnan(fVal);
    coeff = polyfit(timeLineRep(~idx),fVal(~idx),1);
    xFit  = linspace(min(timeLineRep),max(timeLineRep),1000);
    yFit  = polyval(coeff,xFit); mdl = fitlm(timeLineRep(~idx),fVal(~idx));

    plot(xFit,yFit,'-k','LineWidth',1);

    text(600, 0.2,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
    text(600,0.15,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);

    xlabel('Days'); ylabel('Correlations'); 
    title(bandLabels{iBand}); box off; ylim([-0.3 1]);sgtitle('FOV');
end

% Temporal correlations
figure;
for iBand = 1:5
    clear coeff xFit yFit mdl fVal idx
    subplot(2,3,iBand);
    fVal = reshape(roiGroupsNew(:,:,iBand)',[98 1]); 
    plot(timeLineRep,fVal,'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on;
    
    idx   = isnan(fVal);
    coeff = polyfit(timeLineRep(~idx),fVal(~idx),1);
    xFit  = linspace(min(timeLineRep),max(timeLineRep),1000);
    yFit  = polyval(coeff,xFit); mdl = fitlm(timeLineRep(~idx),fVal(~idx));
  
    plot(xFit,yFit,'-k','LineWidth',1);
   
    text(600, -0.6,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
    text(600,-0.65,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);
    
    xlabel('Days'); ylabel('Correlations'); 
    title(bandLabels{iBand}); box off; ylim([-0.7 0.1]);sgtitle('ROI');
end

% Show distribution of correlations as a boxplot
figure; % Spatial
for iBand = 1:5
    subplot(2,3,iBand);
    boxplot(squeeze(fovGroupsNew(:,:,iBand))',timelineNew);
    xlabel('Days'); ylabel('Correlations');
    title(bandLabels{iBand}); box off; ylim([-1 0.3]); sgtitle('FOV');
end

figure; % Temporal
for iBand = 1:5
    subplot(2,3,iBand);
    boxplot(squeeze(roiGroupsNew(:,:,iBand))',timelineNew);
    xlabel('Days'); ylabel('Correlations');
    title(bandLabels{iBand}); box off; ylim([-0.7 0.1]); sgtitle('ROI');
end
