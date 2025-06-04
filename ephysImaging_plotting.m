% Plot all figures - ephysImaging_plotting.m
% October 9, 2024 - KM
% ISOI vs LFP paper figures

clc; clear;
commonDir = 'C:\Users\kem294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\Imaging']));
addpath(genpath([commonDir '\Codes\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\fly_track\videoIO']));

%% Load all compiled variables
clear allMonkeyVars
hemisphere = 'Left'; spatialBin = 3;
monkeys = {'CharlieSheen';'Whiskey'};
dateVals1 = char(caldiff(datetime({'11/29/2021' ; '01/11/2022'; '08/07/2023'; '11/20/2023'},'InputFormat','MM/dd/yyyy'),'Days')); 
dateVals2 = char(caldiff(datetime({'08/14/2023';'10/16/2023';'12/04/2023';'02/20/2024';'04/29/2024'},'InputFormat','MM/dd/yyyy'),'Days'));

dateVals1 = cumsum([1; str2num(dateVals1(:,1:end-1))]); %#ok<*ST2NM> 
dateVals2 = cumsum([1; str2num(dateVals2(:,1:end-1))]);

% Load monkey variables
for iM = 1:2
    allMonkeyVars(iM) =  load(['D:\Data\' monkeys{iM} '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat']); 
end
bandLabels = {'Theta';'Alpha';'Beta';'Gamma';'Spiking'}; 

%% Sort based on dates
dateValsCh = ([dateVals1(1) dateVals1(1) NaN NaN NaN;... % Date x run
    dateVals1(2) dateVals1(2) NaN NaN NaN; ...
    dateVals1(3) dateVals1(3) dateVals1(3) 0 dateVals1(3); ...
    dateVals1(4) dateVals1(4) dateVals1(4) dateVals1(4) NaN]);

dateValsCh = reshape(dateValsCh,[size(dateValsCh,1)*size(dateValsCh,2) 1]);
dateValsCh(isnan(dateValsCh)) = [];
dateValsCh(dateValsCh==0) = [];

dateValsW = ([dateVals2(1) dateVals2(1) dateVals2(1) NaN NaN NaN NaN; ... % Date x run
    dateVals2(2) dateVals2(2) dateVals2(2) dateVals2(2) dateVals2(2) dateVals2(2) dateVals2(2) ; ...
    dateVals2(3) 0 dateVals2(3) dateVals2(3) NaN NaN NaN; ...
    dateVals2(4) dateVals2(4) 0 dateVals2(4) dateVals2(4) dateVals2(4) dateVals2(4); ...
    dateVals2(5) 0 0 dateVals2(5) dateVals2(5) dateVals2(5) 0]);

dateValsW = reshape(dateValsW,[size(dateValsW,1)*size(dateValsW,2) 1]);
dateValsW(isnan(dateValsW)) = [];
dateValsW(dateValsW==0) = [];

dateValsAll = [dateValsCh; dateValsW];
fovCorr = [allMonkeyVars(1).peakNegValsAllT(:,:,2) ; allMonkeyVars(2).peakNegValsAllT(:,:,2)];     
roiCorr = [allMonkeyVars(1).allChCorr ; allMonkeyVars(2).allChCorr];
smFlag  = ['M';'M';'S';'S';'S';'M';'M';'M';'M';'M';'M';'M';'S';'S';'M';'S';'M';'M';'S';'M';'S';'S';'M'];

smFlagFOV = [allMonkeyVars(1).smFlagFOV; smFlag];
smFlagROI = [allMonkeyVars(1).smFlagROI; allMonkeyVars(2).smFlagROI];

[tAll,loc] = sort(dateValsAll,'ascend');
fovCorrS = fovCorr(loc,:);
roiCorrS = roiCorr(loc,:);

smFlagFOVS = smFlagFOV(loc); smFlagFOVT = NaN(size(smFlagFOVS));
smFlagROIS = smFlagROI(loc); smFlagROIT = NaN(size(smFlagROIS));

smFlagFOVT(smFlagFOVS=='S') = 1; smFlagFOVT(smFlagFOVS=='M') = 2;
smFlagROIT(smFlagROIS=='S') = 1; smFlagROIT(smFlagROIS=='M') = 2;

comb = cumsum([0; groupcounts(tAll)]);

fovGroups = NaN(length(comb)-1,7,5);
roiGroups = NaN(length(comb)-1,7,5);
smFOVGroups = NaN(7,length(comb)-1) ;
smROIGroups = NaN(7,length(comb)-1) ;

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

timeLineRep = reshape(repmat(timelineNew,[7 1]),[98 1]);

figure;
for iBand = 1:5
    clear coeff xFit yFit mdl fVal idx
    subplot(2,3,iBand);
    fVal = reshape(-fovGroupsNew(:,:,iBand)',[98 1]); 
    plot(timeLineRep,fVal,'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on;
    idx = isnan(fVal);
    coeff = polyfit(timeLineRep(~idx),fVal(~idx),1);
    xFit = linspace(min(timeLineRep),max(timeLineRep),1000);
    yFit = polyval(coeff,xFit); mdl = fitlm(timeLineRep(~idx),fVal(~idx));
    plot(xFit,yFit,'-k','LineWidth',1);
    text(600, 0.2,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
    text(600,0.15,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);
    xlabel('Days'); ylabel('Correlations'); 
    title(bandLabels{iBand}); box off; ylim([-0.3 1]);sgtitle('FOV');
end

figure;
for iBand = 1:5
    clear coeff xFit yFit mdl fVal idx
    subplot(2,3,iBand);
    fVal = reshape(roiGroupsNew(:,:,iBand)',[98 1]); 
    plot(timeLineRep,fVal,'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on;
    idx = isnan(fVal);
    coeff = polyfit(timeLineRep(~idx),fVal(~idx),1);
    xFit = linspace(min(timeLineRep),max(timeLineRep),1000);
    yFit = polyval(coeff,xFit); mdl = fitlm(timeLineRep(~idx),fVal(~idx));
    plot(xFit,yFit,'-k','LineWidth',1);
    text(600, -0.6,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
    text(600,-0.65,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);
    xlabel('Days'); ylabel('Correlations'); 
    title(bandLabels{iBand}); box off; ylim([-0.7 0.1]);sgtitle('ROI');
end

figure;
for iBand = 1:5
    subplot(2,3,iBand);
    boxplot(squeeze(fovGroupsNew(:,:,iBand))',timelineNew); 
     xlabel('Days'); ylabel('Correlations');
    title(bandLabels{iBand}); box off; ylim([-1 0.3]); sgtitle('FOV');
end

figure;
for iBand = 1:5
    subplot(2,3,iBand);
    boxplot(squeeze(roiGroupsNew(:,:,iBand))',timelineNew); 
    xlabel('Days'); ylabel('Correlations');
    title(bandLabels{iBand}); box off; ylim([-0.7 0.1]); sgtitle('ROI');
end


%% Correlations vs iso %
isoLevelGoodRunsW = [1;1.2;0.8;0.9;0.7;1.1;1;1;1.25;1;1.75;1.1;1.75;1;1.1;1.1;1;1.3;1.1;1;1.3;1.1;1.3];
isoLevelSpatialW  = [1;1.2;0.8;0.9;0.7;1.1;1;1;1.25;1;1.75;1.1;1.75;1;1.1;1.1;1;1.3;1.1;1;1.3;1.1;1];

isoLevelGoodRunsC = [0.75;0.75;0.9;0.75;0.75;0.75;0.9;0.8;0.9;0.8;0.9;0.9];
isoLevelSpatialC  = [0.75;0.75;0.9;0.75;0.75;0.75;0.9;0.8;0.9;0.8;0.9;0.9];

combinedIsoGoodRuns = [isoLevelGoodRunsC; isoLevelGoodRunsW];
combinedIsoSpatial  = [isoLevelSpatialC; isoLevelSpatialW];

fovCorr = [allMonkeyVars(1).peakNegValsAllT(:,:,2) ; allMonkeyVars(2).peakNegValsAllT(:,:,2)];     
roiCorr = [allMonkeyVars(1).allChCorr ; allMonkeyVars(2).allChCorr];

figure;
for iBand = 1:5
    clear coeff xFit yFit mdl
   subplot(2,3,iBand);
   plot(combinedIsoSpatial,fovCorr(:,iBand),'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on;
   coeff = polyfit(combinedIsoSpatial,fovCorr(:,iBand),1);
   xFit = linspace(min(combinedIsoSpatial),max(combinedIsoSpatial),1000); 
   yFit = polyval(coeff,xFit); mdl = fitlm(combinedIsoSpatial,fovCorr(:,iBand));
   plot(xFit,yFit,'-k','LineWidth',1); 
   text(1.5, 0.2,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
   text(1.5,0.15,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);
   xlabel('iso %'); ylabel('Correlations'); xlim([0.6 1.8]); xticks(0.6:0.1:1.8);
   title(bandLabels{iBand}); box off; ylim([-1 0.3]);sgtitle('FOV');
end

figure;
for iBand = 1:5
    clear coeff xFit yFit mdl
   subplot(2,3,iBand);
   plot(combinedIsoGoodRuns,roiCorr(:,iBand),'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on;
   coeff = polyfit(combinedIsoGoodRuns,roiCorr(:,iBand),1);
   xFit = linspace(min(combinedIsoGoodRuns),max(combinedIsoGoodRuns),1000);
   yFit = polyval(coeff,xFit); mdl = fitlm(combinedIsoGoodRuns,roiCorr(:,iBand));
   plot(xFit,yFit,'-k','LineWidth',1);
   text(1.5, 0,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
   text(1.5,-0.05,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);
   xlabel('iso %'); ylabel('Correlations'); xlim([0.6 1.8]); xticks(0.6:0.1:1.8);
    title(bandLabels{iBand}); box off; ylim([-0.7 0.1]); sgtitle('ROI');
end

%% Show temporal profiles for ROI
clear tempProfiles
for iM = 1:3
    if iM~=3
        tempProfiles = allMonkeyVars(iM).tempProfilesAll;
        monkeyName   = monkeys{iM};
    else
        tempProfiles = [allMonkeyVars(:).tempProfilesAll];
        monkeyName   = 'Combined data'; 
    end
    figure; 
    for iBand = 1:5
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

%% Show correlation distributions for ROI 
clear roiCorr
for iM = 1:3
clear roiCorr
    if iM~=3
        roiCorr = allMonkeyVars(iM).allChCorr;
        monkeyName   = monkeys{iM};
    else
        roiCorr = [allMonkeyVars(1).allChCorr ; allMonkeyVars(2).allChCorr];
        monkeyName   = 'Combined data'; 
    end
    figure; 
    boxplot(roiCorr,bandLabels); ylim([-0.6 0.15]); box off;   
    title(monkeyName);
end
[p,t,s] = anova1(roiCorr,bandLabels);
m = multcompare(s,'Alpha',0.005); % 10 comparisons

%% Show correlation distributions for ROI conditioned on layer and frequency
for iM = 1:3
clear super mid deep
    if iM~=3
        super      = allMonkeyVars(iM).allChCorrSuper;
        mid        = allMonkeyVars(iM).allChCorrMid;
        deep       = allMonkeyVars(iM).allChCorrDeep;
        smFlagTemp = allMonkeyVars(iM).smFlagROI;
        monkeyName = monkeys{iM};
    else
        super      = [allMonkeyVars(1).allChCorrSuper; allMonkeyVars(2).allChCorrSuper];
        mid        = [allMonkeyVars(1).allChCorrMid;   allMonkeyVars(2).allChCorrMid];
        deep       = [allMonkeyVars(1).allChCorrDeep;  allMonkeyVars(2).allChCorrDeep];
        smFlagTemp = [allMonkeyVars(1).smFlagROI;      allMonkeyVars(2).smFlagROI];
        monkeyName   = 'Combined data'; 
    end

    nanRows = find(isnan(mid(:,1)));
    super(nanRows,:) = [];
    mid(nanRows,:) = [];
    deep(nanRows,:) = [];
    smFlagTemp(nanRows) = [];

 figure;
    for iType = 1:3
        switch iType
            case 1
                superTemp = super(smFlagTemp =='S',:);
                midTemp   = mid(smFlagTemp =='S',:);
                deepTemp  = deep(smFlagTemp=='S',:);
                typeLabel = 'Sensory';
            case 2
                superTemp = super(smFlagTemp =='M',:);
                midTemp   = mid(smFlagTemp =='M',:);
                deepTemp  = deep(smFlagTemp=='M',:);
                typeLabel = 'Motor';
            case 3
                superTemp = super;
                midTemp   = mid;
                deepTemp  = deep;
                typeLabel = 'All sites';
        end
        figPlot = [median(superTemp,1,'omitnan');median(midTemp,1,'omitnan' ); median(deepTemp,1,'omitnan')];
        subplot(1,3,iType); imagesc(figPlot); xticks(1:5); xticklabels(bandLabels); yticks(1:3);
        yticklabels({'Superficial';'Middle';'Deep'}); colorbar;colormap(flipud(jet)); 
        caxis([-0.35 0]); title(typeLabel); sgtitle(monkeyName); axis square;
%         figure;
%         subplot(131); boxplot(superTemp,bandLabels); ylim([-0.5 0.5]); box off; title('Superficial');
%         subplot(132); boxplot(midTemp,bandLabels); ylim([-0.5 0.5]); box off; title('Middle');
%         subplot(133); boxplot(deepTemp,bandLabels); ylim([-0.5 0.5]); box off; title('Deep');
%         sgtitle([monkeyName ' - ' typeLabel]);
    end
end

%% Show correlation with FC map for FOV
clear fovProfile
iType = 2;
for iM = 1:3
clear roiCorr
    if iM~=3
        fovProfile = -squeeze(allMonkeyVars(iM).corrFCHybridT(:,:,iType,:));
        monkeyName   = monkeys{iM};
    else
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

%% Show distributions of peak negative correlations with FC map for FOV
clear fovCorr iType
iType = 3;
for iM = 1:3
clear fovCorr
    if iM~=3
        fovCorr = allMonkeyVars(iM).peakNegValsAllT(:,:,iType);
        monkeyName   = monkeys{iM};
    else
        fovCorr = [allMonkeyVars(1).peakNegValsAllT(:,:,iType) ; allMonkeyVars(2).peakNegValsAllT(:,:,iType)];
        monkeyName   = 'Combined data'; 
    end
    figure; 
    boxplot(fovCorr,bandLabels); ylim([-1 1]); box off;   
    title(monkeyName);
end
[p,t,s] = anova1(fovCorr,bandLabels);
m = multcompare(s,'Alpha',0.005); % 10 comparisons

%% Show correlation with FC map for superficial/middle/deep layers
% Boxplots of frequencies for each layer
clear super mid deep
iType =2;
for iM = 1:3
    clear super mid deep
    if iM~=3
        super      = allMonkeyVars(iM).corrFCSuperT(:,:,iType);
        mid        = allMonkeyVars(iM).corrFCMidT(:,:,iType);
        deep       = allMonkeyVars(iM).corrFCDeepT(:,:,iType);
        smFlagTemp = allMonkeyVars(iM).smFlagFOV;
        monkeyName   = monkeys{iM};
    else
        super      = [allMonkeyVars(1).corrFCSuperT(:,:,iType); allMonkeyVars(2).corrFCSuperT(:,:,iType)];
        mid        = [allMonkeyVars(1).corrFCMidT(:,:,iType); allMonkeyVars(2).corrFCMidT(:,:,iType)];
        deep       = [allMonkeyVars(1).corrFCDeepT(:,:,iType); allMonkeyVars(2).corrFCDeepT(:,:,iType)];
        smFlagTemp = [allMonkeyVars(1).smFlagFOV; allMonkeyVars(2).smFlagFOV];
        monkeyName = 'Combined data';
    end

    for iArea = 1:3
        switch iArea
            case 1
                superTemp = super(smFlagTemp =='S',:);
                midTemp   = mid(smFlagTemp =='S',:);
                deepTemp  = deep(smFlagTemp=='S',:);
                typeLabel = 'Sensory';
            case 2
                superTemp = super(smFlagTemp =='M',:);
                midTemp   = mid(smFlagTemp =='M',:);
                deepTemp  = deep(smFlagTemp=='M',:);
                typeLabel = 'Motor';
            case 3
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
    colormap((jet));caxis([0 0.65 ]); colorbar;title('Sensory'); axis image square; 

    subplot(132); imagesc(-motorAll);xticks(1:5); yticks(1:3);
    xticklabels(bandLabels); yticklabels({'Superficial';'Middle';'Deep'});
    colormap((jet));caxis([0 0.65]); colorbar;title('Motor'); axis image square;

    subplot(133); imagesc(-medAll); xticks(1:5); yticks(1:3);
    xticklabels(bandLabels); yticklabels({'Superficial';'Middle';'Deep'});
    colormap((jet));caxis([0 0.65]); colorbar;title('All sites');axis image square;
    sgtitle(monkeyName);

    % Boxplot of correlations with FC for a frequency
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
    subplot(131); boxplot(sensoryOrg,labels); box off; ylim([-1 0.8]);  title('Sensory');
    subplot(132); boxplot(motorOrg,labels); box off; ylim([-1 0.8]);  title('Motor');
    subplot(133); boxplot(allOrg,labels); box off;ylim([-1 0.8]);  title('All areas');


%     for iBand = 1:5
%         subplot(2,3,iBand); boxplot([super(:,iBand) mid(:,iBand) deep(:,iBand)],{'Superficial';'Middle';'Deep'});
%         title(bandLabels{iBand}); box off; ylim([-1 1]);
%     end 
%      sgtitle(monkeyName);
end

%% Comparison between superfical, middle and deep layers
clear superDeep superMid midDeep 
iType = 2;
for iM = 1:3
    if iM~=3
        superDeep  = allMonkeyVars(iM).super_DeepAvgFramesT(:,:,iType);
        superMid   = allMonkeyVars(iM).super_MidAvgFramesT(:,:,iType);
        midDeep    = allMonkeyVars(iM).deep_MidAvgFramesT(:,:,iType);
        controls   = allMonkeyVars(iM). superDeepPosCheckT;
        monkeyName = monkeys{iM};
    else
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

 [pSpCorr,tblSpCorr,statsSpCorr] = anova1([superMid(:,iBand) midDeep(:,iBand) superDeep(:,iBand) controls(:,iBand-1)],{'S/M';'M/D';'S/D';'Controls'});
 [rSpCorr,mSpCorr,~,gnamesSpCorr] = multcompare(statsSpCorr,"Alpha",0.01,"CriticalValueType","bonferroni");


%% Show correlations between hybrids across frequencies for a layer compartment 
clear super mid deep
iType = 2;
for iM = 1:3
clear super mid deep
    if iM~=3
        super      = squeeze(allMonkeyVars(iM).superHybridAllBandsT(:,iType,:,:));
        mid        = squeeze(allMonkeyVars(iM).midHybridAllBandsT(:,iType,:,:));
        deep       = squeeze(allMonkeyVars(iM).deepHybridAllBandsT(:,iType,:,:));
        smFlagTemp = allMonkeyVars(iM).smFlagFOV;
        monkeyName = monkeys{iM};
    else
        super      = [squeeze(allMonkeyVars(1).superHybridAllBandsT(:,iType,:,:)); squeeze(allMonkeyVars(2).superHybridAllBandsT(:,iType,:,:))];
        mid        = [squeeze(allMonkeyVars(1).midHybridAllBandsT(:,iType,:,:)); squeeze(allMonkeyVars(2).midHybridAllBandsT(:,iType,:,:))];
        deep       = [squeeze(allMonkeyVars(1).deepHybridAllBandsT(:,iType,:,:)); squeeze(allMonkeyVars(2).deepHybridAllBandsT(:,iType,:,:))];
        smFlagTemp = [allMonkeyVars(1).smFlagFOV; allMonkeyVars(2).smFlagFOV];
        monkeyName = 'Combined data';
    end

    figure;
    for iArea = 1:3
        switch iArea
            case 1
                superTemp = squeeze(median(super(smFlagTemp =='S',:,:),1,'omitnan'));
                midTemp   = squeeze(median(mid(smFlagTemp =='S',:,:),1,'omitnan'));
                deepTemp  = squeeze(median(deep(smFlagTemp=='S',:,:),1,'omitnan'));
                typeLabel = 'Sensory';
            case 2
                superTemp = squeeze(median(super(smFlagTemp =='M',:,:),1,'omitnan'));
                midTemp   = squeeze(median(mid(smFlagTemp =='M',:,:),1,'omitnan'));
                deepTemp  = squeeze(median(deep(smFlagTemp=='M',:,:),1,'omitnan'));
                typeLabel = 'Motor';
            case 3
                superTemp = squeeze(median(super,1,'omitnan'));
                midTemp   = squeeze(median(mid,1,'omitnan'));
                deepTemp  = squeeze(median(deep,1,'omitnan'));
                typeLabel = 'All sites';
        end

        
        subplot(3,3,(3*(iArea-1)+1)); imagesc(superTemp); xticks(1:5); yticks(1:5); xticklabels(bandLabels);
        yticklabels(bandLabels); colormap jet; axis image; colorbar; caxis([0 1]); box off; title([typeLabel ' - Superficial']);
        subplot(3,3,(3*(iArea-1)+2)); imagesc(midTemp);   xticks(1:5); yticks(1:5); xticklabels(bandLabels);
        yticklabels(bandLabels); colormap jet;axis image; colorbar;  caxis([0 1]); box off; title([typeLabel ' - Middle']);
        subplot(3,3,(3*(iArea-1)+3)); imagesc(deepTemp);  xticks(1:5); yticks(1:5); xticklabels(bandLabels);
        yticklabels(bandLabels); colormap jet; axis image; colorbar; caxis([0 1]); box off; title([typeLabel ' - Deep']);
        sgtitle(monkeyName);
    end

end

%% Show correlations between infraslow ephys across frequencies for a layer compartment 
clear super mid deep
for iM = 1:3
clear super mid deep
    if iM~=3
        super      = squeeze(allMonkeyVars(iM).superInfraAllBandsT);
        mid        = squeeze(allMonkeyVars(iM).midInfraAllBandsT);
        deep       = squeeze(allMonkeyVars(iM).deepInfraAllBandsT);
        smFlagTemp = allMonkeyVars(iM).smFlagFOV;
        monkeyName = monkeys{iM};
    else
        super      = [squeeze(allMonkeyVars(1).superInfraAllBandsT); squeeze(allMonkeyVars(2).superInfraAllBandsT)];
        mid        = [squeeze(allMonkeyVars(1).midInfraAllBandsT); squeeze(allMonkeyVars(2).midInfraAllBandsT)];
        deep       = [squeeze(allMonkeyVars(1).deepInfraAllBandsT); squeeze(allMonkeyVars(2).deepInfraAllBandsT)];
        smFlagTemp = [allMonkeyVars(1).smFlagFOV; allMonkeyVars(2).smFlagFOV];
        monkeyName = 'Combined data';
    end

    smMat = logical(diag([1 1 1 1 1 1])); smMat = reshape(smMat,[36 1]);
    dMat  = logical(diag([1 1 1 1 1 1 1 1 1])); dMat = reshape(dMat,[81 1]);



    figure;
    for iType = 1:3
        clear superTemp midTemp deepTemp
%         superTemp = squeeze(reshape(super,[size(super,1) 5 5 36]));
%         midTemp   = squeeze(reshape(mid,[size(mid,1) 5 5 36]));
%         deepTemp  = squeeze(reshape(deep,[size(deep,1) 5 5 81]));

        switch iType
            case 1
                superTemp = squeeze(median(super(smFlagTemp =='S',:,:,:,:),[1 4 5],'omitnan'));
                midTemp   = squeeze(median(mid(smFlagTemp =='S',:,:,:,:),[1 4 5],'omitnan'));
                deepTemp  = squeeze(median(deep(smFlagTemp=='S',:,:,:,:),[1 4 5],'omitnan'));

%                 superTemp = squeeze(median(superTemp(smFlagTemp =='S',:,:,smMat),[1 4],'omitnan'));
%                 midTemp   = squeeze(median(midTemp(smFlagTemp =='S',:,:,smMat),[1 4],'omitnan'));
%                 deepTemp  = squeeze(median(deepTemp(smFlagTemp =='S',:,:,dMat),[1 4],'omitnan'));  
                  typeLabel = 'Sensory';
                  %
            case 2
                %                 superTemp = squeeze(median(superTemp(smFlagTemp =='M',:,:,smMat),[1 4],'omitnan'));
                %                 midTemp   = squeeze(median(midTemp(smFlagTemp =='M',:,:,smMat),[1 4],'omitnan'));
                %                 deepTemp  = squeeze(median(deepTemp(smFlagTemp =='M',:,:,dMat),[1 4],'omitnan'));
                superTemp = squeeze(median(super(smFlagTemp =='M',:,:,:,:),[1 4 5],'omitnan'));
                midTemp   = squeeze(median(mid(smFlagTemp =='M',:,:,:,:),[1 4 5],'omitnan'));
                deepTemp  = squeeze(median(deep(smFlagTemp=='M',:,:,:,:),[1 4 5],'omitnan'));
                typeLabel = 'Motor';

            case 3
                superTemp = squeeze(median(super,[1 4 5],'omitnan'));
                midTemp   = squeeze(median(mid,[1 4 5],'omitnan'));
                deepTemp  = squeeze(median(deep,[1 4 5],'omitnan')); 

%                 superTemp = squeeze(median(superTemp(:,:,:,smMat),[1 4],'omitnan'));            
%                 midTemp   = squeeze(median(midTemp(:,:,:,smMat),[1 4],'omitnan'));
%                 deepTemp  = squeeze(median(deepTemp(:,:,:,dMat),[1 4],'omitnan'));
                typeLabel = 'All sites';
        end

        
        subplot(3,3,(3*(iType-1)+1)); imagesc(superTemp); xticks(1:5); yticks(1:5); xticklabels(bandLabels);
        yticklabels(bandLabels); colormap jet; axis image; colorbar; caxis([0 1]); box off; title([typeLabel ' - Superficial']);
        subplot(3,3,(3*(iType-1)+2)); imagesc(midTemp);   xticks(1:5); yticks(1:5); xticklabels(bandLabels);
        yticklabels(bandLabels); colormap jet;axis image; colorbar;  caxis([0 1]); box off; title([typeLabel ' - Middle']);
        subplot(3,3,(3*(iType-1)+3)); imagesc(deepTemp);  xticks(1:5); yticks(1:5); xticklabels(bandLabels);
        yticklabels(bandLabels); colormap jet; axis image; colorbar; caxis([0 1]); box off; title([typeLabel ' - Deep']);
        sgtitle(monkeyName);
    end
end

%% Spatial controls, but correlating FC map at reference with FC map at test locations 
% Get Whiskey/Charlie's run-wise data before you run this.
% Codes have to be organized. Big time
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
        serverDir = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\kem294\Data\' monkeyName '_SqM\' ...
            hemisphere ' Hemisphere\' expDate '\' runName];
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

        allSpatialVars = load([dataDir '\spatialControlVarsFOV.mat']);
        locAll                         = allSpatialVars.locAll;
        seedRad =  round(roiSize{iDate}(iRun)/(spatialBin));

         clc; disp(['Analyzing data for ' monkeyName ' '  expDate ' run: ' runName]);

        for iShift = 1:5
            clear loc
            loc = locAll{iShift};
            for iPoint = 1:size(loc,1)
                % Get Gaussian weighted seed signal
                if ~isnan(loc(iPoint,1))
                    seedSigT = calculateSeedSignal(greenFig,clipMaskCortex,loc(iPoint,:),seedRad,pDatTemp);
                    corrMapT = reshape(plotCorrMap(seedSigT,pDatTemp,0),[imSize(1)*imSize(2) 1]);
                    corrMapT(~corrMaskT) = NaN;
                    fc_TestCorr{iRun,iDate}(iShift,iPoint) = corr(corrMapT, fcMap,'rows','complete');
                else
                     fc_TestCorr{iRun,iDate}(iShift,iPoint) = NaN;
                end
            end
        end
    end
end
save(['X:\Data\' monkeyName '_SqM\Left Hemisphere\fcTest_RefCorr.mat'],'fc_TestCorr');

% Grouping the data from the controls
fc_TestCorr = load(['X:\Data\' monkeyName '_SqM\Left Hemisphere\fcTest_RefCorr.mat']);
fc_TestCorr = fc_TestCorr.fc_TestCorr';
fc_TestCorr = reshape(fc_TestCorr,[size(fc_TestCorr,1)*size(fc_TestCorr,2) 1]);
fc_TestCorr(cellfun(@isempty,fc_TestCorr)) = [];
fc_TestCorr(~goodRunsSpatial) = [];

gammaCorr = fovCorr(:,4);
fcTestCorr = cellfun(@(x, y) -(x - y), fc_TestCorr, num2cell(gammaCorr), 'UniformOutput', false);

%%
spCorrControlR = reshape(spCorrControl,[size(spCorrControl,1)*size(spCorrControl,2) 1]);
zeroInd   = cell2mat(cellfun(@(x) isempty(x),spCorrControlR,'un',0));
spCorrControlR(zeroInd) = [];
spCorrControlR(~goodRunsSpatial) = [];

clear spCorrTest
for iRow = 1:size(gammaCorr,1)
    spCorrTest{iRow,1} = (spCorrControlR{iRow,1}-gammaCorr(iRow));
end

spCorrAvg = cellfun(@(x) median(x,2,'omitnan'),spCorrTest,'un',0);
spCorrAvg = (cat(2,spCorrAvg{:}))'; 
figure; boxplot(spCorrAvg,{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'}); hold on;
xlabel('Distance from probe (mm)'); ylabel('Relative correlation between FC and cross-modal map');
ylim([-1 1]); yticks(-1:0.2:1); box off;

figure; violin(spCorrAvg,'facecolor','b','edgecolor','none','bw',0.1);
box off; legend off; ylim([-1 1.2]);yticks(-1:0.2:1); xticklabels([0.5 1 2 3 4]);

spCorrAvgT = cat(2,spCorrTest{:});
figure; boxplot(spCorrAvgT',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'}); hold on; 
ylim([-1 1]);yticks(-1:0.2:1); 
xlabel('Distance from probe (mm)'); ylabel('Correlation between FC map and peak negative map');

figure; violin(spCorrAvgT','facecolor','b','edgecolor','none','bw',0.1);
box off; legend off; ylim([-1 1.2]);yticks(-1:0.2:1); xticklabels([0.5 1 2 3 4]);

[pSpCorrT,tblSpCorrT,statsSpCorrT] = anova1(spCorrAvgT',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'},'off');
[rSpCorrT,mSpCorrT,~,gnamesSpCorrT] = multcompare(statsSpCorrT,"CriticalValueType","bonferroni","Alpha", 0.008);

tblSpCorrMT = array2table(rSpCorrT,"VariableNames",["Group","Control Group","Lower Limit",...
    "Difference","Upper limit","p-val"]);
tblSpCorrMT.("Group") = gnamesSpCorrT(tblSpCorrMT.("Group"));
tblSpCorrMT.("Control Group") = gnamesSpCorrT(tblSpCorrMT.("Control Group"));


%%
clear spCorrMinFC
spCorrMinFC = cellfun(@(x) x./abs(max(x,[],'all','omitnan')),fc_TestCorr,'un',0);
spCorrMinFC = reshape(spCorrMinFC,[size(fc_TestCorr,1)*size(fc_TestCorr,2) 1]);
zeroInd   = cell2mat(cellfun(@(x) isempty(x),spCorrMinFC,'un',0));
spCorrMinFC(zeroInd) = [];
spCorrMinFC(~goodRunsSpatial) = [];
spCorrMinFC = (cat(2,spCorrMinFC{:}));

figure; violin(spCorrMinFC','facecolor','r','edgecolor','none','bw',0.1);
box off; legend off; hold on; 

x1 = (reshape(repmat(1:5,[934 1]),[934*5 1]));
y1 = reshape(spCorrMinFC',[934*5 1]);
s = swarmchart(x1,y1,5,'r','filled');
s.XJitterWidth = 0.5;
ylim([-1 1.3]);yticks(-1:0.1:1.3);box off;
xticks(1:5);xticklabels({'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});

figure; boxplot(spCorrMinFC',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'}); hold on;

% Plot individual penetration data
clear testFCCorr
testFCCorr =  fc_TestCorr{2,6};
figure; violin(testFCCorr','facecolor','r','edgecolor','none','bw',0.1);
box off; legend off; ylim([-1 1.3]); hold on; 

x2 = (reshape(repmat(1:5,[41 1]),[41*5 1]));
y2 = reshape(testFCCorr',[41*5 1]);
s = swarmchart(x2,y2,5,'r','filled');
s.XJitterWidth = 0.5;
ylim([-1 1.3]);yticks(-1:0.1:1.3);box off;
xticks(1:5);xticklabels({'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});


%% Time series shuffling
% Get Whiskey/Charlie's run-wise data before you run this.
% Codes have to be organized. Big time
clear runWiseCorrAllShuffled runWiseLagAllShuffled runWiseCorrShuffled corrNegShuffle corrNegTimes
for iDate = 1:size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);
    for iRun = 1: size(allRuns{iDate,1},1)
        clear runName dataDir clipMask elecMask clipMaskCortex corrMask lowIdx ...
            pDatTemp greenFig seedLocIn crossCorrTemp allHybridMaps mapsAll...
            peakNegHybridMap mapsAllTemp probeCh badTimes szLFP skullMask ...
            inDatSize infraEphys allCortexMask
        x = -200:200;
        negIdx = (-100<=x)&(x<=0); negVals = x(negIdx);

        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
        serverDir = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\kem294\Data\' monkeyName '_SqM\' ...
            hemisphere ' Hemisphere\' expDate '\' runName];
       
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
        if ~exist([dataDir '\phaseControlVarsFOV.mat'],'file')
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
            clear  runWiseCorrShuffled runWiseLagShuffled gammaEphys probeCh...
                processedDat10 ch badChannels badTimes szLFP timeStamp...
                badTimeThreshTemp badTimes

            clc; disp(['Obtaining temporal controls for ' monkeyName ' ' expDate ' ' runName]);

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
            winLen  = [0.001 0.01 0.1 1 5 10 50 100 500 timeLen./1e3].*1e3;

            gammaFFT = fft(gammaEphys);
            magGamma = abs(gammaFFT);
            phaseGamma = angle(gammaFFT);

            for iWin = 1: length(winLen)
                for iRep = 1:10
                    repFlag = 1;
                    disp(['Window length: ' num2str(winLen(iWin)) ' Rep: ' num2str(iRep)]);
                    rng('shuffle');
                    comb1 = randperm(round(timeLen/winLen(iWin)));
                    clear newPhase gammaNew gammaNewFFT magGammaNew

                    if iWin == 1
                        newPhase = phaseGamma(comb1,:);
                        magGammaNew = magGamma((1:size(newPhase,1)),:);
                        gammaNewFFT = magGammaNew.*newPhase;
                        gammaNew = real(ifft(gammaNewFFT));

                    elseif iWin == 10
                        gammaNew = gammaEphys;
                        if iRep>1 
                            repFlag = 0; 
                        end 

                    else
                        newPhase = ones(size(gammaEphys));
                        rowIdx = 1;
                        for iL = 1:length(comb1)
                            clear win1
                            win1 = ((comb1(iL)-1)*winLen(iWin)+1 : (comb1(iL)-1)*winLen(iWin)+winLen(iWin));
                            win1(win1>timeLen) = [];
                            numWin1 = length(win1);
                            newPhase(rowIdx:rowIdx+numWin1-1, :) = phaseGamma(win1, :);
                            rowIdx = rowIdx + numWin1;        
                        end
                        magGammaNew = magGamma((1:size(newPhase,1)),:);
                        gammaNewFFT = magGammaNew.*newPhase;
                        gammaNew = real(ifft(gammaNewFFT));
                    end

                    if repFlag
                        % Get infraslow ephys
                        envelopeDat = envelope(abs(gammaNew),5);

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
                    else
                        runWiseCorrShuffled(iWin,iRep) = runWiseCorrShuffled(iWin,1);
                        runWiseLagShuffled(iWin,iRep) = runWiseLagShuffled(iWin,1);
                    end
                end
            end

            toc;

            save([dataDir '\phaseControlVarsFOV.mat'],'runWiseCorrShuffled','runWiseLagShuffled');

            runWiseCorrAllShuffled{iDate,iRun} = runWiseCorrShuffled;
            runWiseLagAllShuffled{iDate,iRun}  = runWiseLagShuffled;
            corrNegShuffle(iDate,iRun,:)       = median(runWiseCorrShuffled,2,'omitnan');%runWiseCorrShuffled;%
            corrNegTimes(iDate,iRun,:)         = median(runWiseLagShuffled,2,'omitnan');%runWiseLagShuffled;%
       
        else
            clear allVars runWiseCorrShuffled
            allVars = load([dataDir '\phaseControlVarsFOV.mat']);
            runWiseCorrAllShuffled{iDate,iRun} = allVars.runWiseCorrShuffled;
            runWiseLagAllShuffled{iDate,iRun}  = allVars.runWiseLagShuffled;
            runWiseCorrShuffled                = allVars.runWiseCorrShuffled;
            corrNegShuffle(iDate,iRun,:)       = median(runWiseCorrShuffled,2,'omitnan');
            corrNegTimes(iDate,iRun,:)         = median(allVars.runWiseLagShuffled,2,'omitnan');
        end
    end
end

%% Grouping temporal controls 
winLen  = [0.001 0.01 0.1 1 5 10 50 100 500 900];
runWiseCorrAllShuffledT = reshape(runWiseCorrAllShuffled,[size(runWiseCorrAllShuffled,1)*size(runWiseCorrAllShuffled,2) 1]);
zeroInd   = cell2mat(cellfun(@(x) isempty(x),runWiseCorrAllShuffledT,'un',0));
runWiseCorrAllShuffledT(zeroInd) = [];
runWiseCorrAllShuffledT(~goodRunsSpatial) = [];
runWiseCorrAllShuffledT = cellfun(@(x) x./(min(x(10,:))),runWiseCorrAllShuffledT,'un',0)';
allTimePoints = (cat(2,runWiseCorrAllShuffledT{:}))';
runWiseCorrAllShuffledT = cell2mat(cellfun(@(x) median(x,2,'omitnan'),runWiseCorrAllShuffledT,'un',0));

% figure; plot(smoothdata(runWiseCorrAllShuffledT,2,'movmean',10),'Color',[0.65 0.65 0.65],'LineWidth',1); hold on;
figure; plot(smoothdata(runWiseCorrAllShuffledT,2,'movmean',10),'Color',[0.65 0.65 0.65],'LineWidth',1); hold on;
plot(median(smoothdata(runWiseCorrAllShuffledT,2,'movmean',2),2),'k','LineWidth',2);
xticklabels(winLen);  hold on;  ylim([-0.5 1]); yticks(-1:0.1:1);box off; ylim([-0.5 1]);
xlabel('Length of window for shuffling (s)'); ylabel('Cross correlation'); grid off;


stdAll = std(runWiseCorrAllShuffledT,[],2)/sqrt(size(runWiseCorrAllShuffledT,2));
c95 = tinv([0.025 0.975],size(runWiseCorrAllShuffledT,2)-1);
y95 = bsxfun(@times,stdAll', c95(:));

meanAll = mean(runWiseCorrAllShuffledT,2,'omitnan');
posValAll = meanAll+y95';

figure; plot(median(smoothdata(runWiseCorrAllShuffledT,2,'movmean',2),2),'k','LineWidth',2);
hold on
%  xVar = [(1:length(winLen)) fliplr((1:length(winLen)))];
%  yVar = [(squeeze(corrNegShuffle(iDate,iRun,:))-2.*semAll)' ...
%  flipud((squeeze(corrNegShuffle(iDate,iRun,:))+2.*semAll))'];
%  patch(xVar,yVar,'blue','FaceAlpha',0.3,'EdgeColor','none');
xVar = [1:size(posValAll,1) fliplr((1:size(posValAll,1)))];
patch(xVar,[posValAll(:,1)' fliplr(posValAll(:,2)')],[0.65 0.65 0.65],'FaceAlpha',0.3,'EdgeColor','none')
xticklabels(winLen);  hold on;  ylim([-0.5 1]); yticks(-1:0.1:1);box off; ylim([-0.5 1]);

% figure; plot(smoothdata(allTimePoints',1,'movmean',1),'Color',[0.65 0.65 0.65],'LineWidth',1); hold on;
figure;plot(median(smoothdata(allTimePoints,1,'movmean',2),1),'k','LineWidth',2);hold on; 
meanAll = mean(allTimePoints,1,'omitnan');
stdAll = std(allTimePoints,[],1)/sqrt(size(allTimePoints,1));
c95 = tinv([0.025 0.975],size(allTimePoints,1)-1);
y95 = bsxfun(@times,stdAll, c95(:));

posValAll = (meanAll+y95)';
xVar = [1:size(posValAll,1) fliplr((1:size(posValAll,1)))];
patch(xVar,[posValAll(:,1)' fliplr(posValAll(:,2)')],[0.65 0.65 0.65],'FaceAlpha',0.3,'EdgeColor','none')
xticklabels(winLen);  hold on;  ylim([-0.5 1]); yticks(-1:0.1:1);box off; ylim([-0.5 1]);
xlabel('Length of window for shuffling (s)'); ylabel('Cross correlation'); grid off;

%% Plotting data from both monkeys
figure; plot(median(smoothdata(runWiseShuffledAll,2,'movmean',2),2),'k','LineWidth',2); hold on
stdErr = std(runWiseShuffledAll,[],2)/sqrt(size(runWiseShuffledAll,2)); 
xVar = [1:size(stdErr,1) fliplr((1:size(stdErr,1)))];
meanAll = mean(runWiseShuffledAll,2,'omitnan');
posVal(:,1) = meanAll-2.*stdErr;
posVal(:,2) = meanAll+2.*stdErr;
patch(xVar,[posVal(:,1)' fliplr(posValAll(:,2)')],[0.65 0.65 0.65],'FaceAlpha',0.3,'EdgeColor','none')
box off; ylim([-0.1 1]); yticks(-0.1:0.1:1);
xticklabels(winLen(2:10)); 
xlabel('Length of window for shuffling (s)'); ylabel('Cross correlation'); grid off;

%% Combining data from both monkeys
clear runWiseShuffledAll
runWiseShuffledAll = [runWiseCorrAllShuffledTWhiskey(2:10,:) runWiseCorrAllShuffledT(2:10,:)];
stdAll = std(runWiseShuffledAll,[],2)/sqrt(size(runWiseShuffledAll,2)); 
c95 = tinv([0.025 0.975],size(runWiseShuffledAll,2)-1);
y95 = bsxfun(@times,stdAll', c95(:));

meanAll = mean(runWiseShuffledAll,2,'omitnan');
posValAll = meanAll+y95';

figure; plot(median(smoothdata(runWiseShuffledAll,2,'movmean',2),2),'k','LineWidth',2); hold on
xVar = [1:size(posValAll,1) fliplr((1:size(posValAll,1)))];
patch(xVar,[posValAll(:,1)' fliplr(posValAll(:,2)')],[0.65 0.65 0.65],'FaceAlpha',0.3,'EdgeColor','none')
xticklabels(winLen(2:10));  hold on;  ylim([-0.5 1]); yticks(-1:0.1:1);box off; ylim([-0.5 1]);
xlabel('Length of window for shuffling (s)'); ylabel('Cross correlation'); grid off;

% Variance from the shuffles
clear allTimePointsAll
allTimePointsAll = [allTimePointsWhiskey(:,2:10); allTimePoints(:,2:10)];
figure;plot(median(smoothdata(allTimePointsAll,1,'movmean',2),1),'k','LineWidth',2);hold on; 
meanAll = mean(allTimePointsAll,1,'omitnan');
stdAll = std(allTimePointsAll,[],1)/sqrt(size(allTimePointsAll,1));
c95 = tinv([0.025 0.975],size(allTimePointsAll,1)-1);
y95 = bsxfun(@times,stdAll, c95(:));

posValAll = (meanAll+y95)';
xVar = [1:size(posValAll,1) fliplr((1:size(posValAll,1)))];
patch(xVar,[posValAll(:,1)' fliplr(posValAll(:,2)')],[0.65 0.65 0.65],'FaceAlpha',0.3,'EdgeColor','none')
xticklabels(winLen(2:10));  hold on;  ylim([-0.5 1]); yticks(-1:0.1:1);box off; ylim([-0.5 1]);
xlabel('Length of window for shuffling (s)'); ylabel('Cross correlation'); grid off;


figure; plot(median(smoothdata(runWiseShuffledAll,2,'movmean',2),2),'k','LineWidth',2);hold on;
plot(median(smoothdata(allTimePointsAll,1,'movmean',2),1),'r','LineWidth',2);hold on;
legend('Variance from sessions','Variance from shuffles');box off;
xticklabels(winLen(2:10));  hold on;  ylim([-0.5 1]); yticks(-1:0.1:1);box off; ylim([-0.5 1]);

%% Plotting time courses for verification
iDate = 2; iRun = 6;
probeCh                = probe{iRun,iDate}.probeCh;
ch                     = estChInCortex{iDate}(iRun,:);
badChannels            = badCh{iDate,iRun};
badTimes               = badTimesLFP{iDate,iRun};
probeCh(:,badChannels) = [];
szLFP                  = size(probeCh,1);

probeCh        = probeCh(1:szLFP,:);
probeCh(badTimes,:) = [];

gammaBand  = [30 90]; [bG,aG] = butter(3,gammaBand./(1e3/2),'bandpass');% Gamma band filtering parameters
gammaEphys = single(filtfilt(bG,aG,double(probeCh(:,ch(1):ch(2)))));

timeLen = size(gammaEphys,1);
winLen  = [0.001 0.01 0.1 1 5 10 50 100 500 timeLen./1e3].*1e3;

gammaFFT = fft(gammaEphys);
magGamma = abs(gammaFFT);
phaseGamma = angle(gammaFFT);

%%%%% 1 ms shuffle %%%%%
clear comb1 newPhase magGammaNew gammaNewFFT gammaNew
rng('shuffle');
comb1 = randperm(round(timeLen));
newPhase = phaseGamma(comb1,:);
magGammaNew = magGamma((1:size(newPhase,1)),:);
gammaNewFFT = magGammaNew.*newPhase;
gammaNew = real(ifft(gammaNewFFT));

subplot(511);
plot(gammaNew(:,10)); title('1 ms shuffle');box off;
xlim([30e3 31e3]); ylim([-15 15]); yticks(-15:5:15);

%%%%% 10 ms shuffle %%%%
clear comb1 newPhase magGammaNew gammaNewFFT gammaNew
rng('shuffle');
comb1 = randperm(round(timeLen/10));

newPhase = ones(size(gammaEphys));
     rowIdx = 1;
     for iL = 1:length(comb1)
         clear win1
         win1 = ((comb1(iL)-1)*winLen(2)+1 : (comb1(iL)-1)*winLen(2)+winLen(2));
         win1(win1>timeLen) = [];
         numWin1 = length(win1);
         newPhase(rowIdx:rowIdx+numWin1-1, :) = phaseGamma(win1, :);
         rowIdx = rowIdx + numWin1;
     end
     magGammaNew = magGamma((1:size(newPhase,1)),:);
     gammaNewFFT = magGammaNew.*newPhase;
     gammaNew = real(ifft(gammaNewFFT));

subplot(512);
plot(gammaNew(:,10)); title('10 ms shuffle'); box off;
xlim([30e3 31e3]); ylim([-15 15]); yticks(-15:5:15);

%%%%% 5 s shuffle %%%%%

clear comb1 newPhase magGammaNew gammaNewFFT gammaNew
rng('shuffle');
comb1 = randperm(round(timeLen/5000));

newPhase = ones(size(gammaEphys));
     rowIdx = 1;
     for iL = 1:length(comb1)
         clear win1
         win1 = ((comb1(iL)-1)*winLen(5)+1 : (comb1(iL)-1)*winLen(5)+winLen(5));
         win1(win1>timeLen) = [];
         numWin1 = length(win1);
         newPhase(rowIdx:rowIdx+numWin1-1, :) = phaseGamma(win1, :);
         rowIdx = rowIdx + numWin1;
     end
     magGammaNew = magGamma((1:size(newPhase,1)),:);
     gammaNewFFT = magGammaNew.*newPhase;
     gammaNew = real(ifft(gammaNewFFT));

subplot(513);
plot(gammaNew(:,10)); title('5 s shuffle');box off;
xlim([30e3 31e3]); ylim([-15 15]); yticks(-15:5:15);


%%%%% 500 s shuffle %%%%%
clear comb1 newPhase magGammaNew gammaNewFFT gammaNew
rng('shuffle');
comb1 = randperm(round(timeLen/winLen(9)));
newPhase = ones(size(gammaEphys));
     rowIdx = 1;
 
     for iL = 1:length(comb1)
         clear win1
         win1 = ((comb1(iL)-1)*winLen(9)+1 : (comb1(iL)-1)*winLen(9)+winLen(9));
         win1(win1>timeLen) = [];
         numWin1 = length(win1);
         newPhase(rowIdx:rowIdx+numWin1-1, :) = phaseGamma(win1, :);
         rowIdx = rowIdx + numWin1;
     end
     magGammaNew = magGamma((1:size(newPhase,1)),:);
     gammaNewFFT = magGammaNew.*newPhase;
     gammaNew = real(ifft(gammaNewFFT));

subplot(514);
plot(gammaNew(:,10)); title('500 s shuffle');box off;
xlim([30e3 31e3]); ylim([-15 15]); yticks(-15:5:15);

%%%%% Unshuffled %%%%%
 subplot(515);
 plot(gammaEphys(:,10)); title('unshuffled');box off;
 xlim([30e3 31e3]); ylim([-15 15]); yticks(-15:5:15);