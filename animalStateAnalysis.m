% Developing a parameter to measure animal state for ISOI_Ephys recordings
% November 13,2024 - KM
clc; clear;
commonDir = 'C:\Users\KEM294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\Imaging']));
addpath(genpath([commonDir '\Codes\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\fly_track\videoIO']));
rmpath('C:\Users\KEM294\Documents\Data\Codes\chronux_2_12\spectral_analysis\continuous\dupes')

%% Load ISOI_Ephys data
clear allMonkeyVars
hemisphere = 'Left'; 
monkeys = {'CharlieSheen';'Whiskey'};
dateVals1 = char(caldiff(datetime({'11/29/2021' ; '01/11/2022'; '08/07/2023'; '11/20/2023'},'InputFormat','MM/dd/yyyy'),'Days')); 
dateVals2 = char(caldiff(datetime({'08/14/2023';'10/16/2023';'12/04/2023';'02/20/2024';'04/29/2024'},'InputFormat','MM/dd/yyyy'),'Days'));

dateVals1 = cumsum([1; str2num(dateVals1(:,1:end-1))]); %#ok<*ST2NM> 
dateVals2 = cumsum([1; str2num(dateVals2(:,1:end-1))]);

% Load monkey variables
for iM = 1:2
    allMonkeyVars(iM) =  load(['X:\Data\' monkeys{iM} '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat']); %#ok<SAGROW> 
end
bandLabels = {'Theta';'Alpha';'Beta';'Gamma';'Spiking'}; 

%% Load the iso levels for each run 
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

%% Look at power spectra for EEG for Whiskey - 10/17/22 - runs 5-8
expDate = '10_17_2022';
fileNum = 5:8;
monkeyName  = 'Whiskey'; 
hemisphere = 'Left';
saveFolder = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Processed Data'];

gammaBand   = [30 90]; [bG,aG] = butter(3,gammaBand./(1e3/2),'bandpass'); % Gamma band filtering parameters
params.Fs       = 1e3; % PSD parameter initiation
params.fpass    = [2 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;
paramsNew = params;
% paramsNew.fpass = [2 120];

eegPSD = NaN(120,length(fileNum));
animalStateIso = [1.2 2.25 3.25 1];

for iFile = 1:length(fileNum)
    clear eegCh eegPSDT
    var =  load([saveFolder '\datafile000'  num2str(fileNum(iFile)) '_lfp.mat']);
    eegCh = var.eegCh;

    clear spec powMeansS badTimeThresh badTimeIndOld badTimeInd
    [spec,timeValsSpec,~] = mtspecgramc(eegCh,[5 3],params); 
    
    powMeanS = squeeze(sum(spec,2));
    badTimeThresh   = (median(powMeanS,1)+4.5*mad(powMeanS,1,1));
    badTimeIndOld = floor(timeValsSpec(powMeanS>badTimeThresh))*1e3;
    badTimes = [];
    if ~isempty(badTimeIndOld)
        badTimes = [];
        badTimeInd =[(badTimeIndOld-1e3)'  (badTimeIndOld+1e3)']; % Taking one second before and one second after bad time segments

        for iL = 1:size(badTimeInd,1)
            badTimes = [badTimes badTimeInd(iL,1): badTimeInd(iL,2)];
        end
         badTimes = unique(badTimes);
    end

    eegCh(badTimes,:) = [];
%     [spec, freqValsSpec,t] = spectrogram(eegCh,128,120,0:120,1000,'yaxis');
%     cumPowFreq = cumsum(abs(spec),1);
%     mat = cumPowFreq>= 0.9.*cumPowFreq(end,:);
%     [~,highPowerIdx] = max(mat,[],1);

    [spec,timeValsSpec,freqValsSpec] = mtspecgramc(eegCh,[2 1],paramsNew);
    cumPowFreq = cumsum(spec,2);
    mat = cumPowFreq>= 0.9.*cumPowFreq(:,end);
    [~,highPowerIdx] = max(mat,[],2);

%     figure; imagesc(timeValsSpec,freqValsSpec,20.*log10(spec'));
%     set(gca,'YDir','normal'); colorbar; colormap jet; caxis([0 100]);
%     title(['Anesthesia: ' num2str(animalStateIso(iFile))]); hold on;
%     plot(timeValsSpec,freqValsSpec(highPowerIdx),'w','LineWidth',2);    

    meanFreqW(iFile) = median(freqValsSpec(highPowerIdx)); 
    % Get power for frequencies from 17-27 Hz
    idx = freqValsSpec>=20 & freqValsSpec<=60;

    animalStateCheckW(iFile) = median(20.*log10(spec(:,idx)),'all','omitnan');
end

figure; scatter(animalStateIso,meanFreqW,45,[1 2 3 1],'filled'); xlabel('Iso %'); ylabel('Edge frequency');
figure; scatter(animalStateIso,animalStateCheckW,'filled');

%% Get data from Bordeaux and Rambo
fs        = 1e3;
gammaBand = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); % Gamma band filtering parameters

params.Fs       = 1e3; % PSD parameter initiation
params.fpass    = [2 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;

for iM = 1:2
    clear monkeyName allDates allRuns rippleName meanFreq isoLevels eegPSD specAll animalStateCheck powEnvelope
    switch iM
        case 1
            monkeyName = '15-18_Bordeaux';
            allDates   = ['01_06_2020'; '12_03_2019';'07_09_2019'];
            allRuns    = {['run01';'run02';'run03';'run04'];'run00';['run00';'run01';'run02']};
            rippleName  = {[];[];['0001';'0002';'0003']};
            isoLevels  = [1.1 2.25 1.25 1.75; 0.6 NaN NaN NaN; 1 2 1 NaN];
        case 2
            monkeyName = '16-18_Rambo';
            allDates   = ['07_15_2019'; '11_19_2019'];
            allRuns    = {['run00'; 'run01';'run02';'run03'];['run00'; 'run01';'run02';'run03']};
            rippleName  = {['0001';'0002';'0003';'0004'];['0001';'0002';'0003';'0004']};
            isoLevels  = [0.9 2 2 1; 1.5 2 1.4 1.9];
    end
    serverDataPath = '\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\Euthanized Animal Data\' ;
    
    for iDate = 1:size(allDates,1)
        expDate  = allDates(iDate,:);
        clear eegChT eegCh

        for iRun = 1:size(allRuns{iDate},1)
            runName = allRuns{iDate}(iRun,:);

            if strcmp(expDate,'01_06_2020')
                dat = load([serverDataPath monkeyName '_SqM\Left Hemisphere\' expDate '\' runName ...
                    '\Bordeaux_01062020_' runName '_EEG_EKG_data.mat']); 
                eegChT = dat.storedData(:,1).*1e3;        

            elseif strcmp(expDate,'12_03_2019')
                dat = load([serverDataPath monkeyName '_SqM\Left Hemisphere\' expDate '\' runName ...
                    '\Bordeaux_1232019_' runName '_EEG_EKG_data.mat']);
                eegChT = dat.storedData(:,1).*1e3;
            else
                 datName = ([serverDataPath monkeyName '_SqM\Left Hemisphere\' expDate  ...
                    '\RippleFiles\' rippleName{iDate}(iRun,:)]);
                      
                 [nsResult,hFile] = ns_OpenFile(datName);
                 [nsResult2, fileInfoVar] = ns_GetFileInfo(hFile);
                
                 clear entityInfo
                 for iEntity = 1: fileInfoVar.EntityCount
                     [~,entityInfo(iEntity,1)] = ns_GetEntityInfo(hFile,iEntity);
                 end

                 clear elecList elecID
                 elecList  = find([entityInfo.EntityType] == 2);
                 elecLabel = {entityInfo(elecList).EntityLabel}; % Gets the label of all the channels in lfpList
                 elecIdx   = cell2mat(cellfun(@(x) strcmp(x(1:end-2),'analog'),elecLabel,'un',0));
                 elecList(~elecIdx)  = [];
                 elecLabel(~elecIdx) = [];
                 
                 clear listVal
                 if length(elecList)==2
                     listVal = 2;
                 else
                     listVal = 1;
                 end
       
                 [~, analogInfo2] = ns_GetAnalogInfo(hFile, elecList(listVal));
                 elecCount        = entityInfo(elecList(listVal)).ItemCount;
                 fs_ns5           = analogInfo2.SampleRate;
                 [~, ~, eegChT]    = ns_GetAnalogData(hFile,elecList(listVal),1,elecCount);
                 if fs_ns5 == 30e3
                     eegChT = downsample(eegChT,fs_ns5/1e3);                
                 end
                 
            end

            clear spec powMeansS badTimeThresh badTimeIndOld badTimeInd
            [spec,timeValsSpec,freqValsSpec] = mtspecgramc(eegChT,[5 2],params);

            powMeanS = squeeze(sum(spec,2));
            badTimeThresh   = (median(powMeanS,1)+4.5*mad(powMeanS,1,1));
            badTimeIndOld = floor(timeValsSpec(powMeanS>badTimeThresh))*1e3;
            badTimes = [];
            
            if ~isempty(badTimeIndOld)
                badTimes = [];
                badTimeInd =[(badTimeIndOld-1e3)'  (badTimeIndOld+1e3)']; % Taking one second before and one second after bad time segments

                for iL = 1:size(badTimeInd,1)
                    badTimes = [badTimes badTimeInd(iL,1): badTimeInd(iL,2)];
                end
                badTimes = unique(badTimes);
            end

            eegChT(badTimes,:) = [];
            eegCh(:,iRun) = eegChT(1:500e3);
        end           

        eegChNorm = eegCh;

        for iRun = 1:size(allRuns{iDate},1)
            clear eegChT 
            eegChT = (eegChNorm(:,iRun));
            [eegPSDT,psdFreq]  = mtspectrumc(buffer(eegChT,1000),params);
            eegPSD(:,iDate,iRun)= median(eegPSDT,2,'omitnan');
            
            paramsNew = params;%paramsNew.fpass = [2 120];
            [spec,timeValsSpec,freqValsSpec] = mtspecgramc(eegChT,[2 1],paramsNew);
            cumPowFreq = cumsum(spec,2);
            mat = cumPowFreq>= 0.9.*cumPowFreq(:,end);
            [~,highPowerIdx] = max(mat,[],2);

%             [s,freqValsSpec,~,~] = spectrogram(eegChT,128,120,256,1000,'yaxis');           
%             cumPowFreq = cumsum(abs(s'),2);
%             figure; imagesc(timeValsSpec,freqValsSpec,20.*log10(spec'));
%             set(gca,'YDir','normal'); colorbar; colormap jet; caxis([0 100]);
% figure; histogram(freqValsSpec(highPowerIdx));
% title(strrep([monkeyName ': ' expDate ' run:' num2str(iRun)],'_','\_')); hold on;
%             plot(timeValsSpec,freqValsSpec(highPowerIdx),'w','LineWidth',2);

            meanFreq(iDate,iRun) = median(freqValsSpec(highPowerIdx));
            specAll(:,:,iDate,iRun)= spec; %specAll(iDate,iRun).timeVals = timeValsSpec;
% 
            animalStateCheck(iDate,iRun) = median(20.*log10(eegPSD(20:60,iDate,iRun)),'all','omitnan');        
        end

    end

    switch iM
        case 1
            eegB  = eegPSD;
            specB = specAll; 
            isoB  = [1.1 2.25 1.25 1.75; 0.6 NaN NaN NaN; 1 2 1 NaN];%reshape([1.1 2.25 1.25 1.75; 0.6 NaN NaN NaN; 1 2 1 NaN],[size(eegPSD,2)*size(eegPSD,3) 1]);
            animalB = animalStateCheck; 
            meanFreqB = meanFreq;
        case 2
            eegR  = eegPSD; 
            specR = specAll; 
            isoR  = [0.9 2 2 1; 1.5 2 1.4 1.9];
            animalR = animalStateCheck;
            meanFreqR = meanFreq; 
    end
end 


bordeauxLight = [squeeze(eegB(:,1,[1 3])) squeeze(eegB(:,2,1)) squeeze(eegB(:,3,[1 3]))];
bordeauxSemi  = [squeeze(eegB(:,1,4)) squeeze(eegB(:,3,2))];
bordeauxDeep  = [squeeze(eegB(:,1,2))];

medBordeauxLight = 20.*log10(median(bordeauxLight(:,4:5),2,'omitnan')); 
medBordeauxSemi  = 20.*log10(bordeauxSemi(:,2));
d = medBordeauxLight - medBordeauxSemi; 

medL = 20.*log10(median(bordeauxLight(:,1:3),2,'omitnan'));
medD  = 20.*log10(bordeauxDeep);
dB = medL - medD;

bordeauxLightSpec = cat(3,squeeze(specB(:,:,1,[1 3])),squeeze(specB(:,:,2,1)),squeeze(specB(:,:,3,[1 3])));
bordeauxSemiSpec  = cat(3, squeeze(specB(:,:,1,4)),squeeze(specB(:,:,3,2)));
bordeauxDeepSpec  = squeeze(specB(:,:,1,2));

medBordeauxLightSpec = 20.*log10(median(bordeauxLightSpec(:,:,1:3),3,'omitnan'));
medBordeauxSemiSpec  = 20.*log10(bordeauxSemiSpec(:,:,1));
medBordeauxDeepSpec  = 20.*log10(bordeauxDeepSpec); 

lightPowersB = median(medBordeauxLightSpec,1,'omitnan');
semiPowersB  = median(medBordeauxSemiSpec,1,'omitnan');
deepPowersB  = median(medBordeauxDeepSpec,1,'omitnan'); 

lightDeepB = lightPowersB-deepPowersB; 

ramboLight = [squeeze(eegR(:,1,[1 4])) squeeze(eegR(:,2,[1 3]))] ; 
ramboSemi  = [squeeze(eegR(:,1,2)) squeeze(eegR(:,2,4))] ; 
ramboDeep  = [squeeze(eegR(:,1,3)) squeeze(eegR(:,2,2))] ; 

medRamboLight = 20.*log10(median(ramboLight,2,'omitnan'));
medRamboDeep  = 20.*log10(median(ramboDeep,2,'omitnan'));
dR = medRamboLight - medRamboDeep; 

ramboLightSpec = cat(3,squeeze(specR(:,:,1,[1 4])),squeeze(specR(:,:,2,[1 3]))) ; 
ramboSemiSpec  = cat(3,squeeze(specR(:,:,1,2)),squeeze(specR(:,:,2,4))) ; 
ramboDeepSpec  = cat(3,squeeze(specR(:,:,1,3)),squeeze(specR(:,:,2,2))) ; 

medRamboLightSpec = 20.*log10(median(ramboLightSpec,3,'omitnan'));
medRamboSemiSpec  = 20.*log10(median(ramboSemiSpec,3,'omitnan'));
medRamboDeepSpec  = 20.*log10(median(ramboDeepSpec,3,'omitnan')); 

lightPowersR = median(medRamboLightSpec,1,'omitnan');
semiPowersR  = median(medRamboSemiSpec,1,'omitnan');
deepPowersR  = median(medRamboDeepSpec,1,'omitnan'); 

lightDeepR = lightPowersR-deepPowersR; 

colorR = [1 2 3 1; 1 3 1 2];
colorB = [1 3 1 2; 1 NaN NaN NaN; 1 2 1 NaN];

animalB = reshape(animalB,[12 1]);
animalR = reshape(animalR,[8 1]);

colorB = reshape(colorB,[12 1]);
colorR = reshape(colorR,[8 1]);

isoR    = reshape(isoR,[8 1]);
isoB    = reshape(isoB,[12 1]);
nanVals = isnan(isoB);

meanFreqB = reshape(meanFreqB,[12 1]);
meanFreqR = reshape(meanFreqR,[8 1]); 

freqIdx = freqValsSpec>=20 & freqValsSpec<=60;
specBPow = reshape(20.*log10(squeeze(median(specB(:,freqIdx,:,:),[1 2],'omitnan'))),[12 1]);
specRPow = reshape(20.*log10(squeeze(median(specR(:,freqIdx,:,:),[1 2],'omitnan'))),[8 1]);

isoB(nanVals)         = [];
animalB(nanVals)      = [];
colorB(nanVals)       = [];
specBPow(nanVals)     = [];
meanFreqB(nanVals)    = [];

figure; scatter(isoB,animalB,35,colorB,'filled');ylim([-20 70]);
figure; scatter(isoR,animalR,35,colorR,'filled');ylim([-20 70]);
figure; scatter(isoB,specBPow,35,colorB,'filled'); ylim([-20 70]);
figure; scatter(isoR,specRPow,35,colorR,'filled');ylim([-20 70]);

figure; scatter(isoB,meanFreqB,35,colorB,'filled'); %ylim([-20 70]);
figure; scatter(isoR,meanFreqR,35,colorR,'filled');%ylim([-20 70]);

%% Combining all three monkeys
isoAllMonkeys   = [isoB ;isoR; animalStateIso'];
colorAllMonkeys = [colorB; colorR; [1;2;3;1]];
specPowAll      = [specBPow; specRPow; animalStateCheckW'];
edgeFreqAll     = [meanFreqB; meanFreqR; meanFreqW']; 

figure; scatter(isoAllMonkeys,specPowAll,60,colorAllMonkeys,'filled');
ylim([-25 45]);ylabel('Power (dB)'); xlabel('Iso (%)');xlim([0.5 3.5]);

figure; scatter(isoAllMonkeys,edgeFreqAll,60,colorAllMonkeys,'filled');
ylim([0 20]);ylabel('Edge Frequency (Hz)'); xlabel('Iso (%)');xlim([0.5 3.5]);

% K means clustering to see if the edge frequencies can be separable
[idx,c] = kmeans(edgeFreqAll,3);
figure; plot(isoAllMonkeys(idx==1),edgeFreqAll(idx==1),'r.','MarkerSize',25); hold on; 
plot(isoAllMonkeys(idx==2),edgeFreqAll(idx==2),'m.','MarkerSize',25);
plot(isoAllMonkeys(idx==3),edgeFreqAll(idx==3),'y.','MarkerSize',25);

%% Load Ephys+imaging data and get EEG powers from 20-60 Hz
hemisphere = 'Left'; spatialBin = 3;
for iM = 1:2
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

            isoLevel = ([0.75 0.75 NaN NaN NaN;... % Date x run
                0.75 0.75 NaN NaN NaN; ...
                0.9 0.9 0.9 0.9 0.9; ...
                0.75 0.8 0.8 0.9 NaN]);

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

            isoLevel = ([1 1.1 1.25 NaN NaN NaN NaN; ... % Date x run
                1.2 1 1 1.1 1.1 1.1 1.1 ; ...
                0.8 1.05 1.75 1.75 NaN NaN NaN; ...
                0.9 1 1 1 1 1 1; ...
                0.7 0.8 1 1.1 1.3 1.3 1.3]);
    end
    goodRuns = reshape(goodRuns,[size(goodRuns,1)*size(goodRuns,2) 1]);
    goodRuns(isnan(goodRuns)) = []; goodRuns = logical(goodRuns);

    singleChFlag = reshape(singleChFlag,[size(singleChFlag,1)*size(singleChFlag,2) 1]);
    singleChFlag(isnan(singleChFlag)) = []; singleChFlag = logical(singleChFlag);

    isoLevel = reshape(isoLevel,[size(isoLevel,1)*size(isoLevel,2) 1]);
    isoLevel(isnan(isoLevel)) = [];
    isoLevelGoodRuns = isoLevel; isoLevelGoodRuns(~goodRuns) = [];

    % Get monkey data and parameters
    [allDates,allRuns, refDate, refDir,lensCombo, roiSize, ephysFileNameAll, serverPath,probeLabel,...
        chInCortexNotes, greenMapRef] = getMonkeyParams_Imaging_Ephys(monkeyName,commonDir, hemisphere);

    % Get monkey data....
    [~,~,probe,badCh,badTimesLFP,~,estChInCortex] = ...
        getAllData_Imaging_Ephys(monkeyName,hemisphere,commonDir,serverPath,allDates,allRuns,...
        ephysFileNameAll,greenMapRef,chInCortexNotes,probeLabel,spatialBin);

    clc; disp(['All physiology and imaging data for ' monkeyName ' loaded']);

    % Calculating EEG powers
    % Get filter parameters...
    clear fovCorr roiCorr gIntraProbeCorr gPowerCorr gInfraSlowCorr eegMaxPow meanFreq
    fs = 1e3;
    gammaBand   = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); % Gamma band filtering parameters
    params.Fs       = fs;
    params.fpass    = [2 120];
    params.pad      = -1;
    params.tapers   = [3 5];
    params.trialave = 0;

    gIntraProbeCorr = NaN(size(probe,2), size(probe,1));
    gPowerCorr      = NaN(size(probe,2), size(probe,1));
    gInfraSlowCorr  = NaN(size(probe,2), size(probe,1));
    eegMaxPow       = NaN(size(probe,2), size(probe,1));
    meanFreq        = NaN(size(probe,2), size(probe,1));
    
    clear fovCorr roiCorr

    fovCorr = [allMonkeyVars(iM).peakNegValsAllT(:,:,2) ];
    roiCorr = [allMonkeyVars(iM).allChCorr];


    for iDate = 1:size(allDates,1)
        clear expDate eegCh
        expDate = allDates(iDate,:);

        for iRun = 1:size(allRuns{iDate,1})
            clear runName dataDir probeCh rawCh lfpBadTimes lfpBadCh chInCortex

            runName = allRuns{iDate,1}(iRun,:);
            dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

            clc; disp(['Calculating powers and correlations for  ' monkeyName ' '  expDate ' run: ' runName]);

            % Get run related variables
            probeCh     = probe{iRun,iDate}.probeCh;
            eegChT       = probe{iRun,iDate}.eegCh;
            lfpBadTimes = badTimesLFP{iDate,iRun};
            badChDat    = badCh{iDate,iRun};

            probeCh(:,badChDat)    = [];
            probeCh(lfpBadTimes,:) = [];

            chInCortex = estChInCortex{1,iDate}(iRun,:);
            gammaCh    = filtfilt(bG,aG,double(probeCh(:,chInCortex(1):chInCortex(2))));
            gPower     = envelope(gammaCh,5);
            gInfraSlow = getInfraSlowPowerLFP(probeCh,bG,aG,chInCortex);

            gIntraProbeCorr(iDate,iRun) = median(corr(gammaCh),'all','omitnan');
            gPowerCorr(iDate,iRun)      = median(corr(gPower),'all','omitnan');
            gInfraSlowCorr(iDate,iRun)  = median(corr(gInfraSlow),'all','omitnan');
            
            eegChT(lfpBadTimes,:)   = []; % Remove bad times from LFP
          
%             clear spec powMeansS badTimeThresh badTimeIndOld badTimeInd
%             [spec,timeValsSpec,~] = mtspecgramc(eegChT,[5 2],params);
% 
%             powMeanS = squeeze(sum(spec,2));
%             badTimeThresh   = (median(powMeanS,1)+4.5*mad(powMeanS,1,1));
%             badTimeIndOld = floor(timeValsSpec(powMeanS>badTimeThresh))*1e3;
%             badTimes = [];
%             
%             if ~isempty(badTimeIndOld)
%                 badTimes = [];
%                 badTimeInd =[(badTimeIndOld-1e3)'  (badTimeIndOld+1e3)']; % Taking one second before and one second after bad time segments
% 
%                 for iL = 1:size(badTimeInd,1)
%                     badTimes = [badTimes badTimeInd(iL,1): badTimeInd(iL,2)];
%                 end
%                 badTimes = unique(badTimes);
%             end
% 
%              eegChT(badTimes,:) = [];
             %eegCh(:,iRun) = eegChT(1:700e3,:); 

             paramsNew = params;%paramsNew.fpass = [2 120];
             [spec,timeValsSpec,freqValsSpec] = mtspecgramc(eegChT,[2 1],paramsNew);
             cumPowFreq = cumsum(spec,2);
             mat = cumPowFreq>= 0.9.*cumPowFreq(:,end);
             [~,highPowerIdx] = max(mat,[],2);
             meanFreq(iDate,iRun) = median(freqValsSpec(highPowerIdx));
            
             if iM==2
                 specW{iDate,iRun} = spec;
             end

%         end
% 
%         % Normalize it within a day
% %         clear eegChNorm
% %         if size(eegCh,2)~=1
% %             eegChNorm = (eegCh - median(eegCh,2,'omitnan'))./mad(eegCh,1,2);
% %         else
%             eegChNorm = eegCh;
% %         end
% 
%         for iRun = 1:size(allRuns{iDate},1)
            clear spec timeValsSpec freqValsSpec
            [spec,timeValsSpec,freqValsSpec] = mtspecgramc(eegChT,[5 2],params);          
            freqIdx = freqValsSpec>=20 & freqValsSpec<=60;
            eegMaxPow(iDate,iRun) = 20.*log10(squeeze(median(spec(:,freqIdx),[1 2],'omitnan')));
        end 
    end

    eegMaxPow = reshape(eegMaxPow,[size(eegMaxPow,1)*size(eegMaxPow,2) 1]);
    eegMaxPow(isnan(eegMaxPow)) = [];  eegMaxPow(~goodRuns) = [];

    gIntraProbeCorr = reshape(gIntraProbeCorr,[size(gIntraProbeCorr,1)*size(gIntraProbeCorr,2) 1]);
    gIntraProbeCorr(isnan(gIntraProbeCorr)) = [];  gIntraProbeCorr(~goodRuns) = [];

    gPowerCorr = reshape(gPowerCorr,[size(gPowerCorr,1)*size(gPowerCorr,2) 1]);
    gPowerCorr(isnan(gPowerCorr)) = [];  gPowerCorr(~goodRuns) = [];

    gInfraSlowCorr = reshape(gInfraSlowCorr,[size(gInfraSlowCorr,1)*size(gInfraSlowCorr,2) 1]);
    gInfraSlowCorr(isnan(gInfraSlowCorr)) = [];  gInfraSlowCorr(~goodRuns) = [];

    meanFreq = reshape(meanFreq,[size(meanFreq,1)*size(meanFreq,2) 1]);
    meanFreq(isnan(meanFreq)) = [];  meanFreq(~goodRuns) = [];

    % Plotting for individual monkeys
    figure; subplot(131); showLinearFit(isoLevelGoodRuns,gIntraProbeCorr,1.5,0.8,0.7); % Iso vs Ephys params
    xlabel('Iso %'); ylabel('Correlations');
    axis square;title('iso vs gamma time series'); xlim([0.6 1.8]); ylim([-1 1]);

    subplot(132); showLinearFit(isoLevelGoodRuns,gPowerCorr,1.5,0.8,0.7);
    xlabel('Iso %'); ylabel('Correlations');
    axis square; title('iso vs gamma power'); xlim([0.6 1.8]);ylim([-1 1]);

    subplot(133);showLinearFit(isoLevelGoodRuns,gInfraSlowCorr,1.5,0.8,0.7); 
    xlabel('Iso %'); ylabel('Correlations');
    axis square;title('iso vs gamma infraslow power');xlim([0.6 1.8]);ylim([-1 1]);

    sgtitle(['iso % vs ephys parameters for ' monkeyName]);

    maxLim = round(max(eegMaxPow)); minLim = floor(min(eegMaxPow)); 
    figure; subplot(131); showLinearFit(eegMaxPow,gIntraProbeCorr,maxLim-10,0.9,0.8); % EEG power vs Ephys params
    xlabel('EEG metric'); ylabel('Correlations');
    axis square; title('EEG power vs gamma time series');xlim([minLim-10 maxLim+10]);ylim([-1 1]);

    subplot(132); showLinearFit(eegMaxPow,gPowerCorr,maxLim-10,0.9,0.8);
    xlabel('EEG metric'); ylabel('Correlations');
    axis square; title('EEG power vs gamma power');xlim([minLim-10 maxLim+10]);ylim([-1 1]);
    
    subplot(133); showLinearFit(eegMaxPow,gInfraSlowCorr,maxLim-10,0.9,0.8);
    axis square;xlabel('EEG metric'); ylabel('Correlations'); 
    title('EEG power vs gamma infraslow power');xlim([minLim-10 maxLim+10]);ylim([-1 1]);
    
    sgtitle(['EEG Power (au) vs ephys parameters for ' monkeyName]);

    figure; % EEG Power vs FOV
    for iBand = 1:5
        clear coeff xFit yFit mdl
        subplot(2,3,iBand);
        showLinearFit(eegMaxPow,fovCorr(:,iBand),maxLim-10,0.2,0.1);
        xlabel('EEG metric'); ylabel('Correlations'); xlim([minLim-10 maxLim+10]);
        title(bandLabels{iBand}); box off; ylim([-1 0.3]);
    end
     sgtitle(['FOV - ' monkeyName]);

    figure; % EEG Power vs ROI 
    for iBand = 1:5
        clear coeff xFit yFit mdl
        subplot(2,3,iBand);
        showLinearFit(eegMaxPow,roiCorr(:,iBand),maxLim-10,0,-0.1);
        xlabel('EEG metric'); ylabel('Correlations'); xlim([minLim-10 maxLim+10]);
        title(bandLabels{iBand}); box off; ylim([-0.7 0.1]); 
    end
  sgtitle(['ROI - ' monkeyName]);

    figure; % Iso vs FOV 
    for iBand = 1:5
        clear coeff xFit yFit mdl
        subplot(2,3,iBand);
        showLinearFit(isoLevelGoodRuns,fovCorr(:,iBand),1.5,0.2,0.1);
        xlabel('Iso %'); ylabel('Correlations'); xlim([0.6 1.8]);
        title(bandLabels{iBand}); box off; ylim([-1 0.3]);
    end
   sgtitle(['FOV - ' monkeyName]);

    figure; % Iso vs ROI 
    for iBand = 1:5
        clear coeff xFit yFit mdl
        subplot(2,3,iBand);
        showLinearFit(isoLevelGoodRuns,roiCorr(:,iBand),1.5,0,-0.1);
        xlabel('Iso %'); ylabel('Correlations'); xlim([0.6 1.8]);
        title(bandLabels{iBand}); box off; ylim([-0.7 0.1]);
    end
    sgtitle(['ROI - ' monkeyName]);
  
    if iM==1
        eegMaxPowC       = eegMaxPow;
        gIntraProbeCorrC = gIntraProbeCorr;
        gInfraSlowCorrC = gInfraSlowCorr;
        gPowerCorrC      = gPowerCorr;
        isoC             = isoLevelGoodRuns;
        meanFreqC        = meanFreq;
    else
        eegMaxAll          = [eegMaxPowC;eegMaxPow];
        gIntraProbeCorrAll = [gIntraProbeCorrC; gIntraProbeCorr];
        gInfraSlowCorrAll = [gInfraSlowCorrC; gInfraSlowCorr];
        gPowerCorrAll      = [gPowerCorrC; gPowerCorr];
        isoAll             = [isoC; isoLevelGoodRuns];
        meanFreqAll        = [meanFreqC; meanFreq];

    end
end

% figure; plot(isoLevelGoodRuns,eegMaxPow,'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]);
% figure; plot(isoLevelGoodRuns,gIntraProbeCorr,'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]);

figure; scatter(gIntraProbeCorrAll,gPowerCorrAll,'filled');
figure; scatter(gIntraProbeCorrAll,gInfraSlowCorrAll,'filled');

figure; scatter(isoAll, meanFreqAll,'filled');

%% Plot all data from Whiskey, Charlie with Bordeaux and Rambo
% isoAllMonkeys   = [isoB ;isoR; animalStateIso'];
% colorAllMonkeys = [colorB; colorR; [1;2;3;1]];
% specPowAll      = [specBPow; specRPow; animalStateCheckW'];
% edgeFreqAll     = [meanFreqB; meanFreqR; meanFreqW']; 

isoHigh = isoAll>=2; % to ensure that all data points are included
iso = [isoAll(~isoHigh); isoAllMonkeys];
mFreq = [meanFreqAll(~isoHigh); edgeFreqAll];
eegPowAll = [eegMaxAll(~isoHigh); specPowAll];

figure; scatter(iso,mFreq,70,[ones(size(meanFreqAll(~isoHigh))); colorAllMonkeys],'filled');
ylabel('Edge Frequency (Hz)'); xlabel('Iso (%)'); ylim([0 25]); xlim([0.5 3.5]);

modelfun = @(b,x) b(1) * exp(-b(2).*x);  
beta0 = [10 2]; 
mdl = fitnlm(iso,mFreq, modelfun, beta0);
X = 0.5:0.1:3.5;
coefficients = mdl.Coefficients{:, 'Estimate'};

yFitted = coefficients(1) * exp(-coefficients(2).*X) ;
hold on;
plot(X, yFitted, 'r-', 'LineWidth', 2);
text(2.5, 22,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
text(2.5, 20,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);

figure; scatter(iso,eegPowAll,60,[ones(size(meanFreqAll(~isoHigh))); colorAllMonkeys],'filled');
xlabel('Iso (%)'); ylabel('EEG metric (dB)'); ylim([-20 50]); xlim([0.5 3.5]);

%%
figure; subplot(141);
scatter(iso,mFreq,70,[ones(size(meanFreqAll(~isoHigh))); colorAllMonkeys],'filled');
ylabel('Edge Frequency (Hz)'); xlabel('Iso (%)'); ylim([0 25]); xlim([0.5 3.5]);

modelfun = @(b,x) b(1) * exp(-b(2).*x);  
beta0 = [10 2]; 
mdl = fitnlm(iso,mFreq, modelfun, beta0);
X = 0.5:0.1:3.5;
coefficients = mdl.Coefficients{:, 'Estimate'};

yFitted = coefficients(1) * exp(-coefficients(2).*X) ;
hold on;
plot(X, yFitted, 'r-', 'LineWidth', 2);
text(2.5, 22,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
text(2.5, 20,['p-val: ' num2str(mdl.Coefficients.pValue(2))]); axis square;
xlim([0.5 3.5]); xticks(0.5:0.5:3.5); ylim([2 22]); yticks(2:2:22);

fovCorr = [allMonkeyVars(1).peakNegValsAllT(:,:,2) ; allMonkeyVars(2).peakNegValsAllT(:,:,2)];     
roiCorr = [allMonkeyVars(1).allChCorr ; allMonkeyVars(2).allChCorr];

oneIdx = (gInfraSlowCorrAll== 1); 
infraPow = gInfraSlowCorrAll; infraPow(oneIdx) =[]; 
meanFreqAllNew = meanFreqAll; meanFreqAllNew(oneIdx) = [];
isoAllNew = isoAll; isoAllNew(oneIdx) = [];
eegPowNew  = eegMaxAll; eegPowNew(oneIdx) = [];

subplot(142); 
lowID = infraPow<0.1; 
showLinearFit(meanFreqAllNew(~lowID),infraPow(~lowID),20,0.4,0.35);
xlabel('Edge Frequency'); ylabel('Correlations');
axis square;%title('Edge freq vs infraslow power correlations'); 
xlim([5 22]); ylim([0.4 1]); xticks(0:2:22); yticks(0:0.1:1);

subplot(143);
showLinearFit(meanFreqAll, roiCorr(:,4),20,0,-0.1)
xlabel('Edge Frequency'); ylabel('Correlations');
axis square;%title('Edge freq vs ROI correlations - gamma');
xlim([5 22]); ylim([-0.5 0]); xticks(0:2:22); yticks(-1:0.1:0);

subplot(144); 
showLinearFit(meanFreqAll,fovCorr(:,4),20,0,-0.1)
xlabel('Edge Frequency'); ylabel('Correlations');
axis square;%title('Edge freq vs FOV correlations - gamma'); 
xlim([5 22]); ylim([-1 -0.2]); xticks(0:2:22); yticks(-1:0.1:0);

%% EEG Metric

figure; subplot(141);
scatter(iso,eegPowAll,70,[ones(size(eegMaxAll(~isoHigh))); colorAllMonkeys],'filled');
ylabel('EEG Metric'); xlabel('Iso (%)'); ylim([-20 40]); xlim([0.5 3.5]);

modelfun = @(b,x) b(1) * exp(-b(2).*x);  
beta0 = [10 2]; 
mdl = fitnlm(iso,eegPowAll, modelfun, beta0);
X = 0.5:0.1:3.5;
coefficients = mdl.Coefficients{:, 'Estimate'};

yFitted = coefficients(1) * exp(-coefficients(2).*X) ;
hold on;
plot(X, yFitted, 'r-', 'LineWidth', 2);
text(2.5, 22,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
text(2.5, 20,['p-val: ' num2str(mdl.Coefficients.pValue(2))]); axis square;
xlim([0.5 3.5]); xticks(0.5:0.5:3.5); ylim([-20 40]); yticks(-20:5:40);

subplot(142);
lowID = infraPow<0.1; 
showLinearFit(eegPowNew(~lowID),infraPow(~lowID),35,0.4,0.35);
xlabel('EEG Metric'); ylabel('Correlations');
axis square; xlim([20 40]); xticks(-20:2:38); ylim([0.4 1]);

subplot(143);
showLinearFit(eegMaxAll, roiCorr(:,4),35,0,-0.1)
xlabel('EEG Metric'); ylabel('Correlations');
axis square; xlim([20 38]);xticks(-20:2:38);ylim([-0.6 0]);

subplot(144); 
showLinearFit(eegMaxAll,fovCorr(:,4),35,0,-0.1)
xlabel('EEG Metric'); ylabel('Correlations');
axis square; xlim([20 38]);xticks(-20:2:38);ylim([-1 -0.1]);


% figure; % EEG Power vs FOV
% for iBand = 4%1:5
%     clear coeff xFit yFit mdl
%      %subplot(2,3,iBand);
%     showLinearFit(eegMaxAll,fovCorr(:,iBand),35,0,-0.1); xlim([20 38]);ylim([-1 0]) 
%     xlabel('EEG metric'); ylabel('Correlations');  xlim([minLim-10 maxLim+10]);
%     title(bandLabels{iBand}); box off; ylim([-1 0.3]);
% end
% sgtitle('FOV - Combined');
% 
% figure; % EEG Power vs ROI
% for iBand = 4%1:5
%     clear coeff xFit yFit mdl
% %     subplot(2,3,iBand);
%     showLinearFit(eegMaxAll,roiCorr(:,iBand),35,0,-0.1);xlim([20 38]);ylim([-0.6 0]) 
%     xlabel('EEG metric'); ylabel('Correlations');  xlim([minLim-10 maxLim+10]);
%     title(bandLabels{iBand}); box off; ylim([-0.7 0.1]); 
% end
% sgtitle('ROI - Combined');

%% Plot
fovCorr = [allMonkeyVars(1).peakNegValsAllT(:,:,2) ; allMonkeyVars(2).peakNegValsAllT(:,:,2)];     
roiCorr = [allMonkeyVars(1).allChCorr ; allMonkeyVars(2).allChCorr];

figure; subplot(131); showLinearFit(isoAll,gIntraProbeCorrAll,1.5,0.8,0.7); % Iso vs Ephys params
xlabel('Iso %'); ylabel('Correlations');
axis square;title('iso vs gamma time series'); xlim([0.6 1.8]); ylim([-1 1]);

subplot(132); showLinearFit(isoAll,gPowerCorrAll,1.5,0.8,0.7);
xlabel('Iso %'); ylabel('Correlations');
axis square; title('iso vs gamma power'); xlim([0.6 1.8]);ylim([-1 1]);

subplot(133);showLinearFit(isoAll,gInfraSlowCorrAll,1.5,0.8,0.7);
xlabel('Iso %'); ylabel('Correlations');
axis square;title('iso vs gamma infraslow power');xlim([0.6 1.8]);ylim([-1 1]);

sgtitle('iso % vs ephys parameters for all monkeys');



%% 
isoNew = isoAll; 
highIso = isoNew== 1.75;
isoNew(highIso) = [];
eegNew = eegMaxAll; eegNew(highIso) = [];
figure; showLinearFit(isoNew,eegNew,1.3,38,37);


%%
maxLim = round(max(eegMaxAll)); minLim = floor(min(eegMaxAll)); 
figure; subplot(131); showLinearFit(eegMaxAll,gIntraProbeCorrAll,maxLim-10,0.9,0.8); % EEG power vs Ephys params
xlabel('EEG metric'); ylabel('Correlations');
axis square; title('EEG power vs gamma time series'); xlim([minLim-10 maxLim+10]);ylim([-1 1]);

subplot(132); showLinearFit(eegMaxAll,gPowerCorrAll,maxLim-10,0.9,0.8);
xlabel('EEG metric'); ylabel('Correlations');
axis square; title('EEG power vs gamma power'); xlim([minLim-10 maxLim+10]);ylim([-1 1]);

subplot(133); 
%%
oneIdx = (gInfraSlowCorrAll== 1); 
infraPow = gInfraSlowCorrAll; infraPow(oneIdx) =[]; 
eegMaxNew = eegMaxAll; eegMaxNew(oneIdx) = [];
isoAllNew = isoAll; isoAllNew(oneIdx) = [];

%%
figure;
showLinearFit(eegMaxNew,infraPow,40,0.9,0.8);
xlabel('EEG metric'); ylabel('Correlations');
title('EEG power vs gamma infraslow power'); xlim([minLim-10 maxLim+10]);ylim([-1 1]);
%%
highIsoNew = isoAllNew== 1.75;
figure; showLinearFit(isoAllNew(~highIsoNew),infraPow(~highIsoNew),1.3,0.9,0.8);

% sgtitle('EEG Power (dB) vs ephys parameters for all monkeys');
%%
figure; % EEG Power vs FOV
for iBand = 4%1:5
    clear coeff xFit yFit mdl
     %subplot(2,3,iBand);
    showLinearFit(eegMaxAll,fovCorr(:,iBand),35,0,-0.1); xlim([20 38]);ylim([-1 0]) 
    xlabel('EEG metric'); ylabel('Correlations');  xlim([minLim-10 maxLim+10]);
    title(bandLabels{iBand}); box off; ylim([-1 0.3]);
end
sgtitle('FOV - Combined');
%%
figure; % EEG Power vs ROI
for iBand = 4%1:5
    clear coeff xFit yFit mdl
%     subplot(2,3,iBand);
    showLinearFit(eegMaxAll,roiCorr(:,iBand),35,0,-0.1);xlim([20 38]);ylim([-0.6 0]) 
    xlabel('EEG metric'); ylabel('Correlations');  xlim([minLim-10 maxLim+10]);
    title(bandLabels{iBand}); box off; ylim([-0.7 0.1]); 
end
sgtitle('ROI - Combined');

%%
figure; % Iso vs FOV
for iBand = 4%1:5
    clear coeff xFit yFit mdl
%     subplot(2,3,iBand);
    showLinearFit(isoAll(~highIso),fovCorr(~highIso,iBand),1.3,0,-0.1); %ylim([20 38]);xlim([0.6 1.8]);
    xlabel('Iso %'); ylabel('Correlations'); 
    title(bandLabels{iBand}); box off; ylim([-1 0.3]);
end
sgtitle('FOV - Combined');
%%
figure; % Iso vs ROI
for iBand = 4%1:5
    clear coeff xFit yFit mdl
%     subplot(2,3,iBand);
    showLinearFit(isoAll(~highIso),roiCorr(~highIso,iBand),1.3,0,-0.1);
    xlabel('Iso %'); ylabel('Correlations'); xlim([0.6 1.8]);
    title(bandLabels{iBand}); box off; ylim([-0.6 0]); 
end
sgtitle('ROI - Combined');

%%
figure; showLinearFit(isoAll,eegMaxAll,1.5,38,37);
ylim([20 38]);xlim([0.6 1.8]); xlabel('Iso %'); ylabel('EEG metric');

%%  Function to fit a line
function showLinearFit(xVal,yVal,textLocX,textLocY1,textLocY2)
    plot(xVal,yVal,'o','MarkerSize',7,'MarkerFaceColor',[0 0.4470 0.7410]); hold on; box off; 
    coeff = polyfit(xVal,yVal,1);
    xFit = linspace(min(xVal),max(xVal),1000); 
    yFit = polyval(coeff,xFit); mdl = fitlm(xVal,yVal);
    plot(xFit,yFit,'-k','LineWidth',1);
    text(textLocX, textLocY1,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
    text(textLocX, textLocY2,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);
end
