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

params.Fs       = 1e3; % PSD parameter initiation
params.fpass    = [1 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;

eegPSD = NaN(120,length(fileNum));
for iFile = 1:length(fileNum)
    clear eegCh eegPSDT
    var =  load([saveFolder '\datafile000'  num2str(fileNum(iFile)) '_lfp.mat']);
    eegCh = var.eegCh;

    clear spec powMeansS badTimeThresh badTimeIndOld badTimeInd
    [spec,timeValsSpec,freqValsSpec] = mtspecgramc(eegCh,[5 2],params);

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
    [eegPSDT,psdFreq]  = mtspectrumc(buffer(eegCh,1000),params);

    eegPSDCW(:,iFile)= median(eegPSDT,2);

    [spec,timeValsSpec,freqValsSpec] = mtspecgramc(eegCh,[5 2],params);
    specAllCW(iFile).spec = spec; specAllCW(iFile).timeVals = timeValsSpec; specAllCW(iFile).freqVals = freqValsSpec;

    % Get power for frequencies from 17-27 Hz
    idx = freqValsSpec>=30 & freqValsSpec<=80;

    % Z-score the eeg relative to the first 100s of the recording
    specVals = 20.*log10(spec(:,idx)); %specVals = (specVals- median(specVals(1:50,:),1))./mad(specVals(1:50,:),1);
    animalStateCheckW(iFile) = median(20.*log10(eegPSDCW(100:120,iFile)),'all','omitnan');%median(specVals,'all'); %#ok<SAGROW>
end 

% figure;plot(1:120,20.*log10(eegPSDCW(:,1)))
% hold on; box off;plot(1:120,20.*log10(eegPSDCW(:,3)));
% ylim([-50 100]); yticks(-50:10:100);xticks(0:10:120);
% ylabel('Power (dB)'); xlabel('Frequency (Hz)');

d = 20.*log10(eegPSDCW(:,1)) - 20.*log10(eegPSDCW(:,3));
% Maximum difference occurs at 22 Hz. 
% Hence choosing frequencies between 17-27 Hz - 22+/- 5Hz

animalStateIso = [1.2 2.25 3.25 1];
figure; scatter(animalStateIso,animalStateCheckW,'filled');

specOnly = {specAllCW(:).spec}; 
specMedian = 20.*log10(cell2mat(cellfun(@(x) median(x,1,'omitnan'),specOnly,'un',0)'));

diffLightDeep = (specMedian(1,:)- specMedian(3,:))';
diffLightSemi = (specMedian(1,:)- specMedian(2,:))';

%% Getting the envelope of gamma band and plotting the psd - Whiskey
gammaBand   = [30 90]; [bG,aG] = butter(3,gammaBand./(1e3/2),'bandpass'); % Gamma band filtering parameters

var =  load([saveFolder '\datafile000'  num2str(fileNum(1)) '_lfp.mat']);
eegLight = var.eegCh;eegLight = envelope(filtfilt(bG,aG,eegLight),100,'peak');

var =  load([saveFolder '\datafile000'  num2str(fileNum(2)) '_lfp.mat']);
eegSemi = var.eegCh; eegSemi = envelope(filtfilt(bG,aG,eegSemi),100,'peak');

var =  load([saveFolder '\datafile000'  num2str(fileNum(3)) '_lfp.mat']);
eegDeep = var.eegCh;eegDeep = envelope(filtfilt(bG,aG,eegDeep),100,'peak');

params.Fs       = 1e3; % PSD parameter initiation
params.fpass    = [0 30];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;

[lightPow,~] = mtspectrumc(buffer(eegLight,1000),params); lightPow = median(lightPow,2);
[semiPow,~]  = mtspectrumc(buffer(eegSemi,1000),params); semiPow = median(semiPow,2);
[deepPow,psdFreq]  = mtspectrumc(buffer(eegDeep,1000),params); deepPow= median(deepPow,2);


%% Get data from Bordeaux and Rambo
for iM = 1:2
    clear monkeyName allDates allRuns rippleName isoLevels eegPSD specAll animalStateCheck
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
          
        
%         clear eegChNorm
%         if size(eegCh,2)~=1
%             eegChNorm = (eegCh - median(eegCh,2,'omitnan'))./mad(eegCh,1,2);
%         else
            eegChNorm = eegCh;
%         end

        for iRun = 1:size(allRuns{iDate},1)
            [eegPSDT,psdFreq]  = mtspectrumc(buffer(eegChNorm(:,iRun),1000),params);
            eegPSD(:,iDate,iRun)= median(eegPSDT,2,'omitnan');

            [spec,timeValsSpec,freqValsSpec] = mtspecgramc(eegChNorm(:,iRun),[5 2],params);
            specAll(:,:,iDate,iRun)= spec; %specAll(iDate,iRun).timeVals = timeValsSpec;
            %specAll(iDate,iRun).freqVals = freqValsSpec;

            % Get power for frequencies from 17-27 Hz
            idx = freqValsSpec>=80 & freqValsSpec<=100;

            % Z-score the eeg relative to the first 100s of the recording
            specVals = 20.*log10(spec(:,idx)); %specVals = (specVals- median(specVals(1:50,:),1))./mad(specVals(1:50,:),1);
%             animalStateCheck(iDate,iRun) = median(specVals,'all','omitnan'); %#ok<SAGROW>
            animalStateCheck(iDate,iRun) = median(20.*log10(eegPSD(80:120,iDate,iRun)),'all','omitnan');
        end   

    end

    switch iM
        case 1
            eegB  = eegPSD;%reshape(eegPSD,[120 size(eegPSD,2)*size(eegPSD,3)]);
            specB = specAll; %reshape(specAll,[size(eegPSD,2)*size(eegPSD,3) 1]);
            isoB  = [1.1 2.25 1.25 1.75; 0.6 NaN NaN NaN; 1 2 1 NaN];%reshape([1.1 2.25 1.25 1.75; 0.6 NaN NaN NaN; 1 2 1 NaN],[size(eegPSD,2)*size(eegPSD,3) 1]);
            animalB = animalStateCheck; 
        case 2
            eegR  = eegPSD; %reshape(eegPSD,[120 size(eegPSD,2)*size(eegPSD,3)]);
            specR = specAll; %reshape(specAll,[size(eegPSD,2)*size(eegPSD,3) 1]);
            isoR  = [0.9 2 2 1; 1.5 2 1.4 1.9];%reshape([0.9 2 2 1; 1.5 2 1.4 1.9],[size(eegPSD,2)*size(eegPSD,3) 1]);
            animalR = animalStateCheck;
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

animalB = reshape(animalB,[12 1]);
animalR = reshape(animalR,[8 1]);
isoR  = reshape(isoR,[8 1]);
isoB  = reshape(isoB,[12 1]);
nanVals = isnan(isoB);
isoB(nanVals) = [];
animalB(nanVals) = [];
figure; scatter(isoB,animalB,'filled')
figure; scatter(isoR,animalR,'filled')

%% Load Ephys+imaging data and get EEG powers from 17-27 Hz
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
    clear fovCorr roiCorr gIntraProbeCorr gPowerCorr gInfraSlowCorr eegMaxPow
    fs = 1e3;
    gammaBand   = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); % Gamma band filtering parameters
    params.Fs       = fs;
    params.fpass    = [1 120];
    params.pad      = -1;
    params.tapers   = [3 5];
    params.trialave = 0;

    gIntraProbeCorr = NaN(size(probe,2), size(probe,1));
    gPowerCorr      = NaN(size(probe,2), size(probe,1));
    gInfraSlowCorr  = NaN(size(probe,2), size(probe,1));
    eegMaxPow       = NaN(size(probe,2), size(probe,1));
    
    clear fovCorr roiCorr

    fovCorr = [allMonkeyVars(iM).peakNegValsAllT(:,:,2) ];
    roiCorr = [allMonkeyVars(iM).allChCorr];


    for iDate = 1:size(allDates,1)
        clear expDate
        expDate = allDates(iDate,:);

        for iRun = 1:size(allRuns{iDate,1})
            clear runName dataDir probeCh rawCh lfpBadTimes lfpBadCh chInCortex

            runName = allRuns{iDate,1}(iRun,:);
            dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

            clc; disp(['Calculating powers and correlations for  ' monkeyName ' '  expDate ' run: ' runName]);

            % Get run related variables
            probeCh     = probe{iRun,iDate}.probeCh;
            eegCh       = probe{iRun,iDate}.eegCh;
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

            eegCh(lfpBadTimes,:)   = []; % Remove bad times from LFP
%             eegCh = (eegCh-mean(eegCh))./std(eegCh); % Z-score
%             [spec,timeValsSpec,freqValsSpec] = mtspecgramc(eegCh,[5 2],params);
% 
%             % Get power for frequencies from 17-27 Hz
%             idx = freqValsSpec>=17 & freqValsSpec<=27;
% 
%             % Z-score the eeg relative to the first 100s of the recording
%             specVals = spec(:,idx); %specVals = (specVals- median(specVals(1:50,:),1))./mad(specVals(1:50,:),1);
            eegMaxPow(iDate,iRun) = real(20.*log10(median(eegCh(100:120,:),'all','omitnan'))); 
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
%%
    figure; subplot(131); showLinearFit(eegMaxPow,gIntraProbeCorr,0.8,0.9,0.8); % EEG power vs Ephys params
    xlabel('EEG metric'); ylabel('Correlations');
    axis square; title('EEG power vs gamma time series');%xlim([-2 2]);ylim([-1 1]);

    subplot(132); showLinearFit(eegMaxPow,gPowerCorr,0.8,0.9,0.8);
    xlabel('EEG metric'); ylabel('Correlations');
    axis square; title('EEG power vs gamma power');%xlim([-2 2]);ylim([-1 1]);
    
    subplot(133); showLinearFit(eegMaxPow,gInfraSlowCorr,0.8,0.9,0.8);
    axis square;xlabel('EEG metric'); ylabel('Correlations'); 
    title('EEG power vs gamma infraslow power');%xlim([-2 2]);ylim([-1 1]);
    
    sgtitle(['EEG Power (au) vs ephys parameters for ' monkeyName]);
%%
    figure; % EEG Power vs FOV
    for iBand = 1:5
        clear coeff xFit yFit mdl
        subplot(2,3,iBand);
        showLinearFit(eegMaxPow,fovCorr(:,iBand),0.8,0.2,0.1);
        xlabel('EEG metric'); ylabel('Correlations'); xlim([-2 2]);
        title(bandLabels{iBand}); box off; ylim([-1 0.3]);
    end
     sgtitle(['FOV - ' monkeyName]);

    figure; % EEG Power vs ROI 
    for iBand = 1:5
        clear coeff xFit yFit mdl
        subplot(2,3,iBand);
        showLinearFit(eegMaxPow,roiCorr(:,iBand),0.8,0,-0.1);
        xlabel('EEG metric'); ylabel('Correlations'); xlim([-2 2]);
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
    else
        eegMaxAll          = [eegMaxPowC;eegMaxPow];
        gIntraProbeCorrAll = [gIntraProbeCorrC; gIntraProbeCorr];
        gInfraSlowCorrAll = [gInfraSlowCorrC; gInfraSlowCorr];
        gPowerCorrAll      = [gPowerCorrC; gPowerCorr];
        isoAll             = [isoC; isoLevelGoodRuns];

    end
end

% figure; plot(isoLevelGoodRuns,eegMaxPow,'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]);
% figure; plot(isoLevelGoodRuns,gIntraProbeCorr,'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]);

figure; scatter(gIntraProbeCorrAll,gPowerCorrAll,'filled');
figure; scatter(gIntraProbeCorrAll,gInfraSlowCorrAll,'filled');


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
figure; subplot(131); showLinearFit(eegMaxAll,gIntraProbeCorrAll,300,0.9,0.8); % EEG power vs Ephys params
xlabel('EEG metric'); ylabel('Correlations');
axis square; title('EEG power vs gamma time series'); xlim([0 500]);ylim([-1 1]);

subplot(132); showLinearFit(eegMaxAll,gPowerCorrAll,300,0.9,0.8);
xlabel('EEG metric'); ylabel('Correlations');
axis square; title('EEG power vs gamma power'); xlim([0 500]);ylim([-1 1]);

subplot(133); showLinearFit(eegMaxAll,gInfraSlowCorrAll,300,0.9,0.8);
axis square;xlabel('EEG metric'); ylabel('Correlations');
title('EEG power vs gamma infraslow power'); xlim([0 500]);ylim([-1 1]);

sgtitle('EEG Power (dB) vs ephys parameters for all monkeys');

figure; % EEG Power vs FOV
for iBand = 1:5
    clear coeff xFit yFit mdl
    subplot(2,3,iBand);
    showLinearFit(eegMaxAll,fovCorr(:,iBand),300,0.2,0.1);
    xlabel('EEG metric'); ylabel('Correlations');  xlim([0 500]);
    title(bandLabels{iBand}); box off; ylim([-1 0.3]);
end
sgtitle('FOV - Combined');

figure; % EEG Power vs ROI
for iBand = 1:5
    clear coeff xFit yFit mdl
    subplot(2,3,iBand);
    showLinearFit(eegMaxAll,roiCorr(:,iBand),300,0,-0.1);
    xlabel('EEG metric'); ylabel('Correlations');  xlim([0 500]);
    title(bandLabels{iBand}); box off; ylim([-0.7 0.1]); 
end
sgtitle('ROI - Combined');

figure; % Iso vs FOV
for iBand = 1:5
    clear coeff xFit yFit mdl
    subplot(2,3,iBand);
    showLinearFit(isoAll,fovCorr(:,iBand),1.5,0.2,0.1);
    xlabel('Iso %'); ylabel('Correlations'); xlim([0.6 1.8]);
    title(bandLabels{iBand}); box off; ylim([-1 0.3]);
end
sgtitle('FOV - Combined');

figure; % Iso vs ROI
for iBand = 1:5
    clear coeff xFit yFit mdl
    subplot(2,3,iBand);
    showLinearFit(isoAll,roiCorr(:,iBand),1.5,0,-0.1);
    xlabel('Iso %'); ylabel('Correlations'); xlim([0.6 1.8]);
    title(bandLabels{iBand}); box off; ylim([-0.7 0.1]); 
end
sgtitle('ROI - Combined');
%
figure; showLinearFit(isoAll,eegMaxAll,1.5,301,300);
ylim([-2 2]);xlim([0 500]); xlabel('Iso %'); ylabel('EEG metric');

%%  Function to fit a line
function showLinearFit(xVal,yVal,textLocX,textLocY1,textLocY2)
    plot(xVal,yVal,'o','MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410]); hold on; box off; 
    coeff = polyfit(xVal,yVal,1);
    xFit = linspace(min(xVal),max(xVal),1000); 
    yFit = polyval(coeff,xFit); mdl = fitlm(xVal,yVal);
    plot(xFit,yFit,'-k','LineWidth',1);
    text(textLocX, textLocY1,['R^2 : ' num2str(mdl.Rsquared.Ordinary*100) '%']);
    text(textLocX, textLocY2,['p-val: ' num2str(mdl.Coefficients.pValue(2))]);
end
