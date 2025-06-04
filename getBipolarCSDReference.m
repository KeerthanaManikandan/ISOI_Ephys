clc; clear;
commonDir = 'C:\Users\KEM294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\Imaging']));
addpath(genpath([commonDir '\Codes\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\fly_track\videoIO']));


hemisphere = 'Left'; spatialBin = 3;
iM = 1;
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

goodRunsSpatial = reshape(goodRunsSpatial,[size(goodRunsSpatial,1)*size(goodRunsSpatial,2) 1]);
goodRunsSpatial(isnan(goodRunsSpatial)) = []; goodRunsSpatial = logical(goodRunsSpatial);

smFlag = reshape(smFlag,[size(smFlag,1)*size(smFlag,2) 1]);
smFlag((smFlag == '#')) = [];
sensoryGoodRuns = (smFlag =='S') & goodRuns; 
sensoryGoodSpatialRuns = (smFlag == 'S') & goodRunsSpatial;
motorGoodRuns = (smFlag == 'M') & goodRuns; 
motorGoodSpatialRuns = (smFlag == 'M') & goodRuns; 

isoLevel = reshape(isoLevel,[size(isoLevel,1)*size(isoLevel,2) 1]);
isoLevel(isnan(isoLevel)) = []; 
isoLevelGoodRuns = isoLevel; isoLevelGoodRuns(~goodRuns) = [];
isoLevelSpatial = isoLevel; isoLevelSpatial(~goodRunsSpatial) = [];

% Get monkey data and parameters
[allDates,allRuns, refDate, refDir,lensCombo, roiSize, ephysFileNameAll, serverPath,probeLabel,...
    chInCortexNotes, greenMapRef] = getMonkeyParams_Imaging_Ephys(monkeyName,commonDir, hemisphere);

[processedDat,greenIm,probe,badCh,badTimesLFP,badTimeThresh,estChInCortex] = ...
    getAllData_Imaging_Ephys(monkeyName,hemisphere,commonDir,serverPath,allDates,allRuns,...
    ephysFileNameAll,greenMapRef,chInCortexNotes,probeLabel,spatialBin);

%% Get raw signal data
clc;
tic;
for iDate = 1:size(allDates,1)
    clear expDate
    expDate    = allDates(iDate,:);
    for iRun = 1:size(allRuns{iDate,1},1)
        clear entityInfo datFileName timeStamp
        % Load all run information
        runName    = allRuns{iDate,1}(iRun,:);
        saveFolder = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Processed Data'];
        datFileNum = ephysFileNameAll{iDate,1};
        fileNum    = str2double(datFileNum(iRun,end));
        
        disp(['Obtaining raw data for ' monkeyName ' ' expDate ' ' runName]); 

        % Get the name of stored file
        datFileName = ephysFileNameAll{iDate,1}(iRun,:);
        datFileName = datFileName(1:end-1);

     
        if (fileNum>=10)
            datFileName = datFileName(1:end-1);
        end

        % Check if LFP is already stored and if camera timestamps are also
        % stored...
        if ~exist([saveFolder '\' datFileName num2str(fileNum) '_lfp.mat'],'file')
            vInfo = [];
        else
            vInfo = who('-file',[saveFolder '\' datFileName num2str(fileNum) '_lfp.mat']);
        end

        % Get the raw data
        clear datName

        if exist([serverPath expDate '\Electrophysiology\' runName '\' datFileName num2str(fileNum) '.nev'],'file')
            datName = [serverPath expDate '\Electrophysiology\' runName '\' datFileName num2str(fileNum)];

        elseif exist([serverPath expDate '\Electrophysiology\' datFileName num2str(fileNum) '.nev'],'file')
            datName = [serverPath expDate '\Electrophysiology\' datFileName num2str(fileNum)];
        end

        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName];

        if ~exist([dataDir '\raw1k.mat'],'file') 
            try
                [nsResult,hFile] = ns_OpenFile(datName);
                if ~strcmp(nsResult,'ns_OK')
                    disp('Data file did not open! - going to the next datafile');
                end

                % Get file information
                [nsResult2, fileInfoVar] = ns_GetFileInfo(hFile);

                if ~strcmp(nsResult2,'ns_OK')
                    disp('Data file information did not load!');
                end

                % Get entity information
                for iEntity = 1: fileInfoVar.EntityCount
                    [~,entityInfo(iEntity,1)] = ns_GetEntityInfo(hFile,iEntity);
                end

                % Sort the entities whether they contain events, neural data, segment data
                % Get the label  and list of all the channels storing LFP data
                allList   = find([entityInfo.EntityType] == 2);
                allLabel  = {entityInfo(allList).EntityLabel};

                % Get indices of the LFP only
                if strcmp(allLabel{1}(1:3),'lfp')
                    lfpIdx = cell2mat(cellfun(@(x) strcmp(x(1:3),'lfp'),allLabel,'un',0)); % 32-channel electrode

                elseif strcmp(allLabel{1}(end-2:end),'lfp')
                    lfpIdx = cell2mat(cellfun(@(x) strcmp(x(end-2:end),'lfp'),allLabel,'un',0)); % 32-channel electrode

                elseif strcmp(allLabel{1}(1:6),'analog')
                    lfpIdx = cell2mat(cellfun(@(x) strcmp(x,'analog 2'),allLabel,'un',0)); % Single electrode
                end

                lfpList = allList(lfpIdx);
                lfpLabel = allLabel(lfpIdx);

                if strcmp(lfpLabel{1},lfpLabel{2}) % To remove 30kHz sampled data
                    lfpList(2) = [];
                    lfpLabel(2) = [];
                end

                if length(lfpList)~= 1
                    rawList = allList(~lfpIdx);
                    rawLabel = allLabel(~lfpIdx);
                else
                    rawList = allList(end); % Single channel electrode condition
                    rawLabel = allLabel(end);
                end

                % Remove EEG label from raw data list
                if find(cell2mat((cellfun(@(x)(strcmp(x(1:5),'analo')),rawLabel,'un',0)))) && length(rawList)~= 1 % To include for single electrode
                    loc = find(cell2mat((cellfun(@(x)(strcmp(x(1:5),'analo')),rawLabel,'un',0))));
                    rawList(loc)  = [];
                    rawLabel(loc) = [];
                end

                % Get channel order
                if strcmp(lfpLabel{1}(1:3),'lfp')
                    chNum = str2num(cell2mat(cellfun(@(x) x(end-1:end),lfpLabel','un',0))); %#ok<ST2NM>
                else
                    chNum = str2num(cell2mat(cellfun(@(x) x(2:3),lfpLabel','un',0))); %#ok<ST2NM>
                end
                [~,elecID] = sort(chNum);

                % Get the LFP for all the channels in sorted order
                %             clear b a bS aS bH aH  rawCh
                %             [b,a]   = butter(3,[1 250]./(fs/2),'bandpass'); % Bandpass filtering parameters across 1-250 Hz
                %             [bS,aS] = butter(3,[57 62]./(fs/2),'stop'); % Bandstop filtering between 57-62 Hz
                %             [bH,aH] = butter(3,250./(30e3/2),'high'); % High pass filtering >250 Hz

                clear rawChTemp         
                for iElec = 1:length(lfpLabel) % Get LFP and raw data
                    clear elecEntityID lfpEntityID lfpCount 
                    if ~isempty(elecID)
                        lfpEntityID   = lfpList(elecID(iElec));
                        rawEntityID   = rawList(elecID(iElec));
                    else
                        lfpEntityID   = lfpList(iElec);
                        rawEntityID   = rawList(iElec);
                    end
                    lfpCount      = entityInfo(lfpEntityID).ItemCount;
                    rawCount      = entityInfo(rawEntityID).ItemCount;

                    % Get raw data
                    [~, ~, rawChTemp(:,iElec)] = ns_GetAnalogData(hFile,rawEntityID,1,rawCount);
                end
                

                rawChTemp = downsample(rawChTemp,30); % Downsample raw signals
                [b,a] = butter(3,1./(1e3/2),'high');
                rawChTemp = single(filtfilt(b,a,rawChTemp));  % Highpass frequencies >1 Hz

                % Obtain camera frame information
                eventList   = find([entityInfo.EntityType] == 1);

                % Grab the times where the camera frame occurred
                timeStamp = []; itemIdx  = []; nsFrame ='';

                % Check if the number of camera timestamps are >=9000
                [~,itemIdx] = max([entityInfo(eventList(1:end)).ItemCount]);

                if ~isempty(itemIdx)
                    for iT = 1:entityInfo(eventList(itemIdx)).ItemCount
                        [nsFrame, timeStamp(iT), ~, ~] = ns_GetEventData(hFile, eventList(itemIdx),iT);
                    end
                end

                if ~strcmp(nsFrame,'ns_OK')
                    disp('Camera frames not recorded...');
                end

                % Keep data between first and last camera timestamp
                t1k = int32(floor(timeStamp).*1e3);
                rawChTemp  = single(rawChTemp(t1k(1):t1k(end),:));

            catch
                disp(['Data did not load for : ' expDate ' ' runName]);
            end
            save([dataDir '\raw1k.mat'],'rawChTemp','-v7.3');
            rawCh{iRun,iDate} = matfile([dataDir '\raw1k.mat']);

        else
            rawCh{iRun,iDate} = matfile([dataDir '\raw1k.mat']);

        end

        if exist([dataDir '\raw30k.mat'],'file')
            delete([dataDir '\raw30k.mat']);
        end
    end
end
toc;

%% Get re-referenced data for ROI 

clear tempProfileBipRef tempProfileSuperBipRef tempProfileDeepBipRef tempProfileMidBipRef...
    tempProfileCSDRef tempProfileSuperCSDRef tempProfileDeepCSDRef tempProfileMidCSDRef
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

        if ~exist([dataDir '\ROIAllVars_v2.mat'],'file')
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

            clear tempProfileBipRefRun tempProfileSuperBipRefRun tempProfileDeepBipRefRun tempProfileMidBipRefRun...
                tempProfileCSDRefRun tempProfileSuperCSDRefRun tempProfileDeepCSDRefRun tempProfileMidCSDRefRun

            [tempProfileBipRefRun,tempProfileSuperBipRefRun,... % 6-channel split Bipolar
                tempProfileDeepBipRefRun,tempProfileMidBipRefRun,~] = ...
                getCrossCorrROI_v2(dataDir,monkeyName,expDate,runName,...
                processedDat{iDate,iRun}.tempBandPass,probe{iRun,iDate}.probeCh,rawCh{iRun,iDate}.rawChTemp,...
                badTimesLFP{iDate,iRun},badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),seedLocIn,...
                seedRad,clipMaskROI,6,'BipolarRef');

            [tempProfileCSDRefRun,tempProfileSuperCSDRefRun,... % 6-channel split No-Ref
                tempProfileDeepCSDRefRun,tempProfileMidCSDRefRun,~] = ...
                getCrossCorrROI_v2(dataDir,monkeyName,expDate,runName,...
                processedDat{iDate,iRun}.tempBandPass,probe{iRun,iDate}.probeCh,rawCh{iRun,iDate}.rawChTemp,...
                badTimesLFP{iDate,iRun},badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),seedLocIn,...
                seedRad,clipMaskROI,6,'CSDRef');

            save([dataDir '\ROIAllVars_v2.mat'],'tempProfileBipRefRun','tempProfileSuperBipRefRun',...
                'tempProfileDeepBipRefRun','tempProfileCSDRefRun','tempProfileDeepCSDRefRun',...
                'tempProfileSuperCSDRefRun','tempProfileMidBipRefRun','tempProfileMidCSDRefRun');

            tempProfileBipRef(iDate,iRun)           = tempProfileBipRefRun; %#ok<*SAGROW> 
            tempProfileSuperBipRef(iDate,iRun)      = tempProfileSuperBipRefRun;
            tempProfileMidBipRef(iDate,iRun)        = tempProfileMidBipRefRun;
            tempProfileDeepBipRef(iDate,iRun)       = tempProfileDeepBipRefRun;
            tempProfileCSDRef(iDate,iRun)      = tempProfileCSDRefRun;
            tempProfileSuperCSDRef(iDate,iRun) = tempProfileSuperCSDRefRun;
            tempProfileMidCSDRef(iDate,iRun)   = tempProfileMidCSDRefRun;
            tempProfileDeepCSDRef(iDate,iRun)  = tempProfileDeepCSDRefRun;

        else
            clear allVars
            allVars                             = load([dataDir '\ROIAllVars_v2.mat']);
            tempProfileBipRef(iDate,iRun)       = allVars.tempProfileBipRefRun;
            tempProfileSuperBipRef(iDate,iRun)  = allVars.tempProfileSuperBipRefRun;
            tempProfileMidBipRef(iDate,iRun)    = allVars.tempProfileMidBipRefRun;
            tempProfileDeepBipRef(iDate,iRun)   = allVars.tempProfileDeepBipRefRun;
            tempProfileCSDRef(iDate,iRun)       = allVars.tempProfileCSDRefRun;
            tempProfileSuperCSDRef(iDate,iRun)  = allVars.tempProfileSuperCSDRefRun;
            tempProfileMidCSDRef(iDate,iRun)    = allVars.tempProfileMidCSDRefRun;
            tempProfileDeepCSDRef(iDate,iRun)   = allVars.tempProfileDeepCSDRefRun;          
        end
    end
end

%% Compile the ROI level variables 
clear bandLabels tempProfilesAll medTempProfileAll allChCorr allChCorrCSD ...
    allChLagVal allChLagValCSD  tempProfilesAllCSD allChCorrSuperCSD allChCorrMidCSD...
    allChCorrDeepCSD allChCorrSuper allChCorrMid allChCorrDeep
bandLabels = {'Theta'; 'Alpha'; 'Beta'; 'Gamma'; 'Spiking'};

% Temporal profiles for all frequency bands
tempProfilesAll(:,:,1)  = [tempProfileBipRef.profileTheta]; 
tempProfilesAll(:,:,2)  = [tempProfileBipRef.profileAlpha]; 
tempProfilesAll(:,:,3)  = [tempProfileBipRef.profileBeta]; 
tempProfilesAll(:,:,4)  = [tempProfileBipRef.profile];      
tempProfilesAll(:,:,5)  = [tempProfileBipRef.profileRaw];   
tempProfilesAll(:,~goodRuns,:) = [];

tempProfilesAllCSD(:,:,1)  = [tempProfileCSDRef.profileTheta]; 
tempProfilesAllCSD(:,:,2)  = [tempProfileCSDRef.profileAlpha]; 
tempProfilesAllCSD(:,:,3)  = [tempProfileCSDRef.profileBeta]; 
tempProfilesAllCSD(:,:,4)  = [tempProfileCSDRef.profile];      
tempProfilesAllCSD(:,:,5)  = [tempProfileCSDRef.profileRaw];   
tempProfilesAllCSD(:,~goodRuns,:) = [];

% Show the distribution of lags and peak negative correlations 
allChCorr(:,1) = [tempProfileBipRef.magLowTheta]'; 
allChCorr(:,2) = [tempProfileBipRef.magLowAlpha]'; 
allChCorr(:,3) = [tempProfileBipRef.magLowBeta]'; 
allChCorr(:,4) = [tempProfileBipRef.magLow]'; 
allChCorr(:,5) = [tempProfileBipRef.magLowRaw]'; 
allChCorr(~goodRuns,:) = []; 

allChCorrCSD(:,1) = [tempProfileCSDRef.magLowTheta]'; 
allChCorrCSD(:,2) = [tempProfileCSDRef.magLowAlpha]'; 
allChCorrCSD(:,3) = [tempProfileCSDRef.magLowBeta]'; 
allChCorrCSD(:,4) = [tempProfileCSDRef.magLow]'; 
allChCorrCSD(:,5) = [tempProfileCSDRef.magLowRaw]'; 
allChCorrCSD(~goodRuns,:) = [];

allChLagVal(:,1)= [tempProfileBipRef.lagLowTheta]'; 
allChLagVal(:,2)= [tempProfileBipRef.lagLowAlpha]'; 
allChLagVal(:,3)= [tempProfileBipRef.lagLowBeta]'; 
allChLagVal(:,4) = [tempProfileBipRef.lagLow]'; 
allChLagVal(:,5) = [tempProfileBipRef.lagLowRaw]'; 
allChLagVal(~goodRuns,:) = []; 

allChLagValCSD(:,1)= [tempProfileCSDRef.lagLowTheta]'; 
allChLagValCSD(:,2)= [tempProfileCSDRef.lagLowAlpha]'; 
allChLagValCSD(:,3)= [tempProfileCSDRef.lagLowBeta]'; 
allChLagValCSD(:,4) = [tempProfileCSDRef.lagLow]'; 
allChLagValCSD(:,5) = [tempProfileCSDRef.lagLowRaw]'; 
allChLagValCSD(~goodRuns,:) = []; 

% Show the distribution of lags and peak negative correlations - Superficial
allChCorrSuper(:,1) = [tempProfileSuperBipRef.magLowTheta]'; 
allChCorrSuper(:,2) = [tempProfileSuperBipRef.magLowAlpha]'; 
allChCorrSuper(:,3) = [tempProfileSuperBipRef.magLowBeta]'; 
allChCorrSuper(:,4) = [tempProfileSuperBipRef.magLow]'; 
allChCorrSuper(:,5) = [tempProfileSuperBipRef.magLowRaw]'; 
allChCorrSuper(~goodRuns,:) = []; 

allChCorrSuperCSD(:,1) = [tempProfileSuperCSDRef.magLowTheta]'; 
allChCorrSuperCSD(:,2) = [tempProfileSuperCSDRef.magLowAlpha]'; 
allChCorrSuperCSD(:,3) = [tempProfileSuperCSDRef.magLowBeta]'; 
allChCorrSuperCSD(:,4) = [tempProfileSuperCSDRef.magLow]'; 
allChCorrSuperCSD(:,5) = [tempProfileSuperCSDRef.magLowRaw]'; 
allChCorrSuperCSD(~goodRuns,:) = []; 

% Show the distribution of lags and peak negative correlations - Middle
allChCorrMid(:,1) = [tempProfileMidBipRef.magLowTheta]'; 
allChCorrMid(:,2) = [tempProfileMidBipRef.magLowAlpha]'; 
allChCorrMid(:,3) = [tempProfileMidBipRef.magLowBeta]'; 
allChCorrMid(:,4) = [tempProfileMidBipRef.magLow]'; 
allChCorrMid(:,5) = [tempProfileMidBipRef.magLowRaw]'; 
allChCorrMid(~goodRuns,:) = []; 

allChCorrMidCSD(:,1) = [tempProfileMidCSDRef.magLowTheta]'; 
allChCorrMidCSD(:,2) = [tempProfileMidCSDRef.magLowAlpha]'; 
allChCorrMidCSD(:,3) = [tempProfileMidCSDRef.magLowBeta]'; 
allChCorrMidCSD(:,4) = [tempProfileMidCSDRef.magLow]'; 
allChCorrMidCSD(:,5) = [tempProfileMidCSDRef.magLowRaw]'; 
allChCorrMidCSD(~goodRuns,:) = []; 

% Show the distribution of lags and peak negative correlations - Deep
allChCorrDeep(:,1) = [tempProfileDeepBipRef.magLowTheta]'; 
allChCorrDeep(:,2) = [tempProfileDeepBipRef.magLowAlpha]'; 
allChCorrDeep(:,3) = [tempProfileDeepBipRef.magLowBeta]'; 
allChCorrDeep(:,4) = [tempProfileDeepBipRef.magLow]'; 
allChCorrDeep(:,5) = [tempProfileDeepBipRef.magLowRaw]'; 
allChCorrDeep(~goodRuns,:) = []; 

allChCorrDeepCSD(:,1) = [tempProfileDeepCSDRef.magLowTheta]'; 
allChCorrDeepCSD(:,2) = [tempProfileDeepCSDRef.magLowAlpha]'; 
allChCorrDeepCSD(:,3) = [tempProfileDeepCSDRef.magLowBeta]'; 
allChCorrDeepCSD(:,4) = [tempProfileDeepCSDRef.magLow]'; 
allChCorrDeepCSD(:,5) = [tempProfileDeepCSDRef.magLowRaw]'; 
allChCorrDeepCSD(~goodRuns,:) = [];


% Show the distributions of lags and correlations for 10/10 and 6/6 split
% for superficial and deep channels
superCorr_Bip_Gamma = [tempProfileSuperBipRef.magLow]'; 
deepCorr_Bip_Gamma  = [tempProfileDeepBipRef.magLow]';
superCorr_CSD_Gamma  = [tempProfileSuperCSDRef.magLow]'; 
deepCorr_CSD_Gamma   = [tempProfileDeepCSDRef.magLow]'; 

superCorr_Bip_Gamma(~(goodRuns & ~singleChFlag)) = []; 
deepCorr_Bip_Gamma(~(goodRuns & ~singleChFlag))  = []; 
superCorr_CSD_Gamma(~(goodRuns & ~singleChFlag))  = []; 
deepCorr_CSD_Gamma(~(goodRuns & ~singleChFlag))  = [];  

superLag_Bip_Gamma = [tempProfileSuperBipRef.lagLow]'./10; 
deepLag_Bip_Gamma  = [tempProfileDeepBipRef.lagLow]'./10;
superLag_CSD_Gamma  = [tempProfileSuperCSDRef.lagLow]'./10; 
deepLag_CSD_Gamma   = [tempProfileDeepCSDRef.lagLow]'./10;

superLag_Bip_Gamma(~(goodRuns & ~singleChFlag)) = []; 
deepLag_Bip_Gamma(~(goodRuns & ~singleChFlag))  = []; 
superLag_CSD_Gamma(~(goodRuns & ~singleChFlag))  = []; 
deepLag_CSD_Gamma(~(goodRuns & ~singleChFlag))   = [];  

% Get the median +/- mad lag 
x = -200:200;
allChLag = [tempProfileBipRef.lagLow]'; allChLag(~(goodRuns & ~singleChFlag)) = []; 
lagFrameRange = find(x == round(median(allChLag)- mad(allChLag))) : find(x == round(median(allChLag)+ mad(allChLag)));
%169:189; % Median +/- MAD across two monkeys

smFlagFOV = smFlag; smFlagFOV(~(goodRunsSpatial & ~singleChFlag)) = [];
smFlagROI = smFlag; smFlagROI(~(goodRuns)) = [];

if ~exist(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars_v2.mat'],'file') 
    save(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars_v2.mat'],...
        'smFlagFOV','tempProfilesAll','allChCorr','allChLagVal','superCorr_Bip_Gamma','deepCorr_Bip_Gamma',...
        'superCorr_CSD_Gamma','deepCorr_CSD_Gamma','superLag_Bip_Gamma','deepLag_Bip_Gamma',...
        'superLag_CSD_Gamma','deepLag_CSD_Gamma','lagFrameRange','allChLag','allChCorrSuper',...
        'allChCorrMid','allChCorrDeep','smFlagROI', 'tempProfilesAllCSD','allChCorrCSD','allChLagValCSD',...
        'allChCorrSuperCSD','allChCorrMidCSD','allChCorrDeepCSD','-append');
end

%% Compiling and plotting
allMonkeyVars(2) =  load(['X:\Data\Whiskey_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars_v2.mat']);
allMonkeyVars(1) =  load(['X:\Data\CharlieSheen_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars_v2.mat']);
monkeys = {'Charlie Sheen'; 'Whiskey';'Combined'};

for iM = 1:3
    clear super mid deep superCSD midCSD deepCSD
    if iM~=3
        super      = allMonkeyVars(iM).allChCorrSuper;
        mid        = allMonkeyVars(iM).allChCorrMid;
        deep       = allMonkeyVars(iM).allChCorrDeep;

        superCSD = allMonkeyVars(iM).allChCorrSuperCSD;
        midCSD   = allMonkeyVars(iM).allChCorrMidCSD;
        deepCSD  = allMonkeyVars(iM).allChCorrDeepCSD;

        smFlagTemp = allMonkeyVars(iM).smFlagROI;
        monkeyName = monkeys{iM};
    else
        super      = [allMonkeyVars(1).allChCorrSuper; allMonkeyVars(2).allChCorrSuper];
        mid        = [allMonkeyVars(1).allChCorrMid;   allMonkeyVars(2).allChCorrMid];
        deep       = [allMonkeyVars(1).allChCorrDeep;  allMonkeyVars(2).allChCorrDeep];

        superCSD      = [allMonkeyVars(1).allChCorrSuperCSD; allMonkeyVars(2).allChCorrSuperCSD];
        midCSD        = [allMonkeyVars(1).allChCorrMidCSD;   allMonkeyVars(2).allChCorrMidCSD];
        deepCSD       = [allMonkeyVars(1).allChCorrDeepCSD;  allMonkeyVars(2).allChCorrDeepCSD];

        smFlagTemp = [allMonkeyVars(1).smFlagROI;      allMonkeyVars(2).smFlagROI];
        monkeyName   = 'Combined data'; 
    end

    nanRows = find(isnan(mid(:,1)));
    super(nanRows,:) = [];
    mid(nanRows,:) = [];
    deep(nanRows,:) = [];

    superCSD(nanRows,:) = [];
    midCSD(nanRows,:) = [];
    deepCSD(nanRows,:) = [];

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
        caxis([-0.35 0]); title(typeLabel); sgtitle([monkeyName ': Bipolar']); axis square;

    end

    figure;
    for iType = 1:3
        switch iType
            case 1
                superTempCSD = superCSD(smFlagTemp =='S',:);
                midTempCSD   = midCSD(smFlagTemp =='S',:);
                deepTempCSD  = deepCSD(smFlagTemp=='S',:);

                typeLabel = 'Sensory';
            case 2

                superTempCSD = superCSD(smFlagTemp =='M',:);
                midTempCSD   = midCSD(smFlagTemp =='M',:);
                deepTempCSD  = deepCSD(smFlagTemp=='M',:);

                typeLabel = 'Motor';
            case 3

                superTempCSD = superCSD;
                midTempCSD   = midCSD;
                deepTempCSD  = deepCSD;
                typeLabel = 'All sites';
        end


        figPlotCSD = [median(superTempCSD,1,'omitnan');median(midTempCSD,1,'omitnan' ); median(deepTempCSD,1,'omitnan')];
        subplot(1,3,iType); imagesc(figPlotCSD); xticks(1:5); xticklabels(bandLabels); yticks(1:3);
        yticklabels({'Superficial';'Middle';'Deep'}); colorbar;colormap(flipud(jet));
        caxis([-0.35 0]); title(typeLabel); sgtitle([monkeyName ': CSD']); axis square;
    end 

end

%% Show correlation with FC map for superficial/middle/deep layers
% Boxplots of frequencies for each layer
clear super mid deep
% iType =2;
for iM = 1:3
    clear super mid deep
    if iM~=3
        super      = allMonkeyVars(iM).allChCorrSuper;
        mid        = allMonkeyVars(iM).allChCorrMid;
        deep       = allMonkeyVars(iM).allChCorrDeep;

        superCSD = allMonkeyVars(iM).allChCorrSuperCSD;
        midCSD   = allMonkeyVars(iM).allChCorrMidCSD;
        deepCSD  = allMonkeyVars(iM).allChCorrDeepCSD;
        
        smFlagTemp = allMonkeyVars(iM).smFlagFOV;
        monkeyName   = monkeys{iM};

    else
        super      = [allMonkeyVars(1).allChCorrSuper; allMonkeyVars(2).allChCorrSuper];
        mid        = [allMonkeyVars(1).allChCorrMid; allMonkeyVars(2).allChCorrMid];
        deep       = [allMonkeyVars(1).allChCorrDeep; allMonkeyVars(2).allChCorrDeep];

        superCSD      = [allMonkeyVars(1).allChCorrSuperCSD; allMonkeyVars(2).allChCorrSuperCSD];
        midCSD        = [allMonkeyVars(1).allChCorrMidCSD;   allMonkeyVars(2).allChCorrMidCSD];
        deepCSD       = [allMonkeyVars(1).allChCorrDeepCSD;  allMonkeyVars(2).allChCorrDeepCSD];

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
        subplot(131); boxplot(superTemp,bandLabels); ylim([-0.5 0.5]); box off; title('Superficial');
        subplot(132); boxplot(midTemp,bandLabels); ylim([-0.5 0.5]); box off; title('Middle');
        subplot(133); boxplot(deepTemp,bandLabels); ylim([-0.5 0.5]); box off; title('Deep');
        sgtitle([monkeyName ' Bipolar: ' typeLabel]);
    end

    

    for iType = 1:3
        switch iType
            case 1
                superTempCSD = superCSD(smFlagTemp =='S',:);
                midTempCSD   = midCSD(smFlagTemp =='S',:);
                deepTempCSD  = deepCSD(smFlagTemp=='S',:);

                typeLabel = 'Sensory';
            case 2

                superTempCSD = superCSD(smFlagTemp =='M',:);
                midTempCSD   = midCSD(smFlagTemp =='M',:);
                deepTempCSD  = deepCSD(smFlagTemp=='M',:);

                typeLabel = 'Motor';
            case 3

                superTempCSD = superCSD;
                midTempCSD   = midCSD;
                deepTempCSD  = deepCSD;
                typeLabel = 'All sites';
        end

        figure;
        subplot(131); boxplot(superTempCSD,bandLabels); ylim([-0.5 0.5]); box off; title('Superficial');
        subplot(132); boxplot(midTempCSD,bandLabels); ylim([-0.5 0.5]); box off; title('Middle');
        subplot(133); boxplot(deepTempCSD,bandLabels); ylim([-0.5 0.5]); box off; title('Deep');
        sgtitle([monkeyName ' CSD: ' typeLabel]);

    end
end

%% Get re-referenced data for FOV
for iDate = 1:size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);
    for iRun = 1: size(allRuns{iDate,1},1)
        tic;
        clear runName dataDir clipMask elecMask clipMaskCortex corrMask...
            x negIdx lowIdx serverDir
        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName];
        serverDir = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\kem294\Data\' ...
            monkeyName '_SqM\' hemisphere ' Hemisphere\'  expDate '\' runName ];

        clc; disp(['Analyzing data for ' monkeyName ' '  expDate ' run: ' runName]);

        if ~exist([serverDir '\crossCorrFOV_6_BipolarRef_BL.mat'],'file')
            disp('Bipolar_ref_BL');
            getCrossCorrFOV_v2(monkeyName,expDate,runName,serverDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,rawCh{iRun,iDate}.rawChTemp,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,6,'BipolarRef_BL');
        end

        if ~exist([serverDir '\crossCorrFOV_6_CSDRef_BL.mat'],'file') 
             disp('CSD_Ref_BL');
            getCrossCorrFOV_v2(monkeyName,expDate,runName,serverDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,rawCh{iRun,iDate}.rawChTemp,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp,6,'CSDRef_BL');
        end


        if ~exist([dataDir '\refProcessedDat.mat'],'file')
           
            disp(['Compiling data for ' monkeyName ' '  expDate ' ' runName]);
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
            corrFCHybridT = NaN(5,2,401); corrFCSuperT    = NaN(5,2,401);
            corrFCDeepT   = NaN(5,2,401); super_DeepCorrT = NaN(5,2,401);

            corrFCMidT          = NaN(5,2,401); super_MidCorrT       = NaN(5,2,401);
            mid_DeepCorrT       = NaN(5,2,401); peakNegValsAllT      = NaN(5,2);
            peakNegTimesAllT    = NaN(5,2);     super_DeepAvgFramesT = NaN(5,2);
            super_MidAvgFramesT = NaN(5,2);     deep_MidAvgFramesT   = NaN(5,2);

            % Correlations between frequencies
            superHybridAllBandsT = NaN(5,5,5);  deepHybridAllBandsT = NaN(5,5,5);
            midHybridAllBandsT   = NaN(5,5,5); crossFreqCrossLayerHybridT = NaN(5,5,5,3,3);

            for iType = 1:2 % 10/10 split or 6/6 split of superficial/deep channels
                clear crossCorrFOV allXCorr superXCorr deepXCorr allLags fileName
                switch iType
                    case 1
                        crossCorrFOV = matfile([serverDir '\crossCorrFOV_6_BipolarRef_BL.mat']);
                        fileName     = '6_BipolarRef_BL';

                    case 2
                        crossCorrFOV = matfile([serverDir '\crossCorrFOV_6_CSDRef_BL.mat']);
                        fileName     = '6_CSDRef_BL';

                end

                allXcorr     = crossCorrFOV.spatialProfile;
                superXcorr   = crossCorrFOV.spatialProfileSuper;
                deepXcorr    = crossCorrFOV.spatialProfileDeep;

                if chLen~=0
                    midXCorr     = crossCorrFOV.spatialProfileMid;
                end
                allLags      = crossCorrFOV.lagFull;

                x = allLags;
                negIdx = x<0 & x>=-150; negVals = x(negIdx);
                lowIdx = x<0 & x>= -80; xLow = x(lowIdx);

                % Show and save the full field of view for the hybrid map
                % Within frequency comparison

                % Theta
                crossCorrTheta      = reshape(allXcorr.ccFullTheta,[401 imSize(1)*imSize(2)]);   crossCorrTheta(:,~corrMaskT)      = NaN;
                crossCorrSuperTheta = reshape(superXcorr.ccFullTheta,[401 imSize(1)*imSize(2)]); crossCorrSuperTheta(:,~corrMaskT) = NaN;
                crossCorrDeepTheta  = reshape(deepXcorr.ccFullTheta,[401 imSize(1)*imSize(2)]);  crossCorrDeepTheta(:,~corrMaskT)  = NaN;
                lagLowTheta         = tempProfileBipRef(iDate,iRun).lagLowTheta;

                % Alpha
                crossCorrAlpha      = reshape(allXcorr.ccFullAlpha,[401 imSize(1)*imSize(2)]);   crossCorrAlpha(:,~corrMaskT)      = NaN;
                crossCorrSuperAlpha = reshape(superXcorr.ccFullAlpha,[401 imSize(1)*imSize(2)]); crossCorrSuperAlpha(:,~corrMaskT) = NaN;
                crossCorrDeepAlpha  = reshape(deepXcorr.ccFullAlpha,[401 imSize(1)*imSize(2)]);  crossCorrDeepAlpha(:,~corrMaskT)  = NaN;
                lagLowAlpha         = tempProfileBipRef(iDate,iRun).lagLowAlpha;

                % Beta
                crossCorrBeta      = reshape(allXcorr.ccFullBeta,[401 imSize(1)*imSize(2)]);   crossCorrBeta(:,~corrMaskT)      = NaN;
                crossCorrSuperBeta = reshape(superXcorr.ccFullBeta,[401 imSize(1)*imSize(2)]); crossCorrSuperBeta(:,~corrMaskT) = NaN;
                crossCorrDeepBeta  = reshape(deepXcorr.ccFullBeta,[401 imSize(1)*imSize(2)]);  crossCorrDeepBeta(:,~corrMaskT)  = NaN;
                lagLowBeta         = tempProfileBipRef(iDate,iRun).lagLowBeta;

                % Gamma
                crossCorrGamma      = reshape(allXcorr.ccFull,[401 imSize(1)*imSize(2)]);   crossCorrGamma(:,~corrMaskT)      = NaN;
                crossCorrSuperGamma = reshape(superXcorr.ccFull,[401 imSize(1)*imSize(2)]); crossCorrSuperGamma(:,~corrMaskT) = NaN;
                crossCorrDeepGamma  = reshape(deepXcorr.ccFull,[401 imSize(1)*imSize(2)]);  crossCorrDeepGamma(:,~corrMaskT)  = NaN;
                lagLowGamma         = tempProfileBipRef(iDate,iRun).lagLow;

                % Spiking
                crossCorrSpiking      = reshape(allXcorr.ccFullRaw,[401 imSize(1)*imSize(2)]);   crossCorrSpiking(:,~corrMaskT)      = NaN;
                crossCorrSuperSpiking = reshape(superXcorr.ccFullRaw,[401 imSize(1)*imSize(2)]); crossCorrSuperSpiking(:,~corrMaskT) = NaN;
                crossCorrDeepSpiking  = reshape(deepXcorr.ccFullRaw,[401 imSize(1)*imSize(2)]);  crossCorrDeepSpiking(:,~corrMaskT)  = NaN;
                lagLowSpiking         = tempProfileBipRef(iDate,iRun).lagLowRaw;

                if  chLen~=0
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

                            if chLen~=0
                                crossCorrMidR   = crossCorrMidTheta;
                            end

                        case 2
                            bandName        = 'Alpha';
                            mapsAll         = crossCorrAlpha;
                            crossCorrSuperR = crossCorrSuperAlpha;
                            crossCorrDeepR  = crossCorrDeepAlpha;
                            lagLow          = lagLowAlpha;

                            if  chLen~=0
                                crossCorrMidR   = crossCorrMidAlpha;
                            end

                        case 3
                            bandName        = 'Beta';
                            mapsAll         = crossCorrBeta;
                            crossCorrSuperR = crossCorrSuperBeta;
                            crossCorrDeepR  = crossCorrDeepBeta;
                            lagLow          = lagLowBeta;

                            if  chLen~=0
                                crossCorrMidR   = crossCorrMidBeta;
                            end

                        case 4
                            bandName        = 'Gamma';
                            mapsAll         = crossCorrGamma;
                            crossCorrSuperR = crossCorrSuperGamma;
                            crossCorrDeepR  = crossCorrDeepGamma;
                            lagLow          = lagLowGamma;

                            if chLen~=0
                                crossCorrMidR   = crossCorrMidGamma;
                            end

                        case 5
                            bandName        = 'Spiking';
                            mapsAll         = crossCorrSpiking;
                            crossCorrSuperR = crossCorrSuperSpiking;
                            crossCorrDeepR  = crossCorrDeepSpiking;
                            lagLow          = lagLowSpiking;

                            if chLen~=0
                                crossCorrMidR   = crossCorrMidSpiking;
                            end
                    end

                    for iMap = 1:size(crossCorrGamma,1)
                        corrFCHybridT(iBand,iType,iMap)   = corr(fcMap,mapsAll(iMap,:)','rows','complete');
                        corrFCSuperT(iBand,iType,iMap)    = corr(fcMap,crossCorrSuperR(iMap,:)','rows','complete');
                        corrFCDeepT(iBand,iType,iMap)     = corr(fcMap,crossCorrDeepR(iMap,:)','rows','complete');
                        super_DeepCorrT(iBand,iType,iMap) = corr(crossCorrSuperR(iMap,:)',crossCorrDeepR(iMap,:)','rows','complete');

                        if  chLen~=0
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

                    if chLen~=0
                        midHybridAvg = mean(crossCorrMidR(lagFrameRange,:),1,'omitnan');
                        super_MidAvgFramesT(iBand,iType) = corr(superHybridAvg',midHybridAvg','rows','complete');
                        deep_MidAvgFramesT(iBand,iType)  = corr(deepHybridAvg',midHybridAvg','rows','complete');
                    end

                    if ~exist([dataDir '\HybridMapFOV_' fileName '_' bandName '.png'],'file') || 1
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

                    if ~exist([dataDir '\HybridMapFOV_Super_' fileName '_' bandName '.png'],'file')  || 1
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

                    if ~exist([dataDir '\HybridMapFOV_Deep_' fileName '_' bandName '.png'],'file')  || 1
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

                    if ((~exist([dataDir '\HybridMapFOV_Mid_' fileName '_' bandName '.png'],'file'))|| 1) && chLen~=0 
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
                            if  chLen~=0
                                mid1 = mean(crossCorrMidTheta(lagFrameRange,:),1,'omitnan');
                            end

                        case 2
                            super1 = mean(crossCorrSuperAlpha(lagFrameRange,:),1,'omitnan');
                            deep1  = mean(crossCorrDeepAlpha(lagFrameRange,:),1,'omitnan');
                            if chLen~=0
                                mid1 = mean(crossCorrMidAlpha(lagFrameRange,:),1,'omitnan');
                            end

                        case 3
                            super1 = mean(crossCorrSuperBeta(lagFrameRange,:),1,'omitnan');
                            deep1  = mean(crossCorrDeepBeta(lagFrameRange,:),1,'omitnan');
                            if  chLen~=0
                                mid1 = mean(crossCorrMidBeta(lagFrameRange,:),1,'omitnan');
                            end

                        case 4
                            super1 = mean(crossCorrSuperGamma(lagFrameRange,:),1,'omitnan');
                            deep1  = mean(crossCorrDeepGamma(lagFrameRange,:),1,'omitnan');
                            if chLen~=0
                                mid1 = mean(crossCorrMidGamma(lagFrameRange,:),1,'omitnan');
                            end

                        case 5
                            super1 = mean(crossCorrSuperSpiking(lagFrameRange,:),1,'omitnan');
                            deep1  = mean(crossCorrDeepSpiking(lagFrameRange,:),1,'omitnan');
                            if  chLen~=0
                                mid1 = mean(crossCorrMidSpiking(lagFrameRange,:),1,'omitnan');
                            end
                    end

                    for iBand2 = 1:5
                        clear super2 deep2 mid2
                        switch iBand2
                            case 1
                                super2 = mean(crossCorrSuperTheta(lagFrameRange,:),1,'omitnan');
                                deep2  = mean(crossCorrDeepTheta(lagFrameRange,:),1,'omitnan');
                                if  chLen~=0
                                    mid2 = mean(crossCorrMidTheta(lagFrameRange,:),1,'omitnan');
                                end

                            case 2
                                super2 = mean(crossCorrSuperAlpha(lagFrameRange,:),1,'omitnan');
                                deep2  = mean(crossCorrDeepAlpha(lagFrameRange,:),1,'omitnan');
                                if chLen~=0
                                    mid2 = mean(crossCorrMidAlpha(lagFrameRange,:),1,'omitnan');
                                end

                            case 3
                                super2 = mean(crossCorrSuperBeta(lagFrameRange,:),1,'omitnan');
                                deep2  = mean(crossCorrDeepBeta(lagFrameRange,:),1,'omitnan');
                                if chLen~=0
                                    mid2 = mean(crossCorrMidBeta(lagFrameRange,:),1,'omitnan');
                                end

                            case 4
                                super2 = mean(crossCorrSuperGamma(lagFrameRange,:),1,'omitnan');
                                deep2  = mean(crossCorrDeepGamma(lagFrameRange,:),1,'omitnan');
                                if  chLen~=0
                                    mid2 = mean(crossCorrMidGamma(lagFrameRange,:),1,'omitnan');
                                end

                            case 5
                                super2 = mean(crossCorrSuperSpiking(lagFrameRange,:),1,'omitnan');
                                deep2  = mean(crossCorrDeepSpiking(lagFrameRange,:),1,'omitnan');
                                if  chLen~=0
                                    mid2 = mean(crossCorrMidSpiking(lagFrameRange,:),1,'omitnan');
                                end
                        end

                        superHybridAllBandsT(iType,iBand1,iBand2) = corr(super1',super2','rows','complete');
                        deepHybridAllBandsT(iType,iBand1,iBand2)  = corr(deep1',deep2','rows','complete');
                        if  chLen~=0
                            midHybridAllBandsT(iType,iBand1,iBand2) = corr(mid1',mid2','rows','complete');
                        end

                        if iBand1~=iBand2
                            crossFreqCrossLayerHybridT(iType,iBand1,iBand2,1,1) = corr(super1',super2','rows','complete');
                            crossFreqCrossLayerHybridT(iType,iBand1,iBand2,1,3) = corr(super1',deep2','rows','complete');
                            crossFreqCrossLayerHybridT(iType,iBand1,iBand2,3,1) = corr(deep1',super2','rows','complete');
                            crossFreqCrossLayerHybridT(iType,iBand1,iBand2,3,3) = corr(deep1',deep2','rows','complete');

                            if  chLen~=0
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

            save([dataDir '\refProcessedDat.mat'],'corrFCHybridT','corrFCSuperT',...
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
            allHybridVars                    = matfile([dataDir '\refProcessedDat.mat']);
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

%% Save all FOV level variables

corrFCHybridT = reshape(corrFCHybrid,[size(corrFCHybrid,1)*size(corrFCHybrid,2) size(corrFCHybrid,3) size(corrFCHybrid,4) size(corrFCHybrid,5)]);
nanRow        = (corrFCHybridT(:,1,1,1)==0);
corrFCHybridT(nanRow,:,:,:) = [];
corrFCHybridT(~goodRunsSpatial,:,:,:) = [];

% Compile correlations between layer compartments and FC map
corrFCSuperT = mean(corrFCSuper(:,:,:,:,lagFrameRange),5,'omitnan'); 
corrFCSuperT = reshape(corrFCSuperT,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]); 
corrFCSuperT(nanRow,:,:) = []; 
corrFCSuperT(~(goodRunsSpatial & ~singleChFlag),:,:) = []; 

corrFCMidT = mean(corrFCMid(:,:,:,:,lagFrameRange),5,'omitnan'); 
corrFCMidT = reshape(corrFCMidT,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]); 
corrFCMidT(nanRow,:,:) = []; 
corrFCMidT(~(goodRunsSpatial & ~singleChFlag),:,:) = []; 

corrFCDeepT = mean(corrFCDeep(:,:,:,:,lagFrameRange),5,'omitnan'); 
corrFCDeepT = reshape(corrFCDeepT,[size(super_DeepCorr,1)*size(super_DeepCorr,2) size(super_DeepCorr,3) size(super_DeepCorr,4)]); 
corrFCDeepT(nanRow,:,:) = []; 
corrFCDeepT(~(goodRunsSpatial & ~singleChFlag),:,:) = []; 

% Compile the peak negative correlations and lags for good runs
peakNegValsAllT = reshape(peakNegValsAll,[size(peakNegValsAll,1)*size(peakNegValsAll,2) size(peakNegValsAll,3) size(peakNegValsAll,4)]);
peakNegValsAllT(nanRow,:,:) = [];
peakNegValsAllT(~goodRunsSpatial,:,:) = [];

% Compile the across layers hybrid map correlations ie; S/D, S/M and M/D 
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

% clear varInfo;
% varFlag = 0;
% try matfile(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars_v2.mat']);
%     varInfo = who('-file',['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars_v2.mat']);
%     if sum(ismember(varInfo,'corrFCHybridT'))==0 || sum(ismember(varInfo,'superHybridAllBandsT'))==0
%         varFlag = 1;
%     end
% catch
%     varFlag = 1;
% end
        
if 1%varFlag 
    save(['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars_v2.mat'],...
       'corrFCHybridT','corrFCSuperT','corrFCMidT','corrFCDeepT','peakNegValsAllT',...
       'super_DeepAvgFramesT','super_MidAvgFramesT','deep_MidAvgFramesT',...
       'superHybridAllBandsT','midHybridAllBandsT','deepHybridAllBandsT');
end

%% Plot
bandLabels = {'Theta';'Alpha';'Beta';'Gamma';'Spiking'}; 

for iType = 2%:2
    clear super mid deep
    super = squeeze(corrFCSuperT(:,:,iType));
    mid   = squeeze(corrFCMidT(:,:,iType));
    deep  = squeeze(corrFCDeepT(:,:,iType));

    figure; 
    subplot(131); boxplot(super,bandLabels); ylim([-1 1]); box off; title('Superficial');
    subplot(132); boxplot(mid,bandLabels); ylim([-1 1]); box off; title('Middle');
    subplot(133); boxplot(deep,bandLabels); ylim([-1 1]); box off; title('Deep');

end


%% Plot the data
allMonkeyVars(2) =  load(['X:\Data\Whiskey_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars_v2.mat']);
allMonkeyVars(1) =  load(['X:\Data\CharlieSheen_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars_v2.mat']);
monkeys = {'Charlie Sheen'; 'Whiskey';'Combined'};

% Get gamma Bipolar and CSD data for superficial, middle and deep layers
iType = 1; iBand = 4;
superBip  = [allMonkeyVars(1).corrFCSuperT(:,iBand,iType); allMonkeyVars(2).corrFCSuperT(:,iBand,iType)];
midBip    = [allMonkeyVars(1).corrFCMidT(:,iBand,iType); allMonkeyVars(2).corrFCMidT(:,iBand,iType)];
deepBip   = [allMonkeyVars(1).corrFCDeepT(:,iBand,iType); allMonkeyVars(2).corrFCDeepT(:,iBand,iType)];

iType = 2;
superCSD = [allMonkeyVars(1).corrFCSuperT(:,iBand,iType); allMonkeyVars(2).corrFCSuperT(:,iBand,iType)];
midCSD   = [allMonkeyVars(1).corrFCMidT(:,iBand,iType); allMonkeyVars(2).corrFCMidT(:,iBand,iType)];
deepCSD  = [allMonkeyVars(1).corrFCDeepT(:,iBand,iType); allMonkeyVars(2).corrFCDeepT(:,iBand,iType)];

iType = 2;
super      = [allMonkeyVars(1).corrFCSuperT(:,iBand,iType); allMonkeyVars(2).corrFCSuperT(:,iBand,iType)];
mid        = [allMonkeyVars(1).corrFCMidT(:,iBand,iType); allMonkeyVars(2).corrFCMidT(:,iBand,iType)];
deep       = [allMonkeyVars(1).corrFCDeepT(:,iBand,iType); allMonkeyVars(2).corrFCDeepT(:,iBand,iType)];

% Show gamma powers for the three referencing schemes
bip = [superBip; midBip; deepBip];
csd = [superCSD; midCSD; deepCSD];
noRef = [super; mid; deep];

figure;boxplot([noRef bip csd],{'No reference';'Bipolar';'CSD'});
box off; ylim([-1 1]); title('Gamma powers between referencing schemes');

% Show different compartments for the referencing schemes 
figure; 
subplot(131); boxplot([super superBip superCSD],{'No reference';'Bipolar';'CSD'});
box off; ylim([-1 1]); title('Superficial');

subplot(132); boxplot([mid midBip midCSD],{'No reference';'Bipolar';'CSD'});
box off; ylim([-1 1]); title('Middle');

subplot(133); boxplot([deep deepBip deepCSD],{'No reference';'Bipolar';'CSD'});
box off; ylim([-1 1]); title('Deep');

% Combine across references to show lack of laminar difference
figure; 
superAll = [super; superBip; superCSD];
midAll = [mid; midBip; midCSD];
deepAll = [deep; deepBip; deepCSD]; 
figure; boxplot([superAll midAll deepAll],{'Superficial';'Middle';'Deep'});
box off; ylim([-1 1]); title('Gamma correlations between compartments');

% Show compartments for different ref schemes
figure; 
subplot(131); boxplot([super mid deep],{'Superficial';'Middle';'Deep'});
ylim([-1 1]); box off; title('No reference');
subplot(132); boxplot([superBip midBip deepBip],{'Superficial';'Middle';'Deep'});
ylim([-1 1]); box off; title('Bipolar reference');
subplot(133); boxplot([superCSD midCSD deepCSD],{'Superficial';'Middle';'Deep'}); 
ylim([-1 1]); box off; title('CSD reference');

% Show correlation with FC map for FOV
clear fovProfile
iType = 2;
for iM = 1:3
clear roiCorr super mid deep
    if iM~=3
        fovProfile = squeeze(allMonkeyVars(iM).corrFCHybridT(:,:,iType,:));
        super = squeeze(allMonkeyVars(iM).corrFCSuperT(:,:,iType));
        mid   = squeeze(allMonkeyVars(iM).corrFCMidT(:,:,iType));
        deep  = squeeze(allMonkeyVars(iM).corrFCDeepT(:,:,iType));

    else
        fovProfile = [squeeze(allMonkeyVars(1).corrFCHybridT(:,:,iType,:)) ; squeeze(allMonkeyVars(2).corrFCHybridT(:,:,iType,:))];
        super      = [allMonkeyVars(1).corrFCSuperT(:,:,iType); allMonkeyVars(2).corrFCSuperT(:,:,iType)];
        mid        = [allMonkeyVars(1).corrFCMidT(:,:,iType); allMonkeyVars(2).corrFCMidT(:,:,iType)];
        deep       = [allMonkeyVars(1).corrFCDeepT(:,:,iType); allMonkeyVars(2).corrFCDeepT(:,:,iType)];
        
    end

    figure;
    subplot(131); boxplot(super,bandLabels); ylim([-1 1]); box off; title('Superficial');
    subplot(132); boxplot(mid,bandLabels); ylim([-1 1]); box off; title('Middle');
    subplot(133); boxplot(deep,bandLabels); ylim([-1 1]); box off; title('Deep');
    sgtitle(monkeys{iM});

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
     sgtitle(monkeys{iM});
end

%% Show distributions of peak negative correlations with FC map for FOV
clear fovCorr iType


iType = 3;
for iM = 1%:3
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
allMonkeyVarsOld(2) =  load(['X:\Data\Whiskey_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat']);
allMonkeyVarsOld(1) =  load(['X:\Data\CharlieSheen_SqM\' hemisphere ' Hemisphere\ISOI_Ephys_allVars.mat']);

iType =2;
for iM = 1:3
    clear super mid deep
    if iM~=3
        super      = allMonkeyVars(iM).corrFCSuperT(:,:,iType);
        mid        = allMonkeyVars(iM).corrFCMidT(:,:,iType);
        deep       = allMonkeyVars(iM).corrFCDeepT(:,:,iType);
        smFlagTemp = allMonkeyVarsOld(iM).smFlagFOV;
%         monkeyName   = monkeys{iM};
    else
        super      = [allMonkeyVars(1).corrFCSuperT(:,:,iType); allMonkeyVars(2).corrFCSuperT(:,:,iType)];
        mid        = [allMonkeyVars(1).corrFCMidT(:,:,iType); allMonkeyVars(2).corrFCMidT(:,:,iType)];
        deep       = [allMonkeyVars(1).corrFCDeepT(:,:,iType); allMonkeyVars(2).corrFCDeepT(:,:,iType)];
        smFlagTemp = [allMonkeyVarsOld(1).smFlagFOV; allMonkeyVarsOld(2).smFlagFOV];
%         monkeyName = 'Combined data';
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
        sgtitle([monkeys{iM} ' - ' typeLabel]);
    end
end


