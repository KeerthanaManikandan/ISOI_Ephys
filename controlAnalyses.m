%% Control analyses
% Script to perform spatial and temporal control analyses for one monkey
% Make sure cross-modal maps are saved before running this script
% June 13, 2025 - KM
% Set paths
clc; clear;
commonDir = 'C:\Users\kem294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\ISOI_Ephys\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12\fly_track\videoIO']));

%% Initialize variables and get monkey data
hemisphere = 'Left'; spatialBin = 3;
iM = 2; % 1 - Charlie Sheen, 2 - Whiskey

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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Spatial controls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correlate FC maps obtained from seeds placed at varying distances from 
% the eectrode with the cross-modal map with maximum correspondence to FC
% Initialize variables
x = -200:200;
negIdx = (-100<=x)&(x<=0); negVals = x(negIdx);

% Perform spatial controls for each recording
for iDate = 1:size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:);

    for iRun = 1: size(allRuns{iDate,1},1)
        clear runName dataDir clipMask elecMask clipMaskCortex corrMask lowIdx ...
            pDatTemp greenFig seedLocIn crossCorrTemp allHybridMaps mapsAll...
            peakNegHybridMap mapsAllTemp probeCh badTimes szLFP skullMask ...
            inDatSize infraEphys allCortexMask

        runName = allRuns{iDate,1}(iRun,:);
        dataDir = ['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
        serverDir = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\kem294\Data\' monkeyName '_SqM\' ...
            hemisphere ' Hemisphere\' expDate '\' runName];
       
        % IMAGING: Load the appropriate masks for the imaging data
        [clipMaskCortex, corrMask] = getMasks(dataDir,runName);

        % Get rs-ISOI data for the run
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

        fcMap             = plotCorrMap(seedSigT,pDatTemp,0); % Get FC map
        corrMaskT         = reshape(corrMask,[imSize(1)*imSize(2) 1]);
        fcMap             = reshape(fcMap,[361*438 1]);
        fcMap(~corrMaskT) = NaN;

        % Get the cross-modal maps and determine lag  for the run
        try
            allHybridMaps  = matfile([serverDir '\crossCorrFOV_6_NoRef.mat']);
        
        catch % Get maps if they aren't saved already
            getCrossCorrFOV(monkeyName,expDate,runName,serverDir,processedDat{iDate,iRun}.tempBandPass,...
                probe{iRun,iDate}.probeCh,probe{iRun,iDate}.rawCh,badTimesLFP{iDate,iRun},badTimeThresh{iDate,iRun},...
                badCh{iDate,iRun},estChInCortex{1,iDate}(iRun,:),probe{iRun,iDate}.timeStamp);

            allHybridMaps = matfile([serverDir '\crossCorrFOV_6_NoRef.mat']);
        end

        mapsAllTemp           = allHybridMaps.spatialProfile; % Get cross-modal maps
        mapsAll               = mapsAllTemp.ccFull;
        mapsAll               = reshape(mapsAll,[401 361*438]);
        mapsAll(:,~corrMaskT) = NaN;

        disp('Obtaining spatial controls...');
        pDatTemp = reshape(pDatTemp,[imSize(1) imSize(2) imSize(3)]);
       
        tic;
        if ~exist([dataDir '\spatialControlVarsFOV.mat'],'file') 
            clear  lagLow  hybridLowMap runWiseSpatialCorr runWiseSpatialTimes locAll corrHybridMap
            % Get 40 seeds at 0.5, 1, 2, 3,4 mm away from the electrode each
   
            minShift  = round(roiSize{iDate}(iRun)/spatialBin);
            theta     = linspace(0,2*pi, round(pi*minShift)); % number of angles
            maxPoints = length(theta);
            distShift = {'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'};

            for iShift = 1:5 
                clear loc locShift row col numPoints pixelLoc
                if iShift == 1
                    locShift = round(roiSize{iDate}(iRun)/spatialBin); % Get the distance from the electrode (radius)
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
                    
                    % Check if any of the seed locations are on a vessel or
                    % suli or out of the image. 
                    if any((fliplr(loc(iPoint,:))+seedRad)>size(greenFig)) || any(loc(iPoint,:)-seedRad<= 0)
                        runWiseSpatialCorr(iShift,iPoint)  = NaN;
                        runWiseSpatialTimes(iShift,iPoint) = NaN;
                        loc(iPoint,:) = NaN;
                    
                    else % Get FC maps
                        seed = loc(iPoint,:);
                        clipMask_seed = ~corrMask(seed(2)-seedRad:seed(2)+seedRad,seed(1)-seedRad:seed(1)+seedRad);
                        allCortex_seed = ~allCortexMask(seed(2)-seedRad:seed(2)+seedRad,seed(1)-seedRad:seed(1)+seedRad);
                        maskSize = size(clipMask_seed,1) *size(clipMask_seed,2);

                        % Check if the masks occupy more than 50% of seed
                        % or if the seeds are close to the edge (the edges
                        % of the image should at max occupy 25% of seed)
                        if ((sum(clipMask_seed,'all')/maskSize)>0.5) || ((sum(allCortex_seed,'all')/maskSize)>0.25)
                            runWiseSpatialCorr(iShift,iPoint)      = NaN;
                            runWiseSpatialTimes(iShift,iPoint)     = NaN;
                            loc(iPoint,:)                          = NaN;
                            fc_TestCorr{iRun,iDate}(iShift,iPoint) = NaN;
                            continue;

                        else
                            % Get Gaussian weighted seed signal
                            seedSigT = calculateSeedSignal(greenFig,clipMaskCortex,loc(iPoint,:),seedRad,pDatTemp);
                            corrMapT = reshape(plotCorrMap(seedSigT,pDatTemp,0),[imSize(1)*imSize(2) 1]);

                            % Correlate FC map at reference to test maps
                            fc_TestCorr{iRun,iDate}(iShift,iPoint) = corr(corrMapT, fcMap,'rows','complete'); 

                            % Correlate FC map with all hybrid maps. This
                            % correlation makes sure that we remain agnostic
                            % about the peak negative lag and this gives a
                            % distribution of lags and/or correlations. The
                            % variables that change are lags, distance between
                            % probe and seed locations.
                            for iMap = 1:401; mapVals(iMap) = corr(corrMapT,mapsAll(iMap,:)','rows','complete'); end

                            [runWiseSpatialCorr(iShift,iPoint),runWiseSpatialTimes(iShift,iPoint)] = min(mapVals(negIdx));
                            runWiseSpatialTimes(iShift,iPoint) = negVals(runWiseSpatialTimes(iShift,iPoint))./10;
                            

                        end
                    end
                    locAll{iShift} = loc;
                end
            end
            % Save the location and the correlation between cross-modal and
            % FC maps
            save([dataDir '\spatialControlVarsFOV.mat'],'runWiseSpatialCorr','runWiseSpatialTimes','locAll','corrHybridMap');
            spCorrControl{iRun,iDate}      = runWiseSpatialCorr;
            spCorrControlTimes{iRun,iDate} = runWiseSpatialTimes;

            toc;
        else % Load all correlations 
            clear allSpatialVars
            allSpatialVars = load([dataDir '\spatialControlVarsFOV.mat']);
            spCorrControl{iRun,iDate}      = allSpatialVars.runWiseSpatialCorr;
            spCorrControlTimes{iRun,iDate} = allSpatialVars.runWiseSpatialTimes;
            locAll                         = allSpatialVars.locAll;
        end

        % Show the seeds, peak correlations
        if ~exist([dataDir '\spatialControlFOV_v2.png'],'file')
            cVals = {'w','k','b','g','m'};

            % Show the seeds on the blood vessel map
            figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(131); imagesc(greenFig); hold on; colormap gray; axis image off;
            plot(seedLocProbe(1),seedLocProbe(2),'Marker','pentagram','MarkerSize',15,...
                'MarkerFaceColor','r','MarkerEdgeColor','none');

            for iShift = 1:5 % Show the seeds
                if ~isempty(locAll{iShift})
                    plot(locAll{iShift}(:,1),locAll{iShift}(:,2),'.','Color',cVals{iShift},'MarkerSize',10);
                end
            end
            
            % Show the peak correlations vs distance
            subplot(132);boxplot(spCorrControl{iRun,iDate}',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});
            xlabel('Distance from probe (mm)'); ylabel('Correlation between FC map and peak negative map');
            ylim([-1 0.3]);
            
            % Show the lag values where peak correlations were observed
            subplot(133);boxplot(spCorrControlTimes{iRun,iDate}',{'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});
            xlabel('Distance from probe (mm)'); ylabel('Lag at peak negative correlation'); ylim([-11 3]);

            sgtitle(strrep(['FC maps vs Hybrid maps (varying lag) for ' monkeyName ' ' expDate ' ' runName],'_','\_'));
            f = gcf; exportgraphics(f,[dataDir '\spatialControlFOV_v3.png'],'Resolution',300); close gcf;
        end
    end
end

% Save correlations between reference FC map to test map
if ~exist(['D:\Data\' monkeyName '_SqM\Left Hemisphere\fcTest_RefCorr.mat'],'file')
    save(['D:\Data\' monkeyName '_SqM\Left Hemisphere\fcTest_RefCorr.mat'],'fc_TestCorr');
else
    fc_TestCorr = load(['D:\Data\' monkeyName '_SqM\Left Hemisphere\fcTest_RefCorr.mat'],'fc_TestCorr');
end

% Compiling spatial control data for all recordings
clear spCorrNorm medSPCorrControl medSpCorrAll spCorrMinT

% Normalize correlations within a recording relative to the minima and pool
% all correlations (variance comes from the seeds placed)
spCorrMinT = reshape(cellfun(@(x) x./abs(min(x,[],'all','omitnan')),spCorrControl,'un',0),[size(spCorrControl,1)*size(spCorrControl,2) 1]);
zeroInd    = cell2mat(cellfun(@(x) isempty(x),spCorrMinT,'un',0));

spCorrMinT(zeroInd)          = [];
spCorrMinT(~goodRunsSpatial) = []; % Remove bad recordings

spCorrMinT = -(cat(2,spCorrMinT{:})); % Pool the data

% Show the distributions
figure; violinplot(spCorrMinT);
ylim([-1.2 1]);yticks(-1:0.2:1);xticklabels([0.5 1 2 3 4]); hold on; 

% Show the data points
pointSize = size(spCorrMinT);
x1 = (reshape(repmat(1:5,[pointSize(2) 1]),[pointSize(2)*5 1]));
y1 = reshape(spCorrMinT',[pointSize(2)*5 1]);
s = swarmchart(x1,y1,5,'b','filled');
s.XJitterWidth = 0.5; box off;
xticks(1:5);xticklabels({'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});

% Show how FC at test sites relate to the reference
clear spCorrMinFC zeroInd
spCorrMinFC = cellfun(@(x) x./abs(max(x,[],'all','omitnan')),fc_TestCorr,'un',0); % Normalize relative to the minima
spCorrMinFC = reshape(spCorrMinFC,[size(fc_TestCorr,1)*size(fc_TestCorr,2) 1]);
zeroInd     = cell2mat(cellfun(@(x) isempty(x),spCorrMinFC,'un',0));

spCorrMinFC(zeroInd)          = [];
spCorrMinFC(~goodRunsSpatial) = []; 
spCorrMinFC                   = (cat(2,spCorrMinFC{:}));

% Show the distributions
figure; violinplot(spCorrMinFC);
ylim([-1.2 1]);yticks(-1:0.2:1);xticklabels([0.5 1 2 3 4]); hold on; 

% Show the data points
pointSize = size(spCorrMinFC);
x1 = (reshape(repmat(1:5,[pointSize(2) 1]),[pointSize(2)*5 1]));
y1 = reshape(spCorrMinFC',[pointSize(2)*5 1]);
s = swarmchart(x1,y1,5,'b','filled');
s.XJitterWidth = 0.5; box off;
xticks(1:5);xticklabels({'0.5 mm' ; '1 mm'; '2 mm' ; '3 mm'; '4 mm'});

%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Temporal controls %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shuffle Ephys using variable windows and determine the peak correlation 
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

        dataDir = ['D:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
        serverDir = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\kem294\Data\' monkeyName '_SqM\' ...
            hemisphere ' Hemisphere\' expDate '\' runName];
       
        % IMAGING: Load the appropriate masks for the imaging data
        % Load clipmask
        [clipMaskCortex, corrMask] = getMasks(dataDir,runName);

        if ~exist([dataDir '\phaseControlVarsFOV.mat'],'file')
            % Get rs-ISOI data for the run
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

            fcMap             = plotCorrMap(seedSigT,pDatTemp,0); % Get FC map
            corrMaskT         = reshape(corrMask,[imSize(1)*imSize(2) 1]);
            fcMap             = reshape(fcMap,[361*438 1]);
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

            % Get LFP data
            probeCh                = probe{iRun,iDate}.probeCh;
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
            
            % Get timestams
            timeStamp       = probe{iRun,iDate}.timeStamp;
            timeStampSorted = timeStamp- timeStamp(1);
            badTimes10Hz    = unique(badTimeThreshTemp./1000);
            badTimeIm       = [];

            % Identifying frames to be removed from rs-ISOI
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

            % Set the time windows to shuffle LFP
            tic;
            timeLen = min([size(gammaEphys,1) size(processedDat10,2)*1e2]);
            winLen  = [0.001 0.01 0.1 1 5 10 50 100 500 timeLen./1e3].*1e3;
            
            % Calculate phase of the signal
            gammaFFT = fft(gammaEphys); % Fourier transform
            magGamma = abs(gammaFFT); % Magnitude
            phaseGamma = angle(gammaFFT); % Phase
            
            % Shuffle phases in different window length for 10 times each
            for iWin = 1: length(winLen)
                for iRep = 1:10
                    repFlag = 1;
                    disp(['Window length: ' num2str(winLen(iWin)) ' Rep: ' num2str(iRep)]);
                    rng('shuffle');
                    comb1 = randperm(round(timeLen/winLen(iWin)));
                    clear newPhase gammaNew gammaNewFFT magGammaNew

                    if iWin == 1 % Shuffle all samples
                        newPhase = phaseGamma(comb1,:);
                        magGammaNew = magGamma((1:size(newPhase,1)),:);
                        gammaNewFFT = magGammaNew.*newPhase;
                        gammaNew = real(ifft(gammaNewFFT));

                    elseif iWin == 10 % No shuffle
                        gammaNew = gammaEphys;
                        if iRep>1 
                            repFlag = 0; 
                        end 

                    else % Shuffle phase in windows
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
                        % Get infraslow powers
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
                        
                        % Cross-correlate rs-ISOI with infraslow power and
                        % get cross-modal maps
                        parfor iP = 1:size(processedDat10R,1)
                            if iP == 1
                                [ccFull(:,iP),lagFull(:,iP)]  = xcorr(infraEphysS',processedDat10R(iP,:),200,'normalized');
                            else
                                [ccFull(:,iP),~]      = xcorr(infraEphysS',processedDat10R(iP,:),200,'normalized');
                            end
                        end

                        ccFull(:,~corrMaskT) = NaN;
                        
                        % Correlate cross-modal maps to FC
                        for iMap = 1:401; mapVals(iMap) = corr(fcMap,ccFull(iMap,:)','rows','complete'); end

                        % Get peak correlations 
                        [runWiseCorrShuffled(iWin,iRep),runWiseLagShuffled(iWin,iRep)] = min(mapVals(negIdx));
                        runWiseLagShuffled(iWin,iRep) = negVals(runWiseLagShuffled(iWin,iRep))./10;
                    else
                        runWiseCorrShuffled(iWin,iRep) = runWiseCorrShuffled(iWin,1);
                        runWiseLagShuffled(iWin,iRep) = runWiseLagShuffled(iWin,1);
                    end
                end
            end

            toc;
            % Save the correlations
            save([dataDir '\phaseControlVarsFOV.mat'],'runWiseCorrShuffled','runWiseLagShuffled');

            runWiseCorrAllShuffled{iDate,iRun} = runWiseCorrShuffled;
            runWiseLagAllShuffled{iDate,iRun}  = runWiseLagShuffled;
            corrNegShuffle(iDate,iRun,:)       = median(runWiseCorrShuffled,2,'omitnan');%runWiseCorrShuffled;%
            corrNegTimes(iDate,iRun,:)         = median(runWiseLagShuffled,2,'omitnan');%runWiseLagShuffled;%
       
        else % Load the correlations
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


% Grouping temporal controls 
winLen                  = [0.001 0.01 0.1 1 5 10 50 100 500 900]; 
runWiseCorrAllShuffledT = reshape(runWiseCorrAllShuffled,[size(runWiseCorrAllShuffled,1)*size(runWiseCorrAllShuffled,2) 1]);
zeroInd                 = cell2mat(cellfun(@(x) isempty(x),runWiseCorrAllShuffledT,'un',0));

runWiseCorrAllShuffledT(zeroInd) = [];
runWiseCorrAllShuffledT(~goodRunsSpatial) = [];

% Normalize relative to the unshuffled data
runWiseCorrAllShuffledT = cellfun(@(x) x./(min(x(10,:))),runWiseCorrAllShuffledT,'un',0)'; 
allTimePoints           = (cat(2,runWiseCorrAllShuffledT{:}))';
runWiseCorrAllShuffledT = cell2mat(cellfun(@(x) median(x,2,'omitnan'),runWiseCorrAllShuffledT,'un',0));

% Plot the medians
runWiseCorrAllShuffledT = runWiseCorrAllShuffledT(2:end,:);
plot(median(smoothdata(runWiseCorrAllShuffledT,2,'movmean',2),2),'k','LineWidth',2);

% Calculate 95% confidence interval or mean+/- 2 SEM
stdAll = std(runWiseCorrAllShuffledT,[],2)/sqrt(size(runWiseCorrAllShuffledT,2)); % SEM
c95    = tinv([0.025 0.975],size(runWiseCorrAllShuffledT,2)-1);
y95    = bsxfun(@times,stdAll', c95(:));

meanAll   = mean(runWiseCorrAllShuffledT,2,'omitnan');
posValAll = meanAll+y95';

% Show the shaded area (median +/- 2 sem)
xVar = [1:size(posValAll,1) fliplr((1:size(posValAll,1)))];
patch(xVar,[posValAll(:,1)' fliplr(posValAll(:,2)')],[0.65 0.65 0.65],'FaceAlpha',0.3,'EdgeColor','none')
xticklabels(winLen);  hold on;  ylim([-0.5 1]); yticks(-1:0.1:1);box off; ylim([-0.5 1]);
xlabel('Length of window for shuffling (s)'); ylabel('Cross correlation'); grid off; hold on;
