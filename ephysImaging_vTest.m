clear;clc
commonDir = 'C:\Users\KEM294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\neuroshare']));
addpath(genpath([commonDir '\Codes\Ephys']));
addpath(genpath([commonDir '\Codes\Imaging']));
addpath(genpath([commonDir '\Codes\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\chronux_2_12\fly_track\videoIO']));


%% Initializing all the relevant variables
% Imaging parameters
monkeyName   = 'Whiskey';
spatialBin   = 3;
hemisphere   = 'Left';
saveFlag     = 0; 
checkBadRuns = 0; 
fs           = 1e3; 

% Ephys parameters
gammaBand      = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); % Gamma band filtering parameters
alphaBand      = [8 12];  [bA,aA] = butter(3,alphaBand./(fs/2),'bandpass'); % Alpha band filtering parameters
betaBand       = [13 30]; [bB,aB] = butter(3,betaBand./(fs/2),'bandpass');
probeBList     = [1:13 15:21 23:32];
chOutCortex    = 1:3;
chDeep         = 30:32;

% Load all the runs, experiment dates, green reference images, reference directory for the corresponding monkey
if strcmp(monkeyName,'CharlieSheen')
    allDates     = ['11_29_2021'; '01_11_2022'; '08_07_2023'];
    allRuns      = {['run00'; 'run01']; ['run03'; 'run04']; ['run01';'run02'; 'run03'; 'run04'; 'run05';'run06'; 'run07';'run08']};

    % Get imaging parameters
    refDate      = '08_31_2021';
    refDir       = [commonDir '\' monkeyName '_SqM\' hemisphere ' Hemisphere\' refDate '\Master Green Images\'];
    refImageName = 'Charlie Sheen Combined Green 08_31_2021';


    % Get ephys parameters
    ephysFileNameAll = {['run00'; 'run01'];['run0003';'run0004'];['datafile0001';'datafile0002';'datafile0003'; ...
        'datafile0003';'datafile0004';'datafile0005';'datafile0006' ;'datafile0007';'datafile0008']};
    serverPath       = '\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\302-19_CharlieSheen_SqM\Left Hemisphere\';
    chInCortex       = {[9 9];[10 10];[7 7 10 8 8 6 4 7]};
    anesthesiaLevels = {[0.75 0.75]' ;[0.75 0.75]';[0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9]'};
    heartRate        = {[230 210]'; [204 190]';[301 301 300 294 294 294 294 276 ]'};

elseif strcmp(monkeyName,'Whiskey')
    allDates                = '08_14_2023';% '09_19_2022'; '10_17_2022' ;'02_21_2023';'04_11_2023'];
    allRuns                 = {['run02'; 'run03'; 'run04']};
    datFileNumAll           = {2:4};
    serverPath              = '\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\303-19_Whiskey_SqM\Left Hemisphere\';
    refDate                 = '05_09_2022';
    refImageName            = 'Combined Green';
    refDir                  = [commonDir '\' monkeyName '_SqM\' hemisphere ' Hemisphere\' refDate '\Master Green Images\'];
    datFileNameAll          = {'datafile_000'};% 'datafile000'; 'datafile000'; 'datafile000'; 'datafile000'};
    eegChName               = 'analog1';
    ephysFileNameAll        = {['datafile0002'; 'datafile0003' ; 'datafile0004']};
%     chInCortexProbeA        = {[10 10 10 10]; [6 6 6 6 6 9 10 10 10 10 10]; [6 1 10 10 5 5 5 5 10 10 6 7 9 8 12]; [11 10 11 10 10 10 10 7 7 7 6 11 10 14 7 10 ]; [7 7 6 6 6 11 11 11 11 7 7 7 7 13 13]};
%     chInCortexProbeB        = {[1 1 1 1];     [9 9 10 10 6 9 8 8 11 9 9 ];  [8 10 10 10 6 6 6 6 6 6 6 6 6 6 6 ]; [10 4 4 4 4 12 11 10 10 10 10 6 6 6 6 6]; [10 10 10 8 6 7 11 8 10 10 10 10 10 5 12]};
%     chInCortexProbeAUpdated = {[10 13 8 8];[12 12 11 11 13 13 13 12 11 13 11];[4 1 9 10 8  8 8 8 10 1 4 10 9 14 8]; [17 17 18 16 13 17 16 12 12 11 12 15 16 19 15 14];[9 6 5 4 7 14 10 13 9 3 3 3 5 7 10 ]};
%     chInCortexProbeBUpdated = {[1 1 1 1]; [17 11 13 13 13 1 11 14 16 13 9];[16 13 7 16 8 8 8 8 14 13 5 13 13 14 10]; [17 11 12 12 18 18 13 16 16 16 16 12 11 12 12 12];[13 10 6 8 8 10 13 9 11 4 4 10 11 17 16]};
%     probeLabelA             = {[]; ['U592' ;'U592'; 'U592' ;'U592'; 'U592'; 'U592' ;'U592'; 'U592' ;'U592'; 'U592'; 'U592'];['D553'; 'D553'; 'D553'; 'D553';'D553'; 'D553'; 'D553'; 'D553'; 'D553'; 'D553'; 'D553'; 'D553'; 'D553'; 'D553'; 'D553']; ['BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'];['BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29'; 'BD29']};
%     probeLabelB             = {[]; ['D89D'; 'D89D'; 'D89D'; 'D89D'; '967C'; '967C';'967C'; '967C'; '967C'; '967C'; '967C'];['U592'; 'U592'; 'U592'; 'U592'; 'U592';'U592'; 'U592';'U592'; 'U592'; 'U592';'U592'; 'U592'; 'U592'; 'U592'; 'U592']; ['CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'];['CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1'; 'CDE1']};
%     anesthesiaLevels        = {[ 1.1 1.1 1.1 1.1]'; [1 1 1 1 1.25 1.25 1.25 1.1 1.25 1.25 1.25]'; [1.2 1.2 1.2 1.2 1.2 2.25 3.25 1 1 1.1 1.1 1.1 1.1 1.1 1.1]';[0.8 0.8 0.9 1 1 1 1 1 1 1 1 1 1 1 1 1]';[1.3 1.3 1.3 1.3 1.3 1.5 1.75 1.75 1.75 1.6 1.5 1.3 1.3 1.3 1.3]'};
%     heartRate               = {[NaN NaN NaN NaN]';[280 262 285 285 287 274 274 286 242 242 285 ]

end

%%

if exist([refDir refImageName '.bmp'],'file') == 0 % Make sure to get what the reference run is outside the function
    greenMapRef = imread([refDir refImageName '.png']);
else
    greenMapRef = imread([refDir refImageName '.bmp']);
end
greenMapRef  = greenMapRef(:,:,1);
greenMapSize = size(greenMapRef);
greenMap     = double(greenMapRef);
greenMap_RGB = ind2rgb(greenMap,gray(256));

%%
iDate = 1; iRun = 2;
expDate    = allDates(iDate,:);
runName    = allRuns{iDate,1}(iRun,:);
saveFolder = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\Electrophysiology\Processed Data'];
datFileNum = ephysFileNameAll{iDate,1};
fileNum    = str2double(datFileNum(iRun,end));

dataDir     = [commonDir '\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];
numFiles    = length(dir([dataDir '\Spatial Downsample SS3' '/*.mat'])); % loading the downsampled data only
datName     = 'Data_RS_10Hz_SS3_';
templateDir = [commonDir '\' monkeyName '_SqM\Left Hemisphere\'];

% Load the vessel mask
if exist([dataDir '\clipMask0' runName(end) '.BMP'],'file') == 0
    clipMask = imread([dataDir '\clipMask0' runName(end) '.png']);
else
    clipMask = imread([dataDir '\clipMask0' runName(end) '.bmp']);
end

if exist([dataDir '\skullMask.bmp'],'file') == 0
    allCortexMask = imread([dataDir '\skullMask.png']); % Has vessels in this mask
else
    allCortexMask = imread([dataDir '\skullMask.bmp']); % Has vessels in this mask
end

if exist([dataDir '\maskSkull0' runName(end) '.bmp'],'file') == 0
    elecMask = imread([dataDir '\maskSkull0' runName(end) '.png']); % Has vessels in this mask
else
    elecMask = imread([dataDir '\maskSkull0' runName(end) '.bmp']); % Has vessels in this mask
end


clipMask       = imresize(clipMask,1/3); % Resize mask
clipMask       = clipMask(:,:,1)>0; % Converting to 0s and 1s
allCortexMask  = imresize(allCortexMask,1/3); % Resizing cortex mask with vessels
allCortexMask  = allCortexMask(:,:,1)>0;
elecMask       = imresize(elecMask,1/3);
elecMask       = elecMask(:,:,1)>0; 
vesselMask     = ~clipMask & allCortexMask; % Vessel mask
tempVesselMask = double(vesselMask);

clipMaskCortex = clipMask & ~elecMask; 


if exist([templateDir  expDate '\' allRuns{iDate,1}(iRun,:) '\green0' allRuns{iDate,1}(iRun,5) '_Edited.bmp'],'file') == 0
    greenTemp = imread([dataDir '\green0' allRuns{iDate,1}(iRun,5) '_Edited.png']);
else
    greenTemp = imread([dataDir '\green0' allRuns{iDate,1}(iRun,5) '_Edited.bmp']);
end
greenIm{iDate,iRun}  = greenTemp(:,:,1);

% Store/Retrieve processed imaging data
if strcmp(expDate,'06_22_2021')
    serverDataPath    = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\302-19_CharlieSheen_SqM\Left Hemisphere\06_22_2021 Craniotomy\' runName];
else
    serverDataPath = [serverPath expDate '\' runName];
end

if ~exist([dataDir '\processedFrames.mat'],'file')
    disp(['Processing imaging data for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)]);
    [~,~,~,~,tempBandPass] = getPreProcessedDataRestingState(serverDataPath,dataDir,runName,numFiles,spatialBin,datName);
    save([dataDir '\processedFrames.mat'],'tempBandPass');
    processedDat(:,:,:,iDate,iRun) = tempBandPass;

else
    disp(['Loading imaging data for ' monkeyName ' ' expDate ' File: ' num2str(fileNum)]);
    clear processedDatTemp
    processedDatTemp = load([dataDir '\processedFrames.mat']);
    processedDatTemp = processedDatTemp.tempBandPass;
    processedDat(:,:,:,iDate,iRun) = processedDatTemp;

end
cd(commonDir);

datFileName = ephysFileNameAll{iDate,1}(iRun,:);
datFileName = datFileName(1:end-1);
% datName = 'X:\Data\Whiskey_SqM\Left Hemisphere\08_14_2023\Electrophysiology\run02';
datName = [serverPath expDate '\Electrophysiology\' datFileName num2str(fileNum)];

[nsResult,hFile] = ns_OpenFile(datName);
[nsResult2, fileInfo] = ns_GetFileInfo(hFile);

for iEntity = 1: fileInfo.EntityCount
    [~,entityInfo(iEntity,1)] = ns_GetEntityInfo(hFile,iEntity);
end

% Obtain camera frames
eventList   = find([entityInfo.EntityType] == 1);

% Grab the times where the camera frame occurred....
for iT = 1:entityInfo(eventList(2)).ItemCount
    [~, timeStamp(iT), ~, ~] = ns_GetEventData(hFile, eventList(2),iT);
end

%% Sort the entities whether they contain events, neural data, segment data
lfpList   = find([entityInfo.EntityType] == 2);
lfpLabel  = {entityInfo(lfpList).EntityLabel}; % Gets the label of all the channels in lfpList

% Remove the raw signal channels and get the indices of the LFP only
if strcmp(lfpLabel{1}(1:3),'lfp')
    lfpIdx = cell2mat(cellfun(@(x) strcmp(x(1:3),'lfp'),lfpLabel,'un',0));
elseif strcmp(lfpLabel{1}(end-2:end),'lfp')
    lfpIdx = cell2mat(cellfun(@(x) strcmp(x(end-2:end),'lfp'),lfpLabel,'un',0));
elseif strcmp(lfpLabel{1}(1:6),'analog')
    lfpIdx = cell2mat(cellfun(@(x) strcmp(x,'analog 2'),lfpLabel,'un',0));
end
lfpList(~lfpIdx) = [];
lfpLabel(~lfpIdx) = [];

if strcmp(lfpLabel{1},lfpLabel{2}) % To remove 30kHz sampled data
    lfpList(2) = [];
    lfpLabel(2) = [];
end 

if strcmp(lfpLabel{1}(1:3),'lfp')
    chNum = str2num(cell2mat(cellfun(@(x) x(end-1:end),lfpLabel','un',0))); %#ok<ST2NM>
elseif strcmp(lfpLabel{1}(end-2:end),'lfp')
    chNum = str2num(cell2mat(cellfun(@(x) x(2:3),lfpLabel','un',0))); %#ok<ST2NM>
elseif strcmp(lfpLabel{1},'analog 2')
    chNum = 1; 
end

if chNum>1
    [~,elecID] = sort(chNum);
else
    elecID = 1; 
end 

% Get the LFP for all the channels in sorted order
clear b a bS aS probeCh
[b,a]   = butter(3,[1 250]./(fs/2),'bandpass'); % Bandpass filtering parameters across 1-250 Hz
[bS,aS] = butter(3,[57 62]./(fs/2),'stop'); % Bandstop filtering between 57-62 Hz

disp('Getting the LFP for the electrode and filtering it ... ');
for iElec = 1:length(lfpLabel)
    clear elecEntityID lfpEntityID lfpCount
    lfpEntityID   = lfpList(elecID(iElec));
    lfpCount      = entityInfo(lfpEntityID).ItemCount;

    if ~exist('analogInfo','var')
        [~, analogInfo] = ns_GetAnalogInfo(hFile, lfpEntityID);
        fs  = analogInfo.SampleRate;
    end

    % Get LFP data
    [~, ~, probeCh(:,iElec)] = ns_GetAnalogData(hFile,lfpEntityID,1,lfpCount);
end
probeCh = filtfilt(b,a,probeCh); % Bandpass filtering across 1-250 Hz
probeCh = single(filtfilt(bS,aS,probeCh)); % Bandstop filtering between 57-62 Hz

%% Truncate LFP data to match the imaging time points. 
t1k = unique(floor(timeStamp.*1000));

probeCh = single(probeCh(t1k(1):t1k(end),:));
probeCh = probeCh(1:905000,:);

%% Remove noisy data from both LFP and imaging dataset
fs              = 1000;
params.Fs       = fs;
params.fpass    = [1 120];
params.pad      = -1;
params.tapers   = [3 5];
params.trialave = 0;

channels  = 1:size(probeCh,2);
[spec,timeValsSpec,freqValsSpec] = mtspecgramc(probeCh(:,channels),[5 2],params);
powTimeBin = squeeze(sum(10.*log10(abs(spec)),2));

if length(channels)~= 1
    badChVal = [];
    badElecThresh = (median(powTimeBin,2)+5*mad(powTimeBin,1,2));
    badChVal = [badChVal channels(sum((powTimeBin>badElecThresh),1) >= floor(0.75*size(powTimeBin,1)))];
    if ~isempty(badChVal); channels(ismember(channels,badChVal)) = []; end
end

badTimeInd = []; badTimes = [];
clear spec 
if length(channels)~= 1
    [spec,~,~]    = mtspecgramc(probeCh(:,channels(channels>=15)),[5 2],params);
    meanS         = mean(10.*log10(abs(spec)),3,'omitnan'); % Mean across channels 15-32

else
    [spec,~,~]    = mtspecgramc(probeCh,[5 2],params);
    meanS         = spec; 
end
powMeanS      = squeeze(sum(meanS,2));
badTimeThresh = (median(powMeanS,1)+4.5*mad(powMeanS,1,1));
badTimeIndOld = floor((timeValsSpec(powMeanS>badTimeThresh))*fs);

if ~isempty(badTimeIndOld)
    badTimeInd =[(badTimeIndOld-fs/20)'  (badTimeIndOld+fs/20)']; % Taking 50 ms before and one second after bad time segments

    for iL = 1:size(badTimeInd,1)
        badTimes = [badTimes badTimeInd(iL,1): badTimeInd(iL,2)];
    end
%     probeCh(unique(badTimes),:) = [];
end
badTimes = unique(badTimes);

% Picking pixels that are in cortex
clipMaskR = reshape(clipMaskCortex,[size(clipMaskCortex,1)*size(clipMaskCortex,2) 1]); 
processedDatR = reshape(processedDat,[size(clipMaskCortex,1)*size(clipMaskCortex,2) size(processedDat,3)]);
% proc10 = processedDatR;


% Upsample the imaging dataset
parfor iP = 1:size(processedDatR,1)
    processedDat10(iP,:) = interp(processedDatR(iP,:),5);
end 

% Remove bad times from ephys and imaging dataset 
probeCh(badTimes,:) = []; 
badTimes10Hz = unique(floor(badTimeIndOld./100));
processedDat10(:,badTimes10Hz) = [];

clipMaskCortexR = reshape(clipMaskCortex,[size(clipMaskCortex,1)*size(clipMaskCortex,2) 1]);
cortexOnlyMat = processedDat10; 
cortexOnlyMat(~clipMaskCortexR,:) = NaN;

%% Find bad frames

% Use only pixels in cortex

medDatTime = squeeze(median(cortexOnlyMat,1,'omitnan')); 
highT = mean(medDatTime) + 3*std(medDatTime);
lowT  = mean(medDatTime) - 3*std(medDatTime); 
crossing =  find(medDatTime>=highT | medDatTime<=lowT);

% Remove bad frames identified from both ephys and imaging dataset
% To try doing the opposite - imaging first followed by ephys 
processedDat10(:,crossing) =[]; 
cortexOnlyMat(:,crossing) = [];

crossingEphys = crossing.*100; 
diffDuration = [1 diff(crossingEphys)];
diffIdx = find(diffDuration~=100);
ephysBadTimes = [];
% Taking 50 ms before and after bad frame 
for iL = 1:length(crossingEphys)
    ephysBadTimes = [ephysBadTimes crossingEphys(iL)-50: crossingEphys(iL)+50];
%     if iL~= length(diffIdx)
%         ephysBadTimes = [ephysBadTimes crossingEphys(diffIdx(iL)): crossingEphys(diffIdx(iL+1)-1)+100];
%     else
%         ephysBadTimes = [ephysBadTimes crossingEphys(diffIdx(iL)): crossingEphys(end)+100];
%     end
end 
ephysBadTimes = unique(ephysBadTimes);
probeCh(ephysBadTimes,:)=[];

% madDatTime = mad(meanDatTime,1);

% Finding the median frame
medianFrame = median(processedDat10,2,'omitnan');
madFrame = mad(processedDat10,1,2);

threshFrame = medianFrame+3.*madFrame;
threshFrameLow = medianFrame - 3.*madFrame;
numPixels = numel(threshFrame); 
% threshFrameR = reshape(threshFrame,[numPixels 1]);

badFrameIdx = zeros(size(processedDat10,2),1);
for iL = 1: size(processedDat10,2)
    if (numel(find(processedDat10(:,iL)>threshFrame))> numPixels/10)|| (numel(find(processedDat10(:,iL)<threshFrameLow))> numPixels/10)
        badFrameIdx(iL) = 1; 
    end 
end 
badFrameVals = find(badFrameIdx);

%% Get infraslow changes in gamma band power 
probeBandLimited = filtfilt(bG,aG,double(probeCh));
[z,p,k]  = butter(3,[0.01 0.1]./(fs/2),'bandpass');
[sos,g] = zp2sos(z,p,k);
probeBandLimited = abs(single(probeBandLimited));

% Get the envelope of the signal
envelopeBandLimited = envelope(probeBandLimited);

infraEphys = single(downsample(filtfilt(sos,g,double(envelopeBandLimited)),100));

% probeCh = probeCh(t1k,:);
% probeCh = downsample(probeCh,100);

% Keep data constrained within camera frames 
% probe = probeCh(t(1):t(end),:);
% p1 =  probeCh(t1k,:);
% probe = p1; 

%%

figure; imagesc(greenIm{1,1}); colormap gray; axis image off; hold on;
title('Pick a seed pixel that is inside a patch ');
seedLocIn = ginput(1);
% [seedLocInMaster,~] = cpselect(greenMapRef,greenIm{iDate,1},'Wait',true);
seedLocIn = fliplr(round(seedLocIn./3));
% seedLocInMaster = fliplr(round(seedLocInMaster));

seedRad = 12;
processedDat10Im = reshape(processedDat10,[361 438 size(processedDat10,2)]);
inDat = processedDat10Im(seedLocIn(2)-seedRad:seedLocIn(2)+seedRad,seedLocIn(1)-seedRad:seedLocIn(1)+seedRad,:);
inDatR = median(reshape(inDat,[size(inDat,1)*size(inDat,2) size(inDat,3)]),1,'omitnan');
title('Pick a seed pixel that is outside a patch');

seedLocOut = ginput(1); close gcf;
seedLocOut = fliplr(round(seedLocOut./3));
outDat = processedDat10Im(seedLocOut(2)-seedRad:seedLocOut(2)+seedRad,seedLocOut(1)-seedRad:seedLocOut(1)+seedRad,:);
outDatR = median(reshape(outDat,[size(outDat,1)*size(outDat,2) size(outDat,3)]),1,'omitnan');



%%
iDate = 1; iRun = 1;
greenImRGB = ind2rgb(greenIm{iDate,iRun},gray(256));
figure; imagesc(imresize(greenImRGB,[361 438]));
hold on; axis image;

in = [seedLocIn(2)-seedRad:seedLocIn(2)+seedRad; seedLocIn(1)-seedRad:seedLocIn(1)+seedRad];
out = [seedLocOut(2)-seedRad:seedLocOut(2)+seedRad; seedLocOut(1)-seedRad:seedLocOut(1)+seedRad];
% plot([in(2,1) in(2,1)+(seedRad*2) in(2,1)+(seedRad*2) in(2,1) in(2,1)],[in(1,1) in(1,1) in(1,1)+(seedRad*2) in(1,1)+(seedRad*2) in(1,1)],'LineWidth',4);
% plot( [out(2,1) out(2,1)+(seedRad*2) out(2,1)+(seedRad*2) out(2,1) out(2,1)],[out(1,1) out(1,1) out(1,1)+(seedRad*2) out(1,1)+(seedRad*2) out(1,1)],'LineWidth',4);

plot([in(1,1) in(1,1) in(1,1)+(seedRad*2) in(1,1)+(seedRad*2) in(1,1)],[in(2,1) in(2,1)+(seedRad*2) in(2,1)+(seedRad*2) in(2,1) in(2,1)],'LineWidth',4);
plot( [out(1,1) out(1,1) out(1,1)+(seedRad*2) out(1,1)+(seedRad*2) out(1,1)],[out(2,1) out(2,1)+(seedRad*2) out(2,1)+(seedRad*2) out(2,1) out(2,1)],'LineWidth',4);


%% Check for the quality of imaging data by obtaining FC maps
iDate = 1; iRun = 1; 

spatialBin = 3;
figure; imagesc(greenIm{iDate,iRun}); colormap gray; axis image off; hold on;
disp('Pick a seed pixel for getting a FC map ');
seedLoc = ginput(1);
seedLoc = fliplr(round(seedLoc));
plot(seedLoc(2),seedLoc(1),'w.','MarkerSize',15); close gcf;

if exist([templateDir '\' allDates(iDate,:) '\' allRuns{iDate,1}(iRun,:) '\clipMask0' allRuns{iDate,1}(iRun,end) '.BMP'],'file') == 0
    clipMask = imread([templateDir  allDates(iDate,:) '\' allRuns{iDate,1}(iRun,:) '\clipMask0' allRuns{iDate,1}(iRun,end) '.png']);
else
    clipMask = imread([templateDir allDates(iDate,:) '\' allRuns{iDate,1}(iRun,:) '\clipMask0' allRuns{iDate,1}(iRun,end) '.bmp']);
end
clipMask   = imresize(clipMask(:,:,1),1/spatialBin);
%
seedSig = calculateSeedSignal(imresize(greenIm{iDate,iRun},1/spatialBin), clipMask,fliplr(round(seedLoc./spatialBin)),12,processedDat(:,:,:,iDate,iRun)); % Get Gaussian weighted seed signal
corrMap = plotCorrMap(seedSig,processedDat(:,:,:,iDate,iRun),0);

% Overlay correlation map on top of greens

greenImRGB = ind2rgb(greenIm{iDate,iRun},gray(256));
frameSize = size(greenImRGB); 
figure; imagesc(greenImRGB);
hold on;
corrMapR = imresize(corrMap,spatialBin,'OutputSize',[frameSize(1) frameSize(2)]);
imagesc(corrMapR,'alphadata',corrMapR.*0.8);
colormap jet; colorbar; caxis([0 1]);
axis image off;
plot(seedLoc(2),seedLoc(1),'w.','MarkerSize',15);


%%
clear cc lags ccOut lagsOut
sz = length(infraEphys(:,1));

for iCh = 1:size(probeCh,2)
    [cc(:,iCh),lags(:,iCh)] = xcorr(infraEphys(:,iCh)',inDatR(1:sz),10,'normalized'); 
    [ccOut(:,iCh),lagsOut(:,iCh)] = xcorr(infraEphys(:,iCh)',outDatR(1:sz),10,'normalized');
%     [cc(:,iCh),lags(:,iCh)] = xcorr(infraEphysT(:,iCh),inDat10,1000,'normalized'); 
%     [ccOut(:,iCh),lagsOut(:,iCh)] = xcorr(infraEphysT(:,iCh),outDat10,1000,'normalized');
end 

[maxC,ind] = max(cc,[],1);
maxLag = lags(ind);

[maxCOut,ind] = max(ccOut,[],1);
maxLagOut = lagsOut(ind);

figure; stem(maxLag, maxC, 'filled'); hold on; 
stem(maxLagOut, maxCOut, 'filled');

% figure; stem(maxLag(1:24),maxC(1:24),'filled'); hold on; stem(maxLag(25:end),maxC(25:end),'black','filled');xlim([-200 200]); xticklabels(-20:5:20); ylim([-0.1 0.3]); xlabel(' Lag (s)'); ylabel('Correlation');
% figure; stem(maxLagOut,maxCOut,'filled','Color',[0.85 0.325 0.098]);hold on; stem(maxLag(25:end),maxCOut(25:end),'black','filled');xlim([-200 200]); xticklabels(-20:5:20);ylim([-0.1 0.3]);xlabel(' Lag (s)'); ylabel('Correlation');

%%
% imageSize = size(processedDat(:,:,:,iDate,iRun));
imageDataT = reshape(processedDat(:,:,:,iDate,iRun),[imageSize(1)*imageSize(2) imageSize(3)]);

parfor iDat = 1:size(imageDataT,1)
    imData10(iDat,:) = interp(imageDataT(iDat,:),5);
end

imData10(:,badTimes(badTimes<=length(imData10)))= [];
corrInfra10 = corr(imData10', infraEphys);
corrInfra10Im= reshape(corrInfra10,[361 438 size(corrInfra10,2)]);

 meanCorrInfra10 = max(corrInfra10Im,[],3);

figure; greenImRGB = ind2rgb(greenIm{iDate,iRun},gray(256));
imagesc(imresize(greenImRGB,[361 438]));  hold on; axis image;
imagesc(imgaussfilt(meanCorrInfra10,3),'alphadata',imgaussfilt(meanCorrInfra10.*4,1)); colormap jet; caxis([0 0.5]); colorbar;

%% Calculating the averaged seed signal over a few pixels
function seedSig = calculateSeedSignal(green, clipMask, seed, seedRad,frames)
    % choose seed location (gaussian average of a radius around chosen point)
    % seedRad = 6; % 6 = 250um radius, 12 = 500um radius, 24 = 1000um radius
    clear gaussianDisk green_Seed clipMask_Seed 
    if exist('seedRad','var') == 0; seedRad = 6; end
    seedSig = zeros(size(frames,3),1);
    
    gaussianDisk = fspecial('gaussian',2*seedRad+1,seedRad);
    green_seed = green(seed(2)-seedRad:seed(2)+seedRad,seed(1)-seedRad:seed(1)+seedRad,:);
    clipMask_seed = clipMask(seed(2)-seedRad:seed(2)+seedRad,seed(1)-seedRad:seed(1)+seedRad);
    gaussianDisk = gaussianDisk .* double(clipMask_seed);
    gaussianDisk = gaussianDisk ./ sum(gaussianDisk(:));

    for x = 1:2*seedRad+1     % rows are observations, columns are variables
        for y = 1:2*seedRad+1
            seedSig = seedSig + gaussianDisk(x,y) * squeeze(frames(seed(2)+x-seedRad-1,seed(1)+y-seedRad-1,:));

        end
    end
end
