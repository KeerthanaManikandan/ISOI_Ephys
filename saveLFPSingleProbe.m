function [probeCh,rawCh,eegCh,timeStamp] = saveLFPSingleProbe(runName,datName,fs)
% This function saves/loads LFP data and time series of high frequency
% oscillations from a single electtode recorded during simultaneous imaging
% and physiology. 
% December 13,2023 - KM
% Modified on April 03,2024 to remove 1-5 Hz before saving
% Modified again on April 24 to include time series beyond 250 Hz

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
clear b a bS aS bH aH  rawCh
[b,a]   = butter(3,[1 250]./(fs/2),'bandpass'); % Bandpass filtering parameters across 1-250 Hz
[bS,aS] = butter(3,[57 62]./(fs/2),'stop'); % Bandstop filtering between 57-62 Hz
[bH,aH] = butter(3,250./(30e3/2),'high'); % High pass filtering >250 Hz

disp('Getting the LFP and raw signals for the electrode and filtering them ... ');

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
  
    % Get LFP data
    [~, ~, probeCh(:,iElec)] = ns_GetAnalogData(hFile,lfpEntityID,1,lfpCount);
    [~, ~, rawCh(:,iElec)] = ns_GetAnalogData(hFile,rawEntityID,1,rawCount);
end

probeCh = filtfilt(b,a,probeCh); % Bandpass filtering across 1-250 Hz
probeCh = filtfilt(bS,aS,probeCh); % Bandstop filtering between 57-62 Hz
rawCh   = downsample(filtfilt(bH,aH,rawCh),30); % Highpass frequencies beyond 250 Hz

% Save LFP that is recorded simultaneously with imaging
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
t1k = int32(floor(timeStamp.*1000));
probeCh = single(probeCh(t1k(1):t1k(end),:));
rawCh   = single(rawCh(t1k(1):t1k(end),:));

% if size(probeCh,1)>905000;  probeCh = probeCh(1:905000,:); end

% Getting the EEG
clear elecList elecID
elecList  = find([entityInfo.EntityType] == 2);
elecLabel = {entityInfo(elecList).EntityLabel}; % Gets the label of all the channels in lfpList
elecIdx   = cell2mat(cellfun(@(x) strcmp(x(1:end-2),'analog'),elecLabel,'un',0));
elecList(~elecIdx) = [];
elecLabel(~elecIdx) = [];

if ~isempty(elecList)
    if strcmp(elecLabel{1},'analog 1') % analog 1 has EEG
        clear probe2Temp
        [~, analogInfo2] = ns_GetAnalogInfo(hFile, elecList(1));
        elecCount      = entityInfo(elecList(1)).ItemCount;
        fs_ns5           = analogInfo2.SampleRate;
        [~, ~, eegCh]    = ns_GetAnalogData(hFile,elecList(1),1,elecCount);
        if fs_ns5 == 30e3; eegCh = downsample(eegCh,fs_ns5/fs); end
        eegCh = filtfilt(b,a,eegCh);
        eegCh = eegCh(t1k(1):t1k(end));
    end
end

% Filter the LFP to remove 1-5 Hz
[bL,aL] = butter(3,([6 250]./(fs/2)),'bandpass');
probeCh = single(filtfilt(bL,aL,double(probeCh)));

% Remove 60 Hz power line noise from EEG
eegCh = filtfilt(bS,aS,eegCh); % Bandstop filtering between 57-62 Hz

catch
    disp(['Data did not load for : ' runName]);
end
end