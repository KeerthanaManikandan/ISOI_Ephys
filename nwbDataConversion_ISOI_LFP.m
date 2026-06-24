% Convert ephys and ISOI data to NWB format for DANDI
% This code is only for ISOI-LFP project
% June 23, 2026 - KM (with the help of Claude)

% Refer here: https://docs.dandiarchive.org/getting-started/data-standards/#neurodata-without-borders-nwb
% Install MATNWB by cloning from GitHub : git clone https://github.com/NeurodataWithoutBorders/matnwb.git 
% Then add matnwb to path and run this function generateCore() to
% initialize NWB
% Refer: https://matnwb.readthedocs.io/en/latest/pages/tutorials/ecephys.html

% Set paths
clc; clear;
commonDir = 'C:\Users\kem294\Documents\Data';
cd(commonDir);
addpath(genpath(commonDir)); rmpath(genpath([commonDir '\Codes\nonlinear\functions']));clc;
addpath(genpath([commonDir '\Codes\ISOI_Ephys\neuroshare']));
addpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12']));
rmpath(genpath([commonDir '\Codes\ISOI_Ephys\chronux_2_12\fly_track\videoIO']));

%% Initialize variables and get monkey data
hemisphere = 'Left'; spatialBin = 3;
iM = 1; % 1 - Charlie Sheen, 2 - Whiskey

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

        ageVal = 'P15Y';
        monkeyNameShort = 'CS';

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

         ageVal = 'P6Y';
         monkeyNameShort = 'W';
end


% Get monkey data and parameters
[allDates,allRuns, refDate, refDir,lensCombo, roiSize, ephysFileNameAll, serverPath,probeLabel,...
    chInCortexNotes, greenMapRef] = getMonkeyParams_Imaging_Ephys(monkeyName,commonDir, hemisphere);

% Get monkey data....
[processedDat,greenIm,probe,badCh,badTimesLFP,badTimeThresh,estChInCortex] = ...
    getAllData_Imaging_Ephys(monkeyName,hemisphere,commonDir,serverPath,allDates,allRuns,...
    ephysFileNameAll,greenMapRef,chInCortexNotes,probeLabel,spatialBin);

clc; disp(['All physiology and imaging data for ' monkeyName ' loaded']);

%% Saving data in NWB format to share in DANDI database
clc;
for iDate = 1:size(allDates,1)
    clear expDate;
    expDate = allDates(iDate,:); % Get experiment dates 

    for iRun = 1:size(allRuns{iDate,1},1) % Get individual runs within an experiment 
        clear runName dataDir clipMask elecMask clipMaskCortex corrMask x negIdx ...
            lowIdx recording nwbfile device egroup nchans location_col bad_col ...
            group_col group_name_col electrode_table_region  raw_es lfp_es ...
            eeg_table_region eeg_es lfp_obj bad_ts isoi_series ...


        runName = allRuns{iDate,1}(iRun,:);
        % Replace the 001874 with the DANDI ID you receive when you create
        % a DANDI set.
        dataDir = ['C:\Users\kem294\Documents\Data\DANDI\001874\sub-' monkeyNameShort];

        if ~goodRuns(iDate,iRun); continue; end

        if ~exist(dataDir,'dir'); [~,~] = mkdir(dataDir); end

        if exist(fullfile(dataDir,[monkeyNameShort '_' expDate '_' runName '_nwb.mat']), 'file')
            delete(fullfile(dataDir,[monkeyNameShort '_' expDate '_' runName '_nwb.mat']));
        end

        if exist(fullfile(dataDir,[monkeyNameShort '_' expDate '_' runName '.nwb']), 'file')
            delete(fullfile(dataDir,[monkeyNameShort '_' expDate '_' runName '.nwb']));
        end

        disp(['Saving data for ' monkeyNameShort ' Exp: ' expDate ' run: ' runName(end-1:end)]);

        % Ephys
        recording.lfp       = single(probe{iRun,iDate}.probeCh);  % time x nchan
        recording.raw       = single(probe{iRun,iDate}.rawCh);
        recording.eeg       = single(probe{iRun,iDate}.eegCh);
        recording.fs        = 1e3;
        recording.channels  = single(1:size(probe{iRun,iDate}.probeCh,2));

        recording.bad_channels    =  single(badCh{iDate,iRun});
        recording.bad_times       = single(badTimesLFP{iDate,iRun});
        if isempty(recording.bad_times); recording.bad_times = 0; end
        recording.cortex_channels = single(estChInCortex{iDate}(iRun,1): estChInCortex{iDate}(iRun,2));

        % rs-ISOI
        recording.isoi = processedDat{iDate,iRun}.tempBandPass; % Processed rs-ISOI data
        recording.isoi_fs = single(2); % frame rate
        recording.wavelength = single(630); % nm

        % recording info
        expDateNew = datetime(expDate, 'InputFormat', 'MM_dd_yyyy');
        expDateNew = char(datetime(expDateNew, 'format','yyyyMMdd'));

        recording.subject      = monkeyNameShort;
        recording.session_id   = [expDateNew runName];
        recording.date         = expDateNew;
        recording.brain_region = smFlag(iDate,iRun);

        % Create the NWB file container
        nwbfile = NwbFile(...
            'identifier',recording.session_id,...
            'session_description', 'SimulataneousRS ephys and ISOI recording',...
            'session_start_time', datetime([recording.date ' 00:00:00'], 'InputFormat', 'yyyyMMdd HH:mm:ss', 'TimeZone', 'local'),...
            'general_institution','University of Pittsburgh',...
            'general_lab','Gharbawie Lab',...
            'general_source_script','https://github.com/KeerthanaManikandan/ISOI_Ephys',...
            'general_source_script_file_name','nwbDataConversion_ISOI_LFP.m');

        % Add the animal info
        nwbfile.general_subject = types.core.Subject(...
            'subject_id',recording.subject,...
            'species','http://purl.obolibrary.org/obo/NCBITaxon_27679',...
            'sex','M', ...
            'age',ageVal);

        % Add recording device details
        device = types.core.Device(...
            'description','A1x32-15 mm',...
            'manufacturer','Neuronexus');

        nwbfile.general_devices.set('probe', device);

        % Adding electrode group (one group per shank)
        egroup = types.core.ElectrodeGroup( ...
            'description',  'cortical probe', ...
            'location',     'cortex', ...
            'device',       types.untyped.SoftLink(device));
        nwbfile.general_extracellular_ephys.set('shank1', egroup);

        % Build electrode tables (one row per channel)
        % column descriptions: brain region, bad channel flag,
        % electrode group, group name
        nchans = length(recording.channels);
        location_col   = repmat({recording.brain_region}, nchans, 1);
        bad_col        = ismember(recording.channels, recording.bad_channels)';
        group_col      = repmat(types.untyped.ObjectView(egroup), nchans, 1);
        group_name_col = repmat({'shank1'}, nchans, 1);

        nwbfile.general_extracellular_ephys_electrodes = types.core.ElectrodesTable(...
            'description',  'electrode metadata',...
            'colnames',     {'location', 'bad', 'group', 'group_name'},...
            'id',           types.hdmf_common.ElementIdentifiers('data', int64(0:nchans-1)'),...
            'location',     types.hdmf_common.VectorData('data', location_col,   'description', 'brain region'),...
            'bad',          types.hdmf_common.VectorData('data', bad_col,        'description', 'bad channel flag'),...
            'group',        types.hdmf_common.VectorData('data', group_col,      'description', 'electrode group'),...
            'group_name',   types.hdmf_common.VectorData('data', group_name_col, 'description', 'electrode group name'));

        % Create electrode table region
        electrode_table_region = types.hdmf_common.DynamicTableRegion(...
            'table',        types.untyped.ObjectView(nwbfile.general_extracellular_ephys_electrodes),...
            'data',         int32(0:nchans-1)',...
            'description',  'all electrodes');

        % Converting the ephys into volts from micro volts
        recording.raw = double(recording.raw') / 1e6;
        recording.lfp = double(recording.lfp') / 1e6;
        recording.eeg = double(recording.eeg') / 1e6;


        % Store raw signal
        raw_es = types.core.ElectricalSeries(...
            'data',               recording.raw,...
            'starting_time',        0.0,...
            'starting_time_rate',   double(recording.fs),...
            'electrodes',         electrode_table_region,...
            'description',        'broadband signal downsampled to 1kHz');

        % Store LFP 6-250 Hz
        lfp_es = types.core.ElectricalSeries(...
            'data',               recording.lfp,...
            'starting_time',        0.0,...
            'starting_time_rate',   double(recording.fs),...
            'electrodes',         electrode_table_region,...
            'description',        'LFP 6-250 Hz, bandstop at 60 Hz');

        % Wrap raw signal and LFP in a single module
        lfp_obj = types.core.LFP('ElectricalSeries', lfp_es);
        ecephys = types.core.ProcessingModule(...
            'description',  'ecephys',...
            'LFP',          lfp_obj,...
            'raw',          raw_es);

        nwbfile.processing.set('ecephys', ecephys);

        % Store EEG
        eeg_table_region = types.hdmf_common.DynamicTableRegion(...
            'table',        types.untyped.ObjectView(nwbfile.general_extracellular_ephys_electrodes),...
            'data',         int32(0),...
            'description',  'EEG electrode');

        eeg_es = types.core.ElectricalSeries(...
            'data',         recording.eeg,...
            'starting_time',        0.0,...
            'starting_time_rate',   double(recording.fs),...
            'electrodes',   eeg_table_region,...
            'description',  'Intracranial EEG');

        nwbfile.acquisition.set('EEG', eeg_es);

        % Store bad times
        bad_ts = types.core.TimeSeries(...
            'data',         recording.bad_times,...
            'timestamps',   recording.bad_times / recording.fs,...
            'description',  'indices of bad time samples', ...
            'data_unit',    'samples');

        nwbfile.acquisition.set('bad_times', bad_ts);


        % Store ISOI file
        isoi_series = types.core.ImageSeries(...
            'data',         recording.isoi,...
            'starting_time',        0.0,...
            'starting_time_rate',   double(recording.isoi_fs),...
            'description',  'ISOI frames',...,
            'data_unit' , 'Percent change in reflectance',...
            'format',       'processed');

        nwbfile.acquisition.set('ISOI', isoi_series);
        nwbfile.general_protocol = 'R01NS143962';

        % Export nwb file
        nwbExport(nwbfile, fullfile(dataDir,['sub-' monkeyNameShort '_ses-' expDateNew runName '_ecephys+image.nwb']));
    end
end

%% How to read the NWB data
monkeyName = 'CharlieSheen';
expDate    = '11_20_2023';
runName    = 'run06';
dataDir = ['C:\Users\kem294\Documents\Data\DANDI\001874\sub-' monkeyNameShort];

nwb = nwbRead(ullfile(dataDir,['sub-' monkeyNameShort '_ses-' expDateNew runName '_ecephys+image.nwb']));

% Read LFP
lfpDat = nwb.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get('ElectricalSeries').data;

% Read raw
rawDat = nwb.processing.get('ecephys').nwbdatainterface.get('raw').data;

% Read EEG
eegDat = nwb.acquisition.get('EEG').data;

% Read ISOI
isoiDat = nwb.acquisition.get('ISOI').data;

% Retrieve bad times
badTimesDat = nwb.acquisition.get('bad_times').data;

% Retrieve bad times
badElec = nwb.general_extracellular_ephys_electrodes.vectordata.get('bad').data;

% Example plot
figure; 
subplot(3,1,1); plot(lfpDat(10,:)); box off; 
title('LFP - channel 10'); xlabel('Time (samples)'); ylabel('Voltage (V)');

subplot(3,1,2); plot(rawDat(10,:)); box off;
title('Raw signal - channel 10'); xlabel('Time (samples)'); ylabel('Voltage (V)');

subplot(3,1,3); plot(eegDat(:));box off;
title('EEG'); xlabel('Time (samples)'); ylabel('Voltage (V)');

