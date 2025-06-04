function [allDates,allRuns, refDate, refDir,lensCombo, roiSize, ephysFileNameAll, serverPath, probeLabel,chInCortexNotes, greenMapRef] ...
    = getMonkeyParamsEphys_Imaging(monkeyName, commonDir, hemisphere)
% This function retrieves all the variables that are essential to analyze
% imaging and physiology recorded simultaneously for one monkey 

% Load all the runs, experiment dates, green reference images, reference directory for the corresponding monkey
if strcmp(monkeyName,'CharlieSheen') % Charlie Sheen
    allDates     = ['11_29_2021'; '01_11_2022'; '08_07_2023';'11_20_2023'];
    allRuns      = {['run00'; 'run01']; ['run03'; 'run04']; ['run02'; 'run05';'run06'; 'run07'; 'run08'];...
        ['run02'; 'run05'; 'run06'; 'run07' ]};

    % Get imaging parameters
    refDate      = '08_31_2021';
    refDir       = [commonDir '\' monkeyName '_SqM\' hemisphere ' Hemisphere\' refDate '\Master Green Images\'];
    refImageName = 'Charlie Sheen Combined Green 08_31_2021';
    lensCombo        = {['50/50 MF'; '50/50 MF'];['20/50';'20/50'];['50/50 MF';'50/50 MF';'50/50 MF'; '50/50 MF'; '50/50 MF'];...
        ['50/50 MF';'50/80 MF'; '50/80 MF';'50/80 MF']};
    roiSize          = {[64 64];[24 24];[64 64 64 64 64];[38 38 38 38];}; %  500um x 500um

    % Get ephys parameters
    ephysFileNameAll = {['run00'; 'run01'];['run0003';'run0004'];['datafile0002'; 'datafile0005';'datafile0006';...
        'datafile0007';'datafile0008']; ['datafile0002';'datafile0005'; 'datafile0006'; 'datafile0007']};
    serverPath       = '\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\302-19_CharlieSheen_SqM\Left Hemisphere\';
    probeLabel       = {[]; []; ['CDE1';'CDE1';'CDE1';'CDE1';'CDE1';'CDE1']; ['B25A';'B25A'; 'B25A'; 'B25A']};
    chInCortexNotes  = {[1 1];[1 1];[7 10 8 6 4 7];[3 4 5 7]}; % transition determined during experiment


elseif strcmp(monkeyName,'Whiskey') % Whiskey
    allDates = ['08_14_2023'; '10_16_2023'; '12_04_2023';'02_20_2024'];
    allRuns  = {['run02'; 'run04'; 'run05';];[ 'run03'; 'run04'; 'run05'; ...
        'run06'; 'run07'; 'run08'; 'run09'];['run04'; 'run05';'run06'; 'run07']; ...
        ['run01'; 'run02';'run03'; 'run04'; 'run05'; 'run06'; 'run07']};

    % Get imaging parameters
    refDate      = '05_09_2022';
    refDir       = [commonDir '\' monkeyName '_SqM\' hemisphere ' Hemisphere\' refDate '\Master Green Images\'];
    refImageName = 'Combined Green';
    lensCombo        = {['28/50'; '28/50'; '28/50'];[ '50/80 MF'; '50/80 MF'; '50/80 MF'; ...
        '50/80 MF'; '50/80 MF'; '50/80 MF'; '50/80 MF'];['50/80 MF'; '50/80 MF';'50/80 MF'; '50/80 MF'];...
        ['50/80 MF'; '50/80 MF';'50/80 MF'; '50/80 MF';'50/80 MF' ; '50/80 MF'; '50/80 MF']};
    roiSize          = {[36 36 36];[38 38 38 38 38 38 38];[38 38 38 38];[38 38 38 38 38 38 38]}; % for 500 um x 500 um

    % Get ephys parameters
    ephysFileNameAll = {['datafile0002'; 'datafile0004'; 'datafile0005']; ['datafile0003' ;...
        'datafile0004'; 'datafile0005'; 'datafile0006'; 'datafile0007'; 'datafile0008' ;  'datafile0009']; ...
        ['datafile0004'; 'datafile0005'; 'datafile0006'; 'datafile0007']; ['datafile0002'; 'datafile0003';...
        'datafile0004'; 'datafile0005'; 'datafile0006'; 'datafile0007';'datafile0008' ]};
    serverPath       = '\\smb2.neurobio.pitt.edu\Gharbawie\Lab\Data\303-19_Whiskey_SqM\Left Hemisphere\';
    probeLabel       = {[]; ['CDE1';'B25A'; 'B25A'; 'B25A'; 'B25A'; 'B25A' ;'B25A'];['B25A'; 'B25A'; '0763'; '0763'];...
        ['0763'; '0763'; '0763'; '0763'; '0763'; '0763'; '0763']};
    chInCortexNotes  = {[1 1 1]; [5 5 7 4 5 5 6];[6 6 8 5]; [7 6 6 8 7 7 9]};  % transition determined during experiment
end
% Get the master green image
if exist([refDir refImageName '.bmp'],'file') == 0 % Make sure to get what the reference run is outside the function
    greenMapRef = imread([refDir refImageName '.png']);
else
    greenMapRef = imread([refDir refImageName '.bmp']);
end
greenMapRef  = greenMapRef(:,:,1);



end

