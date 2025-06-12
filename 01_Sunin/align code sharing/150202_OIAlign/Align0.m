% Alignuser.m  Optical Imaging image alignment     HDL 060417 modified from 'suninalign.m'
% This file is not a function. 
% This file will be called by 'sunincore', and use parameters defined in 'suninuser' 
% This file also can be used alone (set manual==1 and provide parameters)
% This file calls 'OIAlignCore.m' which aligns two frames.
% input: 
%      Optical Imaging Data location
%      Align type (see below)
%      Obj function (see below)
%      Base frame  (for certain alignment need a common base frame)
%      Preprocess (see below)
% output: 
%      Shiftlog (dx, dy to shift frame back)
%      goodstim.txt: evaluation of  alignment results (0 for good conditions), note this depends on 'Aligntype'.

% aligntype:
%      1:  General alignment, align all frames to one selected frame (xxx), output 'shiftlog.1.txt'
%      2:  First frame alignment, align all first frames (can be a range) in each stim to a selected frame 
%          (usually a good frame in the middle of the imaging session,e.g. first frame in block 100), result: 'shiftlog.2.txt', 'goodstim-ff.txt'
%      3:  Within condition alignment, align all subsquent frame to the first frame selected in each stim presentation. result: 'shiftlog.3.txt', 'goodstim-within.txt'
%      4:  combination of 2 and 3
% objfun: Objective function
%      1: fast correlation, use 'normxcorr2.m', precision=1 pixel
%      2: slow correlation, a combination of fast correlation (1) and sub-pixel correlation, default precision=0.1 pixel (change 'precision) to modify this default.  
%      3: use stdev of difference image (default precision=1 pixel)
%      4: min-difference score (default precision=1 pixel)
%      5: use mutual info to align (note: only uint8 image can be use here, precision=1 pixel)
% Pre-Process:
%      *  illumination bkground subtraction
%      *  low-pass filtering
%      *  cropping
%      1: masking, need implement
%      2: contrast enhancement, need implement
%      3: dark sub, need implement
% HR 121210 used on L3 awake imaging alignment.
% HR 130317 used on mask making.
% HR 141217 for 141216_L45 data

function Align0(datadriver)
manual=1;   % set manual=0 when it is called by 'suninuser.m'
if manual==1;
    clear              % do not use this if this function is not used alone
    system='v';             % 'v' for VDAQ, 'r' for RedShirt
    datadriver = 'G:\';     % Data disk name
    datafolder = '00_data\';   % Data folder name on data disk, results will be saved in 'expresult'
    expname = '140101_L36\'; % Exp folder name (in both data folder and result folder)
    runname = 'Run02_G8\';      % Run foler name (in both data folder and result folder)
    blkfilename={
    };
    aligntype=4;    % 1: general, 2: fframe, 3:within stim, 4: combine of 2&3, 5: one frame
    objfun = 6;     % 1: fast correlation, 2: slow correlation, 3: stdev of difference, 4: min-difference, 5: mutual info, 6: ratio
    shiftrange = [-6 6 -6 6];  % how many pixels to shift to search for best fit, towards left, right, up, down
    baseframe_shiftframe_not_same = 1;
    if baseframe_shiftframe_not_same == 1
        b_datadriver = 'G:\';     % Data disk name
        b_datafolder = '00_data\';   % Data folder name on data disk, results will be saved in 'expresult'
        b_expname = '140101_L36\'; % Exp folder name (in both data folder and result folder)
        b_runname = 'Run01_OD8\';      % Run foler name (in both data folder and result folder)
    else
        b_datadriver = datadriver;     % Data disk name
        b_datafolder = datafolder;   % Data folder name on data disk, results will be saved in 'expresult'
        b_expname = expname; % Exp folder name (in both data folder and result folder)
        b_runname = runname;      % Run foler name (in both data folder and result folder)
    end
    baseframefile={'OD8_E01B000.BLK', 1, 1};  % filename, stim number, frame number for 'aligntype=2 or 4'
    ffrange=[1];      % the range of first frame to be aligned in 'aligntyp=2 or 4'

    LPKernel=2;       % low-pass kernel, part of pre-process befor alignment, put 0 if no need 
%     LPKernel=0;       % low-pass kernel, part of pre-process befor alignment, put 0 if no need 
    HPKernel=0;
    crop1=[280 162 299 181]; %140101_SponRun00 [x1, y1, x2, y2] of clip region, part of pre-process, put 0 if no need
    flagillubksub=0;  % illum bkground sub, calculated once per block with 50kernel, part of pre-porcess, 0 for none
    percentage2=0.90;	% percent of trials to be included after shift2 process
    percentage3=0.90;	% percent of trials to be included after shift3 process
    percentage4=1;	% percent of trials to be included after shift3 process
end
% end of input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if aligntype~=2 
    if size(ffrange, 2)~=1
        fprintf('Error: Check ffrange, only aligntype==2 can have more than one first frame selected!\n');
    else
        fframe=ffrange(1,1);    % used in aligntype=3
    end    
end
    
blockfolder = strcat(datadriver, datafolder, expname, runname);
b_blockfolder = strcat(b_datadriver, b_datafolder, b_expname, b_runname);
ctime=fix(clock);
% A shift log with date and time as file name will always saved in data folder
if aligntype==1 | aligntype==2 | aligntype==3
    shiftlog = strcat(blockfolder, '_shiftlog', strcat(num2str(ctime(1)), '-', num2str(ctime(2)),'-',num2str(ctime(3)),'-',num2str(ctime(4)),'-',num2str(ctime(5))),'.', num2str(aligntype), '.txt');
    fidshiftlog = fopen(shiftlog, 'w');  %for shiftlog output
end
if aligntype==4
    shiftlog1 = strcat(blockfolder, '_shiftlog', strcat(num2str(ctime(1)), '-', num2str(ctime(2)),'-',num2str(ctime(3)),'-',num2str(ctime(4)),'-',num2str(ctime(5))),'.2.txt');
    fidshiftlog1 = fopen(shiftlog1, 'w');  %for shiftlog output
    shiftlog2 = strcat(blockfolder, '_shiftlog', strcat(num2str(ctime(1)), '-', num2str(ctime(2)),'-',num2str(ctime(3)),'-',num2str(ctime(4)),'-',num2str(ctime(5))),'.3.txt');
    fidshiftlog2 = fopen(shiftlog2, 'w');  %for shiftlog output
    
    fidshiftinput1=fopen(strcat(blockfolder, 'shiftinput.1.txt'), 'w');
    fidshiftinput=fopen(strcat(blockfolder, 'shiftinput.2.txt'), 'w');

end
if aligntype==5
    shiftlog3 = strcat(blockfolder, '_shiftlog', strcat(num2str(ctime(1)), '-', num2str(ctime(2)),'-',num2str(ctime(3)),'-',num2str(ctime(4)),'-',num2str(ctime(5))),'.4.txt');
    fidshiftlog3 = fopen(shiftlog3, 'w');  %for shiftlog output
    fidshiftinput3=fopen(strcat(blockfolder, 'shiftinput.txt'), 'w');
end


if isempty(blkfilename)
    if system=='v'
        tempfilename=struct2cell(dir([blockfolder, '*.blk']));
    elseif system=='r'
        fprintf('Note: you may need delete non-block "*.da" files in data folder\n');
        tempfilename=struct2cell(dir([blockfolder, '*.da']));
    end
    blkfilename=sort(tempfilename(1,:)');
    for i=1:size(blkfilename,1)
        fprintf('''%s''\n', getfield(cell2struct(blkfilename(i), 'junk'), 'junk'));
    end
    fprintf('\nfound %d blk files(sorted, check sequence).\n', size(blkfilename,1));
end
blocknum=size(blkfilename, 1);      % how many blocks

% Read head info
anapar=OIHeadRead(strcat(blockfolder,getfield(cell2struct(blkfilename(1), 'junk'), 'junk')), system);
FrameWidth=anapar.FrameWidth;
FrameHeight=anapar.FrameHeight;
FramesPerStim=anapar.FramesPerStim;
NStim=anapar.NStim;

baseframe=OIReadFrame(strcat(b_blockfolder, getfield(cell2struct(baseframefile(1), 'junk'), 'junk')), 'v', getfield(cell2struct(baseframefile(2), 'junk'), 'junk'), getfield(cell2struct(baseframefile(3), 'junk'), 'junk'));

% Preprocess
if flagillubksub
    baseillubk=OIMeanFilt(baseframe, 50);
    imbase=OIAPreProcess(baseframe, baseillubk, LPKernel, [0 0], crop1,HPKernel);
else
    imbase=OIAPreProcess(baseframe, 0, LPKernel, [0 0], crop1,HPKernel);    
end

if system=='v'
    headerlength = 1716;    % byte
elseif system=='r'
    headerlength = 5120;
end
timertrigger=1;
ctime=clock;
for k=1:blocknum  % Start read/process block by block
    fprintf('\rblock=%d  ',k);
    filename=getfield(cell2struct(blkfilename(k), 'junk'), 'junk');
    if flagillubksub    % calculate illumination bk
        illubk=OIReadFrame(strcat(blockfolder,filename), system, 3,1);  % note: 3 is just for selecting a middle frame to calculate illumination background
        illubk=OIMeanFilt(illubk, 50);              
    else
        illubk=0;
    end    
    for i=1:NStim
        switch aligntype
        case 1  % General alignment
            tDCFrames=OIReadStim(strcat(blockfolder,filename), i, system);
            for j = 1:FramesPerStim            
                if timertrigger==1|timertrigger==2  % estimate time
                    [ctime timertrigger]=OITimer(clock, ctime, timertrigger, blocknum, NStim*FramesPerStim);
                end
                imshift=OIAPreProcess(tDCFrames(:,:,j), illubk, LPKernel, [0 0], crop1,HPKernel);
                [junk, pixoff, bestobj, currentobj] = OIAlignCore(baseframe, imshift, shiftrange, objfun); % find the best alignment
                OIAWriteShiftLog(fidshiftlog, filename, k, i, j, pixoff(1), pixoff(2), bestobj, getfield(cell2struct(baseframefile(1), 'junk'), 'junk'), getfield(cell2struct(baseframefile(2), 'junk'), 'junk'), getfield(cell2struct(baseframefile(3), 'junk'), 'junk'), currentobj);
            end
        case 2        % First frame alignment
            tDCFrames=OIReadStim(strcat(blockfolder,filename), i, system);
            if timertrigger==1|timertrigger==2  % estimate time
                [ctime timertrigger]=OITimer(clock, ctime, timertrigger, blocknum, NStim*size(ffrange, 2));
            end
            for j=ffrange
                imshift=OIAPreProcess(tDCFrames(:,:,j), illubk, LPKernel, [0 0], crop1,HPKernel);
                [imshifted, pixoff, bestobj, currentobj] = OIAlignCore(imbase, imshift, shiftrange, objfun);
                OIAWriteShiftLog(fidshiftlog, filename, k, i, j, pixoff(1), pixoff(2), bestobj, getfield(cell2struct(baseframefile(1), 'junk'), 'junk'), getfield(cell2struct(baseframefile(2), 'junk'), 'junk'), getfield(cell2struct(baseframefile(3), 'junk'), 'junk'), currentobj);
%                imwrite(nc(imshifted-imbase), strcat('tempd\', num2str(k), '-', num2str(i), '.bmp')); 
            end
        case 3        % Within condition alignment
            if timertrigger==1|timertrigger==2  % estimate time
                 [ctime timertrigger]=OITimer(clock, ctime, timertrigger, blocknum, NStim);
            end            
            tDCFrames=OIReadStim(strcat(blockfolder,filename), i, system);
            imfframe=OIAPreProcess(tDCFrames(:,:,fframe), illubk, LPKernel, [0 0], crop1,HPKernel);        % note: a bug here, fframe maybe more than one
            OIAWriteShiftLog(fidshiftlog, filename, k, i, 1, 0, 0, 1, filename, i, 1, 1);          % just for consistency
            for j = 2:FramesPerStim
                if timertrigger==1|timertrigger==2  % estimate time
                     [ctime timertrigger]=OITimer(clock, ctime, timertrigger, blocknum, NStim*(FramesPerStim-1));
                end            
                imshift=OIAPreProcess(tDCFrames(:,:,j), illubk, LPKernel, [0 0], crop1,HPKernel);
                [junk, pixoff, bestobj, currentobj] = OIAlignCore(imfframe, imshift, shiftrange, objfun); 
                OIAWriteShiftLog(fidshiftlog, filename, k, i, j, pixoff(1), pixoff(2), bestobj, filename, i, 1, currentobj);
            end
        case 4        % combine 2 and 3
            tDCFrames=OIReadStim(strcat(blockfolder,filename), i, system);
            if timertrigger==1|timertrigger==2  % estimate time
                 [ctime timertrigger]=OITimer(clock, ctime, timertrigger, blocknum, NStim);
            end            
            for j=ffrange
                imshift=OIAPreProcess(tDCFrames(:,:,j), illubk, LPKernel, [0 0], crop1,HPKernel);
                [imshifted, pixoff, bestobj, currentobj] = OIAlignCore(imbase, imshift, shiftrange, objfun);
                OIAWriteShiftLog(fidshiftlog1, filename, k, i, j, pixoff(1), pixoff(2), bestobj, getfield(cell2struct(baseframefile(1), 'junk'), 'junk'), getfield(cell2struct(baseframefile(2), 'junk'), 'junk'), getfield(cell2struct(baseframefile(3), 'junk'), 'junk'), currentobj);
                fprintf(fidshiftinput, '%s \t%d \t%d \t%d \t%3.1f \t%3.1f \t%6.6f\r\n', filename, k, i, j, pixoff(1), pixoff(2), bestobj);
            end
            imfframe=OIAPreProcess(tDCFrames(:,:,fframe), illubk, LPKernel, [0 0], crop1,HPKernel);        % note: a bug here, fframe maybe more than one
%             OIAWriteShiftLog(fidshiftlog2, filename, k, i, 1, 0, 0, 1, filename, i, 1, 1);          % just for consistency
            for j = 1:FramesPerStim
                imshift=OIAPreProcess(tDCFrames(:,:,j), illubk, LPKernel, [0 0], crop1,HPKernel);
                [junk, pixoff, bestobj, currentobj] = OIAlignCore(imfframe, imshift, shiftrange, objfun);        % Note: here forced 'objfun' to 1
                OIAWriteShiftLog(fidshiftlog2, filename, k, i, j, pixoff(1), pixoff(2), bestobj, filename, i, 1, currentobj);
                fprintf(fidshiftinput1, '%s \t%d \t%d \t%d \t%3.1f \t%3.1f \t%6.6f\r\n', filename, k, i, j, pixoff(1), pixoff(2), bestobj);
           
            end
            case 5        % align to one frame
            tDCFrames=OIReadStim(strcat(blockfolder,filename), i, system);
            if timertrigger==1|timertrigger==2  % estimate time
                 [ctime timertrigger]=OITimer(clock, ctime, timertrigger, blocknum, NStim);
            end            
 %             OIAWriteShiftLog(fidshiftlog2, filename, k, i, 1, 0, 0, 1, filename, i, 1, 1);          % just for consistency
            for j = 1:FramesPerStim
                imshift=OIAPreProcess(tDCFrames(:,:,j), illubk, LPKernel, [0 0], crop1,HPKernel);
                [junk, pixoff, bestobj, currentobj] = OIAlignCore(imbase, imshift, shiftrange, objfun);        % Note: here forced 'objfun' to 1
                OIAWriteShiftLog(fidshiftlog3, filename, k, i, j, pixoff(1), pixoff(2), bestobj, filename, i, 1, currentobj);
                fprintf(fidshiftinput3, '%s \t%d \t%d \t%d \t%3.1f \t%3.1f \t%6.6f\r\n', filename, k, i, j, pixoff(1), pixoff(2), bestobj);
           
            end
        otherwise
            fprintf('error: wrong alignment method\n');
        end        
        % close & reopen files
		switch aligntype
		case 1
            fclose(fidshiftlog);
		case 2
            fclose(fidshiftlog);
		case 3
            fclose(fidshiftlog);
		case 4
            fclose(fidshiftlog1);
            fclose(fidshiftlog2);
            fclose(fidshiftinput);
            fclose(fidshiftinput1);
        case 5
            fclose(fidshiftlog3);
            fclose(fidshiftinput3);            
		end
		switch aligntype
		case 1
            fidshiftlog = fopen(shiftlog, 'a');
		case 2
            fidshiftlog = fopen(shiftlog, 'a'); 
		case 3
            fidshiftlog = fopen(shiftlog, 'a'); 
		case 4
            fidshiftlog1 = fopen(shiftlog1, 'a'); 
            fidshiftlog2 = fopen(shiftlog2, 'a'); 
            fidshiftinput=fopen(strcat(blockfolder, 'shiftinput.2.txt'), 'a');
            fidshiftinput1=fopen(strcat(blockfolder, 'shiftinput.1.txt'), 'a');
        case 5
            fidshiftlog3 = fopen(shiftlog3, 'a'); 
            fidshiftinput3=fopen(strcat(blockfolder, 'shiftinput.txt'), 'a');     
		end
    end % for i=1:NStim
end % for k1:blocknum

switch aligntype
case 1
    fclose(fidshiftlog);
%    OIAshift2goodstim1(shiftlog);  % not implemented yet
case 2
    fclose(fidshiftlog);
    OIAshift2goodstim2(shiftlog1, 0.0001, NStim, percentage2);
case 3
    fclose(fidshiftlog);
    OIAshift2goodstim3(shiftlog2, percentage3, 1, NStim, FramesPerStim, [7:16]);
case 4
    fclose(fidshiftlog1);
% 	OIAshift2goodstim2(shiftfile2, threshold, NStim, percentage)
    OIAshift2goodstim2(shiftlog1, 0.0001, NStim, percentage2);
    fclose(fidshiftlog2);
%   OIAshift2goodstim3(shiftfile3, percentage, shiftthreshold, NStim, NFrame, framerange)    
    OIAshift2goodstim4(shiftlog2, percentage3, 0.45, NStim, FramesPerStim, [1:16]);
    fclose(fidshiftinput);
    fclose(fidshiftinput1);
    
case 5
    fclose(fidshiftlog3);
    OIAshift2goodstim4(shiftlog3, percentage4, 0.3, NStim, FramesPerStim, [1:16]);
    fclose(fidshiftinput3);
    
otherwise
end


return