function BlockSum(expfolder)

% function BlockSum(path)
% Do automatic block average
% Currently only work for VDAQ 3001, integer32 type (data type 13)

% need add:	goodstim, random, shifting, redshirt

clear all;
system='v';     % currently only process 'vdaq'
expfolder='J:\expt\070315BB_\';
flaggoodstim=0; % =1: 
flagshift=0;    % =1: shift all frames in each stim (need provide 'shiftinput.2.txt' in block folder), =0 will be no shift
runname={       % which run need average.
    'Run01'
    'Run27'
    };
ext='.all';     % filename extension for averaged block.
% summed block will have the same name as first block and above extension, e.g. 'Els_E0B00.all'
% if there is a file with this name already exist, this file will be deleted. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nrun=size(runname, 1);

for i=1:Nrun
    runfolder=strcat(expfolder, getfield(cell2struct(runname(i), 'junk'), 'junk'), '\');
    tempfilename=struct2cell(dir([runfolder, '*.blk']));
    filename=sort(tempfilename(1,:)');
    Nblock=size(filename,1);
    fprintf('Found following ''*.blk'' files:\r');
    for j=1:Nblock
        fprintf('''%s''\n', [runfolder, getfield(cell2struct(filename(j), 'junk'), 'junk')]);
    end
    anapar=OIHeadRead([runfolder, getfield(cell2struct(filename(1), 'junk'), 'junk')], system);
    FrameWidth=anapar.FrameWidth;
    FrameHeight=anapar.FrameHeight;
    FramesPerStim=anapar.FramesPerStim;
    NCond=anapar.NStim;
    DataType=anapar.DataType;
    
    if flaggoodstim
        goodstim=textread(strcat(runfolder, 'goodstim.txt'), '%d');
        if isempty(goodstim)|size(goodstim,1)~=Nblock*NCond;
            fprintf('Error, "goodstim.txt" does not contain right number of conditions\r');
        end
        goodstim=reshape(goodstim, [NCond, Nblock]);
        goodstim=goodstim';
    else
        goodstim=ones(Nblock, NCond);	% if no goodstim.txt provided, use all conditions.    
    end

    if flagshift
        [sfname, sfblock, sfstim, sfframe, sfx1, sfy1, sfcoor, junk1, junk2, junk3, junk4]=textread(strcat(runfolder, 'shiftinput.2.txt'), '%s %d %d %d %f %f %f %s %d %d %f');
        if size(sfname,1)~=blockfilenum*NCond
            fprintf('Error: "shiftinput.2.txt" doesnot contain Nblock*NCond number of entries!\r');
        end
        shiftstructname=cell2struct(sfname, 'sfname', 2);   % sfname is cell type, others are arrays
        for k=1:Nblock
            for i=1:NCond
                if getfield(cell2struct(filename(k), 'junk'), 'junk')~=shiftstructname((k-1)*NCond+i).sfname
                    fprintf('\rwrong file name match: %s vs %s', getfield(cell2struct(filename(k), 'junk'), 'junk'), shiftstruct(k*NCond*NFrames).sfname);
                end
            end
        end
    end

    
    for j=1:Nblock  % check if all blocks are the same type
        anapar1=OIHeadRead([runfolder, getfield(cell2struct(filename(j), 'junk'), 'junk')], system);     
        if  FrameWidth~=anapar1.FrameWidth | ...    
            FrameHeight~=anapar1.FrameHeight | ...
            FramesPerStim~=anapar1.FramesPerStim | ...
            NCond~=anapar1.NStim | ...
            DataType~=anapar1.DataType
            fprintf('not all blocks are the same, fail for this run\r');
            break
        end
    end

    fidfirstblock=fopen([runfolder, getfield(cell2struct(filename(1), 'junk'), 'junk')], 'r');   % simply copy an existing file header since they are the same, except data type
    header=fread(fidfirstblock, 1716, 'uint8');
    header(29)=uint8(13);   % force to interger data type
    fclose(fidfirstblock);
    sumfilename=[runfolder, getfield(cell2struct(filename(1), 'junk'), 'junk')];
    sumfilename=[sumfilename(1:end-4), ext];  % summed block file name
    if ~isempty(dir(sumfilename))    % check if this file is already exist
        fprintf('this file already exists: \r\t\t%s,\r press Enter to delete, ctrl+c to abort\r', sumfilename);
        pause;
        delete(sumfilename);
        fprintf('file deleted\r');
    end

    
    fidsumblock=fopen(sumfilename, 'a');
    fwrite(fidsumblock, header, 'uint8');   % write header
    % To avoid memory problem, Read/average stim by stim
    for j=1:NCond
        fprintf('Process stim # %d: ', j);
        AvgStim=zeros(FrameHeight*FrameWidth*FramesPerStim,1);
        count=0;
        for k=1:Nblock
            fprintf('blk%d ', k);
            blkfilename=[runfolder, getfield(cell2struct(filename(k), 'junk'), 'junk')];
            fidblk=fopen(blkfilename,'rb','ieee-le'); % open block file
            % random (determin stimID here if it's random)
                        stimloc=j;  
                        
            fseek(fidblk, 1716+FrameWidth*FrameHeight*FramesPerStim*(j-1)*4, -1);
    	    currentstim=fread(fidblk,[FrameHeight*FrameWidth*FramesPerStim],'uint32=>double'); 
    	    if goodstim(k, j)==1
                if flagshift
                    for i=1:FramePerStim
                        currentstim(:,:,i)=OIShift(currentstim(:,:,i), sfx1((k-1)*NCond+stimloc), sfy1((k-1)*NCond+stimloc));
                    end
                end        	    
                AvgStim=AvgStim+currentstim;
                count=count+1;
            end
            fclose(fidblk);
        end
        fprintf('\r');
        fwrite(fidsumblock, AvgStim./count, 'uint32');
    end
    fclose(fidsumblock);
end
return;
    
    
    
% if (flagrandom==1);
%     file=[blockfolder, '_stimseq.txt'];
%     stimseq=textread([blockfolder, '_stimseq.txt'], '%d');
% 	if size(stimseq,1)~=blockfilenum*NCond
%         fprintf('Error: number of stim in ''_stimseq.txt'' is %d\r', size(stimseq));
%     end
%     stimseq=reshape(stimseq, [NCond, blockfilenum])';
%     % check if it's a valid stim sequence file
%     if max(max(stimseq))>NCond | min(min(stimseq))<1
%         fprintf('Error: stimid in ''_stimseq.txt'' must >0 and <%d\n', NCond);
%     end
%     for i=1:blockfilenum
%         checksum=zeros(1, NCond);
%         checksum(stimseq(i,:))=1;
%         if sum(checksum)~=NCond
%             fprintf('Error: checksum of ''_stimseq.txt'' is not equal to stim number (%d) at line %d\n', NCond, i);
%         end
%     end
% else 
%     stimseq=repmat([1:NCond], blockfilenum, 1);    
% end
% 
% 
