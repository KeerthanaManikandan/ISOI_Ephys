function OIAshift2goodstim5(shiftfile3, percentage, shiftthreshold, NStim, NFrame, framerange)

% For processing within stim shifts (shift3), discard large shift stims. % 
% input: 
%   shiftfile3: output file from OIAlign (text file, 11 columns, Nblock*NStim rows)
%   shiftthreshold:  % pixels: if mean pixel shift for framerange exceeds this value, the trial is marked as '0'
%   percentage: % percentage of trials to be inluded, this percentage is for absolute correlation value only. 
%   NStim: number of stimulus, for output format
%   NFrame: number of frames per stim
%   framerange: frames range for sum frames.
% output a file with 'goodstim' info, default filename "withingoodstim + input filename"

% Method: find out those frames has large shift or poor correlation. 

manual=0;   % for testing
if manual
    shiftfile3='_shiftlog2011-10-13-18-13.3.txt';
    shiftthreshold=1;   
    percentage=0.9;     % percentage of '1's
    NStim=9;
    NFrame=16;
    framerange=[5:16];  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outputfile=strcat(shiftfile3(1:end-4), '_all_gs.txt');

[filename, block, stim, frame, dx, dy, bestobj, reffilename, refstim, refframe, currentobj]=OIAReadShiftLog(shiftfile3);
NBlock=size(filename, 1)/NStim/NFrame;
A=zeros(size(filename,1), 4);
A(:,1)=dx;
A(:,2)=dy;
A(:,3)=bestobj;
A(:,4)=currentobj;
A=reshape(A', [4, NFrame, NStim, NBlock]);

if percentage==1
    index=1;
else
   index=ceil(NBlock*NStim*(1-percentage));
end
tempobj=min(A(4, framerange, :, :), [], 2);
tempobj=reshape(tempobj, [NStim*NBlock, 1]);
currentobj=sort(tempobj);
corrthreshold=currentobj(index);

C=ones(4, NStim, NBlock);      % goodstim
for i=1:NBlock
    for j=1:NStim
        shift1=max(A(1, framerange, j, i),[], 2)-min(A(1, framerange, j, i),[], 2);
        shift2=max(A(2, framerange, j, i),[], 2)-min(A(2, framerange, j, i),[], 2);
        skip1=abs(sum(A(1, framerange, j, i), 2))-sum(abs(A(1, framerange, j, i)), 2);
        skip2=abs(sum(A(2, framerange, j, i), 2))-sum(abs(A(2, framerange, j, i)), 2);
%         mean3=mean((A(3, framerange, j, i)), 2);
%         mean4=mean((A(4, framerange, j, i)), 2);
        if shift1>shiftthreshold || skip1<0
                C(1, j, i)=0;
        end
        if shift2>shiftthreshold || skip2<0
                C(2, j, i)=0;
        end
        if min(A(4, :, j, i))<corrthreshold
            C(4, j, i)=0;
        end
    end
end
fid=fopen(outputfile, 'w');
for i=1:NBlock
    for j=1:NStim
        fprintf(fid, '%d \t', C(1,j,i)*C(2,j,i)*C(4,j,i));
    end
    fprintf(fid, '\r\n');
end
fclose(fid);
fprintf('Goodstim3:\rBased on shifts(threshold=%1.1f), goodstim=%2.1f%%\r', shiftthreshold, sum(sum(C(1,:,:).*C(2,:,:)))*100/NBlock/NStim);
fprintf('Based on correlation(threshold=%0.5f), goodstim=%2.1f%%\r', corrthreshold, sum(sum(C(4,:,:)))*100/NBlock/NStim);
fprintf('Based on both, goodstim=%2.1f%%\r', sum(sum(C(1,:,:).*C(2,:,:).*C(4,:,:)))*100/NBlock/NStim);
return
