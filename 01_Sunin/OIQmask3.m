function values=OIQmask3(data, mask)
% Get average pixel values from matrix 'data' using 'mask' map
% mask map may contain multiple domains, this program will indentify 
% each domain (uses imfill.m, >R13) and average pixels in each domain. 
% Different from OIQmask1: this function returns average value in each domain while OIQmask1 gives pixel value
% This function calls OIQmask1

% input
% data: map to quantify
% mask: a black white mask

% output:
% values: nx1 array, n=domain number
manual=0;
if manual
    data=OIReadIVF('mag.ivf');
    mask=imread('hzvtv1.bmp');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgsize=size(data);
data=double(data);
masktemp1=double(mask);		% for searching for 1's, double() is not necessory
masktemp2=logical(1-mask);	% for filling
Ndomain=0;
repeat=1;
while repeat
    [y, x]=find(masktemp1, 1);
    if isempty(y)
        repeat=0;
        fprintf('total domain number: %d\r', Ndomain)
        break;
    else
        Ndomain=Ndomain+1;
        masktemp2=imfill(masktemp2, [y,x]);
        masktemp3=double(~masktemp2);
        singledomain=masktemp1-masktemp3;
        masktemp1=masktemp3;
%       imwrite(norm_to_uint8(double(singledomain)), strcat(num2str(Ndomain), '.bmp'));
		values(Ndomain)=mean(OIQmask1(data, singledomain)); 
%		values(Ndomain)=OIQmask2bak(data, singledomain);    % a little faster
    end
end

if manual
    fprintf('%3.2f\t', values);
    fprintf('\r%3.2f(N=%d)\r', mean(values), Ndomain);
end
return


    