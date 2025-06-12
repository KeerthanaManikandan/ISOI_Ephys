function Polarmap=OIPolar(singlemap, lut, mask, clipsd, LPMethod, LPFKernel, HPMethod, HPFKernel)

% Correct angle bug 061121 HDL
% Modified 051007 HDL
% Polarmap=OIPolar(singlemap, lut, mask, clipsd, lowpass, highpass)
% Calculates polar map from a set of single/difference maps
% input: singlemap: single condition maps (height, width, anglenum)
%           Note: pixel value high means higher response (so regular OI map should be inversed befor calling this function)            
%        lut: color lookuptable, if not defined, no color for polar
%		 mask: a mask used for normalization and filtering, pixels has to be "1" (used) or "0" unused
%        clipsd: input singlemap is clipped before filtering, use 0 for no clip
%        LPMethod: low-pass filter type, usually 'gaussian', or 'oidisk'(better)
%        LPFKernel: LP kernel size, 1-10;
%        HPMethod: high-pass, 'oidisk' is best, 'ribot' is fastest
%        HPFKernel: HP kernel size, 50-100 depends on the map size and noise size
% output: polarmap.ang: angle map (value: 1-256)
%         polarmap.mag: magnitude map (unit: same as input map)
%         polarmap.sum: a simple summation of all conditions (like allo or alld)
% Modified from Xiangming Xu 04 'Ivf2PolarMapsFIN4HD.m'


% try modify following flags to see the effect, which may vary from case to
% case
flagVectNorm=1; % whether or not use Bosking normalization, default 1, no use for 040113GarRun2
flagSubMap=0;   % whether or not use subtraction map (H V A O -> HV AO), default 0, no use for 040113GarRun2

%--------------------------------------
fprintf('-------------- Note: OI map should be reversed before calling "OIPolar", because of negative response\r\r');

n=size(singlemap, 3);
[r, c] = size(singlemap(:, :, 1));
if nargin==1	
	lut=textread('bwfix.lut'); %% if no lookup table provided, use standard 
	mask=ones([r,c]);	%if no mask provided, use default: all map area. 
	clipsd=0;
    LPMethod='gaussian';
    LPFKernel=0;
    HPMethod='oidisk';
    HPFKernel=0;
end	
maskpixnum=sum(sum(mask));

Polarmap.ang=zeros([r,c]);	% angle
Polarmap.mag=zeros([r,c]);	% magnitude
Polarmap.sum=mean(singlemap, 3);	% simple average of all conditions (like lumo)

% Get difference map, according to Bosking 97, may not necessory
% skip this when n is a odd number
if ~mod(n, 2)&flagSubMap
    tempmap=singlemap;
    for i=1:n
        if i<=n/2
            singlemap(:,:,i)=tempmap(:,:,i)-tempmap(:,:,i+n/2);
        else
            singlemap(:,:,i)=tempmap(:,:,i)-tempmap(:,:,i-n/2);
        end	
    end
    clear tempmap;
end

for k = 1:n;
	% clipping
	if clipsd~=0
		singlemap(:,:,k)=OIClip(singlemap(:,:,k), 2, clipsd, mask);
    end
    fprintf('Filtering %d...\r', k);
    singlemap(:,:,k)=OIEasyFilter(singlemap(:,:,k), LPMethod, LPFKernel, HPMethod, HPFKernel, mask);
end

% normalization	(using Bosking 97 method)
if flagVectNorm==1
    imgmean=zeros(n, 1);	% image mean for each image
    imgMAD=zeros(n, 1);	% Mean Absolute Deviation from the mean for each image
    OADmax=-999;		% OverAll max/min of the deviation, (across all maps)
    OADmin=999;
    for i=1:n
        maskarray=reshape(mask, prod(size(mask)),1);
        [junk, maskindex]=sort(maskarray);
        imgarray=reshape(singlemap(:,:,i), prod(size(singlemap(:,:,i))), 1);
        b=imgarray(maskindex(end-maskpixnum+1:end));
        imgmean(i)=mean(b);
        imgMAD(i)=mad(b);
        cmax=max((b-imgmean(i))/imgMAD(i));
        cmin=min((b-imgmean(i))/imgMAD(i));
        if cmax>OADmax
            OADmax=cmax;
        end
        if cmin<OADmin
            OADmin=cmin;
        end	
    end
    imgscale=63.0/(OADmax-OADmin);
    for i=1:n
        b=(singlemap(:,:,i)-imgmean(i))/imgMAD(i);
        singlemap(:,:,i)=(b-OADmin)*imgscale;
    end
end

% calculate polar map
fprintf('Creating polar image, vector calculating...\n');
for k = 1: n
    angle = (360/n)*(k-1);  % note that angle start at 0 and increase
    xcomponent(k) = cos(angle*pi/180);
    ycomponent(k) = sin(angle*pi/180);
end
maxtotalmag =0;
magbuff = zeros(r,c);
preferor = zeros(r,c);
magor = zeros(r,c);
count=0;	% for online display calculating progress
for i = 1:r
    for j = 1:c
         xmag = 0;
         ymag = 0;
         tmpsum=0;      
         for k = 1:n
             response = singlemap(i, j, k);
             xmag = xmag + response*xcomponent(k);
             ymag = ymag + response*ycomponent(k);
             tmpsum = tmpsum + response;
         end
         totalmag = sqrt((xmag*xmag) + (ymag*ymag));
         Polarmap.mag(i, j) = totalmag;
         angle = atan(ymag/xmag)*180/pi;	% atan returns -pi/2 to pi/2 for -infinity to +infinity
         if (xmag <0)
            angle = angle + 180; 		% e.g. y/x=-1 can be two conditions (y=-1 or x=-1)
         end
         angle = mod(angle + 360, 360);		% make angle range from 0 to 359
         Polarmap.ang(i, j) = floor((angle/360)*256);  %% dvided by 360 to half the angle
    end  %%j
end  %%i
return;

% Note:
%compiling OR maps like, Batshelet 1981; Drgao et al., 2000; 
%Worgotter & Eysel, 1991; cos(2*ang) sin(2*ang);  resultant strength r=
%sqrt((sum(xcos 2*ang))^2 + (sum(ycos 2*ang))^2).  In the end, 
%resultant vector angle: atan((sum(ycos 2*ang))^2/(sum(xcos 2*ang))^2)/2  (doubletheta/2)