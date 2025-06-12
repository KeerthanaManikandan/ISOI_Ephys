function [imresult, pixoff, bestobj, currentobj] = OIAlignCore(imbase, imshift, range, objfun)
% function [imresult, pixoff, bestobj] = OIAlignCore(imbase, imshift, objfun)
% Align 'imshift' with 'imbase'     HDL (060417 modified from OIAlign2.m)
% Doesn't include preprocess.
%
% input
% 	imbase: pedestal image to be compared, usually is the very first frame
% 	imshift: image to be shifted & aligned with 'imbase', size should be the same as 'imbase'
% 	range(4): (x1, x2, y1, y2) for (left, right, top, bottom) max shift range
%   objfun: Objective function
%       1: fast correlation, use 'normxcorr2.m', precision=1 pixel
%       2: slow correlation, a combination of fast correlation (1) and sub-pixel correlation, default precision=0.1 pixel (change 'precision) to modify this default.  
%       3: use stdev of difference image (default precision=1 pixel)
%       4: min-difference score (default precision=1 pixel)
%       5: use mutual info to align (note: only uint8 image can be use here, precision=1 pixel)
% return: 
%   imresult: shifted image aligned with 'imbase' (now block this function to improve speed, 060417)
%   pixoff(2): (dx, dy) values: shift 'imshift' by (dx, dy) will achieve best alignment with 'imbase'
%   bestobj: value of objective function at best alignment. 
%   currentobj: value of objective function at current (non-shift) position
%
% Note: 'imshift' should be the size as 'imbase'
%       Datatype: 'imshift' and 'imbase' should be 2D grayscale image, any data type (but when objfun=5, they will be transformed to 'uint8).
% improvement: use range() to replace croppixel, will improve some speed.
%

%% Initial check
basesize = size (imbase);
imshiftsize = size (imshift);
imresult=imshift;       % note, 'imresult' is not imresult
if size(basesize)>2 | size(imshiftsize)>2
    fprintf('Error in OIAlignCore: only 2D map can be aligned (no color!)\r\r');
    return;
end
if size(basesize)~=size(imshiftsize)
    fprintf('Error in OIAlignCore: sizes do not match (imbase, imshift)!\r\r');
    return;
end
croppixel= ceil(max(max(abs(range)))+1); % how many pixel will corp off along the edges, for obtain a smaller area for comparison, since shifting creates 'blank' margins
if min(basesize)<2*croppixel
    fprintf('\rError in OIAlignCore: min(basesize)<2*croppixel ! (need increase selected size or decrease shift range)\r\r');
    return;
end

if objfun==1 | objfun==2 |  objfun==6
	imshiftcrop=imshift(croppixel+1:end-croppixel, croppixel+1:end-croppixel);
	cc=normxcorr2(imshiftcrop, imbase);	
	[bestobj, maxlocation] = max(cc(:));	%abs(cc)?
	[yoffset, xoffset] = ind2sub(size(cc), maxlocation(1));
	pixoff(1)=-(size(imbase,2)-croppixel-xoffset);	% note: reverse sign to be comparable with score 1-3 (shift back pixel values)
	pixoff(2)=-(size(imbase,1)-croppixel-yoffset);
    currentobj=cc(size(imbase,1)-croppixel, size(imbase,2)-croppixel);
    if currentobj>bestobj
        fprintf('strange!\r');
    end
%	imresult=OIShift(imshift, pixoff(1), pixoff(2));
end
if objfun==1    % fast correlation, no further alignment
    return;
end             % slow correlation, do fine-tuning alignment based on fast correlation results (normcorr2)

if objfun==2    % Do fine-tuning shift 
    precision=0.1;  % step size in shifting an image
    croppixel=croppixel+1;  
    imbasesm=imbase(croppixel+1:end-croppixel, croppixel+1:end-croppixel);
    pixofftemp=pixoff;
    for xstep=pixoff(1)-0.6:precision:pixoff(1)+0.6    % x shifting         
        for ystep=pixoff(2)-0.6:precision:pixoff(2)+0.6   % y shifting
            imtemp=OIShift(imshift, xstep, ystep);     
            imtempsm=imtemp(croppixel+1:end-croppixel, croppixel+1:end-croppixel);
            cortemp=corr2(imbasesm, imtempsm);
            if (cortemp>bestobj)
                bestobj=cortemp;
                pixofftemp=[xstep ystep];
            end
        end
    end
%    currentobj=corr2(imbasesm, imtempsm);  % no more calculation of currentobj, use the old one (in last coarse shift)
    pixoff=pixofftemp;
%	imresult=OIShift(imshift, pixoff(1), pixoff(2));
    return;
end

if objfun==3    % stdev of difference
    bestobj=9999999999;
    precision=1;    % modify here if need higher precision
    imbasesm=imbase(croppixel+1:end-croppixel, croppixel+1:end-croppixel);
    for xstep=range(1):precision:range(2)    % x shifting
        for ystep=range(3):precision:range(4)   % y shifting
            imtemp=OIShift(imshift, xstep, ystep);     
            imtempsm=imtemp(croppixel+1:end-croppixel, croppixel+1:end-croppixel);
            diffimg=imbasesm-imtempsm;
            objtemp=std2(diffimg);
            if (objtemp<bestobj)
                bestobj=objtemp;
                pixoff=[xstep ystep];
            end
            if xstep==0 & ystep==0
                currentobj=cortemp;
            end
        end
    end
%    imresult=OIShift(imshift, pixoff(1), pixoff(2));
    return;
end
         
if objfun==4    % minimum difference 
    precision=1;    % modify here if need improve precision
    maxobj=prod(basesize)*99999;   %estimate a large number
    imbasesm=imbase(croppixel+1:end-croppixel, croppixel+1:end-croppixel);
    for xstep=range(1):precision:range(2)    % x shifting
        for ystep=range(3):precision:range(4)   % y shifting
            imtemp=OIShift(imshift, xstep, ystep);     
            imtempsm=imtemp(croppixel+1:end-croppixel, croppixel+1:end-croppixel);
            diffimg=imbasesm-imtempsm;
            objtemp=sum(sum(abs(diffimg)));
            if (objtemp<bestobj)
                bestobj=objtemp;
                pixoff=[xstep ystep];
            end
            if xstep==0 & ystep==0
                currentobj=cortemp;
            end
        end
    end
%    imresult=OIShift(imshift, pixoff(1), pixoff(2));
    return;
end

% note: following method need check
if objfun==5	% use "image_registr_MI.m"
   %[h,im_matched, theta,I,J]=image_registr_MI(image1, image2, angle, step,crop);   % larger image2(imshift) is register to smaller image1(imbase)
    [entropy,im_matched, theta,yoffset,xoffset]=image_registr_MI(nc(imbasesm(range(4)+1:end+range(3), range(2)+1:end+range(1))), nc(imshift), 0, step1, 0); % angle 0, no crop
    pixoff(1)=xoffset-range(2)-1;
    pixoff(2)=yoffset-range(4)-1;
    bestobj=max(entropy(:));
    currentobj=entropy(range(2), range(4));
%    imresult=OIShift(imshift, pixoff(1), pixoff(2));
	return    
end

if objfun==6    % Do fine-tuning shift with ratio
    precision=0.1;  % step size in shifting an image
    croppixel=croppixel+1;  
    imbasesm=imbase(croppixel+1:end-croppixel, croppixel+1:end-croppixel);
%     pixofftemp=pixoff;
    tempx=pixoff(1)-0.6:precision:pixoff(1)+0.6;
    tempy=pixoff(2)-0.6:precision:pixoff(2)+0.6;
    dr=zeros(size(tempx,2)*size(tempy,2),3);
    t=1;
    for xstep=pixoff(1)-0.6:precision:pixoff(1)+0.6    % x shifting         
        for ystep=pixoff(2)-0.6:precision:pixoff(2)+0.6   % y shifting
            imtemp=OIShift(imshift, xstep, ystep);     
            imtempsm=imtemp(croppixel+1:end-croppixel, croppixel+1:end-croppixel);
            r=imtempsm./imbasesm;
            dr(t,:)=[std(r(:)),xstep,ystep];
            t=t+1;
        end
    end
    drsort=sortrows(dr);
%    currentobj=corr2(imbasesm, imtempsm);  % no more calculation of currentobj, use the old one (in last coarse shift)
    bestobj=drsort(1,1);
    pixoff=drsort(1,2:3);
%	imresult=OIShift(imshift, pixoff(1), pixoff(2));
    return;
end

printf('Please check shift objfun!\r');
return;