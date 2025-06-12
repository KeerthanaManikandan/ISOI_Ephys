function J = HROIGaussianFilter_new2(I, kerneldiameter, mask)
% avoid error caused by margin effect with Cheng Xu's modification
% (there is marker where edited)

% A disk filter that is better than simple convolution
% original disk filter: B = conv2(A, fspecial('disk', floor(Kernel/2)), 'same');
%   which has border problem (padding zeros beyond border lowered the pixel
%   values around border in the matrix because of averaging).
% 'OItDiskFilter' compensate this effect by multiply the result matrix by a coefficient 
% matrix which is derived from convolving a all "1" matrix of the same size. 
% This method is also used for subtracting the effect of specified regions in
% the input matrix (specified in the 'mask' matrix) 
% Note that if a mask is provided, there are two ways to compute the output pixels 
%    within masked region (same as input, or filted), which is controled by variable "maskoutsame"
% Thanks to GC for the inspiration. 		-- HDL 061115

% input:
% I: a 2D matrix
% kerneldiameter, in pixel, if is <=1, return the original matrix
% mask: a 0,1 matrix of same size as I. pixel value "1" means usful pixel
% output:
% J: filtered matrix


% Note J has the same mean as I
% Note: used disk

maskoutsame=0;	 % "1": keep the masked region the same in output matrix as in input matrix,
% "0" for filtered (takes time)

if kerneldiameter<=1
    J=I;
    return;
end
if nargin==2
    mask=ones(size(I));    % a mask to eliminate edge effect
else
    mask=double(mask);
end
Kernel=kerneldiameter;
cfilter=fspecial('gaussian', Kernel,floor(Kernel));
B = I.*mask;    % first, set all masked pixel 0, so their influence is same
B = conv2(B, cfilter, 'same');  % normal convlution
coef = conv2(mask, cfilter, 'same');   % then calculate mask&edge influence
coef=coef+(coef==0);    % make sure there is no 0 pixels (these pixel should be in mask region, doesn't matter, will be replaced later)
B=B./coef;              % subtract out the influence (
if ~maskoutsame
    C=conv2(I, cfilter, 'same');	% this map is used for masked-out region 
   %%%xc added below
   coef2 = conv2(ones(size(I)), cfilter, 'same');   % then calculate mask&edge influence
   C=C./coef2;
   clear coef2
   %%%xc added above
    J=B.*mask+C.*(1-mask);  % so mask region is also filtered    
else
    J=B.*mask+I.*(1-mask);  % so mask region remain the same value. 
end
J=J*sum(I(:))./sum(J(:));    % to keep the mean same


% coef = conv2(ones(size(I)), cfilter, 'same');   % then calculate mask&edge influence
% factor=conv2(mask, cfilter, 'same');
% factor=factor+(factor==0);
% J=conv2(I, cfilter, 'same')./factor-conv2((1-mask).*I, cfilter, 'same');%.*conv2(1-mask, cfilter, 'same');
return;