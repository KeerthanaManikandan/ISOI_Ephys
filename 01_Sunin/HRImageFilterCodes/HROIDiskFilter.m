function J = HROIDiskFilter(I, kerneldiameter, mask)
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
%% Syntax
% input:
% I: a 2D matrix
% kerneldiameter, in pixel, if is <=1, return the original matrix
% mask: a 0,1 matrix of same size as I. pixel value "1" means usful pixel
% output:
% J: filtered matrix

% Note J has different mean wtih I
% Note: used disk
%% Revise
% add    coef2 = conv2(ones(size(I)), cfilter, 'same');  C=C./coef2; clear coef2
% this is to get rid of edge effect
% -- Cheng Xu 100531(above)

%1:change the value of maskoutsame 
% avoid filtering area which has masked out to save time,but if you are interested in that area ,you can change this value to 0
%2: change 'cfilter=fspecial('gaussian', Kernel,floor(Kernel));' to'cfilter=fspecial('gaussian', Kernel,floor(Kernel/2));'
% change the sigma value to increase the difference of weight between the middle and the edge
%3:add  'C = I.*(1-mask); 
% if  no this sentence,then we will filter the whole image, and the vessel part will be influenced by the other part
%4: coef2=coef2+(coef2==0);
% make sure there is no 0 pixels
%5:drop out 'J=J*sum(I(:))./sum(J(:));'
% this sentence is not necessary,and may get strange image if sum(I(:))./sum(J(:)) is not close to 1
% -- Cheng Xu 110212(above)

% HR 140514: 'J=J*sum(I(:))./sum(J(:));' may be not necessary. this version
% provides with same mean output without doing this.

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
Kernel=floor(kerneldiameter/2);
cfilter=fspecial('disk', Kernel);
B = I.*mask;    % first, set all masked pixel 0, so their influence is same
B = conv2(B, cfilter, 'same');  % normal convlution
coef = conv2(mask, cfilter, 'same');   % then calculate mask&edge influence
coef=coef+(coef==0);    % make sure there is no 0 pixels
% (these 0 pixels are in the mask region, and will be replaced later)
B=B./coef;              % subtract out the influence 
if ~maskoutsame
    C = I.*(1-mask);
    C = conv2(C, cfilter, 'same');
    coef2 = conv2((1-mask), cfilter, 'same');
    coef2=coef2+(coef2==0);
    C=C./coef2;
    clear coef2
    J=B.*mask+C.*(1-mask);
else
    J=B.*mask+I.*(1-mask);  % so mask region remain the same value.
end


% coef = conv2(ones(size(I)), cfilter, 'same');   % then calculate mask&edge influence
% factor=conv2(mask, cfilter, 'same');
% factor=factor+(factor==0);
% J=conv2(I, cfilter, 'same')./factor-conv2((1-mask).*I, cfilter, 'same');%.*conv2(1-mask, cfilter, 'same');
return;