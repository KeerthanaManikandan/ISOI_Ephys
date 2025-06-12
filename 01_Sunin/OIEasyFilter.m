function B=OIEasyFilter(A, LPMethod, LPFKernel, HPMethod, HPFKernel, mask)
% function OIEasyFilter(A, LPFName, LPFKernel, HPFName, HPFKernel)
% kernel=0 means no filter
% LPFKernel and HPFKernel are all diameters in pixel unit

if nargin==5
    mask=ones(size(A));
    fprintf('default mask (all 1) used\r');
end
B=A;
if LPFKernel~=0
    switch LPMethod
        case 'oidisk'
            B = OIDiskFilter(A, LPFKernel, mask);
        case 'fastmean'
            B = conv2(A, fspecial('disk', floor(LPFKernel/2)), 'same');
        case 'slowmean'
            B = OIMeanFilt(A, LPFKernel);
        case 'gaussian'
            B = conv2(A, fspecial('gaussian', LPFKernel, floor(LPFKernel)), 'same');
        case 'fastmedian'
            B = medfilt2(A, [LPFKernel, LPFKernel], 'symmetric');
        otherwise
            fprintf('Error: please specify filting method, no filtering performed\r');
    end
end
if HPFKernel~=0
    switch HPMethod
        case 'oidisk'
            B = B-OIDiskFilter(A, HPFKernel, mask);
        case 'fastmean'
            B = B-conv2(A, fspecial('disk', floor(HPFKernel/2)), 'same');
        case 'slowmean'
            B = B-OIMeanFilt(A, HPFKernel);
        case 'gaussian'
            B = B-conv2(A, fspecial('gaussian', HPFKernel, floor(HPFKernel)), 'same');
        case 'fastmedian'
            B = B-medfilt2(A, [HPFKernel, HPFKernel], 'symmetric');
        case 'ribot'
            B = B-RibotFilter(A, HPFKernel);   % note HPFKernel here is actually order of polynomial function, usually 2 or 3
        otherwise
            fprintf('Error: please specify filting method, no filtering performed\r');
    end
end
return