function pixarray=OIQmask1(data, mask)
% Get individual pixel values from matrix 'data' using 'mask'
% retrun 'pixarray' is a 1-D array containing all pixels in the mask
% If the same mask was used, the order of pixels in pixarray should be the same

% If you don't care about the values of each pixel, you can just use "mean2(data.*mask)".
% If your mask contain more than one domain and you care about the value in each domain, use 'OIQmask2'
% If your mask is coordinates + radius, use 'OIQmask3'

% input
%   data: a 2-D matrix
%   mask: 0/1 mask, same size as data

% output:
%   'pixarray', size = sum(sum(mask))

imgsize=size(data);
masksize=sum(sum(double(mask)));

maskarray=reshape(mask, prod(imgsize),1);
[junk, maskindex]=sort(maskarray);		% maskarray contains (0's,1's), and maskindex contain location of 1's
oneindex=maskindex(end-masksize+1:end); % all location of 1's

dataarray=reshape(data, prod(imgsize), 1);
pixarray=dataarray(oneindex);

return;
