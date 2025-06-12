function J = HROIMeanFilt(I, kerneldiameter, skippix);

% J = OIMeanFilt(I, kerneldiameter, skippix);
% Returns same size matrix with average value within a kernel size (smoothing)
% I:  Input image (Float point format), 
% kerneldiameter: integer to decide the radius size (diameter), e.g. 3 X 3, 5 X 5.  
% J: result image (Float point format)
% skippix: for fasten calculation, =2 means skip every other pixel,

% modified from Xiangming Xu 04 'FloatMeanfiltering'

if nargin<2
	fprintf('Too few input arguments\n');
elseif nargin<3
	skippix=1;
end
[r,c] = size(I);
J = zeros(r,c);
radius = floor(kerneldiameter/2)
if radius==0
	J=I;
return;
end

for i = 1:r
%	if (~mod(i, 10))			% for display filting progress
%		fprintf('%d  ', i);
%	end
%	if (~mod(i, 100))
%		fprintf('\r');
%	end
    for j = 1:c
    	if (i-radius)<1 | (j-radius)<1 |(i+radius)>r |(j+radius)>c
	        numpix=0;
    	    sumpixval=0;
        	for padr = (i - radius):skippix:(i + radius)
            	for  padc = (j -radius):skippix:(j + radius)
                	if ~((padr <1) | (padc <1) | (padr > r) | (padc > c) )
                    	numpix = numpix + 1;
	                    sumpixval = sumpixval + I(padr, padc);
    	            end
        	    end
        	end	
	        J(i, j) = sumpixval/numpix;
        else
        	J(i, j) = mean2(I(i-radius:i+radius, j-radius:j+radius));
        end
    end  
%	fprintf('\r');
end
return;