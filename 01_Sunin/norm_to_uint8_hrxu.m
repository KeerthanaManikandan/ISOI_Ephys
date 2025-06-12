% normalize image I to [0,255] then convert it to uint8
% norm to unit8 hrxu
% updated by Haoran Xu
% clip on the range of both map contents within the mask.

function J = norm_to_uint8_hrxu(I1, mask1, I2, mask2)

M = max(max(I1(:)),max(I2(:)));         % find the maximum intensity
m = min(min(I1(mask1~=0)),min(I2(mask2~=0)));         % find the minimum intensity
for i=1:504
    for j=1:504
        if (mask1(i,j))
            I1(i,j) = floor(255*(I1(i,j)-m)/(M - m)); % normalize to [0-255] 
        elseif (mask1(i,j)==0)
            I1(i,j) = I1(i,j);
        end
    end
end

% I = reshape(Iprime, 504, 504);
J = uint8(I1);                 % convert to uint8

return
