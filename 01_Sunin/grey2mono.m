% 8 bit image to 1 bit image

% coded by Haoran Xu
% 2010-08-27

imgname='green0.bmp';
image=imread(imgname);

for i=1:504
    for j=1:504
        mapcombine1(i,j)=sum(image(i,j,:));
    end
end
mapcombine1=(mapcombine1~=0);
imwrite(mapcombine1, strcat(imgname(1:end-4), '-mono.bmp'));
