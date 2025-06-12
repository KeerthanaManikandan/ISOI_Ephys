% Make polar map by calling OIPolar, also output a color legend
% Operations in following order:
% 1) rotate (if specified)
% 2) corp (if specified)
% 3) resize (if specified)
% 4) Polar (include filtering)

clear all;
imgtype='ivf';
mask=ones(504); 	%imread('greenmask2.bmp'); %			% specify a blood vessel mask if available, otherwise use ones(504)
vectorclipsd=1;       % used in OIPolar, default 1
clipsd=1;             % for output map
LPMethod='gaussian';	% low pass method, default 'gaussian'
LPFKernel=10;			% low-pass kernel size, usually 5-10 for 504x504 image, also depends on domain size
HPMethod='oidisk';		% high-pass method, default 'gaussian', also can be 'slowmean', 'ribot', 'oidisk'...
HPFKernel=100;			% high-pass kernel size, usually 80-120

rotateangle=0;     % counter-clockwise rotate the image (in degree)

cropleft=86;   % x1    for crop of lateral domain after 177 rotation
croptop=48;    % y1
cropright=162;  % x2
cropbottom=106; % y2

% cropleft=345;   % x1  for crop of lateral domain from original map
% croptop=375;    % y1
% cropright=455;  % x2
% cropbottom=470; % y2

% cropleft=1;   % x1    for no crop
% croptop=1;    % y1
% cropright=504;  % x2
% cropbottom=504; % y2

magnifysize=6;  % enlarge or reduce map size, use imresize(A, magnifysize, 'bilinear', 0);
saveivf=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempfilename=struct2cell(dir(['*.', imgtype]));
filename=sort(tempfilename(1,:)')
size(filename);
Nvect=size(filename,1);
lut=textread('bwfix.lut');  % this color table should be in sunin folder
mask=mask(croptop:cropbottom, cropleft:cropright);
mask=double(imresize(mask, magnifysize, 'bilinear', 0)>0.5);
for j=1:Nvect
	imgname=getfield(cell2struct(filename(j), 'junk'), 'junk')
    maptemp=-OIReadIVF(imgname);   % note: signal reversed because of negative response
    maptemp=imrotate(maptemp, rotateangle, 'bilinear');    
    maptemp=maptemp(croptop:cropbottom, cropleft:cropright);
    maptemp=imresize(maptemp, magnifysize, 'bilinear', 0);
    singlemap(:,:,j)=maptemp;
end

Polarmap=OIPolar(singlemap, lut, mask, vectorclipsd, LPMethod, LPFKernel, HPMethod, HPFKernel);
if saveivf
    OIWriteIVF(Polarmap.ang, 'ang.ivf');
    OIWriteIVF(Polarmap.mag, 'mag.ivf');
end
% output
imwrite(norm_to_uint8(Polarmap.ang), lut, 'ang.tiff');
imwrite(norm_to_uint8(OIClip(Polarmap.mag, 1,1)), 'mag.tiff');
mag=norm_to_01(Polarmap.mag);
mag=OIClip(mag, 1, clipsd);   % to adjust map darkness    
ang=double(OIColorMap(norm_to_uint8(Polarmap.ang), lut));
polarmap=ang;
polarmap(:,:,1)=ang(:,:,1).*mag;
polarmap(:,:,2)=ang(:,:,2).*mag;
polarmap(:,:,3)=ang(:,:,3).*mag;
imwrite(norm_to_uint8(polarmap), lut, 'polar.tiff');

% output two color tables (one contineous, one with bars)
OIDrawPolarIndex(Nvect, '', lut);

return;

