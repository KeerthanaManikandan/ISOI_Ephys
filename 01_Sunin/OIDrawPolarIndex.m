function OIDrawPolarIndex(NN, foldername, lut)
% draw color index for polar maps, containing orientation bars, and color spectrum 
% input: NN: how many input vectors
%        foldername: where to output maps
%        lut: color lookup table, usually is the 'bwfix.lut' in sunin folder

% define barsize, if size are odd number, will be change to the next even number
% and the final size is 1 pixel larger than the changed ones
barwidth=10;
barlength=35;
backgroundgray=255;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
manual=0;
if manual
    NN=8;
    foldername='';
    lut=textread('bwfix.lut');  % this color table should be in sunin folder
end

if mod(barwidth, 2)
    barwidth=barwidth+1;
end
if mod(barlength, 2)
    barlength=barlength+1;
end

canvas=ones(barlength+3, 1, 3)*backgroundgray;
space=ones(barlength+3,1,3)*backgroundgray;
for i=1:NN
    templet=ones(barlength+3, barlength+3, 3)*backgroundgray;
    center=barlength/2+2;
    x=floor((i-1)*256/NN)+1;
    templet(center-(barwidth/2):center+(barwidth/2), center-(barlength/2):center+(barlength/2), 1)=256*ones(barwidth+1, barlength+1)*lut(x,1);
    templet(center-(barwidth/2):center+(barwidth/2), center-(barlength/2):center+(barlength/2), 2)=256*ones(barwidth+1, barlength+1)*lut(x,2);
    templet(center-(barwidth/2):center+(barwidth/2), center-(barlength/2):center+(barlength/2), 3)=256*ones(barwidth+1, barlength+1)*lut(x,3);
    templet=imrotate(templet, (i-1)*180/NN, 'bilinear', 'crop');
    canvas=[canvas, templet, space];
end
% check if there is any black edge induced by rotation
[ch, cw, junk]=size(canvas);
if backgroundgray ~=0
    for yy=1:ch
        for xx=1:cw
            if canvas(yy, xx, 1)==canvas(yy, xx, 2) & canvas(yy, xx, 3)==canvas(yy, xx, 1)
                canvas(yy, xx, 1)=backgroundgray;
                canvas(yy, xx, 2)=backgroundgray;
                canvas(yy, xx, 3)=backgroundgray;
            end
        end
    end
end
colorspectrum=ones(barlength, cw, 3)*255;
for i=1:cw
    x=floor((i-1)*256/cw)+1;
    colorspectrum(:, i, 1)=256.*lut(x,1);
    colorspectrum(:, i, 2)=256.*lut(x,2);
    colorspectrum(:, i, 3)=256.*lut(x,3);
end
canvas=[canvas; colorspectrum];
imwrite(uint8(canvas), strcat(foldername, 'colortable.tiff'), 'tiff'); 

