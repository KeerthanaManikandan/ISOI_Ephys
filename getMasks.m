function [clipMaskCortex,corrMask] = getMasks(dataDir,runName)
% This function loads the appropriate masks for a run 
% June 12, 2025 - KM

% Load clipmask
if exist([dataDir '\clipMask0' runName(end) '.BMP'],'file') == 0
    clipMask = imread([dataDir '\clipMask0' runName(end) '.png']);
else
    clipMask = imread([dataDir '\clipMask0' runName(end) '.bmp']);
end

% This mask includes vessels and electrode(s)
if exist([dataDir '\maskSkull0' runName(end) '.bmp'],'file') == 0
    elecMask = imread([dataDir '\maskSkull0' runName(end) '.png']);
else
    elecMask = imread([dataDir '\maskSkull0' runName(end) '.bmp']);
end

clipMask       = imresize(clipMask,1/3); % Resize mask
clipMask       = clipMask(:,:,1)>0; % Converting to 0s and 1s
elecMask       = imresize(elecMask-100,1/3); % Binarize image
elecMask       = elecMask(:,:,1)>0;
clipMaskCortex = clipMask & ~elecMask; % Mask of cortex sans vessels

% Mask for comparing the imaging and hybrid map (after removing
% areas beyond lateral sulcus)
if exist([dataDir '\corrMask0' runName(end) '.BMP'],'file')
    corrMask = imread([dataDir '\corrMask0' runName(end) '.bmp']);
    corrMaskFlag = 1;
elseif exist([dataDir '\corrMask0' runName(end) '.PNG'],'file')
    corrMask = imread([dataDir '\corrMask0' runName(end) '.png']);
    corrMaskFlag = 1;
else
    corrMaskFlag = 0;
end

if corrMaskFlag
    corrMask = imresize(corrMask,1/3); % Resize the mask
    corrMask = corrMask(:,:,1)>0;
    corrMask = corrMask & ~elecMask;
else
    corrMask = clipMaskCortex;
end

end