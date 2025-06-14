%% Get CSD for an example recording, plot the spectrogram and get the cross-modal maps
% run ephysImaging_compilation_v5 in order to get all data for the monkey

iDate = 2; iRun = 7; 
clear green expDate runName clipMask elecMask green dataDir probeCh rawCh ...
    imSize corrMask bThresh imData imData10 infraEphys crossCorrFOV allXcorr ...
    superXcorr deepXcorr midXCorr

fs = 1e3;
gammaBand = [30 90]; [bG,aG] = butter(3,gammaBand./(fs/2),'bandpass'); % Gamma band filtering parameters
alphaBand = [8 12];  [bA,aA] = butter(3,alphaBand./(fs/2),'bandpass'); % Alpha band filtering parameters
betaBand  = [13 30]; [bB,aB] = butter(3,betaBand./(fs/2),'bandpass');  % Beta band filtering parameters
thetaBand = [6 8];   [bT,aT] = butter(3,thetaBand./(fs/2),'bandpass'); % Theta band filtering parameters

params.Fs     = fs;
params.fpass  = [1 120];
params.pad    = -1;
params.tapers = [3 5];

expDate = allDates(iDate,:);
runName = allRuns{iDate,1}(iRun,:);
dataDir = ['X:\Data\' monkeyName '_SqM\' hemisphere ' Hemisphere\' expDate '\' runName ];

if exist([dataDir '\clipMask0' runName(end) '.BMP'],'file') == 0
    clipMask = imread([dataDir '\clipMask0' runName(end) '.png']);
else
    clipMask = imread([dataDir '\clipMask0' runName(end) '.bmp']);
end

if exist([dataDir '\clipMask0' runName(end) '.BMP'],'file') == 0
    green = imread([dataDir '\green0' runName(end) '_Edited.png']);
else
    green = imread([dataDir '\green0' runName(end) '_Edited.bmp']);
end

green       = imresize(green(:,:,1),1/3); % Resize mask

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
imSize         = size(clipMask);

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

probeCh    = probe{iRun,iDate}.probeCh;
rawCh      = probe{iRun,iDate}.rawCh;
timeStamp  = probe{iRun,iDate}.timeStamp;
badTimes   = badTimesLFP{iDate,iRun};
chInCortex = estChInCortex{1,iDate}(iRun,:);
szLFP      = size(probeCh,1);

imData     = processedDat{iDate,iRun}.tempBandPass;
imSize     = size(imData);
imData     = reshape(imData,[imSize(1)*imSize(2) imSize(3)]);
bThresh    = badTimeThresh{iDate,iRun};

% Upsample data 
clear imData10
parfor iP = 1:size(imData,1)
    imData10(iP,:) = interp(imData(iP,:),5);
end

szIm = size(imData10,2)*100;
if ~(szLFP == szIm)
    szMin    = min([szLFP, szIm]);
    probeCh  = probeCh(1:szMin,:);
    rawCh    = rawCh(1:szMin,:);
    imData10 = imData10(:,1:floor(szMin/100));

    badTimes(badTimes>szMin) = [];
    bThresh(bThresh>szMin) = [];
end

% Removing bad channels
probeCh(:,badCh{iDate,iRun}) = [];
rawCh(:,badCh{iDate,iRun})   = [];

% Remove bad times
timeStampSorted = timeStamp- timeStamp(1);
badTimes10Hz    = unique(bThresh./1000);
badTimeIm       = [];

% Identifying frames to be removed from RS-ISOI
for iT = 1: length(badTimes10Hz)
    badTimeIm(iT) = find((floor(abs(timeStampSorted - badTimes10Hz(iT))*100)./100)<=0.05,1,'first'); %#ok<*AGROW> 
end

badTimeIm = unique(badTimeIm);
badTimeIm(badTimeIm>size(imData10,2)) = [];
imData10(:,badTimeIm) = [];
probeCh(badTimes,:) = [];
rawCh(badTimes,:)   = [];

% Remove channels out of cortex

if chInCortex(1)-chInCortex(2) ~= 0
    clear probeChCortex spec aLim bLim gLim spec timeValsSpec freqValsSpec
    probeChCortex = probeCh(:,chInCortex(1):chInCortex(2));
    
    % Bipolar reference
    probeRef = -diff(probeChCortex')';

    % CSD reference
    probeCSD = diff(diff(probeChCortex'))';
    figure;
    subplot(131); imagesc(mean(corr(probeChCortex),1,'omitnan')'); colormap jet; title('No ref');caxis([0 1]); colorbar;
    subplot(133); imagesc(mean(corr(probeRef),1,'omitnan')'); colormap jet; title('Bipolar ref');caxis([0 0.5]);colorbar;
    subplot(132); imagesc(mean(corr(probeCSD),1,'omitnan')'); colormap jet; title('CSD');caxis([0 0.5]); colorbar;
    sgtitle('Wideband');

    for iBand = 1:4
        switch iBand
            case 1
                b = bT; a = aT;  
                bandLabel = 'Theta';
            case 2
                b = bA; a = aA;
                bandLabel = 'Alpha';
            case 3
                b = bB; a = aB;
                bandLabel = 'Beta';
            case 4
                b = bG; a = aG;
                bandLabel = 'Gamma';
        end

        figure;
        subplot(131); imagesc(mean(corr(single(filtfilt(b,a,double(probeChCortex)))),1,'omitnan')'); colormap jet; title('No ref');caxis([0 1]); colorbar;
        subplot(133); imagesc(mean(corr(single(filtfilt(b,a,double(probeRef)))),1,'omitnan')'); colormap jet; title('Bipolar ref');caxis([0 0.5]);colorbar;
        subplot(132); imagesc(mean(corr(single(filtfilt(b,a,double(probeCSD)))),1,'omitnan')'); colormap jet; title('CSD');caxis([0 0.5]); colorbar;
        sgtitle(bandLabel);
    end

    % Get cross-modal maps for the LFP
    serverDir = ['\\smb2.neurobio.pitt.edu\Gharbawie\Lab\kem294\Data\' ...
        monkeyName '_SqM\' hemisphere ' Hemisphere\'  expDate '\' runName ];

    crossCorrFOV = matfile([serverDir '\crossCorrFOV_6_NoRef.mat']);
    allXcorr     = crossCorrFOV.spatialProfile;
    superXcorr   = crossCorrFOV.spatialProfileSuper;
    deepXcorr    = crossCorrFOV.spatialProfileDeep;
    midXCorr     = crossCorrFOV.spatialProfileMid;

    % Compute cross-corr maps for CSD and bipolar referencing
    clear infraEphys infraSuper infraMid infraDeep infraSuperCSD infraDeepCSD infraMidCSD...
        ccFullSuper ccFullSuperCSD ccFullMid ccFullMidCSD ccFullDeep ccFullDeepCSD
    for iBand = 1:4
        switch iBand
            case 1
                b = bT; a = aT;
                super = superXcorr.ccFullTheta;
                mid   = midXCorr.ccFullTheta;
                deep  = deepXcorr.ccFullTheta;
            case 2
                b = bA; a = aA;
                super = superXcorr.ccFullAlpha;
                mid   = midXCorr.ccFullAlpha;
                deep  = deepXcorr.ccFullAlpha;
            case 3
                b = bB; a = aB;
                super = superXcorr.ccFullBeta;
                mid   = midXCorr.ccFullBeta;
                deep  = deepXcorr.ccFullBeta;
            case 4
                b = bG; a = aG;
                super = superXcorr.ccFull;
                mid   = midXCorr.ccFull;
                deep  = deepXcorr.ccFull;
        end
        % Get infraslow bandlimited powers for different frequencies
        infraEphys = getInfraSlowPowerLFP(probeRef,b,a,1:19);
        infraSuper = mean(infraEphys(:,1:6),2,'omitnan');
        infraDeep  = mean(infraEphys(:,13:end),2,'omitnan'); 
        infraMid   = mean(infraEphys(:,7:12),2,'omitnan');

        minSize    = min(size(imData10,2),size(infraEphys,1));
        infraSuper = infraSuper(1:minSize);
        infraDeep  = infraDeep(1:minSize);
        infraMid   = infraMid(1:minSize);
        imData10   = imData10(:,1:minSize); 

        infraEphysCSD = getInfraSlowPowerLFP(probeCSD,b,a,1:19);
        infraSuperCSD = mean(infraEphysCSD(:,1:6),2,'omitnan'); infraSuperCSD = infraSuperCSD(1:minSize);
        infraDeepCSD  = mean(infraEphysCSD(:,13:end),2,'omitnan'); infraDeepCSD = infraDeepCSD(1:minSize);
        infraMidCSD   = mean(infraEphysCSD(:,7:12),2,'omitnan'); infraMidCSD = infraMidCSD(1:minSize);

        parfor iP = 1:size(imData10,1)
            [ccFullSuper(:,iP,iBand),~] = xcorr(infraSuper',imData10(iP,:),200,'normalized');
            [ccFullMid(:,iP,iBand),~]   = xcorr(infraMid',imData10(iP,:),200,'normalized');
            [ccFullDeep(:,iP,iBand),~]  = xcorr(infraDeep',imData10(iP,:),200,'normalized');

            [ccFullSuperCSD(:,iP,iBand),~] = xcorr(infraSuperCSD',imData10(iP,:),200,'normalized');
            [ccFullMidCSD(:,iP,iBand),~]   = xcorr(infraMidCSD',imData10(iP,:),200,'normalized');
            [ccFullDeepCSD(:,iP,iBand),~]  = xcorr(infraDeepCSD',imData10(iP,:),200,'normalized');
        end
    end

    ccFullSuper = reshape(ccFullSuper,[401 361 438 4]);
    ccFullMid = reshape(ccFullMid,[401 361 438 4]);
    ccFullDeep = reshape(ccFullDeep,[401 361 438 4]);

    ccFullSuperCSD = reshape(ccFullSuperCSD,[401 361 438 4]);
    ccFullMidCSD = reshape(ccFullMidCSD,[401 361 438 4]);
    ccFullDeepCSD = reshape(ccFullDeepCSD,[401 361 438 4]);


    for iBand = 1:4
        switch iBand
            case 1
                super = superXcorr.ccFullTheta;
                mid   = midXCorr.ccFullTheta;
                deep  = deepXcorr.ccFullTheta;
                bandLabel = 'Theta';

            case 2
                super = superXcorr.ccFullAlpha;
                mid   = midXCorr.ccFullAlpha;
                deep  = deepXcorr.ccFullAlpha;
                bandLabel = 'Alpha';

            case 3
                b = bB; a = aB;
                super = superXcorr.ccFullBeta;
                mid   = midXCorr.ccFullBeta;
                deep  = deepXcorr.ccFullBeta;
                bandLabel = 'Beta';

            case 4
                super = superXcorr.ccFull;
                mid   = midXCorr.ccFull;
                deep  = deepXcorr.ccFull;
                bandLabel = 'Gamma';
        end

        x = -200:200;
        ind = find(x==-9);

        figure;
        subplot(331);
        imagesc(ind2rgb(green,gray(256))); hold on; axis image off;
        imagesc(squeeze(super(ind,:,:)),'AlphaData',squeeze(super(ind,:,:)).*-3);
        caxis([-0.6 0.2]); colormap(flipud(jet)); colorbar; title('Superficial - LFP');

        subplot(334);
        imagesc(ind2rgb(green,gray(256))); hold on; axis image off;
        imagesc(squeeze(mid(ind,:,:)),'AlphaData',squeeze(mid(ind,:,:)).*-3);
        caxis([-0.6 0.2]); colormap(flipud(jet)); colorbar; title('Middle - LFP');

        subplot(337);
        imagesc(ind2rgb(green,gray(256))); hold on; axis image off;
        imagesc(squeeze(deep(ind,:,:)),'AlphaData',squeeze(deep(ind,:,:)).*-3);
        caxis([-0.6 0.2]); colormap(flipud(jet)); colorbar;title('Deep - LFP');

        subplot(332);
        imagesc(ind2rgb(green,gray(256))); hold on; axis image off;
        imagesc(squeeze(ccFullSuperCSD (ind,:,:,iBand)),'AlphaData',squeeze(ccFullSuperCSD (ind,:,:,iBand)).*-3);
        caxis([-0.6 0.2]); colormap(flipud(jet)); colorbar; title('Superficial - CSD');

        subplot(335);
        imagesc(ind2rgb(green,gray(256))); hold on; axis image off;
        imagesc(squeeze(ccFullMidCSD (ind,:,:,iBand)),'AlphaData',squeeze(ccFullMidCSD (ind,:,:,iBand)).*-3);
        caxis([-0.6 0.2]); colormap(flipud(jet)); colorbar; title('Middle - CSD');

        subplot(338);
        imagesc(ind2rgb(green,gray(256))); hold on; axis image off;
        imagesc(squeeze(ccFullDeepCSD (ind,:,:,iBand)),'AlphaData',squeeze(ccFullDeepCSD (ind,:,:,iBand)).*-3);
        caxis([-0.6 0.2]); colormap(flipud(jet)); colorbar;title('Deep - CSD');

        subplot(333);
        imagesc(ind2rgb(green,gray(256))); hold on; axis image off;
        imagesc(squeeze(ccFullSuper(ind,:,:,iBand)),'AlphaData',squeeze(ccFullSuper(ind,:,:,iBand)).*-3);
        caxis([-0.6 0.2]); colormap(flipud(jet)); colorbar; title('Superficial - Bipolar');

        subplot(336);
        imagesc(ind2rgb(green,gray(256))); hold on; axis image off;
        imagesc(squeeze(ccFullMid(ind,:,:,iBand)),'AlphaData',squeeze(ccFullMid(ind,:,:,iBand)).*-3);
        caxis([-0.6 0.2]); colormap(flipud(jet)); colorbar; title('Middle - Bipolar');

        subplot(339);
        imagesc(ind2rgb(green,gray(256))); hold on; axis image off;
        imagesc(squeeze(ccFullDeep(ind,:,:,iBand)),'AlphaData',squeeze(ccFullDeep(ind,:,:,iBand)).*-3);
        caxis([-0.6 0.2]); colormap(flipud(jet)); colorbar;title('Deep - Bipolar');

        sgtitle(bandLabel);
    end
end
