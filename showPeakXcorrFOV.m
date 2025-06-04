function showPeakXcorrFOV(dataDir,bandName,monkeyName,expDate,runName,greenFig,crossCorrFOV,clipMaskCortex,frameNumHighTime,frameNumLowTime,allLags)
% This function plots the following -
% 1. Correlation map at peak positive and negative
% 2. Significant pixels of correlation map on blood vessel map

imSize = size(crossCorrFOV,[2,3]);
grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));

frameNumHigh = find(allLags == frameNumHighTime);
frameNumLow  = find(allLags == frameNumLowTime);

% 1. Correlation map at peak positive and negative
if ~exist([dataDir '\XCorrFOV_' bandName '.png'],'file')
    figure('units','normalized','outerposition',[0 0 1 1]); tiledlayout(1,2);
    clear ax1; ax1 = nexttile;
    imagesc(squeeze(crossCorrFOV(frameNumHigh,:,:))); axis image off; hold on;
    h = imagesc(grayImFull); hold off;
    set(h,'AlphaData',~clipMaskCortex);
    colormap(ax1,'jet'); colorbar; caxis([-0.5 0.5]);
    title(strrep([' Peak Positive: ' num2str(frameNumHigh/10)],'_','\_'));

    ax2 = nexttile; imagesc(squeeze(crossCorrFOV(frameNumLow,:,:))); axis image off; hold on;
    h = imagesc(grayImFull); hold off;
    set(h,'AlphaData',~clipMaskCortex);
    colormap(ax2,flipud(jet)); colorbar; caxis([-0.5 0.5]);
    title(strrep([' Peak Negative: ' num2str(allLags(frameNumLow)/10)],'_','\_'));

    colormap(ax1,'jet');
    sgtitle(strrep([ bandName ' xcorr ' monkeyName ' Date: ' expDate ...
        ' - Run: ' runName(end)],'_','\_'));

    f = gcf; exportgraphics(f,[dataDir,'\XCorrFOV_' bandName '.png'],'Resolution',300); close gcf;
end

% 2. Significant pixels of correlation map on blood vessel map
if ~exist([dataDir '\greenCorrFOV_' bandName '.png'],'file')
    clear figFrameHigh figFrameLow
    figFrameHigh = squeeze(crossCorrFOV(frameNumHigh,:,:));
    figFrameLow  = squeeze(crossCorrFOV(frameNumLow,:,:));

    figure('units','normalized','outerposition',[0 0 1 1]); tiledlayout(1,2);
    clear ax1; ax1 = nexttile;
    imagesc(ind2rgb(greenFig,gray(256))); hold on; axis image off;
    imagesc(figFrameHigh,'alphaData',figFrameHigh.*5);
    colormap(ax1,'jet'); caxis([-0.5 0.5]); colorbar;
    title(strrep([' Peak Positive: ' num2str(allLags(frameNumLow)/10)],'_','\_'));

    ax2 = nexttile;  imagesc(ind2rgb(greenFig,gray(256))); hold on; axis image off;
    imagesc(figFrameLow,'alphaData',figFrameLow.*-5);
    colormap(ax2,flipud(jet)); colorbar; caxis([-0.5 0.5]);
    title(strrep([' Peak Negative: ' num2str(frameNumLow/10)],'_','\_'));
    colormap(ax1,'jet');

    sgtitle(strrep([ bandName ' xcorr ' monkeyName ' Date: ' expDate ...
        ' - Run: ' runName(end)],'_','\_'));

    f = gcf; exportgraphics(f,[dataDir,'\greenCorrFOV_' bandName '.png'],'Resolution',300); close gcf;
end

end
