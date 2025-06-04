function [lagHigh, lagLow] = plotLagProfiles(dataDir,bandName,crossCorr,allLags,monkeyName,expDate,runName,clipMaskROI)
% Functions to plot the following -
% 1.  XCORR: Median xcorr vs lag
% 2.  XCORR: Find the frame with peak positive and peak negative correlations

clear x xNew negIdx vidNew xLow;
x = allLags;
negIdx = x<0 & x>=-150; xNew = x(negIdx);
lowIdx = x<0 & x>= -80; xLow = x(lowIdx);
imSize = size(crossCorr,[2 3]);
grayIm = cat(3, 0.25.*ones((imSize)),0.25.*ones((imSize)), 0.25.*ones((imSize)));

if ~exist(dataDir,'dir'); [~,~] = mkdir(dataDir); end 

% 1. XCORR: Median xcorr vs lag
if ~exist([dataDir '\xcorrVsLag' bandName '.png'],'file') 
    figure('Position',[400 400 475 600]);
    plot(x,movmean(median(crossCorr,[2,3],'omitnan'),3));
    xlabel('Lag (s)'); ylabel('cross correlation');   ylim([-0.4 0.4]);
    legend('Median xcorr','Location','northeast','Autoupdate','off');
    xline(0); xticks(-200:40:200); xticklabels(-20:4:20); grid on; xlim([-200 200]);
    title(strrep([ bandName ' xcorr vs lag for ' monkeyName ' Date: ' expDate ' - Run: ' runName(end)],'_','\_'));
    f = gcf; exportgraphics(f,[dataDir,'\xcorrVsLag' bandName '.png'],'Resolution',300); close gcf;
end

% 2. XCORR: Find the frame with peak positive and peak negative correlations
clear maxMedcorrInd maxMedcorrIndLow frameNumMed frameNumMedLow
[~,maxMedcorrInd] = max(median(crossCorr(negIdx,:,:),[2,3],'omitnan'));
lagHigh           = xNew(maxMedcorrInd);
frameNumMed       = (x == lagHigh);

[~,maxMedcorrIndLow] = min(median(crossCorr(lowIdx,:,:),[2,3],'omitnan'));
lagLow               = xLow(maxMedcorrIndLow);
frameNumMedLow       = (x == lagLow);

% Save the frame where the highest correlation occurred
if ~exist([dataDir '\highestXCorrFrames_' bandName '.png'],'file') 
    figure('units','normalized','outerposition',[0 0 1 1]); tiledlayout(1,2);
    clear ax1; ax1 = nexttile;
    imagesc(squeeze(crossCorr(frameNumMed,:,:))); axis image off; hold on;
    h = imagesc(grayIm); hold off;
    set(h,'AlphaData',~clipMaskROI);
    colormap jet; colorbar; caxis([-0.5 0.5]);
    title(['Peak Positive median correlation at Lag: ' num2str(xNew(maxMedcorrInd)./10) 's ']);

    nexttile; imagesc(squeeze(crossCorr(frameNumMedLow,:,:))); axis image off; hold on;
    h = imagesc(grayIm); hold off; set(h,'AlphaData',~clipMaskROI);
    colormap(flipud(jet)); colorbar; caxis([-0.5 0.5]);
    title(['Peak Negative median correlation at Lag: ' num2str(xLow(maxMedcorrIndLow)./10) 's ']);
    colormap(ax1,'jet');

    sgtitle(strrep([ bandName ' xcorr ' monkeyName ' Date: ' expDate ...
        ' - Run: ' runName(end)],'_','\_'));

    f = gcf; exportgraphics(f,[dataDir '\highestXCorrFrames_' bandName '.png'],'Resolution',300); close gcf;
end

end
