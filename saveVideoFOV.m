function saveVideoFOV(dataDir,bandName,monkeyName,expDate,runName,crossCorrFOV,clipMaskCortex,allLags)
% This function saves the hybrid maps for all lags as a video

if ~exist(fullfile([dataDir,'\XCorrFOV_' bandName 'Video.avi']),'file')
    imSize     = size(crossCorrFOV,[2,3]);
    grayImFull = cat(3, 0.25.*ones(imSize),0.25.*ones(imSize), 0.25.*ones(imSize));

    v = VideoWriter([dataDir,'\XCorrFOV_' bandName 'Video']);
    open(v); % Open and write video files

    lagInd  = find(allLags>=-150 & allLags<=150);
    lagVals = allLags(allLags>=-150 & allLags<=150);

    % Save each hybrid map as a video frame
    for iLag = 1:length(lagVals)
        imagesc(squeeze(crossCorrFOV(lagInd(iLag),:,:)));
        hold on; h = imagesc(grayImFull); hold off; set(h,'AlphaData',~clipMaskCortex);
        axis image off; colormap(flipud(jet)); caxis([-0.5 0.5]); colorbar;

        title(strrep([ bandName ' xcorr ' monkeyName ' Date: ' expDate ...
            ' - Run: ' runName(end) ' Lag at: ' num2str(lagVals(iLag)./10) 's'],'_','\_'));

        pause(0.1); frame = getframe(gcf); writeVideo(v,frame);
    end

    close(v); close gcf; % Close video file after writing
end
end
