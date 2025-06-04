function [probe,rawCh,imDat] = removeBadTimesFromSpec(monkeyName,expDate,runName,probe,rawCh,imDat)
% This function removes the bad times determined from the spectrogram
% This function has to be updated as and when new data is recorded
if isempty(rawCh); rawCh = probe; end % Temporary fix
if isempty(imDat); imDat = zeros(1,round(size(probe,1)/1e2)); end 

if strcmp(monkeyName,'CharlieSheen') % Charlie Sheen
    if strcmp(expDate,'01_11_2022') && strcmp(runName,'run03') % Remove last 150 s of data
        probe(end-150e3+1:end,:) = [];
        rawCh(end-150e3+1:end,:) = [];
        imDat(:,end-1500+1:end) = [];

    elseif strcmp(expDate,'01_11_2022') && strcmp(runName,'run04') % Remove 380-600 s data
        probe(380e3+1:600e3,:) = [];
        rawCh(380e3+1:600e3,:) = [];
        imDat(:,3800+1:6000) = [];
    end

elseif strcmp(monkeyName,'Whiskey') % Whiskey
    if strcmp(expDate,'08_14_2023') && strcmp(runName,'run04') % Remove 320 - 410 s of data
        probe(320e3+1:410e3,:) = [];
        rawCh(320e3+1:410e3,:) = [];
        imDat(:,3200+1:4100) = [];

    elseif strcmp(expDate,'10_16_2023') && strcmp(runName,'run06') % Remove 430 - 510 s of data
        probe(430e3+1:510e3,:) = [];
        rawCh(430e3+1:510e3,:) = [];
        imDat(:,4300+1:5100) = [];

    elseif strcmp(expDate,'12_04_2023') && strcmp(runName,'run04') % Remove 510 - 630 s of data
        probe(510e3+1:630e3,:) = [];
        rawCh(510e3+1:630e3,:) = [];
        imDat(:,5100+1:6300) = [];

    elseif strcmp(expDate,'12_04_2023') && strcmp(runName,'run05') % Remove 750 - 800 s of data
        probe(750e3+1:800e3,:) = [];
        rawCh(750e3+1:800e3,:) = [];
        imDat(:,7500+1:8000) = [];

    elseif strcmp(expDate,'02_20_2024') && strcmp(runName,'run01') % Remove last 100 s of data
        probe(end-100e3+1:end,:) = [];
        rawCh(end-100e3+1:end,:) = [];
        imDat(:,end-1000+1:end)= [];

    elseif strcmp(expDate,'02_20_2024') && strcmp(runName,'run03') % Remove 1-100; 580-650 s of data
        probe([1:100e3 580e3+1:650e3],:) = [];
        rawCh([1:100e3 580e3+1:650e3],:) = [];
        imDat(:,[1:1000 5800+1:6500])  = [];

    elseif strcmp(expDate,'02_20_2024') && strcmp(runName,'run05') % Remove 500 - 700 s of data
        probe(500e3+1:700e3,:) = [];
        rawCh(500e3+1:700e3,:) = [];
        imDat(:,5000+1:7000) = [];

    elseif strcmp(expDate,'04_29_2024') && strcmp(runName,'run01') % Remove 250 - 350 s of data
        probe(680e3+1:end,:)  = [];
        rawCh(680e3+1:end,:)  = [];
        imDat(:,6800+1:end) = [];

    elseif strcmp(expDate,'04_29_2024') && strcmp(runName,'run05') % Remove 250 - 350 s of data
        probe(1:350e3,:) = [];
        rawCh(1:350e3,:) = [];
        imDat(:,1:3500) = [];
    end
end
end