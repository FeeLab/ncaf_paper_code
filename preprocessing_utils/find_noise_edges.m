function [noiseI, noiseF, isNoiseOut] = find_noise_edges(audio, songTonset, trigDelay, daq_fs, isNoise)


    audioLP = bandpass(audio, [10000 20000], daq_fs);

    audioW = .05;
    noiseOnsets = zeros(round(audioW*daq_fs), sum(isNoise));
    songTnoise = songTonset(isNoise);
    for i = 1:numel(songTnoise)
        noiseOnsets(:, i) = audioLP((0:size(noiseOnsets, 1)-1)+round((songTnoise(i)+trigDelay)*daq_fs));
    end

    loudness = zeros(size(noiseOnsets));
    loudW = .001;
    tw = triang(2*round(loudW*daq_fs)+1);
    onsetReflect = cat(1, flipud(noiseOnsets), noiseOnsets, flipud(noiseOnsets));
    Ni = size(loudness, 2);
    Nj = size(loudness, 1);
    parfor i = 1:Ni
        for j = 1:Nj
            loudness(j, i) = sum((tw.*onsetReflect(size(noiseOnsets, 1)+ (j-loudW*daq_fs:j+loudW*daq_fs), i) ).^2);
        end
    end

    [N, edges] = histcounts(loudness(:), 1000);
    noiseL = min(edges)+(max(edges)-min(edges))*otsuthresh(N);

    noiseI = [];
    noiseF = [];
    isNoiseOut = isNoise;
    noiseInd = find(isNoise);
    for i = 1:size(loudness, 2)
        if sum(loudness(:, i)>noiseL)==0
            isNoiseOut(noiseInd(i))=0;
        end
        noiseI = [noiseI; find(loudness(:, i)>noiseL, 1)/daq_fs+trigDelay];
        noiseF = [noiseF; find(loudness(:, i)>noiseL, 1, 'last')/daq_fs+trigDelay];
    end
end