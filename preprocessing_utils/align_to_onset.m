function [onsets, trigTout, songTout] = align_to_onset(audio, trigT, songT, daq_fs)

    audioLP = bandpass(audio, [2000 6000], daq_fs);

    audioW = .1;
    syllableOnsets = zeros(2*audioW*daq_fs+1, numel(trigT));
    for i = 1:numel(trigT)
        syllableOnsets(:, i) = audioLP((-audioW*daq_fs:audioW*daq_fs)+round(trigT(i)*daq_fs));
    end
    

    loudness = zeros(size(syllableOnsets));
    loudW = .005;
    tw = triang(2*round(loudW*daq_fs)+1);
    onsetReflect = cat(1, flipud(syllableOnsets), syllableOnsets, flipud(syllableOnsets));
    Ni = size(loudness, 2);
    Nj = size(loudness, 1);
    parfor i = 1:Ni
        for j = 1:Nj
            loudness(j, i) = sum((tw.*onsetReflect(size(syllableOnsets, 1)+ (j-loudW*daq_fs:j+loudW*daq_fs), i) ).^2);
        end
    end


    imagesc(log(loudness)');
    ax = gca;
    disp("select region for onset detection, then press any key to update");
    roi = drawrectangle(ax);
    pause;
    position = roi.Position;
    xmin = round(position(1));
    xmax = floor(position(1) + position(3));

    loudAlign = log(loudness(round(xmin):round(xmax), :));
    [N, edges] = histcounts(loudAlign(:), 100);
    thresh = min(edges)+(max(edges)-min(edges))*otsuthresh(N);

    onsetDelay = zeros(numel(trigT), 1);
    for i = 1:numel(trigT)
        if sum(diff(loudAlign(:, i)>thresh)>0)>0
            onsetDelay(i) = find(diff(loudAlign(:, i)>thresh)>0, 1)/daq_fs;
        end
    end

    goodI = (onsetDelay>0);
    onsetOffset = (xmin/daq_fs - audioW) + onsetDelay(goodI);
    trigTout = trigT(goodI);
    onsets = trigTout+onsetOffset';
    songTout = songT(goodI);


    %%plot aligned loudness

    syllableOnsetsAligned = zeros(2*audioW*daq_fs+1, numel(onsets));
    for i = 1:numel(onsets)
        syllableOnsetsAligned(:, i) = audioLP((-audioW*daq_fs:audioW*daq_fs)+round(onsets(i)*daq_fs));
    end
    

    loudnessAligned = zeros(size(syllableOnsetsAligned));
    loudW = .005;
    tw = triang(2*round(loudW*daq_fs)+1);
    onsetReflect = cat(1, flipud(syllableOnsetsAligned), syllableOnsetsAligned, flipud(syllableOnsetsAligned));
    Ni = size(loudnessAligned, 2);
    Nj = size(loudnessAligned, 1);
    parfor i = 1:Ni
        for j = 1:Nj
            loudnessAligned(j, i) = sum((tw.*onsetReflect(size(syllableOnsetsAligned, 1)+ (j-loudW*daq_fs:j+loudW*daq_fs), i) ).^2);
        end
    end


    imagesc(log(loudnessAligned)');
    shg;


end