function spikeCount = contingent_spikes(signal, np_fs, trigDelay, cW, normed)

    cT = round([trigDelay-cW trigDelay]*np_fs);

    spikeCount = squeeze(sum(signal(cT(1):cT(2),:,:), 1));
    if size(spikeCount, 1) == 1;
        spikeCount = spikeCount';
    end
    frateW = 0.01;
    binT = 100;

    marginW = 0.025;
    if normed
        for i = 1:size(signal, 3)
            frate = smoothdata(signal(:,:,i), 'gaussian', frateW*np_fs)*np_fs;
            binSig = smoothdata(frate, 2, 'gaussian', binT);
            spikeCount(:, i) = spikeCount(:, i)./median(binSig([1:round((trigDelay-cW-marginW)*np_fs) round((trigDelay+marginW)*np_fs):end], :), 1)';
            spikeCount(:, i) = spikeCount(:, i)*median(binSig([1:round((trigDelay-cW-marginW)*np_fs) round((trigDelay+marginW)*np_fs):end], 1), 1);
        end
    end

end
