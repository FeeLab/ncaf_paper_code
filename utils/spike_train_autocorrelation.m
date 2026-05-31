function [acorr, meanAcorr, lags] = spike_train_autocorrelation(unitSignal, np_fs, unitNum, unitI, maxLag, binW)

    signalDec = downsample_spikes(unitSignal(:,:,unitNum==unitI), binW, 0)*binW;
    acorr = zeros(2*round(maxLag*np_fs/binW)+1, size(signalDec, 2));
    for i = 1:size(signalDec, 2)
        [acorr(:, i), lags] = xcorr(signalDec(:, i), round(maxLag*np_fs/binW), 'normalized');
    end
    meanAcorr = mean(acorr, 2);
    lags = lags/np_fs*binW;

end