function fluc_frac = fr_fluctuations(unitSignal, np_fs, trigDelay, cW)

    buffW = .01; % buffer period around contingent window

    outsideT = [1:round((trigDelay-cW-buffW)*np_fs) round((trigDelay+buffW)*np_fs):size(unitSignal, 1)];

    mean_fr = squeeze(sum(unitSignal(outsideT, :, :), 1)/(size(unitSignal, 1)/np_fs-(2*buffW+cW)));

    fluc_frac = (std(movmean(mean_fr, 100, 1))./mean(mean_fr, 1))';

end