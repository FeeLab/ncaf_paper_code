function [unitx, unity] = fano_plot_count(spikeCount)

    windowW = 100;
    unitx = reshape(movmean(spikeCount, windowW, 'Endpoints', 'discard'), [], 1);
    unity = reshape(movvar(spikeCount, windowW, 'Endpoints', 'discard'), [], 1);

end