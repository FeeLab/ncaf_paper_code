function [unitx, unity] = fano_plot_count_resample(spikeCount)

    windowW = 100;

    unitx = zeros(numel(spikeCount)-windowW+1, 1);
    unity = zeros(numel(spikeCount)-windowW+1, 1);
    for i = 1:numel(unitx)
        thisCount = spikeCount(i:i+windowW-1);
        spikeResample = thisCount(randsample(windowW, windowW, true));
        unitx(i) = mean(spikeResample);
        unity(i) = var(spikeResample);
    end
end