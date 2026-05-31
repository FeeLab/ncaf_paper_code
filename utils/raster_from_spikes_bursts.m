function raster_from_spikes_bursts(signal, np_fs, trigDelay, cW, motifNum, chan, burstTimes)
    set(0, 'DefaultFigureRenderer', 'painters');
    currMotif = 1;
    lineH = 1;
    plotW = .1;
    for i = 1:numel(motifNum)
        burstPoints = [];
        for j = 1:size(burstTimes{motifNum(i)}, 1)
            burstPoints = [burstPoints burstTimes{motifNum(i)}(j, 1):burstTimes{motifNum(i)}(j, 2)];
        end

        for j = 1:size(signal, 1)
            if signal(j, motifNum(i), chan)
                if ismember(j, burstPoints)
                    line([1 1]*j/np_fs*1000, [-0.5 0.5]*lineH+motifNum(i), 'Color', 'red', 'LineWidth', 2);
                else
                    line([1 1]*j/np_fs*1000, [-0.5 0.5]*lineH+motifNum(i), 'Color', 'black', 'LineWidth', 2);
                end
            end
        end
        currMotif = currMotif+1;
    end

    xline([trigDelay-cW trigDelay]*1000, '--r');
    set(gca, 'YDir', 'reverse');
    set(gca, 'Tickdir', 'out');
    ylim([min(motifNum)-lineH/2 max(motifNum)+lineH/2]);

    xlim([0 size(signal, 1)/np_fs*1000]);
    xlabel('time in song (ms)');
    ylabel('motif number');
end