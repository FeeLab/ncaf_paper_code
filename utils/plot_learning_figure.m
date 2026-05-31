function plot_learning_figure(unitSignal, trigDelay, unitI, unitNum, aTemplate, cW, np_fs, noiseI, noiseW)


    frateW = 0.002;
    binT = 100;
    cT = [trigDelay-cW trigDelay]*1000;
    
    signal_smooth = single(smoothdata(unitSignal(:, :, unitNum==unitI), 1, 'gaussian', frateW*np_fs)*np_fs);
    binSig = smoothdata(signal_smooth, 2, "gaussian", binT);
    
    f = figure;
    ax1 = subplot(2, 1, 1);
    imagesc([0 size(binSig, 1)/np_fs*1000], [0 size(binSig, 2)], binSig');
    set(gca, 'TickDir', 'none');
    xlabel('time (ms)');
    ylabel('motif number');
    title(strcat("Unit Number ", num2str(unitI)));
    cb = colorbar;
    cb.Location = 'manual';
    cb.Position = [0.915 0.6 0.02 0.3];
    cb.Label.String = 'firing rate (relative)';
    cb.FontSize = 6;
    colormap(ax1, colorcet('L20'));
    
    ax2 = subplot(2, 1, 2);
    imagesc([0 size(aTemplate, 2)-1], [8000 500], aTemplate);
    xline(cT, '--r');
    xline(1000*(median(noiseI)+[0 noiseW]), '--g');
    set(gca, 'Ydir', 'normal');
    set(gca, 'TickDir', 'none');
    xlabel('time (ms)');
    ylabel('frequency (Hz)');
    colormap(ax2, colorcet('L01'));
end