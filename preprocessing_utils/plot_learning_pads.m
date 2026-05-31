function binSig = plot_learning_pads(signal, np_fs, file_path, aTemplate, padIndices, trigDelay, cW)
    

    binT = 100;
    binSig = zeros(size(signal, 1), size(signal, 2));
    cT = [trigDelay-cW trigDelay]*1000;

    
    if ~isfolder(fullfile(file_path, 'activity_aligned_pads'))
        mkdir(fullfile(file_path, 'activity_aligned_pads'));
    end

    parfor i = 1:size(signal, 3)
        signal_smooth = single(smoothdata(signal(:, :, i), 1, 'gaussian', 0.002*np_fs)*np_fs);
        binSig = smoothdata(signal_smooth, 2, "gaussian", binT);

        f = figure('visible', 'off');
        subplot(2, 1, 1);
        imagesc([0 size(binSig, 1)/np_fs*1000], [0 size(binSig, 2)], binSig');
        set(gca, 'TickDir', 'none');
        xlabel('time (ms)');
        ylabel('motif number');
        title(strcat("Pad Number ", num2str(padIndices(i))));
        cb = colorbar;
        cb.Location = 'manual';
        cb.Position = [0.915 0.6 0.02 0.3];
        cb.Label.String = 'firing rate (Hz)';
        cb.FontSize = 6;
        subplot(2, 1, 2);
        imagesc([0 size(aTemplate, 2)], [8000 500], aTemplate);
        xline(cT, '--r');
        set(gca, 'Ydir', 'normal');
        set(gca, 'TickDir', 'none');
        xlabel('time (ms)');
        ylabel('frequency (Hz)');
        colormap(colorcet('L20'));
        saveas(f, fullfile(file_path, 'activity_aligned_pads', strcat('pad_', num2str(padIndices(i)))), 'svg');
        close(f);
    end

end