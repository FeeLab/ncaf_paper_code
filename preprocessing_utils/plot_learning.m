function plot_learning(signal, unitNum, np_fs, daq_fs, file_path, aTemplate, trigDelay, cW, noiseStart, noiseW)
    
    binT = 100;
    Ntemp = size(signal, 3);
    binSig = zeros(size(signal, 1), size(signal, 2));
    frateW = 0.002;

    cT = [trigDelay-cW trigDelay]*1000;
    
    if ~isfolder(fullfile(file_path, 'activity_aligned_median'))
        mkdir(fullfile(file_path, 'activity_aligned_median'));
    end
    if ~isfolder(fullfile(file_path, 'activity_aligned_raw'))
        mkdir(fullfile(file_path, 'activity_aligned_raw'));
    end


    parfor i = 1:Ntemp
        signal_smooth = single(smoothdata(signal(:, :, i), 1, 'gaussian', frateW*np_fs)*np_fs);
        binSig = smoothdata(signal_smooth, 2, "gaussian", binT);

        f = figure('visible', 'off');
        subplot(2, 1, 1);
        
        plotOut = binSig./repmat(mean(binSig([1:round((trigDelay-cW-.05)*np_fs) round((trigDelay+.05)*np_fs):end], :), 1), size(binSig, 1), 1);
        imagesc([0 size(plotOut, 1)/np_fs*1000], [0 size(binSig, 2)], plotOut');
        
        set(gca, 'TickDir', 'none');
        xlabel('time (ms)');
        ylabel('motif number');
        title(strcat("Unit Number ", num2str(unitNum(i))));
        cb = colorbar;
        cb.Location = 'manual';
        cb.Position = [0.915 0.6 0.02 0.3];
        cb.Label.String = 'firing rate (relative)';
        cb.FontSize = 6;
        subplot(2, 1, 2);
        imagesc([0 size(aTemplate, 2)-1], [8000 500], aTemplate);
        xline(cT, '--r');
        set(gca, 'Ydir', 'normal');
        set(gca, 'TickDir', 'none');
        xlabel('time (ms)');
        ylabel('frequency (Hz)');
        colormap(colorcet('L20'));

        saveas(f, fullfile(file_path, 'activity_aligned_median', strcat('unit_', num2str(unitNum(i)))));
        saveas(f, fullfile(file_path, 'activity_aligned_median', strcat('unit_', num2str(unitNum(i)))), 'svg');
        close(f);

        
        f = figure('visible', 'off');
        subplot(2, 1, 1);
        imagesc([0 size(binSig, 1)/np_fs*1000], [0 size(binSig, 2)], binSig');
        set(gca, 'TickDir', 'none');
        xlabel('time (ms)');
        ylabel('motif number');
        title(strcat("Unit Number ", num2str(unitNum(i))));
        cb = colorbar;
        cb.Location = 'manual';
        cb.Position = [0.915 0.6 0.02 0.3];
        cb.Label.String = 'firing rate (Hz)';
        cb.FontSize = 6;
        subplot(2, 1, 2);
        imagesc([0 size(aTemplate, 2)-1], [8000 500], aTemplate);
        xline(cT, '--r');
        set(gca, 'Ydir', 'normal');
        set(gca, 'TickDir', 'none');
        xlabel('time (ms)');
        ylabel('frequency (Hz)');
        colormap(colorcet('L20'));
        
        saveas(f, fullfile(file_path, 'activity_aligned_raw', strcat('unit_', num2str(unitNum(i)))));
        saveas(f, fullfile(file_path, 'activity_aligned_raw', strcat('unit_', num2str(unitNum(i)))), 'svg');
        close(f);
    end
end