
function plot_learning_trendlines_allTimes(unitSignal, np_fs, trigDelay, thisUnit, unitNum, plotAx)

    %plot mean firing rates in contingent window
    
    cW = .01;
    buffW = .01;
    
    insideTimes = round((trigDelay-cW)*np_fs):round((trigDelay)*np_fs);

    trialCounts = sum(unitSignal(insideTimes, :, unitNum==thisUnit), 1)';
    meanW = 100;
    Nbootstrap = 1000;
    meanF = zeros(size(trialCounts));
    bootstrapF = zeros(size(trialCounts, 1), Nbootstrap, size(trialCounts, 2));
    for i = 1:size(trialCounts, 1)
        meanWindow = trialCounts(max(1, i-meanW):min(size(trialCounts, 1), i+meanW), :);
        meanF(i, :) = mean(meanWindow)/cW;
        parfor j = 1:Nbootstrap
            bootstrapF(i, j, :) = mean(meanWindow(randsample(size(meanWindow, 1), size(meanWindow, 1), true), :))/numel(insideTimes)*np_fs;
        end
    end

    Nt = 0.0001;
    outsideT = [cW:.001:trigDelay-cW-buffW trigDelay+buffW:.001:size(unitSignal, 1)/np_fs];
    outsideTimes = zeros(round(cW*np_fs), floor((trigDelay-cW-buffW)/Nt)+floor((size(unitSignal, 1)/np_fs-trigDelay-buffW-cW)/Nt), 'uint32');
    for i = 1:floor((trigDelay-cW-buffW)/Nt)
        outsideTimes(:, i) = (1:round(cW*np_fs)) + round(Nt*np_fs*(i-1));
    end
    for i = 1:size(outsideTimes, 2)-floor((trigDelay-cW-buffW)/Nt)
        outsideTimes(:, i+floor((trigDelay-cW-buffW)/Nt)) = (1:round(cW*np_fs)) + round(Nt*np_fs*(i-1)) + round((trigDelay+buffW)*np_fs);
    end

    bgCounts = zeros(size(outsideTimes, 2), size(unitSignal, 2));
    for i = 1:size(bgCounts, 1)
        bgCounts(i, :) = sum(unitSignal(outsideTimes(:, i), :, unitNum==thisUnit))/cW;
    end
    bgF = movmean(bgCounts, meanW, 2);
    
    if ~exist("plotAx")
        f = figure;
        plotAx = gca(f);
    end
    lb = prctile(bootstrapF(:, :), 2.5, 2);
    ub = prctile(bootstrapF(:, :), 97.5, 2);
    xconf = [1:size(meanF, 1) size(meanF, 1):-1:1];
    yconf = [lb' ub(end:-1:1)'];
    p=fill(plotAx, xconf, yconf, [0, 0, 1]);
    hold on;
    p.FaceAlpha = .1;
    p.EdgeColor = 'none';
    plot(plotAx, 1:size(meanF, 1), meanF, 'Color', [0, 0, 1]);

    lb = prctile(bgF, 2.5, 1);
    ub = prctile(bgF, 97.5, 1);
    xconf = [1:size(bgF, 2) size(bgF, 2):-1:1];
    yconf = [lb ub(end:-1:1)];
    p=fill(plotAx, xconf, yconf, [1, 0, 0]);
    hold on;
    p.FaceAlpha = .1;
    p.EdgeColor = 'none';
    plot(plotAx, 1:size(bgF, 2), mean(bgF, 1), 'Color', [1, 0, 0]);

    xlim([1, size(meanF, 1)]);
    ylim([0, Inf]);
    xlabel('motif number');
    ylabel('firing rate (Hz)');
    set(gca, 'TickDir', 'out');
end