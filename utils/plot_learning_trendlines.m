
function meanF = plot_learning_trendlines(unitSignal, np_fs, trigDelay, thisUnit, unitNum, correctBaseline, plotAx)

    %plot mean firing rates in contingent window
    
    cW = .01;

    trialCounts = contingent_spikes(unitSignal(:,:,unitNum==thisUnit), np_fs, trigDelay, cW, correctBaseline);

    meanW = 100;
    Nbootstrap = 100;
    meanF = zeros(size(trialCounts));
    stdF = zeros(size(trialCounts));
    for i = 1:size(trialCounts, 1)
        meanWindow = trialCounts(max(1, i-meanW):min(size(trialCounts, 1), i+meanW), :);
        meanF(i) = mean(meanWindow)/cW;
        stdF(i) = std(meanWindow/cW)/sqrt(meanW);
    end
    
    if ~exist("plotAx")
        f = figure;
        plotAx = gca(f);
    end
    lb = meanF-2*stdF;
    ub = meanF+2*stdF;
    xconf = [1:size(meanF, 1) size(meanF, 1):-1:1];
    yconf = [lb' ub(end:-1:1)'];
    pLine = plot(plotAx, 1:size(meanF, 1), meanF);
    hold on;
    p=fill(plotAx, xconf, yconf, pLine.Color);
    p.FaceAlpha = .1;
    p.EdgeColor = 'none';

    xlim([1, size(meanF, 1)]);
    ylim([0, Inf]);
    xlabel('motif number');
    ylabel('firing rate (Hz)');
    set(gca, 'TickDir', 'out');
end