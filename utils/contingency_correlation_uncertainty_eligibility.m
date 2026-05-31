
function [tvals, frateCorr, confC] = contingency_correlation_uncertainty_eligibility(unitSignal, np_fs, trigDelay, cW, unitI, unitNum, isNoise, plotAx, isNormed, ww, binW)
    
    noiseStart = 101;
    
    fetchDelay = 0;
    trigEnd = trigDelay+fetchDelay;
    trigStart = trigEnd-cW;

    medW = 50;
    tOrigin = trigEnd;
    originI = round(tOrigin*np_fs/binW);


    signalDec = downsample_spikes(unitSignal(:,:,unitNum==unitI), binW, mod(round(tOrigin*np_fs)-1, binW))*np_fs;
    iVals = round((trigStart-ww)*np_fs/binW):round((trigEnd+ww)*np_fs/binW);
    
    spikeC = zeros(numel(iVals), size(signalDec, 2));
    spikeC(iVals<=size(signalDec, 1), :) = signalDec(iVals(iVals<=size(signalDec, 1)), :);

    normBuff = .05;
    binT = 100;
    frateW = .01;
    if isNormed
        tOutside = [1:round((trigDelay-cW-normBuff)*np_fs) round((trigDelay+normBuff)*np_fs):size(unitSignal, 1)];
        frate = smoothdata(unitSignal(:,:,unitNum==unitI), 'gaussian', frateW*np_fs)*np_fs;
        binSig = smoothdata(frate, 2, 'gaussian', binT);
        corrF = median(binSig(tOutside, 1))./median(binSig(tOutside, :), 1);
        spikeC = spikeC.*repmat(corrF, size(spikeC, 1), 1);
    end


    spikeC = spikeC - movmedian(spikeC, 50, 2);
    spikeC = spikeC(:, noiseStart:end);

    frateCorr = zeros(size(spikeC, 1), 1);
    confC = zeros(size(spikeC, 1), 2);

    for i = 1:size(spikeC, 1)
        [c, ~, rl, ru] = corrcoef(spikeC(i, :), double(~isNoise(noiseStart:end))');
        frateCorr(i) = c(1, 2);
        confC(i, :) = [rl(1, 2), ru(1, 2)];
    end


    tvals = 1000*((iVals-originI)/np_fs*binW);

    nanI = isnan(frateCorr);
    frateCorr(nanI) = 0;
    confC(nanI, 1) = 0;
    confC(nanI, 2) = 0;

    if ~isempty(plotAx)

        xconf = [tvals tvals(end:-1:1)];
        yconf = [confC(:, 1)' confC(end:-1:1, 2)'];
        p=fill(plotAx, xconf, yconf, 'red');
        hold on;
        p.FaceColor = [1 0.8 0.8];
        p.FaceAlpha = .4;
        p.EdgeColor = 'none';
        xlabel('time relative to contingency midpoint (ms)');
        ylabel('correlation between firing rate and noise')
        title('activity-noise correlation');
        plot(plotAx, tvals, frateCorr, '-r');
    
        xline(([trigStart trigEnd]-tOrigin)*1000, '--k');
        xlim([min(tvals), max(tvals)]);
    
        yline(0);
        set(gca, 'TickDir', 'out');
    end

end