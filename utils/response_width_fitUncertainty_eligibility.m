

function [tvals, learnR, confR] = response_width_fitUncertainty_eligibility(unitSignal, np_fs, trigDelay, cW, unitI, unitNum, plotAx, isNormed, ww, binW)
    

    fetchDelay = 0;
    trigEnd = trigDelay+fetchDelay;
    trigStart = trigEnd-cW;


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


    learnR = zeros(size(spikeC, 1), 1);
    confR = zeros(size(spikeC, 1), 2);

    for i = 1:size(spikeC, 1)
        b = [ones(size(spikeC, 2), 1) (1:size(spikeC, 2))']\spikeC(i,:)';
        learnR(i) = b(2)*np_fs;


        fitobject = fit((1:size(spikeC, 2))', spikeC(i,:)', 'poly1');
        learnR(i) = fitobject.p1;
        c = confint(fitobject);
        confR(i, :) = c(:, 1);

    end

    tvals = 1000*((iVals-originI)/np_fs*binW);

    if ~isempty(plotAx)
        xconf = [tvals tvals(end:-1:1)];
        yconf = [confR(:, 1)' confR(end:-1:1, 2)'];
        p=fill(plotAx, xconf, yconf, 'red');
        hold on;
        p.FaceColor = [1 0.8 0.8];
        p.FaceAlpha = .4;
        p.EdgeColor = 'none';
        xlabel('time relative to contingency midpoint (ms)');
        ylabel('learning rate (Hz/trial)');
        title('learning rate');
        plot(plotAx, tvals, learnR, '-r');
    
        xline(([trigStart trigEnd]-tOrigin)*1000, '--k');
        xlim([min(tvals), max(tvals)]);
        yline(0);
        set(gca, 'TickDir', 'out');
    end
end