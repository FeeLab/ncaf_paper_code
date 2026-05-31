function [learnR, confR] = learning_per_neuron_outsideTimes(unitSignal, np_fs, trigDelay, thisUnit, unitNum, cW)
    

    buffW = .01; % buffer period around contingent window


    Nt = 0.0001;
    outsideT = [cW:.001:trigDelay-cW-buffW trigDelay+buffW:.001:size(unitSignal, 1)/np_fs];
    outsideTimes = zeros(round(cW*np_fs), floor((trigDelay-cW-buffW)/Nt)+floor((size(unitSignal, 1)/np_fs-trigDelay-buffW-cW)/Nt), 'uint32');
    for i = 1:floor((trigDelay-cW-buffW)/Nt)
        outsideTimes(:, i) = (1:round(cW*np_fs)) + round(Nt*np_fs*(i-1));
    end
    for i = 1:size(outsideTimes, 2)-floor((trigDelay-cW-buffW)/Nt)
        outsideTimes(:, i+floor((trigDelay-cW-buffW)/Nt)) = (1:round(cW*np_fs)) + round(Nt*np_fs*(i-1)) + round((trigDelay+buffW)*np_fs);
    end

    bgF = zeros(size(unitSignal, 2), size(outsideTimes, 2));
    for i = 1:size(bgF, 2)
        bgF(:, i) = sum(unitSignal(outsideTimes(:, i), :, unitNum==thisUnit))/cW;
    end

    rTot = zeros(size(bgF, 2), 1);
    for i = 1:size(bgF, 2)
        fitobject = fit((1:size(bgF, 1))', bgF(:, i), 'poly1');
        rTot(i) = fitobject.p1;
    end

    learnR = median(rTot);
    confR = [prctile(rTot, 2.5), prctile(rTot, 97.5)];
end