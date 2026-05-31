function [pSpike, burstFR, burstRate] = parameterize_bursting(unitSignal, burstTimes, np_fs)
        
    %calculate burst count distribution
    burstCounts = [];
    for i = 1:numel(burstTimes)
        for j = 1:size(burstTimes{i}, 1)
            burstCounts = [burstCounts; sum(unitSignal(burstTimes{i}(j, 1):burstTimes{i}(j, 2), i))];
        end
    end
    
    countBins = 3:max(burstCounts);
    countEdges = [countBins-0.5 countBins(end)+0.5];
    burstMag = histcounts(burstCounts, countEdges);
    burstMag = burstMag/sum(burstMag);
    
    f = fit(countBins', burstMag', 'exp1');
    figure;
    plot(countBins, burstMag);
    hold on;
    plot(countBins, feval(f, countBins));
    ylabel('probability');
    xlabel('number of spikes in burst');
    title('distribution of burst sizes outside contingent window');
    legend({'data', 'exponential fit'});
    pSpike = 1-f(3);
    
    %calculate mean firing rate in bursts
    burstISI = [];
    for i = 1:numel(burstTimes)
        for j = 1:size(burstTimes{i}, 1)
            spikeT = find(unitSignal(:, i));
            burstISI = [burstISI; diff(spikeT(spikeT>=burstTimes{i}(j, 1) & spikeT<=burstTimes{i}(j, 2)))/np_fs];
        end
    end
    figure;
    histogram(burstISI);
    burstFR = 1/mean(burstISI);
    
    %calculate burst fraction
    burstTot = 0;
    for i = 1:numel(burstTimes)
        for j = 1:size(burstTimes{i}, 1)
            burstTot = burstTot+1;
        end
    end
    burstRate = burstTot/numel(burstTimes)/size(unitSignal, 1)*np_fs;

end