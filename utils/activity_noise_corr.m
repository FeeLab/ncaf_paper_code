function [phiB, phiConf] = activity_noise_corr(spikeCount, isNoise)
    
    medW = 50;
    movThresh = movmedian(spikeCount, [medW-1 0], 'Endpoints', 'discard');
    movThresh = movThresh(1:end-1, :);
    excessSpikes = spikeCount(medW+1:end, :) - movThresh;
    
    phiB = zeros(size(spikeCount, 2), 1);
    phiConf = zeros(size(spikeCount, 2), 2);
    for i = 1:size(phiB, 1)
        [c, ~, rl, ru] = corrcoef(excessSpikes(:, i), double(~isNoise(medW+1:end))');
        phiB(i) = c(1, 2);
        phiConf(i, :) = [rl(1, 2), ru(1, 2)];
    end
    
end