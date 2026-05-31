function output = fano_to_logLikelihood(Ncounts, binVals, spikeCount)

    ff = .001*max(Ncounts(:));
    offsetDist = (Ncounts+ff)/sum(Ncounts+ff, 'all');
    [datax, datay] = fano_plot_count(spikeCount);
    logLikelihoods = zeros(numel(datax), 1);
    for i = 1:numel(datax)
        [~, binX] = min(abs(binVals-datax(i)));
        [~, binY] = min(abs(binVals-datay(i)));
        logLikelihoods(i) = -log(offsetDist(binY, binX));
    end
    output = sum(logLikelihoods);
end