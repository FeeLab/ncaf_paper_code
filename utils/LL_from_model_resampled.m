function logLikelihood = LL_from_model_resampled(x, fRate, Nsim, np_fs, cW, g, burstFR, datax, datay)


    f_burst = x;
    simSignal = simulate_composite_ramp(fRate, Nsim, np_fs, cW, g, burstFR, f_burst);
    
    [simx, simy] = fano_plot_count(squeeze(sum(simSignal, 1)));
    limVal = ceil(max([simx; simy]));

    binRes = 0.04;
    edgeVals = 0:binRes:limVal;
    binVals = edgeVals(1:end-1)+binRes/2;
    [Ncounts, ~, ~] = histcounts2(simy, simx, edgeVals, edgeVals, 'Normalization', 'probability');

    logLikelihood = fano_to_logLikelihood_resampled(Ncounts, binVals, datax, datay);
end