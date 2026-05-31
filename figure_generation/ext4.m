%% import functions

addpath('../utils')
addpath('../external/SpikeGLX_Datafile_Tools/MATLAB');
addpath('../external/npy-matlab/npy-matlab');
addpath('../external/ndSparse_G4_2021_03_16');
addpath('../external/colorcet');


%% plot ISI distributions

run_name = "nonexistent_run";
unitVec = [256,815,893];
runVec = [repmat("2025-11-12_Yellow33_Post-Advance_LMAN_BOTM_0_aligned.mat", 1, 1);
    repmat("2025-02-18_11208_Post-Advance_LMAN_BOTM_0_aligned.mat", 2, 1)];


figure;
tl = tiledlayout('flow');

for k = 1:numel(unitVec)
    unitI = unitVec(k);
    % load dataset if not already
    if ~contains(runVec(k), run_name)
        load(fullfile('../presorted_data', runVec(k)));
        unitSignal = full(unitSigSparse);
        [isiVals, isiDist] = compute_ISI_distribution_all(unitSignal, np_fs);
    end
    
    % fit relative refractory model to observed spike statistics
    binW = 0.0001;
    edgeVals = (0:binW:1)+binW/2;
    binVals = edgeVals(1:end-1)+binW/2;
    isiHist = histcounts(isiDist{unitNum==unitI}, edgeVals);
    g = fit_relative_model(isiHist, binVals);
    
    % simulate inhomogeneous Poisson process
    binT = 100;
    Nsim = 1;
    simSignal = simulate_poisson_neuron_recovery(unitSignal(:, :, unitNum==unitI), binT, np_fs, g, Nsim);
    
    % detect bursts in data
    sThresh = -log(.05);
    isiThresh = 0.005;
    blankW = 0.01;
    
    [burstTimes, isiBurst, isiPoisson, fracBurst, ISICat] = find_bursts_in_signal(unitSignal(:, :, unitNum==unitI), np_fs, trigDelay, cW, sThresh, isiThresh, blankW);
    [burstTimesSim, isiBurstSim, isiPoissonSim, fracBurstSim, ISICatSim] = find_bursts_in_signal(simSignal, np_fs, trigDelay, cW, sThresh, isiThresh, blankW);
    
    Npoisson = 100000000;
    r = sum(unitSignal(:, :, unitNum==unitI), 'all')/(size(unitSignal, 1)/np_fs)/size(unitSignal, 2);
    isiPoisson = exprnd(1/r, Npoisson, 1);
    
    % generate and plot histograms
    binW = .04;
    edgeVals = 10.^(-3:binW:0);
    binVals = edgeVals(1:end-1)+diff(edgeVals)/2;
    
    histData = histcounts(ISICat, edgeVals);
    histSim = histcounts(ISICatSim, edgeVals);
    histPoisson = histcounts(isiPoisson, edgeVals);
    
    %figure;
    nexttile(tl);
    hold on;
    plot(binVals, histData/sum(histData));
    plot(binVals, histSim/sum(histSim));
    plot(binVals, histPoisson/sum(histPoisson));
    xscale('log');
    xlabel('interspike interval (s)');
    ylabel('probability');
    set(gca, 'TickDir', 'out');
    xline(0.005);
    title("Unit "+unitI);
    shg;
end

%% find composite model fit parameters for range of neurons

figure;
tl = tiledlayout('flow');

unitVec = [75,127,187];
runVec = [repmat("2025-11-12_Yellow33_Post-Advance_LMAN_BOTM_0_aligned.mat", 1, 1);
    repmat("2025-02-21_11208_Post-Advance_LMAN_nCAF_0_aligned.mat", 1, 1);
    repmat("2024-04-27_10872_LMAN-X_nCAF_2_aligned.mat", 1, 1)];

Nstrap = 20;
output = zeros(numel(unitVec), Nstrap);
run_name = "nonexistent_run";
for k = 1:numel(unitVec)
    unitI = unitVec(k);

    % load dataset if not already
    if ~contains(runVec(k), run_name)
        load(fullfile('../presorted_data', runVec(k)));
        unitSignal = full(unitSigSparse);
        [isiVals, isiDist] = compute_ISI_distribution_all(unitSignal, np_fs);
    end
    
    % fit relative refractory model
    binW = 0.0001;
    edgeVals = (0:binW:1)+binW/2;
    binVals = edgeVals(1:end-1)+binW/2;
    isiHist = histcounts(isiDist{unitNum==unitI}, edgeVals);
    g = fit_relative_model(isiHist, binVals);

    % locate bursts
    sThresh = -log(.05);
    isiThresh = 0.005;
    [burstTimes, isiBurst, isiPoisson, fracBurst] = find_bursts_in_signal(unitSignal(:, :, unitNum==unitI), np_fs, songLength, 0, sThresh, isiThresh, 0);
    [pSpike, burstFR, burstRate] = parameterize_bursting(unitSignal(:, :, unitNum==unitI), burstTimes, np_fs);
    burstFRvals(k) = burstFR;
    binT = 100;
    spikeCount = contingent_spikes(unitSignal(:, :, unitNum==unitI), np_fs, trigDelay, cW, false);
    fRate = smoothdata(spikeCount/cW, 'movmean', binT);

    % fit burst fraction to data
    Nsim = 16;
    options = optimset('Display','iter','TolX', 1e-4);
    fun = @(x)LL_from_model(x, fRate, Nsim, np_fs, cW, g, burstFR, spikeCount);
    beta = fminbnd(fun, 0, 1/(fRate(end)-fRate(1)), options);
    alpha = (1-beta*(burstFR-fRate(1)))/(1-beta*(fRate(end)-fRate(1)));
    
    ax1 = nexttile(tl);
    % simulate and plot best fit
    Nsim=64;
    [simSignal, burstTrials] = simulate_composite_ramp(fRate, Nsim, np_fs, cW, g, burstFR, beta);
    [binVals, Ncounts] = plot_fano_comparison(spikeCount, squeeze(sum(simSignal, 1)), 4, ax1);
    c_burst = burstFR*beta;
    c_poisson = alpha*(1-beta*(fRate(end)-fRate(1)));
    title("Unit "+unitI+", c_{burst} = "+c_burst);
    shg;
end