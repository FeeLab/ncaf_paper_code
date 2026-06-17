%% Set random seed
rng(0);

%% import functions

addpath('../utils')
addpath('../external/SpikeGLX_Datafile_Tools/MATLAB');
addpath('../external/npy-matlab/npy-matlab');
addpath('../external/ndSparse_G4_2021_03_16');
addpath('../external/colorcet');

%% generate raster

load('../presorted_data/2025-11-12_Yellow33_Post-Advance_LMAN_BOTM_0_aligned.mat');
unitSignal = full(unitSigSparse);

unitI = 73; %target unit, very well isolated and huge amplitude (~800 uV)
Nstart = 1;
Nend = size(unitSignal, 2);


% detect bursts in data
sThresh = -log(.05);
isiThresh = 0.005;
[burstTimes, isiBurst, isiPoisson, fracBurst] = find_bursts_in_signal(unitSignal(:, :, unitNum==unitI), np_fs, songLength, 0, sThresh, isiThresh, 0);

% plot raster
plotN = 25;
figure;
ax1 = subplot(3, 1, 1);
raster_from_spikes_bursts(unitSignal, np_fs, trigDelay, cW, 1:plotN, find(unitNum==unitI), burstTimes);
ax2 = subplot(3, 1, 2);
raster_from_spikes_bursts(unitSignal, np_fs, trigDelay, cW, Nend-plotN:Nend, find(unitNum==unitI), burstTimes);
linkaxes([ax1, ax2], 'x');
subplot(3, 1, 3);
imagesc(silenceTemplate);
colormap(colorcet('L20'));
%% plot ISI distribution comparison

%calculate isi distribution
[isiVals, isiDist] = compute_ISI_distribution_all(unitSignal, np_fs);    
binW = 0.0001; %histogram bin width of 100 us
edgeVals = (0:binW:1)+binW/2;
binVals = edgeVals(1:end-1)+binW/2;
isiHist = histcounts(isiDist{unitNum==unitI}, edgeVals);
% fit relative refractory model to observed spike statistics
g = fit_relative_model(isiHist, binVals);

% simulate inhomogeneous Poisson process
binT = 100; %calculate mean firing rate over this many motifs
Nsim = 10; %number of times to simulate
simSignal = simulate_poisson_neuron_recovery(unitSignal(:, Nstart:Nend, unitNum==unitI), binT, np_fs, g, Nsim);

% detect bursts in data
sThresh = -log(.05); %surprisal threshold
isiThresh = 0.005; %maximum isi for spike
blankW = 0.01; %ignore data within this width of contingent window

%detect bursts in data
[burstTimes, isiBurst, isiPoisson, fracBurst, ISICat] = find_bursts_in_signal(unitSignal(:, Nstart:Nend, unitNum==unitI), np_fs, trigDelay, cW, sThresh, isiThresh, blankW);
%detect bursts in simulated data
[burstTimesSim, isiBurstSim, isiPoissonSim, fracBurstSim, ISICatSim] = find_bursts_in_signal(simSignal, np_fs, trigDelay, cW, sThresh, isiThresh, blankW);

Npoisson = 1000000000; %size of poisson distribution
%average rate parameter
r = sum(unitSignal(:, :, unitNum==unitI), 'all')/(size(unitSignal, 1)/np_fs)/size(unitSignal, 2);
%generate random waiting times using rate parameter
isiPoisson = exprnd(1/r, Npoisson, 1);

% generate and plot histograms
binW = .04; %bin width to plot
edgeVals = 10.^(-3:binW:0);
binVals = edgeVals(1:end-1)+diff(edgeVals)/2;

%calculate histograms
histData = histcounts(ISICat, edgeVals);
histSim = histcounts(ISICatSim, edgeVals);
histPoisson = histcounts(isiPoisson, edgeVals);

%plot figure
figure;
hold on;
plot(binVals, histData/sum(histData));
plot(binVals, histSim/sum(histSim));
plot(binVals, histPoisson/sum(histPoisson));
xscale('log');
xlabel('interspike interval (s)');
ylabel('probability');
set(gca, 'TickDir', 'out');
xline(0.005); %reference line for definition of burst isi in literature

%% plot pure bursting and pure Poisson comparison

% calculate spike counts and firing rate across experiment
binT = 100; %number of motifs to compute mean firing rate
%find spike count in contingent window
spikeCount = contingent_spikes(unitSignal(:, Nstart:Nend, unitNum==unitI), np_fs, trigDelay, cW, false);
%moving average over motifs
fRate = smoothdata(spikeCount/cW, 'movmean', binT);

% find firing rate within bursts
sThresh = -log(.05); %surprisal threshold
isiThresh = 0.005; %max isi for burst
%detect bursts
[burstTimes, isiBurst, isiPoisson, fracBurst] = find_bursts_in_signal(unitSignal(:, Nstart:Nend, unitNum==unitI), np_fs, songLength, 0, sThresh, isiThresh, 0);
%measure firing rate in bursts and rate of bursting
[pSpike, burstFR, burstRate] = parameterize_bursting(unitSignal(:, Nstart:Nend, unitNum==unitI), burstTimes, np_fs);

% simulate bursting only model
Nsim = 64; %number of times to simulate
simSignal_burstOnly = simulate_burst_ramp(fRate, Nsim, np_fs, cW, g, burstFR);

% simulate poisson only model
Nsim = 64; %number of times to simualte
simSignal_poissonOnly = simulate_poisson_ramp(fRate, Nsim, np_fs, cW, g, burstFR);

% plot both models
plot_fano_twoChannel(spikeCount, squeeze(sum(simSignal_burstOnly, 1)), squeeze(sum(simSignal_poissonOnly, 1)));
%% plot best fit composite model

% fit burst fraction to data
Nsim = 16; %number of times to simulate
options = optimset('Display','iter','TolX', 1e-4);
%function to compute log likelihood
fun = @(x)LL_from_model(x, fRate, Nsim, np_fs, cW, g, burstFR, spikeCount);
%fit to obtain beta parameters
beta = fminbnd(fun, 0, 1/(fRate(end)-fRate(1)), options);
alpha = (1-beta*(burstFR-fRate(1)))/(1-beta*(fRate(end)-fRate(1)));

% simulate and plot best fit
Nsim=64; %number of times to simulate
%simulate composite model using fit value of beta
[simSignal, burstTrials] = simulate_composite_ramp(fRate, Nsim, np_fs, cW, g, burstFR, beta);

%plot results
figure;
ax1 = gca;
[binVals, Ncounts] = plot_fano_comparison(spikeCount, squeeze(sum(simSignal, 1)), 4, ax1);
c_burst = burstFR*beta; %burst fraction
c_poisson = alpha*(1-beta*(fRate(end)-fRate(1)));
title("Unit "+unitI+", c_{burst} = "+c_burst);


%% find burst fractions for range of neurons

%well-isolated single units with substantial learning
unitVec = [73,242,253,256,764,779,785,787,815,858,868,873,893];
runVec = [repmat("2025-11-12_Yellow33_Post-Advance_LMAN_BOTM_0_aligned.mat", 4, 1);
    repmat("2025-02-18_11208_Post-Advance_LMAN_BOTM_0_aligned.mat", 9, 1)];

Nsamp = 100; %bootstrap replicates
blankW = 0.01; %ignore times within this width of contingent window
fracVals = zeros(Nsamp, numel(unitVec));
fracValsSim = zeros(Nsamp, numel(unitVec));
run_name = "nonexistent_run";

for k = 1:numel(unitVec)
    
    % load dataset if not already
    if ~contains(runVec(k), run_name)
        load(fullfile('../presorted_data', runVec(k)));
        unitSignal = full(unitSigSparse);
        %compute isi distribution for this run
        [isiVals, isiDist] = compute_ISI_distribution_all(unitSignal, np_fs);
    end

    % fit relative refractory model to data
    binW = 0.0001; %histogram resolution (100 us)
    edgeVals = (0:binW:1)+binW/2;
    binVals = edgeVals(1:end-1)+binW/2;
    isiHist = histcounts(isiDist{unitNum==unitVec(k)}, edgeVals);
    %fit relative refractory period model
    g = fit_relative_model(isiHist, binVals);
    title("unit "+unitVec(k));

    % simulate Poisson neuron
    binT = 100; %number of motifs to average firing rate over
    Nsim = 10; %number of times to simulate poisson process
    %simulate poisson process with relative refractory period
    simSignal = simulate_poisson_neuron_recovery(unitSignal(:, :, unitNum==unitVec(k)), binT, np_fs, g, Nsim);


    sThresh = -log(.05); %surprisal threshold
    isiThresh = 0.005; %maximum isi for burst
    fracBootstrap = zeros(Nsamp, 1);
    %locate bursts in randomly resampled data
    for i = 1:Nsamp
        [~, ~, ~, fracBootstrap(i)] = find_bursts_in_signal(unitSignal(:, randsample(size(unitSignal, 2), size(unitSignal, 2), true), unitNum==unitVec(k)), np_fs, trigDelay, cW, sThresh, isiThresh, blankW);
    end
    
    fracBootstrapSim = zeros(Nsamp, 1);
    %locate bursts in randomly resampled simulated poisson spiking
    for i = 1:Nsamp
        [~, ~, ~, fracBootstrapSim(i)] = find_bursts_in_signal(simSignal(:, randsample(size(simSignal, 2), size(unitSignal, 2), true)), np_fs, trigDelay, cW, sThresh, isiThresh, blankW);
    end

    fracVals(:, k) = fracBootstrap;
    fracValsSim(:, k) = fracBootstrapSim;

    disp(k);
end

% plot results
figure;
errorbar(1:size(fracVals, 2), mean(fracVals), mean(fracVals)-prctile(fracVals, 2.5), prctile(fracVals, 97.5)-mean(fracVals), 'o', 'MarkerFaceColor', 'auto');
hold on;
bar((size(fracVals, 2)+1)/2, mean(fracVals(:)));
ylim([0 Inf]);
xlim([0 size(fracVals, 2)+1])
ylabel('fraction of spikes in bursts');
set(gca, 'TickDir', 'out');
yline(mean(fracValsSim, 'all'), '-r');
yline(min(mean(fracValsSim)), '--r');
yline(max(mean(fracValsSim)), '--r');

%output values in text
disp("spike fraction = "+mean(fracVals, 'all')+" +/- "+std(mean(fracVals)));
disp("excess factor = "+mean(fracVals, 'all')/mean(fracValsSim, 'all'));

[h, p, ci, stats] = ttest2(mean(fracVals, 1), mean(fracValsSim, 1), 'Tail', 'right');
disp("p-value = "+p);
disp("t-value = "+stats.tstat);

%% find burst rate for range of neurons

%set of very well-isolated units (rp_violation < 0.02)
unitVec = [73,242,253,256,764,779,785,787,815,858,868,873,893];
runVec = [repmat("2025-11-12_Yellow33_Post-Advance_LMAN_BOTM_0_aligned.mat", 4, 1);
    repmat("2025-02-18_11208_Post-Advance_LMAN_BOTM_0_aligned.mat", 9, 1)];
Nsamp = 100; %bootstrap replicates
blankW = 0.01; %ignore data within this width of contingent window
rateVals = zeros(numel(unitVec), 1);

run_name = "nonexistent_run";

burstLengths = [];
for k = 1:numel(unitVec)
    
    % load dataset if not already
    if ~contains(runVec(k), run_name)
        load(fullfile('../presorted_data', runVec(k)));
        unitSignal = full(unitSigSparse);
        %compute isi distribution for this run
        [isiVals, isiDist] = compute_ISI_distribution_all(unitSignal, np_fs);
    end

    sThresh = -log(.05); %surprisal threshold
    isiThresh = 0.005; %maximum isi for burst
    %find bursts in data
    [burstTimes, isiBurst, isiPoisson, fracBurst] = find_bursts_in_signal(unitSignal(:,:, unitNum==unitVec(k)), np_fs, songLength, 0, sThresh, isiThresh, 0);
    %calculate firing rate in bursts and rate of bursting
    [pSpike, burstFR, rateVals(k)] = parameterize_bursting(unitSignal(:,:, unitNum==unitVec(k)), burstTimes, np_fs);

    %concatenate lengths of each burst in ms
    for j = 1:numel(burstTimes)
        burstLengths = [burstLengths; diff(burstTimes{j}, 1, 2)/np_fs*1000];
    end

    disp(k);
end

%output values in text
disp("burst rate = "+mean(rateVals)+" +/- "+std(rateVals));
disp("burst length = "+mean(burstLengths)+" +/- "+std(burstLengths));
%% find composite model fit parameters for range of neurons

%set of single neurons with significant learning
unitVec = [71,73,75,113,127,132,133,134,187];
runVec = [repmat("2025-11-12_Yellow33_Post-Advance_LMAN_BOTM_0_aligned.mat", 3, 1);
    repmat("2025-02-21_11208_Post-Advance_LMAN_nCAF_0_aligned.mat", 5, 1);
    repmat("2024-04-27_10872_LMAN-X_nCAF_2_aligned.mat", 1, 1)];

Nstrap = 100; %bootstrap replicates
output = zeros(numel(unitVec), Nstrap);
run_name = "nonexistent_run";
for k = 1:numel(unitVec)
    unitI = unitVec(k);

    % load dataset if not already
    if ~contains(runVec(k), run_name)
        load(fullfile('../presorted_data', runVec(k)));
        unitSignal = full(unitSigSparse);
        %compute isi distribution for run
        [isiVals, isiDist] = compute_ISI_distribution_all(unitSignal, np_fs);
    end
    
    % fit relative refractory model
    binW = 0.0001; %histogram bin width (100 ms)
    edgeVals = (0:binW:1)+binW/2;
    binVals = edgeVals(1:end-1)+binW/2;
    isiHist = histcounts(isiDist{unitNum==unitI}, edgeVals);
    %fit relative refractory model to data
    g = fit_relative_model(isiHist, binVals);

    % locate bursts
    sThresh = -log(.05); %surprisal threshold
    isiThresh = 0.005; %maximum isi fo rburst
    %detect bursts
    [burstTimes, isiBurst, isiPoisson, fracBurst] = find_bursts_in_signal(unitSignal(:, :, unitNum==unitI), np_fs, songLength, 0, sThresh, isiThresh, 0);
    %calculate burst firing rate and rate of bursting
    [pSpike, burstFR, burstRate] = parameterize_bursting(unitSignal(:, :, unitNum==unitI), burstTimes, np_fs);
    burstFRvals(k) = burstFR;
    binT = 100; %number of motifs to calculate firing rate over
    spikeCount = contingent_spikes(unitSignal(:, :, unitNum==unitI), np_fs, trigDelay, cW, false);
    fRate = smoothdata(spikeCount/cW, 'movmean', binT);

    % bootstrap fit distribution
    betaStrap = zeros(Nstrap, 1);
    Nsim = 16; %bootstrap replicates
    for i = 1:Nstrap
        %calculate mean and variance across motifs from resampled data
        [unitx, unity] = fano_plot_count_resample(spikeCount);
        %compute log-likelihood from resample data
        fun = @(x)LL_from_model_resampled(x, fRate, Nsim, np_fs, cW, g, burstFR, unitx, unity);
        %maximize log-likelihood
        output(k, i) = fminbnd(fun, 0, 1/burstFR)*burstFR;
    end
    disp(k);
end

% plot results
figure;
errorbar(1:size(output, 1), mean(output, 2), 2*std(output, 0, 2), 'o', 'MarkerFaceColor', 'auto');
hold on;
bar((size(output, 1)+1)/2, mean(output(:)));
ylim([0 Inf]);
xlim([0 (size(output, 1)+1)])
ylabel('fraction of spikes in bursts');
set(gca, 'TickDir', 'out');
title('fraction of learned spikes in bursts');