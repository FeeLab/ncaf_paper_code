%% Set random seed
rng(0);

%% import functions

addpath('../utils')
addpath('../external/SpikeGLX_Datafile_Tools/MATLAB');
addpath('../external/npy-matlab/npy-matlab');
addpath('../external/ndSparse_G4_2021_03_16');
addpath('../external/colorcet');

%% plot correlation and learning rate curves

%four units that learned
unitVec = [58,68,148,36];
runVec = ["2024-11-18_11085_LMAN_nCAF_aligned.mat";
    "2024-11-24_11085_LMAN_nCAF_aligned.mat";
    "2025-11-18_Yellow33_Post-Advance_LMAN_Precise-Alignment_BOTM_1_aligned.mat";
    "2025-11-19_Yellow33_Post-Advance_LMAN_Precise-Alignment_BOTM_0_aligned.mat"];

ww = 0.02; %plot width
downsampW = .0005; %downsample to this timebin width (0.5 ms)

for k = 1:numel(unitVec)
    unitI = unitVec(k);
    % load dataset if not already
    load(fullfile('../presorted_data', runVec(k)));
    unitSignal = full(unitSigSparse);

    g = figure;
    title(k);
    plotAx = subplot(2, 1, 1);
    hold on;
    [tvals, frateCorr, confC] = contingency_correlation_uncertainty(unitSignal, np_fs, snipDelay, 0.005, unitI, unitNum, isNoise, plotAx, true, ww, round(downsampW*np_fs));
    yline(0);
    
    plotAx = subplot(2, 1, 2);
    [tvals, learnR, confR] = response_width_fitUncertainty(unitSignal, np_fs, snipDelay, .005, unitI, unitNum, plotAx, true, ww, round(downsampW*np_fs));
    yline(0);
    shg;
end


%% fit gaussians and calculate uncertainties

%list of units to analyze
%selected based on sufficient correlation and learning amplitude
unitVec = [58,59,61,62,65,63,64,68,71,75,148,153,38,36,214];
runVec = [repmat("2024-11-18_11085_LMAN_nCAF_aligned.mat", 3, 1);
    repmat("2024-11-21_11085_LMAN_nCAF_aligned.mat", 2, 1);
    repmat("2024-11-24_11085_LMAN_nCAF_aligned.mat", 5, 1);
    repmat("2025-11-18_Yellow33_Post-Advance_LMAN_Precise-Alignment_BOTM_1_aligned.mat", 2, 1);
    repmat("2025-11-19_Yellow33_Post-Advance_LMAN_Precise-Alignment_BOTM_0_aligned.mat", 2, 1);
    repmat("2025-12-05_11384_Post-Advance_LMAN_Precise-Timing_BOTM_1_aligned.mat", 2, 1)];

Nsim = 1000; %bootstrap replicates
ft = fittype('a*exp(-(x-b)^2/(2*c^2))+d'); %gaussian fit
run_name = "nonexistent_run";
ww = 0.02; %time width to consider
downsampW = 3; %downsample to 100 us
maxW = .001; %width to average over for estimate of maximum (s)

learnMu = zeros(numel(unitVec), 1);
corrMu = zeros(numel(unitVec), 1);
learnSigma = zeros(numel(unitVec), 1);
corrSigma = zeros(numel(unitVec), 1);

simLearnMu = zeros(numel(unitVec), Nsim);
simCorrMu = zeros(numel(unitVec), Nsim);
simLearnSigma = zeros(numel(unitVec), Nsim);
simCorrSigma = zeros(numel(unitVec), Nsim);

for k = 1:numel(unitVec)
    % load dataset if not already
    if ~contains(runVec(k), run_name)
        load(fullfile('../presorted_data', runVec(k)));
        unitSignal = full(unitSigSparse);
    end

    %fit gaussian to learning rates
    [tvals, learnR, confR] = response_width_fitUncertainty(unitSignal, np_fs, snipDelay, .005, unitVec(k), unitNum, [], true, ww, downsampW);
    f=fit(tvals', learnR, ft, 'StartPoint', [max(movmean(learnR, round(maxW*np_fs)/downsampW)), 3, 4.5, 0], 'Lower', [0, -Inf, 0, -Inf]);
    learnParams = coeffvalues(f);
    learnMu(k) = learnParams(2);
    learnSigma(k) = learnParams(3);

    %fit gaussian to correlations
    [tvals, frateCorr, confC] = contingency_correlation_uncertainty(unitSignal, np_fs, snipDelay, 0.005, unitVec(k), unitNum, isNoise, [], true, ww, downsampW);
    f=fit(tvals', frateCorr, ft, 'StartPoint', [max(movmean(frateCorr, round(maxW*np_fs)/downsampW)), 0, 3, 0], 'Lower', [0, -Inf, 0, -Inf]);
    corrParams = coeffvalues(f);
    corrMu(k) = corrParams(2);
    corrSigma(k) = corrParams(3);

    %bootstrap fits
    simData = zeros(size(learnR));
    simLearnParams = zeros(4, Nsim);
    simCorrParams = zeros(4, Nsim);
    for i = 1:Nsim
        for j = 1:numel(simData)
            simData(j) = normrnd(learnR(j), (confR(j, 2)-confR(j, 1))/4);
        end
        f=fit(tvals', simData, ft, 'StartPoint', [max(movmean(learnR, round(maxW*np_fs)/downsampW)), 3, 4.5, 0], 'Lower', [0, -Inf, 0, -Inf]);
        simLearnParams(:, i) = coeffvalues(f);
    
        for j = 1:numel(simData)
            simData(j) = normrnd(frateCorr(j), (confC(j, 2)-confC(j, 1))/4);
        end
        f=fit(tvals', simData, ft, 'StartPoint', [max(movmean(frateCorr, round(maxW*np_fs)/downsampW)), 0, 3, 0], 'Lower', [0, -Inf, 0, -Inf]);
        simCorrParams(:, i) = coeffvalues(f);
    end

    simLearnMu(k, :) = simLearnParams(2, :);
    simCorrMu(k, :) = simCorrParams(2, :);
    simLearnSigma(k, :) = simLearnParams(3, :);
    simCorrSigma(k, :) = simCorrParams(3, :);

    disp(k);
end


%% scatterplots of correlation and learning fits

learnMuErr = prctile(simLearnMu, [2.5,97.5], 2);
corrMuErr = prctile(simCorrMu, [2.5,97.5], 2);
learnSigmaErr = prctile(simLearnSigma, [2.5,97.5], 2);
corrSigmaErr = prctile(simCorrSigma, [2.5,97.5], 2);

figure;
subplot(1, 2, 1);
errorbar(corrMu, corrSigma, corrSigma-corrSigmaErr(:, 1), corrSigmaErr(:, 2)-corrSigma, corrMu-corrMuErr(:, 1), corrMuErr(:, 2)-corrMu, 'o', 'MarkerFaceColor', 'auto', 'CapSize', 0);
xline(0, '--k');
yline(0, '--k');
axis square;
xlabel('correlation center (ms)');
ylabel('correlation width (ms)');
title('Correlation Parameters');
set(gca, 'TickDir', 'out');
axis square;

subplot(1, 2, 2)
errorbar(learnMu, learnSigma, learnSigma-learnSigmaErr(:, 1), learnSigmaErr(:, 2)-learnSigma, learnMu-learnMuErr(:, 1), learnMuErr(:, 2)-learnMu, 'o', 'MarkerFaceColor', 'auto', 'CapSize', 0);
xline(0, '--k');
yline(0, '--k');
axis square;
xlabel('learning center (ms)');
ylabel('learning width (ms)');
title('Learning Parameters');
set(gca, 'TickDir', 'out');
axis square;