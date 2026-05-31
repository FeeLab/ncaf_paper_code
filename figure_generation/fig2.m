

%% add functions to path

addpath('../utils')
addpath('../external/SpikeGLX_Datafile_Tools/MATLAB');
addpath('../external/npy-matlab/npy-matlab');
addpath('../external/ndSparse_G4_2021_03_16');
addpath('../external/colorcet');

%% plot correlation and learning rate curves

%run name
load('../presorted_data/2024-11-21_11085_LMAN_nCAF_aligned.mat');
unitSignal = full(unitSigSparse);

unitI = 65; %example unit
ww = 0.02; %width of plot
downsampW = .0005; %downsample data to this time resolution (ms)
g = figure;

%plot correlation with DAF
plotAx = subplot(2, 1, 1);
hold on;
[tvals, frateCorr, confC] = contingency_correlation_uncertainty(unitSignal, np_fs, snipDelay, cW, unitI, unitNum, isNoise, plotAx, true, ww, round(downsampW*np_fs));
yline(0);

%plot learning rate
plotAx = subplot(2, 1, 2);
[tvals, learnR, confR] = response_width_fitUncertainty(unitSignal, np_fs, snipDelay, cW, unitI, unitNum, plotAx, true, ww, round(downsampW*np_fs));
yline(0);

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

Nsim = 1000; %number of bootstrap iterations
ft = fittype('a*exp(-(x-b)^2/(2*c^2))+d'); %gaussian function to fit
run_name = "nonexistent_run"; %initial value
ww = 0.02; %width of times to consider
downsampW = 3; %downsample to 100 microseconds
maxW = .001; % window width used to smooth for initial amplitude estimate

learnMu = zeros(numel(unitVec), 1);
corrMu = zeros(numel(unitVec), 1);
learnSigma = zeros(numel(unitVec), 1);
corrSigma = zeros(numel(unitVec), 1);

simLearnMu = zeros(numel(unitVec), Nsim);
simCorrMu = zeros(numel(unitVec), Nsim);
simLearnSigma = zeros(numel(unitVec), Nsim);
simCorrSigma = zeros(numel(unitVec), Nsim);

%iterate through units
for k = 1:numel(unitVec)
    % load dataset if not already
    if ~contains(runVec(k), run_name)
        load(fullfile('../presorted_data', runVec(k)));
        unitSignal = full(unitSigSparse);
    end
    
    %calculate learning rates
    [tvals, learnR, confR] = response_width_fitUncertainty(unitSignal, np_fs, snipDelay, cW, unitVec(k), unitNum, [], true, ww, downsampW);
    %fit gaussian model to learning rates
    %use approximate values as starting points
    f=fit(tvals', learnR, ft, 'StartPoint', [max(movmean(learnR, round(maxW*np_fs)/downsampW)), 3, 4.5, 0], 'Lower', [0, -Inf, 0, -Inf]);
    learnParams = coeffvalues(f);
    learnMu(k) = learnParams(2);
    learnSigma(k) = learnParams(3);

    %calculate correlation with DAF
    [tvals, frateCorr, confC] = contingency_correlation_uncertainty(unitSignal, np_fs, snipDelay, cW, unitVec(k), unitNum, isNoise, [], true, ww, downsampW);
    %fit gaussian model to learning rates
    %use approximate values as starting points
    f=fit(tvals', frateCorr, ft, 'StartPoint', [max(movmean(frateCorr, round(maxW*np_fs)/downsampW)), 0, 3, 0], 'Lower', [0, -Inf, 0, -Inf]);
    corrParams = coeffvalues(f);
    corrMu(k) = corrParams(2);
    corrSigma(k) = corrParams(3);

    %bootstrap
    simData = zeros(size(learnR));
    simLearnParams = zeros(4, Nsim);
    simCorrParams = zeros(4, Nsim);
    %randomly jitter curves and refit gaussians
    for i = 1:Nsim
        for j = 1:numel(simData)
            %add random jitter to learning rate
            simData(j) = normrnd(learnR(j), (confR(j, 2)-confR(j, 1))/4);
        end
        %fit to jittered data
        f=fit(tvals', simData, ft, 'StartPoint', [max(movmean(learnR, round(maxW*np_fs)/downsampW)), 0, 5, 0], 'Lower', [0, -Inf, 0, -Inf]);
        simLearnParams(:, i) = coeffvalues(f);
    
        for j = 1:numel(simData)
            %add random jitter to correlation
            simData(j) = normrnd(frateCorr(j), (confC(j, 2)-confC(j, 1))/4);
        end
        %fit to jittered data
        f=fit(tvals', simData, ft, 'StartPoint', [max(movmean(frateCorr, round(maxW*np_fs)/downsampW)), 0, 5, 0], 'Lower', [0, -Inf, 0, -Inf]);
        simCorrParams(:, i) = coeffvalues(f);
    end

    simLearnMu(k, :) = simLearnParams(2, :);
    simCorrMu(k, :) = simCorrParams(2, :);
    simLearnSigma(k, :) = simLearnParams(3, :);
    simCorrSigma(k, :) = simCorrParams(3, :);

    disp(k);
end

%infer kernel parameters
kernelMu = learnMu-corrMu;
kernelSigma = real((sqrt(learnSigma.^2-corrSigma.^2)));
simKernelMu = simLearnMu-simCorrMu;
simKernelSigma = real((sqrt(simLearnSigma.^2-simCorrSigma.^2)));

%fit line to bootstrapped offsets
fMu = fittype('a+x');
muVals = zeros(Nsim, 1);
for i = 1:Nsim
    f=fit(simCorrMu(:, i), simLearnMu(:, i), fMu, 'StartPoint', mean(kernelMu), 'Lower', 0);
    muVals(i) = f.a;
end

%fit hyperbola to bootstrapped widths
fSigma = fittype('sqrt(a^2+x^2)');
sigmaVals = zeros(Nsim, 1);
for i = 1:Nsim
    f=fit(simCorrSigma(:, i), simLearnSigma(:, i), fSigma, 'StartPoint', mean(kernelSigma(real(kernelSigma)>0)), 'Lower', 0);
    sigmaVals(i) = f.a;
end

%% scatterplots of correlation and learning fits

%calculate confidence intervals from bootstrap results
learnMuErr = prctile(simLearnMu, [2.5,97.5], 2);
corrMuErr = prctile(simCorrMu, [2.5,97.5], 2);
learnSigmaErr = prctile(simLearnSigma, [2.5,97.5], 2);
corrSigmaErr = prctile(simCorrSigma, [2.5,97.5], 2);

figure;
%scatterplot of offsets
subplot(1, 2, 1);
errorbar(corrMu, learnMu, learnMu-learnMuErr(:, 1), learnMuErr(:, 2)-learnMu, corrMu-corrMuErr(:, 1), corrMuErr(:, 2)-corrMu, 'o', 'MarkerFaceColor', 'auto', 'CapSize', 0);
xline(0, '--k');
yline(0, '--k');
axis square;
xlabel('correlation center (ms)');
ylabel('learning center (ms)');
title('Gaussian \mu');
set(gca, 'TickDir', 'out');
axis square;

%scatterplot of widths
subplot(1, 2, 2)
errorbar(corrSigma, learnSigma, learnSigma-learnSigmaErr(:, 1), learnSigmaErr(:, 2)-learnSigma, corrSigma-corrSigmaErr(:, 1), corrSigmaErr(:, 2)-corrSigma, 'o', 'MarkerFaceColor', 'auto', 'CapSize', 0);
axis square;
hold on;
xVals = 0:.1:10;
%plot hyperbola of best fit
plot(xVals, sqrt(median(sigmaVals)^2+xVals.^2));
plot(xVals, sqrt(prctile(sigmaVals, 2.5)^2+xVals.^2), '--k');
plot(xVals, sqrt(prctile(sigmaVals, 97.5)^2+xVals.^2), '--k');
line([0 100], [0 100]);
xlabel('correlation width (ms)');
ylabel('learning width (ms)');
title('Gaussian \sigma');
set(gca, 'TickDir', 'out');
xlim([0 max([corrSigmaErr(:, 2); learnSigmaErr(:, 2)])]);
ylim([0 max([corrSigmaErr(:, 2); learnSigmaErr(:, 2)])]);


%output values in text
disp("mean delay = "+mean(kernelMu)+" +/- "+std(kernelMu));
disp("mean width = "+mean(sigmaVals)+" +/- "+std(sigmaVals));

disp("correlation width = "+mean(corrSigma)+" +/- "+std(corrSigma));
disp("learning width = "+mean(learnSigma)+" +/- "+std(learnSigma));
%% individual kernel parameters

%plotting parameters
groupSpace = 1;
groupW = .3;

%plot inferred kernal parameters for each unit individually
figure;
bar([0 1], [mean(kernelMu) median(sigmaVals)]);
hold on;
errorbar(linspace(-groupW, groupW, numel(kernelMu)), kernelMu, kernelMu-prctile(simKernelMu, 2.5, 2), prctile(simKernelMu, 97.5, 2)-kernelMu, 'o', 'MarkerFaceColor', 'auto', 'CapSize', 0);
errorbar(linspace(-groupW, groupW, numel(kernelSigma))+groupSpace, kernelSigma, kernelSigma-prctile(simKernelSigma, 2.5, 2), prctile(simKernelSigma, 97.5, 2)-kernelSigma, 'o', 'MarkerFaceColor', 'auto', 'CapSize', 0);


