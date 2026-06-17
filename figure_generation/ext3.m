%% import functions

addpath('../utils')
addpath('../external/SpikeGLX_Datafile_Tools/MATLAB');
addpath('../external/npy-matlab/npy-matlab');
addpath('../external/ndSparse_G4_2021_03_16');
addpath('../external/colorcet');

%% calculate autocorrelations for range of neurons

%set of well isolated single units with rp_violation < 0.02
unitVec = [73,242,253,256,764,779,785,787,815,858,868,873,893];
runVec = [repmat("2025-11-12_Yellow33_Post-Advance_LMAN_BOTM_0_aligned.mat", 4, 1);
    repmat("2025-02-18_11208_Post-Advance_LMAN_BOTM_0_aligned.mat", 9, 1)];

run_name = "nonexistent_run";
Nmotif = 100; %calculate over first N motifs
blankW = 0.01; %ignore data within this width of contingent window
maxLag = 0.1; %maximum lag for autocorrelation calculation (s)
binW = 3; %number of timepoints to bin, downsampling to 100 us

load(fullfile('../presorted_data', runVec(1)));
unitSignal = full(unitSigSparse);
acorrTot = zeros(2*round(maxLag*np_fs/binW)+1, numel(unitVec));

for k = 1:numel(unitVec)
    
    % load dataset if not already
    if ~contains(runVec(k), run_name)
        load(fullfile('../presorted_data', runVec(k)));
        unitSignal = full(unitSigSparse);
    end
    
    %calculate mean song-locked signal
    meanSig = mean(unitSignal(:, 1:Nmotif, unitNum==unitVec(k)), 2);
    %calculate autocorrelation of variation from song-locked signal
    [acorr, meanAcorr, lags] = spike_train_autocorrelation(unitSignal(:, 1:Nmotif, :)-meanSig, np_fs, unitNum, unitVec(k), maxLag, binW);
    acorrTot(:, k) = meanAcorr;
end

%plot full autocorrelation
windowW = 0.1; %width of window in ms, for calculating density
acorrPlot = mean(acorrTot, 2);
figure;
plot(lags*1000, acorrPlot/windowW);
xlim([-20 20]);
ylim([-Inf, Inf]);
xlabel('time lag (ms)');
ylabel('correlation/time (ms^{-1})');
set(gca, 'TickDir', 'out');

%plot zoomed-in sidelobe
acorrPlot_noDelta = acorrPlot;
acorrPlot_noDelta(abs(lags)<.0005)=0;
figure;
plot(lags*1000, acorrPlot_noDelta/windowW);
xlim([1, 20]);
ylim([-Inf Inf]);
xlabel('time lag (ms)');
ylabel('correlation/time (ms^{-1})');
set(gca, 'TickDir', 'out');
yline(0)

%% compute median correlation profile of actual units
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
run_name = "nonexistent_run"; %initial value
ww = 0.02; %width of times to consider
downsampW = .0005; %downsample data to this time resolution (ms)
maxW = .001; % window width used to smooth for initial amplitude estimate

FRCorrs = cell(numel(unitVec),1);
SEs = cell(numel(unitVec),1);

%iterate through units
for k = 1:numel(unitVec)
    % load dataset if not already
    if ~contains(runVec(k), run_name)
        load(fullfile('../presorted_data', runVec(k)));
        unitSignal = full(unitSigSparse);
    end

    %calculate correlation with DAF
    [tvalsData, frateCorr, confC] = contingency_correlation_uncertainty(unitSignal, np_fs, snipDelay, cW, unitVec(k), unitNum, isNoise, [], true, ww, round(downsampW*np_fs));

    FRCorrs{k} = frateCorr;
    SEs{k} = diff(confC,1,2) ./ 3.92;

end

% Compute mean corr and CI
FRCorrs= cat(2, FRCorrs{:});
SEs = cat(2, SEs{:});

% Calculate estimates
N = numel(unitVec);
medianFRCorr = median(FRCorrs, 2);
medianFRCI = bootci(1000, @(x) median(x, 1), FRCorrs.').';

%% convolve contingent window with autocorrelation

ww = 20;
tvals = -ww:0.1:ww;

%approximate contingent window vector
tW = 2.5;
k_on = 5; %blurring factor
sigmoidCorr = (1./(1+exp(-k_on*(tvals+tW)))).*(1./(1+exp(k_on*(tvals-tW))));

%convolve with autocorrelation
convCorr = conv(sigmoidCorr, acorrPlot, 'same');
convCorr = convCorr/max(convCorr);

%% scale outputs to each other
% We now scale the data to best fit the convCorr
convCorrResampled = interp1(tvals, convCorr, tvalsData, 'linear', 'extrap');
scaleFactor = convCorrResampled(:) \ medianFRCorr(:);

%% plot results
% Colors to use
sigmoidCorrColor = [0,0.447,0.698];    % Blue
convCorrColor = [0.835,0.369,0];       % Orange
measuredCorrColor = [0.8,0.475,0.655]; % Pink

% Make the plot
figure;
hold on;
plot(tvals, sigmoidCorr, "Color", sigmoidCorrColor);
plot(tvals, convCorr, "Color", convCorrColor);
fill([tvalsData, fliplr(tvalsData)], ...
    [medianFRCI(:,1).', fliplr(medianFRCI(:,2).')] / scaleFactor, ...
    measuredCorrColor, ...
    'FaceAlpha', 0.2, ...
    'LineStyle', 'none');
plot(tvalsData, medianFRCorr / scaleFactor, ...
    'Color', measuredCorrColor);
xline([-2.5 2.5], '--k');
xlabel('time (ms)');
ylabel('correlation with DAF (relative)');
set(gca, 'TickDir', 'out');
hold off;

%fit gaussian to convolved correlations
ft = fittype('a*exp(-(x-b)^2/(2*c^2))');
f=fit(tvals', convCorr', ft, 'StartPoint', [1, 0, 2]);
fConf = confint(f);
disp("gaussian fit width = "+f.c+" +/- "+diff(fConf(:, 3))/2);