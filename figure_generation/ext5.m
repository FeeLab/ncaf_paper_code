%% import functions

addpath('../utils')
addpath('../external/SpikeGLX_Datafile_Tools/MATLAB');
addpath('../external/npy-matlab/npy-matlab');
addpath('../external/ndSparse_G4_2021_03_16');
addpath('../external/colorcet');


%% neuron matrix correlations

%dataset with many well-isolated single units
load('../presorted_data/2025-02-18_11208_Post-Advance_LMAN_BOTM_0_aligned.mat');
unitSignal = full(unitSigSparse);

tw = .020; %calculate correlations over time windows with this width (s)
Nmax = 100; %calculate for first N motifs
Nshuff = 10; %shuffles to calculate significant thresholds

%estimate song-locked activity
filtSignal = smoothdata(unitSignal(:, 1:Nmax, ismember(unitNum, lmanNum(lmanSingleUnit))), 1, 'movmean', tw*np_fs);
meanSignal = mean(filtSignal, 2);

%mean-subtract signal
timeseries = zeros(size(filtSignal, 1)*Nmax, size(filtSignal, 3));
shuffTimeseries = zeros(size(filtSignal, 1)*Nmax, Nshuff, size(filtSignal, 3));
for i = 1:Nmax
    timeseries((i-1)*size(filtSignal, 1)+1:i*size(filtSignal, 1), :) = filtSignal(:, i, :)-meanSignal;
end

%calculate pairwise correlations
N = size(timeseries, 2);
corrVals = zeros(N, N);
for i = 1:N
    for j = 1:N
        corrVals(i, j) = corr(timeseries(:, i), timeseries(:, j));
    end
end

%calculate shuffled pairwise correlations
corrValsShuff = zeros(N, N, Nshuff);
for k = 1:Nshuff
    shuffTimeseries = zeros(size(filtSignal, 1)*Nmax, size(filtSignal, 3));
    iShuff = randperm(Nmax);
    %shuffle over motifs
    for i = 1:Nmax
        shuffTimeseries((i-1)*size(filtSignal, 1)+1:i*size(filtSignal, 1), :) = filtSignal(:, iShuff(i), :)-meanSignal;
    end

    %calculate correlations
    for i = 1:N
        for j = 1:N
            corrValsShuff(i, j, k) = corr(timeseries(:, i), shuffTimeseries(:, j));
        end
    end
    disp(k);
end

%unit locations
unitLocs = unitCOM(ismember(unitNum, lmanNum(lmanSingleUnit)), :);

%calculate pairwise distances
distMatrix = zeros(size(corrVals));
for i = 1:N
    for j = 1:N
        distMatrix(i, j) = sqrt(sum((unitLocs(i, :) - unitLocs(j, :)).^2));
    end
end

%remove duplicates
upperVals = logical(triu(ones(N), 1));

%plot results
figure;
Nsamp = sum(upperVals(:));
sigLim = [prctile(corrValsShuff(:), 2.5/Nsamp), prctile(corrValsShuff(:), 100-2.5/Nsamp)];
sigI = corrVals<sigLim(1) | corrVals>sigLim(2);
scatter(distMatrix(~sigI & upperVals), corrVals(~sigI & upperVals), [], 'MarkerEdgeColor', 'none', 'MarkerFaceColor',[1 1 1]*.5, 'MarkerFaceAlpha',0.2);
hold on;
scatter(distMatrix(sigI & upperVals), corrVals(sigI & upperVals), [], 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.835,0.369,0], 'MarkerFaceAlpha',0.4);
xlabel('pairwise distance (microns)');
ylabel('correlation');
yline(0, '--k')
set(gca, 'TickDir', 'out');
    
%% calculate unit correlations with noise

targetI = 868; %target unit in this dataset
unitOrigin = unitCOM(unitNum==targetI, :);
unitDist = sqrt(sum((unitCOM(lmanSingleUnit, :)-repmat(unitOrigin, sum(lmanSingleUnit), 1)).^2, 2));
lmanSingleI = find(ismember(unitNum, lmanNum(lmanSingleUnit)));

spikeCount = contingent_spikes(unitSignal(:,:,lmanSingleI), np_fs, trigDelay, cW, false);
[phi, phiConf] = activity_noise_corr(spikeCount, isNoise);

%plot correlation with DAF vs. correlation with target
figure;
hold on;
scatter(corrVals(unitNum(lmanSingleI)==targetI, unitNum(lmanSingleI)~=targetI), phi(unitNum(lmanSingleI)~=targetI), 'filled', 'MarkerEdgeColor', 'none');
xline(0, '--k');
yline(0, '--k');
xlabel('correlation with target neuron');
ylabel('correlation with DAF');
%fit linear model
f = fit(corrVals(unitNum(lmanSingleI)==targetI, unitNum(lmanSingleI)~=targetI)', phi(unitNum(lmanSingleI)~=targetI), 'poly1');
xLims = xlim;
xvals = linspace(xLims(1), xLims(2), 100);
ci = predint(f, xvals, .95, 'functional', 'on');
confInt = confint(f);
plot(xvals, f(xvals));
plot(xvals, ci, 'k');
set(gca, 'TickDir', 'out');
disp("beta = "+f.p1+" +/- "+diff(confInt(:, 1))/2);

%plot correlation with target vs distance from target
figure;
scatter(distMatrix(unitNum(lmanSingleI)==targetI, unitNum(lmanSingleI)~=targetI), corrVals(unitNum(lmanSingleI)==targetI, unitNum(lmanSingleI)~=targetI), 'filled', 'MarkerEdgeColor', 'none');
yline(0, '--k');
xlabel('distance from target neuron');
ylabel('correlation with target neuron');
set(gca, 'TickDir', 'out');

%plot correlation with noise vs distance from target
figure;
scatter(distMatrix(unitNum(lmanSingleI)==targetI, unitNum(lmanSingleI)~=targetI), phi(unitNum(lmanSingleI)~=targetI), 'filled', 'MarkerEdgeColor', 'none');
yline(0, '--k');
xlabel('distance from target neuron');
ylabel('correlation with DAF');
set(gca, 'TickDir', 'out');

%% load single units

%four targeted single units
runVec = ["2025-02-18_11208_Post-Advance_LMAN_BOTM_0_aligned.mat";
    "2025-02-21_11208_Post-Advance_LMAN_nCAF_0_aligned.mat";
    "2025-11-12_Yellow33_Post-Advance_LMAN_BOTM_0_aligned.mat";
    "2024-04-27_10872_LMAN-X_nCAF_2_aligned.mat"];
targetI = [868,127,73,187];

phiConfTot = [];
confRTot = [];
xdata = [];
learnRTot = [];
unitDistTot=[];
unitNumTot=[];
xdata = [];
isSingleUnit = zeros(0, 'logical');
correctBaseline = true;

for k = 1:numel(targetI)

    load(fullfile('../presorted_data', runVec(k)));
    unitSignal = full(unitSigSparse);

    [learnR, confR] = learning_per_neuron(unitSignal(:,:,lmanUnits), np_fs, trigDelay, cW, correctBaseline);
    spikeCount = contingent_spikes(unitSignal(:,:,lmanUnits), np_fs, trigDelay, cW, correctBaseline);
    [phi, phiConf] = activity_noise_corr(spikeCount, isNoise);
    flucFrac = fr_fluctuations(unitSignal(:, :, lmanUnits), np_fs, trigDelay, cW);
    lman_rpViolation = rp_violation(lmanUnits);
    
    %keep units with sufficiently small baseline fluctuations
    flucLim = .25;
    goodUnit = (flucFrac<flucLim&~isnan(learnR)&~isnan(phi));

    unitOrigin = unitCOM(unitNum==targetI(k), :);
    unitDist = sqrt(sum((unitCOM(lmanUnits, :)-repmat(unitOrigin, sum(lmanUnits), 1)).^2, 2));

    phiConfTot = [phiConfTot; phiConf(goodUnit, :)];
    confRTot = [confRTot; confR(goodUnit, :)];
    learnRTot = [learnRTot; learnR(goodUnit)];
    unitDistTot = [unitDistTot; unitDist(goodUnit)];
    unitNumTot = [unitNumTot; lmanNum(goodUnit)];
    xdata = [xdata; [phi(goodUnit) unitCOM(ismember(unitNum, lmanNum(goodUnit)), :) ones(sum(goodUnit), 1)*k]];
    isSingleUnit = [isSingleUnit; lmanSingleUnit(goodUnit)];
end


phiConfTot = phiConfTot(isSingleUnit, :);
confRTot = confRTot(isSingleUnit, :);
learnRTot = learnRTot(isSingleUnit);
unitDistTot = unitDistTot(isSingleUnit);
xdata = xdata(isSingleUnit, :);
unitNumTot = unitNumTot(isSingleUnit);
%% fit parametric model to correlations

[a,b,c,d,e,mu,sig] = fit_recorded_phi_normed(unitDistTot, xdata);