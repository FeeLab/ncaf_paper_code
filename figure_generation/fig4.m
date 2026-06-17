%% Set random seed
rng(0);

%% import functions

addpath('../utils')
addpath('../external/SpikeGLX_Datafile_Tools/MATLAB');
addpath('../external/npy-matlab/npy-matlab');
addpath('../external/ndSparse_G4_2021_03_16');
addpath('../external/colorcet');

%% plot trendlines

%import data
load('../presorted_data/2025-11-12_Yellow33_Post-Advance_LMAN_BOTM_0_aligned.mat');
unitSignal = full(unitSigSparse);

%correct for baseline fluctuations due to recording drift
correctBaseline = true;

figure;
hold on;
plotAx = gca;
%plot targeted unit
plot_learning_trendlines(unitSignal, np_fs, trigDelay, 73, unitNum, correctBaseline, plotAx);
%plot two untargeted units that learn
plot_learning_trendlines(unitSignal, np_fs, trigDelay, 71, unitNum, correctBaseline, plotAx);
plot_learning_trendlines(unitSignal, np_fs, trigDelay, 76, unitNum, correctBaseline, plotAx);



%% load single and multiunit

%four targeted single units
targetI = [868,127,73,187];
runVec = ["2025-02-18_11208_Post-Advance_LMAN_BOTM_0_aligned.mat";
    "2025-02-21_11208_Post-Advance_LMAN_nCAF_0_aligned.mat";
    "2025-11-12_Yellow33_Post-Advance_LMAN_BOTM_0_aligned.mat";
    "2024-04-27_10872_LMAN-X_nCAF_2_aligned.mat"];

phiConfTot = []; %confidence intervals for correlation
confRTot = []; %confidence itnervals for learning rate
learnRTot = []; %learning rates
unitDistTot=[]; %distances from target neuron
unitNumTot=[]; %unit indices
%xdata is matrix containing correlations, 3D unit coordinates, and runID
xdata = []; %holdover from fitting parametric functions and needing to pass in data
isSingleUnit = zeros(0, 'logical'); %single unit
correctBaseline = true; %correct for baseline fluctuations due to recording drift

for k = 1:numel(targetI)

    load(fullfile('../presorted_data', runVec(k)));
    unitSignal = full(unitSigSparse);

    %calculate learning in contingent window
    [learnR, confR] = learning_per_neuron(unitSignal(:,:,lmanUnits), np_fs, trigDelay, cW, correctBaseline);
    %compute spikes in contingent window
    spikeCount = contingent_spikes(unitSignal(:,:,lmanUnits), np_fs, trigDelay, cW, correctBaseline);
    %compute correlation with noise
    [phi, phiConf] = activity_noise_corr(spikeCount, isNoise);
    
    %comptue background firing rate fluctuations
    flucFrac = fr_fluctuations(unitSignal(:, :, lmanUnits), np_fs, trigDelay, cW);
    %import refractory period violation metric
    lman_rpViolation = rp_violation(lmanUnits);
    
    flucLim = .25; %reject units with baseline fluctuations larger than limit
    rpLim = .3; %reject multiunits with refractory period violations larger than limit
    goodUnit = (flucFrac<flucLim&lman_rpViolation<rpLim&~isnan(learnR)&~isnan(phi));

    %set spatial origin to target neuron
    unitOrigin = unitCOM(unitNum==targetI(k), :);
    %compute euclidian distances to target neuron
    unitDist = sqrt(sum((unitCOM(lmanUnits, :)-repmat(unitOrigin, sum(lmanUnits), 1)).^2, 2));

    phiConfTot = [phiConfTot; phiConf(goodUnit, :)];
    confRTot = [confRTot; confR(goodUnit, :)];
    learnRTot = [learnRTot; learnR(goodUnit)];
    unitDistTot = [unitDistTot; unitDist(goodUnit)];
    unitNumTot = [unitNumTot; lmanNum(goodUnit)];
    xdata = [xdata; [phi(goodUnit) unitCOM(ismember(unitNum, lmanNum(goodUnit)), :) ones(sum(goodUnit), 1)*k]];
    isSingleUnit = [isSingleUnit; lmanSingleUnit(goodUnit)];
end

%% plot learning vs correlation

%sort for plotting order
[~, sortOrder] = sort(unitDistTot, 'descend');
plotSingleUnit = ismember(sortOrder, find(isSingleUnit));

figure;
hold on;
%plot units in trendlines
scatter(xdata(unitNumTot==73, 1), learnRTot(unitNumTot==73), 100, [0.835,0.369,0], 'LineWidth', 2);
scatter(xdata(unitNumTot==71, 1), learnRTot(unitNumTot==71), 100, [0,0.620,0.451], 'LineWidth', 2);
scatter(xdata(unitNumTot==76, 1), learnRTot(unitNumTot==76), 100, [0.8,0.475,0.655], 'LineWidth', 2);
%plot multiunit
scatter(xdata(sortOrder(~plotSingleUnit), 1), learnRTot(sortOrder(~plotSingleUnit)), 50, unitDistTot(sortOrder(~plotSingleUnit)), 'filled', 'LineWidth', 2, 'MarkerEdgeColor', 'flat', 'MarkerFaceColor','none');
%plot single unit
scatter(xdata(sortOrder(plotSingleUnit), 1), learnRTot(sortOrder(plotSingleUnit)), 60, unitDistTot(sortOrder(plotSingleUnit)), 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor','flat');
xlabel('correlation between unit and noise');
ylabel('learning rate (Hz/trial)');
xline(0, '--k');
yline(0, '--k');
a = flipud(colorcet('L04'));
colormap(a(end/4:end, :));
cb = colorbar;
cb.Label.String = "distance from target neuron";

%fit and plot line
f = fit(xdata(:, 1), learnRTot, 'poly1');
xvals = min(xdata(:, 1))-.08:.001:max(xdata(:, 1)+.05);
ci = predint(f, xvals, .95, 'functional', 'on');
plot(xvals, f(xvals));
plot(xvals, ci, 'k');
set(gca, 'TickDir', 'out');

%plot limits
lrMin = -0.06;
lrMax = max(learnRTot)+0.01;
corrMin = min(xdata(:, 1))-0.02;
corrMax = max(xdata(:, 1))+0.02;
xlim([corrMin corrMax]);
ylim([lrMin lrMax]);


%output fit parameters
fCI = confint(f);
disp("alpha = "+f.p2+" +/- "+diff(fCI(:, 2))/2);
disp("beta = "+f.p1+" +/- "+diff(fCI(:, 1))/2);
mdl = fitlm(xdata(:, 1), learnRTot);
disp("r-squared = "+mdl.Rsquared.Ordinary);
%% neuron matrix correlations

%dataset with many well-isolated single units
load('../presorted_data/2025-02-18_11208_Post-Advance_LMAN_BOTM_0_aligned.mat');
unitSignal = full(unitSigSparse);

tw = .020; %evaluate pairwise correlations over time windows with this width (20 ms)
Nmax = 100; %compute correlations over first N motifs before feedback

%compute firing rates over sliding window
filtSignal = smoothdata(unitSignal(:, 1:Nmax, ismember(unitNum, lmanNum(lmanSingleUnit))), 1, 'movmean', tw*np_fs);
%song-locked average activity
meanSignal = mean(filtSignal, 2);

%compute trial-by-trial fluctuations from song-locked activity
timeseries = zeros(size(filtSignal, 1)*Nmax, size(filtSignal, 3));
for i = 1:Nmax
    timeseries((i-1)*size(filtSignal, 1)+1:i*size(filtSignal, 1), :) = filtSignal(:, i, :)-meanSignal;
end

%compute pairwise correlations
N = size(timeseries, 2);
corrVals = zeros(N, N);
for i = 1:N
    for j = 1:N
        corrVals(i, j) = corr(timeseries(:, i), timeseries(:, j));
    end
end

%unit locations for single units
unitLocs = unitCOM(ismember(unitNum, lmanNum(lmanSingleUnit)), :);

%compute pairwise distance matrix
distMatrix = zeros(size(corrVals));
for i = 1:N
    for j = 1:N
        distMatrix(i, j) = sqrt(sum((unitLocs(i, :) - unitLocs(j, :)).^2));
    end
end

%remove duplicate values
upperVals = logical(triu(ones(N), 1));
rPair = distMatrix(upperVals);
CPair = corrVals(upperVals);


%% restrict subsequent analysis to single units
phiConfTot = phiConfTot(isSingleUnit, :);
confRTot = confRTot(isSingleUnit, :);
learnRTot = learnRTot(isSingleUnit);
unitDistTot = unitDistTot(isSingleUnit);
xdata = xdata(isSingleUnit, :);
unitNumTot = unitNumTot(isSingleUnit);

%% plot learning vs. distance

phiTot = xdata(:, 1); %correlation values
distC = 400; %maximum plot distance

yMin = min(learnRTot-f(phiTot)); %plot limits

figure;
ax1 = subplot(2, 1, 1);
axis square;
hold on;
%highlight units from trendlines
scatter(unitDistTot(unitNumTot==73, 1), learnRTot(unitNumTot==73), 100, [0.835,0.369,0], 'LineWidth', 2);
scatter(unitDistTot(unitNumTot==71, 1), learnRTot(unitNumTot==71), 100, [0,0.620,0.451], 'LineWidth', 2);
scatter(unitDistTot(unitNumTot==76, 1), learnRTot(unitNumTot==76), 100, [0.8,0.475,0.655], 'LineWidth', 2);
%plot learning rate vs distance
scatter(ax1, unitDistTot, learnRTot, 50, phiTot, "MarkerFaceColor", "flat", "MarkerEdgeColor", [0 0 0]);
cb2 = colorbar(ax1);
colormap(ax1, flipud(colorcet('D02')));
clim(ax1, [-max(abs(phiTot)), max(abs(phiTot))]);
cb2.Label.String = 'correlation with noise';
set(ax1, 'TickDir', 'out');
yline(ax1, 0, '--k');
xlim([0 distC]);
ylim([lrMin lrMax]);
xlabel('distance from target neuron');
ylabel('learning rate');

%fit exponential to data
fExp = fit(unitDistTot, learnRTot, 'exp1');
xVals = 0:distC;
plot(xVals, fExp(xVals));
ci = predint(fExp, xVals, .95, 'functional', 'on');
plot(xVals, ci, 'k');

%fit line to learning rate vs. correlation
f = fit(phiTot, learnRTot, 'poly1');
ax2 = subplot(2, 1, 2);
axis square;
hold on;
%regress out learning rates attributed to correlation with DAF
scatter(ax2, unitDistTot, learnRTot-f(phiTot), 50, phiTot, "MarkerFaceColor", "flat", "MarkerEdgeColor", [0 0 0]);
cb2 = colorbar(ax2);
colormap(ax2, flipud(colorcet('D02')));
clim(ax2, [-max(abs(phiTot)), max(abs(phiTot))]);
cb2.Label.String = 'correlation with noise';
set(ax2, 'TickDir', 'out');
yline(ax2, 0, '--k');
xlim([0 distC]);
ylim([lrMin lrMax]);
xlabel('distance from target neuron');
ylabel('residual learning rate');

%fit exponential to residual learning rates
fExpResid = fit(unitDistTot, learnRTot-f(phiTot), 'exp1', 'Lower', [-Inf fExp.b], 'Upper', [Inf fExp.b]);
plot(xVals, fExpResid(xVals));
ci = predint(fExpResid, xVals, .95, 'functional', 'on');
plot(xVals, ci, 'k');

%compute confidence intervals
fExpCI = confint(fExp);
fExpResidCI = confint(fExpResid);
widthCI = (1/(fExp.b+fExpCI(1, 2))-1/(fExp.b+fExpCI(2, 2)))/2;
sigmaExp = diff(fExpCI(:, 1));
sigmaResid = diff(fExpResidCI(:, 1));
sigmaRatio = -sigmaResid/fExp.a + sigmaExp*fExpResid.a/fExp.a^2; %propagate uncertainty in ratio of amplitudes

%output values in text
disp("ratio = "+(fExpResid.a/fExp.a)+" +/- "+sigmaRatio);
disp("width (1/e) = "+(-1/fExp.b)+" +/- "+widthCI);



%% simulate LMAN neuron grid

%linear fit of learning rate vs correlation
fEff = fit(xdata(:, 1), learnRTot, 'poly1');

gridW = 20; %grid spacing in microns
xW = gridW*ceil(200/gridW); %x extents
yW = gridW*ceil(200/gridW); %y extents
zW = gridW*ceil(400/gridW); %z extents

%create grid of points
[xCoord, yCoord, zCoord] = meshgrid(-xW:gridW:xW, -yW:gridW:yW, -zW:gridW:zW);
unitLocs = [xCoord(:) yCoord(:) zCoord(:)];

%set origin as target neuron
simTargetI = find(unitLocs(:, 1)==0 & unitLocs(:, 2)==0 & unitLocs(:, 3)==0);

recX = 20; %width of recorded region in x (microns)
recY = 40; %width of recorded region in y (microns)
%define recorded neurons
recordI = unitLocs(:, 1)>=0 & unitLocs(:, 1)<=recX & unitLocs(:, 2)>=-recY & unitLocs(:, 2)<=recY;

%fit parametric model to observed correlations as a function of distance
[a,b,c,d,e,mu,sig] = fit_recorded_phi_normed(unitDistTot, xdata);

%estimate probability distribution of pairwise correlations
[rPair, CPair] = pairwiseCorrelations('../presorted_data/2025-02-18_11208_Post-Advance_LMAN_BOTM_0_aligned.mat');
[pairCorrCDF, xBin, yBin] = pairwise_corr_CDF(rPair, CPair);

jitterDist = 5; %add random jitter to neuron positions with this amplitude (microns)

sigma = 2.^(0:.5:8); %vector of values of projection divergence

Nsamp = sum(recordI); %number of simulated recorded neurons per simulation
Nsim = 10000; %number of simulations
phiSim = zeros(Nsamp*Nsim, 1); %simulated correlations
learnRSim = zeros(Nsamp*Nsim, numel(sigma)); %simulated learning rates
unitDistSim = zeros(Nsamp*Nsim, 1); %distances from target neuron
phiTot = xdata(:, 1);

%initialize gpu arrays
unitLocs = gpuArray(single(unitLocs));
phi = gpuArray(zeros(size(unitLocs, 1), 1, 'single'));

pairwiseCorrF = 0.54; %computed relationship between target and DAF correlations (see ED fig 5)

for l = 1:Nsim
    tic;

    simLocs = unitLocs+randn(size(unitLocs))*jitterDist; %add jitter to datapoints

    %set origin as target neuron
    unitOrigin = simLocs(simTargetI, :);
    %compute euclidian distances from target neuron
    unitDist = sqrt(sum((simLocs-repmat(unitOrigin, size(simLocs, 1), 1)).^2, 2));

    %sample target DAF correlation from fit normal distribution
    phiTarget = normrnd(mu,sig);
    %sample recorded neuron DAF correlations from parameterized model
    phi(recordI) = phiTarget*(a*exp(b*unitDist(recordI)) + randn(sum(recordI), 1).*sqrt(c*exp(d*unitDist(recordI))+e));
    phi(simTargetI) = phiTarget;

    %sample unrecorded neuron correlations from estimated probability distribution
    [~, I] = min(abs(repmat(xBin, sum(~recordI), 1)-repmat(unitDist(~recordI), 1, numel(xBin))), [], 2);
    [~, J] = min(abs(pairCorrCDF(:, I)-repmat(rand(1, sum(~recordI)), numel(yBin), 1)));
    phi(~recordI) = yBin(J)*phiTarget*pairwiseCorrF;

    %add values to results
    unitDistSim((l-1)*Nsamp+(1:Nsamp)) = unitDist(recordI);
    phiSim((l-1)*Nsamp+(1:Nsamp)) = phi(recordI);

    %simulate learning rates for each neuron across projection divergence values
    for k = 1:numel(sigma)
        learnR = simulate_learning_GPU_fast(1, sigma(k), phi, simLocs);
        learnR = learnR * phi(simTargetI) * fEff.p1 / learnR(simTargetI);
        learnRSim((l-1)*Nsamp+(1:Nsamp),k) = learnR(recordI);
    end

    tocT = toc;
    disp("iter " + l + " of " + Nsim);
    disp("remaining time = " + (Nsim-l)*tocT/3600+ " hr");
end

%add observed uncertanties in correlation and learning rate
learnRCI = mean(confRTot(:, 2)-confRTot(:, 1))/4;
phiCI = mean(phiConfTot(:, 2)-phiConfTot(:, 1))/4;
phiSim = phiSim + randn(size(phiSim))*phiCI;
learnRSim = learnRSim + randn(size(learnRSim))*learnRCI;


%% plot simulation results

sigmaPlot = sigma(5:2:end); %values to plot

distC = 400; %max distance to plot
plotN = 50; %plot resolution

colorVals = colorcet('R4', 'N', numel(sigmaPlot)); %color palette for sigma values

figure;
%plot learning rate vs correlation
ax1 = subplot(2, 1, 1);
axis square;
xlim(ax1, [-0.1, 0.5]);
xlabel('correlation with noise');
ylabel('learning rate');
xline(0, '--k');
yline(0, '--k');
%plot learning rate vs distance
ax2 = subplot(2, 1, 2);
axis square;
xlim(ax2, [0, distC]);
xlabel('distance from contingent neuron');
ylabel('learning rate');
yline(0, '--k');
hold([ax1 ax2], 'on');
set([ax1 ax2], 'TickDir', 'out');

alphaVal = 0.1;
yl = [0 0];

for kPlot = 1:numel(sigmaPlot)
    k = find(sigma == sigmaPlot(kPlot)); %index for this plotted value of sigma
    
    learnR_phi = zeros(plotN, 1); %learning per correlation values
    learnR_r = zeros(plotN, 1); %learning per distance values
    rVals = linspace(0, max(unitDistSim), plotN); %vector of distance values
    rPitch = mean(diff(rVals)); %distance spacing
    phiVals = linspace(min(phiSim), max(phiSim), plotN); %vector of correlation values
    phiPitch = mean(diff(phiVals)); %correlation spacing
    %calculate mean learning rates across correlation and distance bins
    for i = 1:plotN
        windowI = phiSim>phiVals(i)-phiPitch/2 & phiSim<phiVals(i)+phiPitch/2;
        learnR_phi(i) = mean(learnRSim(windowI, k));
        windowI = unitDistSim>rVals(i)-rPitch/2 & unitDistSim<rVals(i)+rPitch/2;
        learnR_r(i) = mean(learnRSim(windowI, k));
    end
    
    %interpolate NaN values (undersampled regions)
    nVals = 1:plotN;
    nanI = isnan(learnR_phi(:));
    learnR_phi_smooth(nanI) = interp1(nVals(~nanI), learnR_phi(~nanI), nVals(nanI));
    nanI = isnan(learnR_r(:));
    learnR_r_smooth(nanI) = interp1(nVals(~nanI), learnR_r(~nanI), nVals(nanI));
    
    %downsample using Gaussian kernel
    binW = 1; %downsampling off with binW=1
    learnR_phi_smooth = smoothdata(learnR_phi, 1, 'gaussian', 5*binW/2.355);
    learnR_phi_smooth = downsample(learnR_phi_smooth, binW);
    learnR_r_smooth = smoothdata(learnR_r, 1, 'gaussian', 5*binW/2.355);
    learnR_r_smooth = downsample(learnR_r_smooth, binW);
    %downsample x axis
    phiVals_plot = downsample(phiVals, binW);
    rVals_plot = downsample(rVals, binW);

    %plot curves
    plot(ax1, phiVals_plot, learnR_phi_smooth(:), 'Color', colorVals(kPlot, :));
    plot(ax2, rVals_plot, learnR_r_smooth(:), 'Color', colorVals(kPlot, :));

end

%UNCOMMENT BELOW to plot data over sigma curves

%{
phiTot = xdata(:, 1); %correlation from xdata matrix
edgeVals = prctile(phiTot, 0:5:100);

binVals = edgeVals(1:end-1)+diff(edgeVals)/2;
edgeVals(end) = edgeVals(end)+1;
learnMeans = zeros(size(binVals));
learnStd = zeros(size(binVals));
corrMeans = zeros(size(binVals));
for i = 1:numel(binVals)
    learnVals = learnRTot(xdata(:, 1)>=edgeVals(i) & xdata(:, 1)<edgeVals(i+1));
    learnMeans(i) = mean(learnVals);
    learnStd(i) = std(learnVals)/sqrt(numel(learnVals));
    corrMeans(i) = mean(xdata(xdata(:, 1)>edgeVals(i) & xdata(:, 1)<edgeVals(i+1), 1));
end
errorbar(ax1, corrMeans, learnMeans, learnStd, 'o', 'LineStyle', 'none', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'CapSize', 0);

edgeVals = prctile(unitDistTot(unitDistTot>0), linspace(0, 100, 20));
edgeVals = [-1 edgeVals];
edgeVals(end) = edgeVals(end)+1;
binVals = edgeVals(1:end-1)+diff(edgeVals)/2;
binVals(1) = 0;

learnMeans = zeros(size(binVals));
learnStd = zeros(size(binVals));
distMeans = zeros(size(binVals));
for i = 1:numel(binVals)
    learnVals = learnRTot(unitDistTot>=edgeVals(i) & unitDistTot<edgeVals(i+1));
    learnMeans(i) = mean(learnVals);
    learnStd(i) = std(learnVals)/sqrt(numel(learnVals));
    distMeans(i) = mean(unitDistTot(unitDistTot>=edgeVals(i) & unitDistTot<edgeVals(i+1)));
end
errorbar(ax2, distMeans, learnMeans, learnStd, 'o', 'LineStyle', 'none', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'CapSize', 0);
%}



%set y limits
ylim([ax1 ax2], [-0.04 0.14]);


%% relative log likelihood across range of sigma

%confidence intervals for learning rates and correlations
learnRCI = mean(confRTot(:, 2)-confRTot(:, 1))/4;
phiCI = mean(phiConfTot(:, 2)-phiConfTot(:, 1))/4;

%evaluation points for 3D probability distribution
phiVals = -0.3:0.01:0.6; %correlation values
rVals = 0:10:500; %distance values
learnVals = -.05:.01:.15; %learning rate values

%3D probability distribution for range of sigma values
simHist = zeros(numel(phiVals), numel(rVals), numel(learnVals), numel(sigma));

%calculate probability distribution from 
for i = 1:size(learnRSim, 1)
    [~, phiI] = min(abs(phiVals-phiSim(i)));
    [~, rI] = min(abs(rVals-unitDistSim(i)));
    for k = 1:numel(sigma)
        [~, learnI] = min(abs(learnVals-learnRSim(i, k)));
        simHist(phiI, rI, learnI, k) = simHist(phiI, rI, learnI, k)+1;
    end
end

%additive smoothing using Jeffreys prior
alpha = 0.5;
simHist = simHist + alpha*repmat(sum(simHist, 3), 1, 1, size(simHist, 3), 1);
simProb = simHist./repmat(sum(simHist, 3), 1, 1, size(simHist, 3), 1);


Nstrap = 10000; %bootstrap replicates
likelihoodVals = zeros(Nstrap, numel(sigma)); %log-likelihoods

%calculate likelihood values over jittered data
for j = 1:Nstrap
    %jitter correlations and learning rates
    phiStrap = xdata(:, 1)+randn(size(xdata, 1), 1)*phiCI;
    learnStrap = learnRTot+randn(size(learnRTot))*learnRCI;
    %update likelihood distribution
    for k = 1:numel(sigma)
        for i = 1:size(xdata, 1)
            [~, phiI] = min(abs(phiVals-phiStrap(i)));
            [~, rI] = min(abs(rVals-unitDistTot(i)));
            [~, learnI] = min(abs(learnVals-learnStrap(i)));
            if ~isnan(simProb(phiI, rI, learnI, k))
                likelihoodVals(j, k) = likelihoodVals(j, k) + log(simProb(phiI, rI, learnI, k));
            end
        end
    end
end

likelihoodVals = likelihoodVals-max(likelihoodVals(:)); %normalize for relative LL

meanVals = mean(likelihoodVals, 1); %mean values
sigmaInterp = 2.^(linspace(0, 8, 10000)); %interpolated curve
lInterp = interp1(sigma, meanVals, sigmaInterp, 'pchip');

%confidence intervals on LL
confMin = 1;
confMax = max(find(lInterp>max(lInterp)-icdf('chi2', .95, 1))); %find cutoff based on chi-squared

%plot results
figure;
errorbar(sigma, mean(likelihoodVals, 1), mean(likelihoodVals, 1)-prctile(likelihoodVals, 2.5, 1), prctile(likelihoodVals, 97.5, 1)-mean(likelihoodVals, 1), 'o', 'LineStyle', 'none', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'CapSize', 0);
hold on;
plot(sigmaInterp, lInterp);
yline(max(lInterp)-icdf('chi2', .95, 1), '--r'); %chi-squared cutoff
xline(sigmaInterp([confMin confMax]));

xlim([0 max(sigma)]);
set(gca, 'TickDir', 'out');
xlabel('projection divergence (microns)');
ylabel('relative log likelihood');
xscale('log');

%output values in text
disp("95% confidence interval = "+sigmaInterp(confMin)+" to "+sigmaInterp(confMax));

%plot on linear x axis
linMax = 30;
linL = 10.^likelihoodVals;
figure;
scatter(sigma, 10.^mean(likelihoodVals, 1), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b');
hold on;
plot(sigmaInterp, 10.^lInterp);
xlim([1 linMax]);
set(gca, 'TickDir', 'out');
xlabel('projection divergence (microns)');
ylabel('relative log likelihood');
