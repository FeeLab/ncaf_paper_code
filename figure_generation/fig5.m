

%% import functions

addpath('../utils')
addpath('../external/SpikeGLX_Datafile_Tools/MATLAB');
addpath('../external/npy-matlab/npy-matlab');
addpath('../external/ndSparse_G4_2021_03_16');
addpath('../external/colorcet');

%% plot example up

%load data from eligibility trace experiment
load('../presorted_data/2024-05-03_10872_LMAN_nCAF_aligned.mat');
unitSignal = full(unitSigSparse);

unitI = 58; %high SNR unit that learned

%plot raster
plot_learning_figure(unitSignal, trigDelay, unitI, unitNum, silenceTemplate, cW, np_fs, noiseI, median(noiseF)-median(noiseI));

downR = 30; %downsample to 1 ms time bins
figure;
ax1 = subplot(2, 1, 1);
%plot correlation with DAF
[tvals, frateCorr, confC] = contingency_correlation_uncertainty_eligibility(unitSignal, np_fs, trigDelay, cW, unitI, unitNum, isNoise, ax1, true, .04, downR);
ax2 = subplot(2, 1, 2);
%plot learning rates
[tvals, learnRcurve, confR] = response_width_fitUncertainty_eligibility(unitSignal, np_fs, trigDelay, cW, unitI, unitNum, ax2, true, .04, downR);


%% combine across all trials

windowW = 0.3;

%all units with phi>0.2 and learnR>2e-3 across contingent window (see ED Figure 6)
unitVec = [163,165,167,168,170,171,172,...
    54,55,56,58,60,64,65,...
    37,39,42,44,...
    147,151,152,153,156,...
    32,36,38,43,...
    33,38,39,41,...
    141,142,144,...
    153,155,156,159,160,161,164,166,167,...
    151,153,154,155,157,158];
runVec = [repmat("2024-05-02_10872_LMAN_nCAF_aligned.mat", 7, 1);
    repmat("2024-05-03_10872_LMAN_nCAF_aligned.mat", 7, 1);
    repmat("2025-11-15_Yellow33_Post-Advance_LMAN_Eligibility_1_aligned.mat", 4, 1);
    repmat("2025-11-16_Yellow33_Post-Advance_LMAN_Eligibility_0_aligned.mat", 5, 1);
    repmat("2025-11-20_Yellow33_Post-Advance_LMAN_Eligibility_0_aligned.mat", 4, 1);
    repmat("2025-11-21_Yellow33_Post-Advance_LMAN_Eligibility_0_aligned.mat", 4, 1);
    repmat("2025-11-22_Yellow33_Post-Advance_LMAN_Eligibility_0_aligned.mat", 3, 1);
    repmat("2025-12-10_11384_Post-Advance_LMAN_Eligibility_0_aligned.mat", 9, 1);
    repmat("2025-12-11_11384_Post-Advance_LMAN_Eligibility_0_aligned.mat", 6, 1)];
birdID = [1*ones(14, 1); 2*ones(20, 1); 3*ones(15, 1)]; %unique ID for each bird

%load a data file
load(fullfile('../presorted_data', runVec(end)));

downR = 60; %downsample to 2 ms time bins

learnRTot = zeros(numel(unitVec), round(windowW*np_fs/downR+1)); %learning rates
corrTot = zeros(numel(unitVec), round(windowW*np_fs/downR+1)); %correlations with DAF
learnRTotConf = zeros(numel(unitVec), round(windowW*np_fs/downR+1), 2); %CI for learning rats
corrTotConf = zeros(numel(unitVec), round(windowW*np_fs/downR+1), 2); %CI for correlations
noiseStartTot = zeros(size(unitVec)); %times of noise onset
noiseEndTot = zeros(size(unitVec)); %times of noise offset
cwI = zeros(size(unitVec)); %times of contingent window rising edge

nT = 1000; %noise distribution number of timepoints
noiseT = linspace(-50, 50, nT); %noise distribution timepoints
pNoiseTot = zeros(numel(unitVec), numel(noiseT)); %noise distribution

buffF = 0.05; %time after noise burst rising edge to include in analysis (seconds)

for k = 1:numel(unitVec)
    disp(k);

    % load dataset if not already
    if ~contains(runVec(k), run_name)
        load(fullfile('../presorted_data', runVec(k)));
        unitSignal = full(unitSigSparse);
    end

    %use noise burst onset as time zero
    tStart = median(noiseI)-trigDelay; %noise burst onset relative to end of contingent window
    noiseStartTot(k) = tStart; %noise burst onset
    noiseEndTot(k) = median(noiseF)-trigDelay; %noise burst offset relative to end of contingent window
    cwI(k) = -(cW+tStart); %start of contingent window relative to noise burst onset

    %calculate correlation with DAF
    [tvals, frateCorr, confC] = contingency_correlation_uncertainty_eligibility(unitSignal, np_fs, trigDelay+tStart+buffF, windowW, unitVec(k), unitNum, isNoise, [], true, 0, downR);
    %calculate learning rates
    [tvals, learnRcurve, confR] = response_width_fitUncertainty_eligibility(unitSignal, np_fs, trigDelay+tStart+buffF, windowW, unitVec(k), unitNum, [], true, 0, downR);
    tvals = tvals+1000*buffF; %offset timepoints to include times after noise burst onset (t>0)
    learnRTot(k, :) = learnRcurve;
    corrTot(k, :) = frateCorr;
    learnRTotConf(k, :, :) = confR;
    corrTotConf(k, :, :) = confC;

    %update noise burst histogram
    pNoise = zeros(size(noiseT));
    for i = 1:numel(noiseI)
        thisI = noiseT/1000>=noiseI(i)-trigDelay-tStart & noiseT/1000<=noiseF(i)-trigDelay-tStart;
        pNoise(thisI) = pNoise(thisI)+1;
    end
    %normalize to get noise burst distribution
    pNoiseTot(k, :) = pNoise/numel(noiseI);
end

%% elementwise thiel-sen


Nstrap = 1000; %bootstrap replicates

regR = zeros(size(tvals)); %normalized learning rates
regRconf = zeros(size(tvals, 1), 2); %confidence intervals in learning rates
slopeTot = zeros(size(corrTot)); %slope estimates for each unit
birdNums = unique(birdID); %IDs per bird

sigmaEta = mean(var(corrTot(:, tvals<-150), 0, 1)); %estimate of noise in correlation measurement

lambda = zeros(size(regR)); %reliability ratio values

for i = 1:numel(regR)

    slopeTot(:, i) = learnRTot(:, i)./corrTot(:, i); %calculate slopes across units

    simLambda = zeros(Nstrap, 1); %bootstrapped reliability ratio
    simR = zeros(Nstrap, 1); %bootstrapped learning rates
    for j = 1:Nstrap
        birdStrap = randsample(max(birdNums), max(birdNums), true); %resample birds
        x = []; %total set of resampled points, xvals
        y = []; %total set of resampled points, yvals
        for k = 1:numel(birdStrap)
            birdUnits = birdID==birdStrap(k); %select resampled bird
            birdCorr = corrTot(birdUnits, :); %correlations from this bird
            birdLearn = learnRTot(birdUnits, :); %learning rates from this bird

            %resample units from this bird and append
            N = sum(birdUnits);
            sampI = randsample(N, N, true);
            corrStrap = birdCorr(sampI, i);
            learnStrap = birdLearn(sampI, i);
            x = [x; corrStrap(:)];
            y = [y; learnStrap(:)];
        end

        R = median(y./x, 'omitnan'); %median of slopes for Thiel-Sen
        simR(j) = R;

        simLambda(j) = (mean(x.^2)-sigmaEta)/mean(x.^2); %calculate reliability ratio for this resampling
    end

    %update values
    lambdaVals(i) = median(simLambda);
    regRconf(i, :) = prctile(simR, [2.5, 97.5]); 
    regR(i) = mean(simR);
end

%set timepoints with low reliatiblity ratio to have large uncertainty
ff = 0.5; %cutoff factor for reliability ratio
magVal = 10; %scaling factor of large fixed uncertainty
railVal = max(abs(regR)); %uncertainty value relative to max learning rate
regRconf(lambdaVals<ff, 1) = -magVal*railVal;
regRconf(lambdaVals<ff, 2) = magVal*railVal;



%% plot learning per correlation

noiseStart = min(noiseStartTot)*1000; %start of noise burst
noiseL = median(noiseEndTot-noiseStartTot)*1000; %median burst width

%plot normalized learning rates
figure;
hold on;
plot(tvals, regR);
%plot confidence intervals
xconf = [tvals tvals(end:-1:1)];
yconf = [regRconf(:, 1)' regRconf(end:-1:1, 2)'];
p=fill(xconf, yconf, 'blue');
p.FaceAlpha = .1;
p.EdgeColor = 'none';
xline(-noiseStart, '--k');
xline(min(cwI)*1000, '--k');
yline(0);
xlabel('time relative to noise burst onset (ms)');
ylabel('learning rate per unit correlation (Hz/trial)');
set(gca, 'TickDir', 'out');
xlim([-200 50]);
ylim([-0.3 .7]);
leftLim = get(gca, 'ylim');

%plot DAF distribution
yyaxis right;
plot(noiseT, mean(pNoiseTot, 1));
ylim([leftLim(1)/leftLim(2) 1]);
ylabel('noise probability');
%plot FWHM of DAF distribution
[~, distF] = min(abs(mean(pNoiseTot, 1)-0.5));
xline([0 noiseT(distF)], '-r');

%% plot deconvolution

%downsample DAF distribution to same timescale as learning rate
noiseVec = decimate(mean(pNoiseTot, 1), round(10*downR/np_fs*1000));
noiseVec = noiseVec/sum(noiseVec);

%cutoff value for frequency content
sigma = 3.3/sqrt(2)/(downR/np_fs*1000); %set by measured width of learning kernel

%approximate autocorrelation as gaussian
gaussX = -20:1:20;
gaussVec = exp(-gaussX.^2/(2*sigma^2))/(sqrt(2*pi)*sigma);

%noise autocorrelation (assume white)
nsize = 20; %size of autocorrelation vector
ncorr = zeros(1, 2*nsize+1);
ncorr(nsize+1) = 1e0;

%sweep over values of lambda using Morozov's discrepancy principle
lambda = 0:0.1:10;
rmdVec = zeros(size(lambda));
for i = 1:numel(lambda)
    %deconvolve
    J = deconvwnr(regR, noiseVec, ncorr, lambda(i)*gaussVec);
    rconv = conv(J, noiseVec, 'same');
    rmdVec(i) = sum((regR-rconv).^2);
end
%estimate cutoff value using learning rates far from contingent window
rmdCut = sum((regR(tvals<-150)).^2)*numel(tvals)/sum(tvals<-150);

%plot lambda cutoff
figure;
plot(lambda, rmdVec);
yline(rmdCut);

%set lambda and recalculate deconvolution
lambdaCut = lambda(find(rmdVec<rmdCut, 1));
J = deconvwnr(regR, noiseVec, ncorr, lambdaCut*gaussVec);
Jout = J;

%jitter learning rate curve to obtain bootstrapped eligibility trace distribution
Nstrap = 10000; %bootstrap replicates
Jstrap = zeros(numel(J), Nstrap);
for i = 1:Nstrap
    %jitter learning rates
    rStrap = regR+randn(size(regR)).*diff(regRconf, 1, 2)'/4;
    %deconvolve
    Jstrap(:, i) = deconvwnr(rStrap, noiseVec, ncorr, lambdaCut*gaussVec);
end
Jconf = prctile(Jstrap, [2.5, 97.5], 2); %confidence intervals

%plot deconvolved eligibility trace
figure;
plot(-tvals, Jout);
hold on;
xconf = [-tvals -tvals(end:-1:1)];
yconf = [Jconf(:, 1)' Jconf(end:-1:1, 2)'];
p=fill(xconf, yconf, 'blue');
p.FaceAlpha = .1;
p.EdgeColor = 'none';
yline(0);
xline(0, '--k');
xlim([-50 200]);
ylim([-0.2, 0.3]);
xlabel('time relative to LMAN activity (ms)');
ylabel('learning rate per unit correlation (Hz/trial)');
set(gca, 'TickDir', 'out');

%calculate distribution of offset times
cutVals = zeros(size(Jstrap, 2), 1);
for i = 1:numel(cutVals)
    cutVals(i) = tvals(find(Jstrap(:, i)'<0&tvals<-50, 1, 'last'));
end

%output values in text
disp("cutoff = "+mean(-cutVals)+", 95% CI = "+prctile(-cutVals, 2.5)+" to "+prctile(-cutVals, 97.5));

%% crosscorrelate eligibility trace with 100 ms noise burst to simulate pCAF 

%upsample eligiblity trace to 1 ms time bin
tUp = min(round(tvals)):max(round(tvals));
Jup = interp1(tvals, Jout, tUp);
JstrapUp = zeros(numel(Jup), size(Jstrap, 2));
for i = 1:size(JstrapUp, 2)
    JstrapUp(:, i) = interp1(tvals, Jstrap(:, i), tUp);
end

%use prior of zero for eligibility trace values at long time delays
tCut = tUp(find(Jup<0&tUp<-100, 1, 'last')); %times after tcut are set to prior
priorW = 50; %time bin to calculate uncertainty in prior value
%approximate std of eligibility trace values within times of high confidence
priorStd = mean(std(JstrapUp(tUp<(tCut+priorW) & tUp>tCut, :), 0, 2), 1);

%set values beyond tcut to a constant small value centered around zero
priorT = tUp(tUp<tCut);
for i = 1:size(JstrapUp, 2)
    JstrapUp(ismember(tUp, priorT), i) = randn*priorStd;
end
Jup(tUp<tCut) = 0;


nW = 100; %width of noise burst
premotor_d = 22; % calculated from premotor latency estimate, rounded to nearest ms

%zero-centered vector of eligiblity trace and confidence intervals
eTrace = zeros(2*sum(round(tUp)<0)+1, 1); 
eTvals = min(round(tUp)):-min(round(tUp)); %timepoints to sample eligiblity trace
eTrace(end-numel(Jup)+1:end) = fliplr(Jup);
eTraceStrap = zeros(numel(eTrace), size(JstrapUp, 2));
eTraceStrap(end-numel(Jup)+1:end, :) = flipud(JstrapUp);

%define noise burst vector
noiseVec_pcaf = zeros(size(eTvals));
noiseVec_pcaf(eTvals>premotor_d & eTvals<=premotor_d+nW) = 1/nW;

%crosscorrelate to get predicted learning rates
pcaf_learning = xcorr(eTrace, noiseVec_pcaf);

%boostrap crosscorrelations with distribution of eligibility traces
pcafStrap = zeros(numel(pcaf_learning), size(eTraceStrap,2));
for i = 1:size(eTraceStrap, 2)
    pcafStrap(:, i) = xcorr(eTraceStrap(:, i), noiseVec_pcaf);
end

plotT = 0:300; %plot limits
plotLearning = pcaf_learning(ceil(numel(pcaf_learning)/2)+plotT);
plotStrap = pcafStrap(ceil(numel(pcaf_learning)/2)+plotT, :);
pcafConf = prctile(plotStrap, [2.5, 97.5], 2);

%plot predicted pcaf learning rates
figure;
plot(plotT, plotLearning);
hold on;
xconf = [plotT plotT(end:-1:1)];
yconf = [pcafConf(:, 1)' pcafConf(end:-1:1, 2)'];
p=fill(xconf, yconf, 'blue');
p.FaceAlpha = .1;
p.EdgeColor = 'none';
xlim([0 200]);
yline(0);
xlabel('time delay (ms)');
ylabel('learning rate (relative)');


