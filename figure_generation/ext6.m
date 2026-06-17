%% Set random seed
rng(0);

%% import functions
addpath('../utils')
addpath('../external/SpikeGLX_Datafile_Tools/MATLAB');
addpath('../external/npy-matlab/npy-matlab');
addpath('../external/ndSparse_G4_2021_03_16');
addpath('../external/colorcet');

%% plot all units

%all eligibility trace experiments
runVec = ["2024-05-02_10872_LMAN_nCAF_aligned.mat";
    "2024-05-03_10872_LMAN_nCAF_aligned.mat";
    "2025-11-15_Yellow33_Post-Advance_LMAN_Eligibility_1_aligned.mat";
    "2025-11-16_Yellow33_Post-Advance_LMAN_Eligibility_0_aligned.mat";
    "2025-11-20_Yellow33_Post-Advance_LMAN_Eligibility_0_aligned.mat";
    "2025-11-21_Yellow33_Post-Advance_LMAN_Eligibility_0_aligned.mat";
    "2025-11-22_Yellow33_Post-Advance_LMAN_Eligibility_0_aligned.mat";
    "2025-12-10_11384_Post-Advance_LMAN_Eligibility_0_aligned.mat";
    "2025-12-11_11384_Post-Advance_LMAN_Eligibility_0_aligned.mat"];


for k = 1:numel(runVec)
    load(fullfile('./paper_figure_generation/presorted_data', runVec(k)));

    unitSignal = full(unitSigSparse);
        
    isNormed = true;
    [learnR, confR] = learning_per_neuron(unitSignal(:,:,lmanUnits), np_fs, trigDelay, cW, isNormed);
    
    spikeCount = contingent_spikes(unitSignal(:,:,lmanUnits), np_fs, trigDelay, cW, isNormed);
    [phi, phiConf] = activity_noise_corr(spikeCount, isNoise);
    
    ww = 0.05; %plotting width outside contingent window
    
    downR = 30; %downsample to 1 ms
    %retain units that are sufficiently correlated with DAF and learn sufficiently
    learnI = lmanNum(phi>0.2&learnR>2e-3);
    relR = zeros(numel(learnI), round(cW*np_fs/downR)+1);
    learnRTot = zeros(numel(learnI), round((cW+2*ww)*np_fs/downR+1));
    corrTot = zeros(numel(learnI), round((cW+2*ww)*np_fs/downR+1));
    learnRTotConf = zeros(numel(learnI), round((cW+2*ww)*np_fs/downR+1), 2);
    corrTotConf = zeros(numel(learnI), round((cW+2*ww)*np_fs/downR+1), 2);
    
    %calculate correlations and learning rates
    for k = 1:numel(learnI)
        [tvals, frateCorr, confC] = contingency_correlation_uncertainty_eligibility(unitSignal, np_fs, trigDelay, cW, learnI(k), unitNum, isNoise, [], true, ww, downR);
        [tvals, learnRcurve, confR] = response_width_fitUncertainty_eligibility(unitSignal, np_fs, trigDelay, cW, learnI(k), unitNum, [], true, ww, downR);
        contI = tvals>=-cW*1000 & tvals<=0;
        learnRTot(k, :) = learnRcurve;
        corrTot(k, :) = frateCorr;
        learnRTotConf(k, :, :) = confR;
        corrTotConf(k, :, :) = confC;
    end
    
    %DAF distribution
    pNoise = zeros(size(tvals));
    for i = 1:numel(noiseI)
        thisI = tvals/1000>=noiseI(i)-trigDelay & tvals/1000<=noiseF(i)-trigDelay;
        pNoise(thisI) = pNoise(thisI)+1;
    end
    pNoise = pNoise/numel(noiseI);
    Ninterp = 100000;
    tInterp = linspace(min(tvals), max(tvals), Ninterp);
    pInterp = interp1(tvals, pNoise, tInterp);
    noiseStart = tInterp(find(pInterp>max(pNoise)/2, 1));
    noiseEnd = tInterp(find(pInterp>max(pNoise)/2, 1, 'last'));
    
    %plot results
    figure;
    ax1 = subplot(3, 1, 1);
    plot(tvals-noiseStart, pNoise);
    set(gca, 'TickDir', 'out');
    xlabel('time (ms)');
    ylabel('probability of noise');
    ax2 = subplot(3, 1, 2);
    imagesc(tvals-noiseStart, 1:numel(learnI), corrTot);
    yticks(1:numel(learnI));
    yticklabels(string(learnI));
    clim([0 Inf]);
    colormap(colorcet('L20'));
    xlabel('time (ms)');
    set(gca, 'TickDir', 'out');
    title('correlation with noise');
    cb = colorbar;
    cb.Location = 'manual';
    cb.Position = [0.93 0.425 0.02 0.2];
    cb.Label.String = 'correlation';
    ax3 = subplot(3, 1, 3);
    imagesc(tvals-noiseStart, 1:numel(learnI), learnRTot);
    yticks(1:numel(learnI));
    yticklabels(string(learnI));
    clim([0 Inf]);
    colormap(colorcet('L20'));
    xlabel('time (ms)');
    set(gca, 'TickDir', 'out');
    title('learning rate');
    cb = colorbar;
    cb.Location = 'manual';
    cb.Position = [0.93 0.125 0.02 0.2];
    cb.Label.String = 'learning rate (Hz/trial)';
    linkaxes([ax1 ax2 ax3], 'x');
    xlim([min(tvals) max(tvals)]-noiseStart);
    xline(ax2, [0 -cW*1000]-noiseStart, '-r');
    xline(ax3, [0 -cW*1000]-noiseStart, '-r');
end
