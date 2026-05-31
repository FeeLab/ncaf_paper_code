

%% initial parameters

addpath('../utils')
addpath('../external/SpikeGLX_Datafile_Tools/MATLAB');
addpath('../external/npy-matlab/npy-matlab');
addpath('../external/ndSparse_G4_2021_03_16');
addpath('../external/colorcet');

%% plot example up

%load dataset
load('../presorted_data/2025-02-21_11208_Post-Advance_LMAN_nCAF_0_aligned_warped.mat');
unitSignal = full(unitSigSparse);

unitI = 127; %targeted unit from this run

plot_learning_figure(unitSignal, trigDelay, unitI, unitNum, aTemplate, cW, np_fs, noiseI, noiseW);

%% plot up trendlines

plot_learning_trendlines_allTimes(unitSignal, np_fs, trigDelay, unitI, unitNum);

%% plot up learning width
ww = 0.02; %width of window to plot
downsampW = 15; %downsample factor (0.5 ms bin width for 30 kHz sample rate)

%plot learning width
figure;
plotAx = gca;
[tvals, learnR, confR] = response_width_fitUncertainty(unitSignal, np_fs, trigDelay, cW, unitI, unitNum, plotAx, true, ww, downsampW);
yline(0);

%% plot example down

load('../presorted_data/2025-06-11_11301_Post-Advance_LMAN_nCAF_0_aligned_warped.mat');
unitSignal = full(unitSigSparse);

unitI = 187; %targeted unit from this run

plot_learning_figure(unitSignal, trigDelay, unitI, unitNum, silenceTemplate, cW, np_fs, noiseI, noiseW);

%% plot down trendlines

plot_learning_trendlines_allTimes(unitSignal, np_fs, trigDelay, unitI, unitNum);

%% plot down learning width
ww = 0.02; %width of window to plot
downsampW = 15; %downsample factor (0.5 ms bin width for 30 kHz sample rate)

%plot learning width
figure;
plotAx = gca;
[tvals, learnR, confR] = response_width_fitUncertainty(unitSignal, np_fs, trigDelay, cW, unitI, unitNum, plotAx, true, ww, downsampW);
yline(0);

%% LMAN summary

%list of targeted units and their respective run files
unitVec = [116,148,187,80,868,127,62,73,168,161,187,168];
runVec = ["2024-01-16_10687_nCAF_LMAN-X_1_aligned.mat";
    "2024-04-23_10847_LMAN_nCAF_aligned.mat";
    "2024-04-27_10872_LMAN-X_nCAF_2_aligned.mat";
    "2024-07-24_10977_LMAN_nCAF_aligned.mat";
    "2025-02-18_11208_Post-Advance_LMAN_BOTM_0_aligned.mat";
    "2025-02-21_11208_Post-Advance_LMAN_nCAF_0_aligned.mat";
    "2025-03-23_11238_LMAN_nCAF_0_aligned.mat";
    "2025-11-12_Yellow33_Post-Advance_LMAN_BOTM_0_aligned.mat";
    "2024-02-28_10762_LMAN_nCAF_0_aligned.mat";
    "2024-02-29_10762_LMAN_nCAF_1_aligned.mat";
    "2025-06-11_11301_Post-Advance_LMAN_nCAF_0_aligned.mat";
    "2025-11-13_Yellow33_Post-Advance_LMAN_BOTM_0_aligned.mat"];
direction = logical([1,1,1,1,1,1,1,1,0,0,0,0]); %1=up, 0=down
singleUnit = [0,0,1,0,1,1,0,1,0,0,0,0]; %1=single unit, 0=multiunit

learnRTot = zeros(numel(unitVec), 1); %learning rates
confRTot = zeros(numel(unitVec), 2); %CI on learning rates
learnROutside = zeros(numel(unitVec), 1); %learning rates outside contingent window
confROutside = zeros(numel(unitVec), 2); %CI on learning rates outside

%compute learning rates for each targeted unit
for k = 1:numel(unitVec)
    unitI = unitVec(k);
    load(fullfile('../presorted_data', runVec(k)));
    unitSignal = full(unitSigSparse);

    %compute learning inside contingent window
    [learnRTot(k), confRTot(k, :)] = learning_per_neuron(unitSignal(:,:,unitNum==unitI), np_fs, trigDelay, cW, false);
    %compute learning outside contingent window
    [learnROutside(k), confROutside(k, :)] = learning_per_neuron_outsideTimes(unitSignal, np_fs, trigDelay, unitI, unitNum, cW);
end

%plot results
figure;
hold on;
xVals = 1:numel(learnRTot);
groupOffset = numel(learnRTot)+1;

%plot single unit up
thisPlot = direction&singleUnit;
offset = 0;
errorbar((1:sum(thisPlot))+offset, learnRTot(thisPlot), learnRTot(thisPlot)-confRTot(thisPlot, 1), confRTot(thisPlot, 2)-learnRTot(thisPlot), '^', 'LineStyle', 'none', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'CapSize', 0);
errorbar((1:sum(thisPlot))+offset+groupOffset, learnROutside(thisPlot), learnROutside(thisPlot)-confROutside(thisPlot, 1), confROutside(thisPlot, 2)-learnROutside(thisPlot), '^', 'LineStyle', 'none', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'CapSize', 0);
offset = offset+sum(thisPlot);

%plot multiunit up
thisPlot = direction&~singleUnit;
errorbar((1:sum(thisPlot))+offset, learnRTot(thisPlot), learnRTot(thisPlot)-confRTot(thisPlot, 1), confRTot(thisPlot, 2)-learnRTot(thisPlot), '^', 'LineStyle', 'none', 'MarkerEdgeColor', 'b', 'CapSize', 0);
errorbar((1:sum(thisPlot))+offset+groupOffset, learnROutside(thisPlot), learnROutside(thisPlot)-confROutside(thisPlot, 1), confROutside(thisPlot, 2)-learnROutside(thisPlot), '^', 'LineStyle', 'none', 'MarkerEdgeColor', 'b', 'CapSize', 0);
offset = offset+sum(thisPlot);

%plot single unit down
thisPlot = ~direction&singleUnit;
errorbar((1:sum(thisPlot))+offset, learnRTot(thisPlot), learnRTot(thisPlot)-confRTot(thisPlot, 1), confRTot(thisPlot, 2)-learnRTot(thisPlot), 'v', 'LineStyle', 'none', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'CapSize', 0);
errorbar((1:sum(thisPlot))+offset+groupOffset, learnROutside(thisPlot), learnROutside(thisPlot)-confROutside(thisPlot, 1), confROutside(thisPlot, 2)-learnROutside(thisPlot), 'v', 'LineStyle', 'none', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'CapSize', 0);
offset = offset+sum(thisPlot);

%plot multiunit down
thisPlot = ~direction&~singleUnit;
errorbar((1:sum(thisPlot))+offset, learnRTot(thisPlot), learnRTot(thisPlot)-confRTot(thisPlot, 1), confRTot(thisPlot, 2)-learnRTot(thisPlot), 'v', 'LineStyle', 'none', 'MarkerEdgeColor', 'b', 'CapSize', 0);
errorbar((1:sum(thisPlot))+offset+groupOffset, learnROutside(thisPlot), learnROutside(thisPlot)-confROutside(thisPlot, 1), confROutside(thisPlot, 2)-learnROutside(thisPlot), 'v', 'LineStyle', 'none', 'MarkerEdgeColor', 'b', 'CapSize', 0);
offset = offset+sum(thisPlot);

%bars for mean values
bar(numel(learnRTot)/2, mean(learnRTot(direction)));
bar(numel(learnRTot)/2, mean(learnRTot(~direction)));
bar(numel(learnRTot)/2+groupOffset, mean(learnROutside));

%output values in text
disp("mean learning rate up = " + mean(learnRTot(direction)) + " Hz/trial");
disp("std = " + std(learnRTot(direction))/sqrt(sum(direction)));
[h, p, ci, stats] = ttest(learnRTot(direction), 0, 'Tail', 'right');
disp("t-value = "+stats.tstat+", p-value = "+p);

disp("mean learning rate down = " + mean(learnRTot(~direction)) + " Hz/trial");
disp("std = " + std(learnRTot(~direction))/sqrt(sum(direction)));
[h, p, ci, stats] = ttest(learnRTot(~direction), 0, 'Tail', 'left');
disp("t-value = "+stats.tstat+", p-value = "+p);

disp("outside learning rate up = " + mean(learnROutside(direction)) + " Hz/trial");
disp("std = " + std(learnROutside(direction))/sqrt(sum(direction)));
[h, p, ci, stats] = ttest(learnROutside(direction));
disp("t-value = "+stats.tstat+", p-value = "+p);

disp("outside learning rate down = " + mean(learnROutside(~direction)) + " Hz/trial");
disp("std = " + std(learnROutside(~direction))/sqrt(sum(direction)));
[h, p, ci, stats] = ttest(learnROutside(~direction));
disp("t-value = "+stats.tstat+", p-value = "+p);
%% nido summary

%list of targeted nido units and their respective run files
unitVec = [110,256,54,11,221];
runVec = ["2025-02-22_11208_Post-Advance_Nido_nCAF_0_aligned.mat";
    "2025-03-22_11238_Nido_nCAF_1_aligned.mat";
    "2025-03-24_11238_Nido_nCAF_0_aligned.mat";
    "2025-06-07_11301_Post-Advance_Nido_nCAF_0_aligned.mat";
    "2025-03-25_11238_Nido_nCAF_0_aligned.mat"];
direction = logical([1,1,1,0,1]); %1=up, 0=down
singleUnit = [1,1,0,1,1]; %1=single unit, 2=multiunit

learnRNido = zeros(numel(unitVec), 1);
confRNido = zeros(numel(unitVec), 2);

%compute learning for each target unit
for k = 1:numel(unitVec)
    unitI = unitVec(k);
    load(fullfile('../presorted_data', runVec(k)));
    unitSignal = full(unitSigSparse);

    [learnRNido(k), confRNido(k, :)] = learning_per_neuron(unitSignal(:,:,unitNum==unitI), np_fs, trigDelay, cW, false);
end


%plot results
%single unit up
thisPlot = direction&singleUnit;
offset = offset+groupOffset + 1;
errorbar((1:sum(thisPlot))+offset, learnRNido(thisPlot), learnRNido(thisPlot)-confRNido(thisPlot, 1), confRNido(thisPlot, 2)-learnRNido(thisPlot), '^', 'LineStyle', 'none', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'CapSize', 0);
offset = offset+sum(thisPlot);

%multiunit up
thisPlot = direction&~singleUnit;
errorbar((1:sum(thisPlot))+offset, learnRNido(thisPlot), learnRNido(thisPlot)-confRNido(thisPlot, 1), confRNido(thisPlot, 2)-learnRNido(thisPlot), '^', 'LineStyle', 'none', 'MarkerEdgeColor', 'b', 'CapSize', 0);
offset = offset+sum(thisPlot);

%single unit down
thisPlot = ~direction&singleUnit;
errorbar((1:sum(thisPlot))+offset, learnRNido(thisPlot), learnRNido(thisPlot)-confRNido(thisPlot, 1), confRNido(thisPlot, 2)-learnRNido(thisPlot), 'v', 'LineStyle', 'none', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'CapSize', 0);
offset = offset+sum(thisPlot);

%multiunit down
thisPlot = ~direction&~singleUnit;
errorbar((1:sum(thisPlot))+offset, learnRNido(thisPlot), learnRNido(thisPlot)-confRNido(thisPlot, 1), confRNido(thisPlot, 2)-learnRNido(thisPlot), 'v', 'LineStyle', 'none', 'MarkerEdgeColor', 'b', 'CapSize', 0);
offset = offset+sum(thisPlot);

%average as bar
bar(numel(learnRNido)/2+2*groupOffset, mean(learnRNido));

%output values in text
disp("outside LMAN rate" + mean(learnRNido) + " Hz/trial");
disp("std = " + std(learnRNido)/sqrt(numel(learnRNido)));
[h, p, ci, stats] = ttest(learnRNido);
disp("t-value = "+stats.tstat+", p-value = "+p);