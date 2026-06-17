%% Set random seed
rng(0);

%% import functions

addpath('../utils')
addpath('../external/SpikeGLX_Datafile_Tools/MATLAB');
addpath('../external/npy-matlab/npy-matlab');
addpath('../external/ndSparse_G4_2021_03_16');
addpath('../external/colorcet');

%% LMAN summary

%targeted single units from Fig 1
unitVec = [187, 868, 73];
runVec = ["2024-04-27_10872_LMAN-X_nCAF_2_aligned.mat";
    "2025-02-18_11208_Post-Advance_LMAN_BOTM_0_aligned.mat";
    "2025-11-12_Yellow33_Post-Advance_LMAN_BOTM_0_aligned.mat"];

for k = 1:numel(unitVec)
    unitI = unitVec(k);
    load(fullfile('../presorted_data', runVec(k)));
    unitSignal = full(unitSigSparse);
    plot_learning_figure(unitSignal, trigDelay, unitI, unitNum, silenceTemplate, cW, np_fs, noiseI, noiseW);
end