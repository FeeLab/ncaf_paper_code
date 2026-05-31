%% import_align_detectLMAN.m
%this script imports data from a spike sorting pipeline, aligns spiking to
%song, and detects the boundaries of LMAN based on that song-aligned
%activity. results are saved to a .mat output file.


%% some initial parameters
close all;
addpath('../utils');
addpath('../preprocessing_utils');
addpath('../external/SpikeGLX_Datafile_Tools/MATLAB');
addpath('../external/npy-matlab/npy-matlab');
addpath('../external/ndSparse_G4_2021_03_16');
addpath('../external/colorcet');
addpath('../external/dkbsl-matlab_smote-41e4193');

%name of folder containing sorted data outputs
run_name = '2025-02-18_11208_Post-Advance_LMAN_BOTM_0';

%change this to location of run folder
file_path = fullfile('/path/to/run/folder', run_name);

catgt_path = fullfile(file_path, 'catgt_outputs');
ks_path = fullfile(file_path, 'phy_outputs');

%where outputs will be saved
out_path = fullfile('../preprocessed_outputs/', run_name);
if ~exist(out_path)
    mkdir(out_path)
end

cW = 0.01; %width of contingent window for this run
%% import and process audio data

%load ADC inputs and sync signals from both streams
[audio, daq_sync, trigT, daq_fs, isCatch] = load_nidaq(catgt_path);
[np_sync, np_fs] = extract_np_sync_fast(catgt_path);

%define start and end of song based on initial motifs
[~, songLength, trigDelay] = define_song_edges(audio, trigT, daq_fs);

%times of song onset
songT = trigT-trigDelay;

%align to onset of target syllable
[onsets, trigT, songT] = align_to_onset(audio, trigT, songT, daq_fs);
songTonset = onsets - median(onsets-songT);
trigTcorr = trigT-songTonset;
trigDelay = median(trigTcorr); %time of contingent window falling edge relative to time in song

%compute spectrograms with 1 ms resolution
[upStack, tspec] = spectrogram_highres(audio, songTonset, daq_fs, songLength);
aTemplate = mean(upStack, 3);

%find noise bursts
isNoise = find_noise_motifs(upStack);
[noiseI, noiseF, isNoise] = find_noise_edges(audio, songTonset, trigDelay, daq_fs, isNoise);

%average song without noise burst
silenceTemplate = mean(upStack(:, :, ~isNoise), 3);

%% plot padwise learning
%calculate threshold crossing signal from pads
[padSignal, padIndices] = align_pads_to_song_spikes(catgt_path, songTonset, daq_sync, np_sync, daq_fs, np_fs, songLength, 4);
%plot rasters
plot_learning_pads(padSignal, np_fs, out_path, aTemplate, padIndices, trigDelay, cW);

%% unitwise learning
%import spike times
[tempNum, spikeT] = load_spike_times(ks_path, daq_fs, np_fs);
%align spikes to song
[unitSignal, unitNum] = align_spikes_to_song(tempNum, songTonset, spikeT, daq_sync, np_sync, daq_fs, np_fs, songLength);
%import unit quality info
isSingleUnit = load_single_units(ks_path, unitNum);
rp_violation = load_rp_violation(ks_path, unitNum);

%% find LMAN bounds using waveform shapes
%import unit locations and waveform shapes
unitCOM = find_unit_locations(file_path, unitNum);
waveforms = load_unit_waveforms(file_path, unitNum);
[waveWidth, posPeak, negPeak] = waveform_width(waveforms, np_fs);

%define LMAN boundaries on each shank independently
shankID = assign_shanks_to_units(catgt_path, unitCOM);
shankInd = unique(shankID);
lmanBounds = zeros(numel(shankInd), 2);
boundMargin = 0.25;
lmanUnits = zeros(size(unitNum), 'logical');
figure;
for i = 1:numel(shankInd)
    plotAx = subplot(1, numel(shankInd), i);
    hold on;
    shankUnits = shankID==shankInd(i);
    %initial bounds from song modulation
    initBounds = lman_bounds_from_pads(padSignal, trigDelay, catgt_path, np_fs, shankInd(i), plotAx);
    defLMAN = initBounds - (initBounds-mean(initBounds))*boundMargin;
    defNotLMAN = initBounds + (initBounds-mean(initBounds))*boundMargin;
    try
        %use SVM to classify waveforms
        lmanBounds(i, :) = find_lman_bounds_waveform(unitCOM(shankUnits, 2), posPeak(shankUnits), negPeak(shankUnits), waveWidth(shankUnits), defLMAN, defNotLMAN);
        yline(lmanBounds(i, :), '-r');
    catch ME
        if strcmp(ME.identifier, "MATLAB:randperm:inputKTooLarge") %%not enough units in putative LMAN
            lmanBounds(i, :) = [0, 0]; %%assume no LMAN on this shank
        else
            rethrow(ME);
        end
    end
    lmanUnits(unitCOM(:, 2)>lmanBounds(i, 1) & unitCOM(:, 2)<lmanBounds(i, 2) & shankUnits) = 1;
    scatter3(plotAx, unitCOM(lmanUnits&shankUnits&isSingleUnit, 1), unitCOM(lmanUnits&shankUnits&isSingleUnit, 2), unitNum(lmanUnits&shankUnits&isSingleUnit), [], 'r', 'filled');
    scatter3(plotAx, unitCOM(lmanUnits&shankUnits&~isSingleUnit, 1), unitCOM(lmanUnits&shankUnits&~isSingleUnit, 2), unitNum(lmanUnits&shankUnits&~isSingleUnit), [], 'r');
end

lmanNum = unitNum(lmanUnits); %unit indices of LMAN units
lmanSingleUnit = isSingleUnit(lmanUnits); %boolean of LMAN unit quality

%% plot unitwise learning
%plot rasters of unitwise activity
plot_learning(unitSignal(:, :, lmanUnits), lmanNum, np_fs, daq_fs, out_path, aTemplate, trigDelay, cW);


%% save alignment output

padSigSparse = ndSparse(padSignal);
unitSigSparse = ndSparse(unitSignal);
save(fullfile(out_path, strcat(run_name, '_aligned.mat')), "run_name",...
    "trigTcorr", "songTonset", "songLength", "aTemplate",...
    "daq_fs", "np_fs", "unitNum", "unitSigSparse",...
    "padSigSparse", "padIndices", "trigDelay",...
    "lmanNum","lmanSingleUnit", "lmanUnits", "isNoise", "unitCOM", "noiseI", ...
    "noiseF", "cW");
