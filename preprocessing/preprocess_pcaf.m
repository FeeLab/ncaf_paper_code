%% Load Data
birdName = "11220"; % Change based on bird

addpath("../classes/"); % pCAF experiment class
addpath("../utils");

filename = sprintf("../presorted_data/pCAF/extracted/%s-extracted.mat", birdName);
load(filename);

filename = sprintf("../presorted_data/pCAF/extracted/%s-bad-dates.mat", birdName);
load(filename);

%% Correct for Circadian Fluctuations
dates = keys(birdExperiments);

rng("default");

% Initialize arrays of times-in-day and pi tches
allTimes = [];
allPitches = [];

% Create figure
figure;
hold on

% For each date, 
for i = 1:numel(dates)

    date = dates(i);
    experiment = birdExperiments(date);

    p = experiment.pitch;     % Pitches
    t = experiment.trigTimes; % Times

    % Remove outliers 
    toUse = ~isoutlier(p);
    p = p(toUse);
    t = t(toUse);

    % Subtract mean pitch
    p = p - mean(p);

    % Add times and pitches to array
    allTimes = [allTimes, t];
    allPitches = [allPitches, p];

    % Plot times and pitches
    scatter(t / 3600, p);

end

% To prevent overfitting, we jitter all times by +-30 mins
% Jitter times by adding random noise within the specified range
jitterAmount = 30 * 60; % 30 minutes in seconds
jitteredTimes = allTimes + (rand(size(allTimes)) * 2 - 1) * jitterAmount;

lengthScale = 6 * 60 * 60; % 6 Hours
signalStd = std(allPitches) / sqrt(2);

% Fit Gaussian Process to pitches
gprMdl = fitrgp( ...
    jitteredTimes.', ...
    allPitches.', ...
    "KernelFunction", "squaredexponential", ...
    "KernelParameters", [lengthScale, signalStd]);

% Predict circadian pitch fluctation
pred = predict(gprMdl, allTimes.').';

% Plot prediction
scatter(allTimes / 3600, pred, 'k');

xlabel("Time of Day (Hours)");
ylabel("Mean-Subtracted Pitch");
hold off;

%% Remove Bad Dates
for i = 1:length(badDates)

    birdExperiments = remove(birdExperiments, badDates{i});

end

% Get updated list of dates
dates = keys(birdExperiments);

%% Perform Preprocessing
nBoot = 1000; % Bootstrap replicates

% Initialize arrays
delays = zeros(1,numel(dates));
LRResidual = zeros(1, numel(dates));
LRResidualUB = zeros(1, numel(dates));
LRResidualLB = zeros(1, numel(dates));
LRResidualStd = zeros(1, numel(dates));

corrs = zeros(1,numel(dates));
numMotifs = zeros(1,numel(dates));
directions = strings(1, numel(dates));
startPitches = zeros(1,numel(dates));
bootstats = zeros(nBoot, numel(dates));

for i = 1:numel(dates)

    date = dates(i);
    data = birdExperiments(date);

    % Get noise latency
    tmpLatency = data.noiseLatency;
    tmpLatency(tmpLatency < 0) = NaN;
    delays(i) = median(tmpLatency, 'omitmissing');

    % For 11325, song features meant that noise burst latency
    % is systematically underestimated on 2025-07-02 and 
    % 2025-07-03 by ~20ms. We correct for that here.
    if strcmp(birdName, "11325") & ...
            (strcmp(date, "2025-07-02") || ...
            (strcmp(date, "2025-07-03")))
        delays(i) = delays(i) + 20e-3;
    end
    
    % Get correlation with noise
    isNoise = data.isNoise;
    pitch = data.pitch;
    toUse = ~isoutlier(pitch) & ~(data.noiseLatency < -0.05);
    pitch = pitch(toUse);
    isNoise = isNoise(toUse);

    % Get moment-to-moment correlations
    corrs(i) = corr(diff(pitch).', diff(isNoise).', "Type", "Spearman");

    % Calculate LRs using Theil-Sen Estimator
    ps = data.pitch;
    toUse = ~isoutlier(ps);
    Noise = data.isNoise(toUse);
    ps = ps(toUse);
    ts = data.trigTimes(toUse);

    res = ps - predict(gprMdl, ts.').';
    startPitches(i) = median(res(1:100));

    x = 1:length(data.pitch);
    x = x(toUse);

    d = [x.', res.']; 
    LRResidual(i) = repeatedMedian(d);

    [ci, bootstat] = bootci(nBoot, {@repeatedMedian, d}, "Type", "per");
    LRResidualLB(i) = ci(1);
    LRResidualUB(i) = ci(2);
    LRResidualStd(i) = std(bootstat);
    bootstats(:,i) = bootstat;

    % Number of motifs
    numMotifs(i) = data.numMotifs;

    % Noise contingency direction
    directions(i) = data.direction;

end

LRCorr = LRResidual;
LRCorrUB = LRResidualUB;
LRCorrLB = LRResidualLB;

%% Save Results
clear bird;

bird.LR = LRCorr;
bird.LRBoot = bootstats;
bird.corrs = corrs;
bird.LRStd = LRResidualStd;
bird.directions = directions;
bird.delays = delays;
bird.startPitches = startPitches;

filename = sprintf("%s-processed.mat", birdName);
filepath = fullfile("../presorted_data/pCAF/processed/", ...
    filename);

if ~isfile(filepath)
    fprintf("Saving processed data for %s\n", birdName);
    save(filepath, "bird");
else
    fprintf("Skipping saving processed data for %s\n", birdName);
end