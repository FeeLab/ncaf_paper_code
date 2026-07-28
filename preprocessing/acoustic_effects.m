%% Set the random seed
rng(42);

%% Get utils
addpath('../utils');
addpath('../utils/acoustic-effects/')
addpath(genpath('../external'));

%% Switch based on which bird we're analysing
% If no bird is selected, return
if ~exist('bird', 'var')
    fprintf("No bird selected. Terminating.");
    return;
end

% Set variables based on selected bird
switch bird

    case "11238"
        runName = '2025-03-23_11238_LMAN_nCAF_0';
        splitMotifs = false;
    
    case "10977"
        runName = "2024-07-24_10977_LMAN_nCAF";
        splitMotifs = false;

    case "Yellow33"
        runName = "2025-11-12_Yellow33_Post-Advance_LMAN_BOTM_0";
        splitMotifs = true;

        if ~exist("motifVariant", "var")
            fprintf("No motif variant given. Defaulting to 1.");
            motifVariant = 1;
        end
    
    otherwise
        fprintf("Invalid bird selected. Terminating.");
        return;

end

%% Load Data
% Load audio
load(fullfile("../presorted_data/audio_files/", ...
    sprintf("%s-audio.mat", bird)));

% Load sorted data
load(sprintf('../presorted_data/%s_aligned.mat', runName));

fprintf("Data loaded.\n");

%% Generate spectrograms
% Compute spectrograms
[upStack, fspec, tspec] = create_motif_spectrograms(audio, ...
    songTonset, daq_fs, songLength);

% Remove truncated syllables from edges
upStackCleaned = clean_spectrogram_edges(upStack);

% Show the new mean spectrogram
figure;
imagesc(tspec*1e3, fspec/1e3, squeeze(mean(upStackCleaned, 3)));
set(gca,'YDir','normal'); 
colormap(colorcet('Gouldian'));
title('Mean Spectrogram');
xlabel('Time (ms)');
ylabel('Frequency (kHz)');
xline(trigDelay*1e3, 'r--', "LineWidth",2);
xline(trigDelay*1e3 - cW*1e3, 'r--', "LineWidth",2);

%% Detect DAF (White Noise Bursts)
loudStack = squeeze(sum(upStackCleaned, 1));
loudVar = var(loudStack, [], 2);

levels = statelevels(loudVar);
noisePeriod = (loudVar - levels(1)) > 0.75 * diff(levels);

noiseStart = tspec(find(diff(noisePeriod) == 1, 1, 'first')) - 5e-3;
noiseEnd = tspec(find(diff(noisePeriod) == -1, 1, 'last')) + 5e-3;

noisePeriod = tspec >= noiseStart & tspec <= noiseEnd;

figure;
plot(1e3*tspec, loudVar);
hold on;
plot(1e3*tspec, noisePeriod * max(levels,[],"all"));
title("Noise Burst Detection");
xlabel("Time (ms)")
ylabel("Loudness variance (a.u.)");

%% Filter out bad spectrograms
% Identify anomalous motifs
loudDtw = detect_anomalous_motifs(upStackCleaned, tspec, true, ...
    noisePeriod);
specGood = upStackCleaned(:, :, ~isoutlier(loudDtw));
unitGood = unitSigSparse(:,~isoutlier(loudDtw), :);
isNoise = isNoise(~isoutlier(loudDtw));

% Plot the loudness traces for the original and reduced 
% groups of motifs.
figure;
subplot(1, 2, 1);
imagesc(tspec*1e3, [], squeeze(sum(upStackCleaned, 1)).');
colormap(colorcet('Gouldian'));
title('Original');
xlabel('Time (ms)');
subplot(1, 2, 2);
imagesc(tspec*1e3, [], squeeze(sum(specGood, 1)).');
colormap(colorcet('Gouldian'));
title('Outliers Removed');
xlabel("Time (ms)");

%% Identify Syllables and Gaps
loudStack = squeeze(mean(specGood,1));

[counts, levels] = histcounts(loudStack(:), ...
    linspace(0,max(loudStack(:)),1e5));
[~, I] = max(counts);
gapL = 0.5 * levels(I) + 0.5 * levels(I+1);

aveSongLevel = prctile(loudStack(~noisePeriod), 75);

loudStack(noisePeriod,:) = aveSongLevel;

gapW = 5;

levelMultipliers = linspace(1,25,100);

gapNum = zeros(size(loudStack, 2), numel(levelMultipliers));

for j = 1:numel(levelMultipliers)

    for i = 1:size(gapNum, 1)
        pks = findpeaks(single(loudStack(:, i)<gapL * levelMultipliers(j)), ...
            'MinPeakWidth', gapW, ...
            'WidthReference', 'halfheight');
        gapNum(i,j) = numel(pks);
    end

end

% Find the threshold value with the most agreement
counts = zeros(numel(levelMultipliers), 1);
for j = 1:numel(levelMultipliers)

    win = mode(gapNum(:,j));

    if win == 0
        counts(j) = 0;
    else
        counts(j) = sum(gapNum(:,j) == win);
    end
end

[~, idx] = max(counts);

threshL = gapL * levelMultipliers(idx);

gapNumOptimal = squeeze(gapNum(:,idx));

figure;
histogram(gapNumOptimal);
title('Gap Number Histogram');
xlabel('Gap Number');
ylabel('Count');

gapNum = gapNumOptimal;

Ngap = mode(gapNum);
specWarp = specGood(:, :, gapNum==Ngap);
loudStack = loudStack(:,gapNum==Ngap);
unitGood = unitGood(:,gapNum==Ngap, :);
isNoise = isNoise(gapNum==Ngap);

pks = zeros(size(loudStack, 2), Ngap);
for i = 1:size(loudStack, 2)
    [~, pks(i, :)] = findpeaks(single(loudStack(:, i)<threshL), ...
        'MinPeakWidth', gapW, 'WidthReference', 'halfheight');
    gapNum(i) = numel(pks);
end

onsets = zeros(size(pks, 1), Ngap+1);
offsets = zeros(size(pks, 1), Ngap+1);
offsets(:, 1:end-1) = pks-1;
for i = 1:size(pks, 1)
    onsets(i, 1) = find(loudStack(:, i)>threshL, 1);
    for j = 1:size(pks, 2)
        onsets(i, j+1) = find(loudStack(pks(i, j):end, i)>threshL, 1)+pks(i, j)-1;
    end
    offsets(i, end) = find(loudStack(:, i)>threshL, 1, 'last');
end

onsetOutliers = any(isoutlier(onsets, "percentiles", [5 95], 1), 2);
toKeep = ~onsetOutliers;

specWarp = specWarp(:, :, toKeep);
loudStack = loudStack(:,toKeep);
offsets = offsets(toKeep, :);
onsets = onsets(toKeep, :);
unitWarp = unitGood(:, toKeep, :);
isNoise = isNoise(toKeep);

%% Split into Motif Variants if Relevant
% As Yellow33 has two motif variants, we need to pick just one
if splitMotifs

    features = [onsets(:,1:2) offsets(:,1:2)];
    variants = kmeans(features, 2);

    toKeep = (variants == motifVariant);
    specWarp = specWarp(:, :, toKeep);
    loudStack = loudStack(:,toKeep);
    offsets = offsets(toKeep, :);
    onsets = onsets(toKeep, :);
    unitWarp = unitWarp(:, toKeep, :);

end

%% Get Template to Warp to
loudTemp = mean(loudStack, 2);
tempOn = zeros(Ngap+1, 1);
tempOff = zeros(Ngap+1, 1);
[~, pkTemp] = findpeaks(single(loudTemp<threshL), 'MinPeakWidth', gapW, 'WidthReference', 'halfheight');
tempOff(1:end-1) = pkTemp-1;
tempOn(1) = find(loudTemp>threshL, 1);
for j = 1:numel(pkTemp)
    tempOn(j+1) = find(loudTemp(pkTemp(j):end)>threshL, 1)+pkTemp(j)-1;
end
tempOff(end) = find(loudTemp>threshL, 1, 'last');

figure;
plot(loudTemp);
xline(tempOn, 'g');
xline(tempOff, 'r');
yline(threshL);
title('Detected Onsets and Offsets in Template');

%% Warp Audio
specCorr = warp_audio(specWarp, onsets, offsets, ...
    tempOn, tempOff);

% Plot warped spectrograms
figure;
subplot(1, 2, 1);
imagesc(tspec*1e3, [], squeeze(sum(specWarp, 1)).');
colormap(colorcet('Gouldian'));
title('Original');
xlabel('Time (ms)');
subplot(1, 2, 2);
imagesc(tspec*1e3, [], squeeze(sum(specCorr, 1)).');
colormap(colorcet('Gouldian'));
title('Warped');
xlabel("Time (ms)");


%% Time warp
onsetT = tspec(onsets);
offsetT = tspec(offsets);

tempOnT = tspec(tempOn);
tempOffT = tspec(tempOff);

onsetSp = round(onsetT * np_fs) + 1;
offsetSp = round(offsetT * np_fs) + 1;
tempOnSp = round(tempOnT * np_fs) + 1;
tempOffSp = round(tempOffT * np_fs) + 1;

unitWarp = unitWarp(:,:,lmanUnits);

% Faster construction of spike time cell array 'st' from logical unitWarp
% unitWarp is [time x motif x unit]
[nT, nMotifs, nUnits] = size(unitWarp);
st = cell(nUnits, nMotifs);

% Find linear indices of true entries and convert once
[idx_r, idx_c, idx_d] = ind2sub(size(unitWarp), find(unitWarp));
% Compute times in seconds
times = (idx_r - 1) / np_fs;

% Group by (unit, motif) using linear indexing into st: position = sub2ind(size(st), unit, motif)
pos = sub2ind([nUnits, nMotifs], idx_d, idx_c);

% Preallocate temporary cell array to accumulate times per position using accumarray
% accumarray requires numeric values; use cell arrays by collecting indices first
grp = accumarray(pos(:), (1:numel(pos)).', [], @(x) {x});

% Fill st using grouped indices
for k = 1:numel(grp)
    if isempty(grp{k})
        continue
    end
    inds = grp{k};
    % Convert linear position back to (unit, motif) to index st
    [u, m] = ind2sub([nUnits, nMotifs], k);
    st{u, m} = times(inds).';
end

spikeWarpOffset = 0; % Can offset spike times to song times. We do not here.
unitCorr = warp_spikes_2(st, onsetT.', offsetT.', tempOnT, tempOffT, 0);

fprintf("Preprocessing done.\n");

%% Bin Activity
% We restrict ourselves to the last syllable before nCAF
% is happening
dangerTime = trigDelay - cW; % Buffer
chopPoint = find(tempOffT < dangerTime, 1, 'last');
chopTime = tempOffT(chopPoint);

binEdgesMotif = 0:1e-3:chopTime;

binCentres = (binEdgesMotif(2:end) + binEdgesMotif(1:end-1)) / 2;

bins = zeros(length(binEdgesMotif)-1,size(unitCorr,2),size(unitCorr,1)); 

for m = 1:size(unitCorr,1)

    for s = 1:size(unitCorr,2)

        bins(:,s,m) = histcounts(cell2mat(unitCorr(m,s)),binEdgesMotif);

    end

end

%% Find Correlated Ensembles
% Sort units by location
unitLocs = unitCOM(lmanUnits,:);

% Identify how many shanks are in the recording
hozRange = range(unitLocs(:,1));
nShanks = ceil(hozRange / 250);

% Split by shank, then sort by depth
shank = kmeans(unitLocs(:,1), nShanks);

fullOrder = zeros(1, size(unitLocs,1));
currentIdx = 0;

for i = 1:nShanks

    shankUnits = find(shank == i);
    [~, idx] = sort(unitLocs(shankUnits,2));

    fullOrder(currentIdx + (1:length(shankUnits))) = shankUnits(idx);
    currentIdx = currentIdx + length(shankUnits);

end

% Reorder bins
bins = bins(:,:,fullOrder);

%% Compute Correlations
binConcat = permute(bins, [3 1 2]);
binConcat = reshape(binConcat, size(binConcat, 1), []);

binConcatSubSample = zeros(size(binConcat, 1), floor(size(binConcat,2) / 10));

for i = 1:size(binConcatSubSample, 2)
    binConcatSubSample(:,i) = sum(binConcat(:,(1:10) + (i-1)*10), 2);
end

corrMat = corr(binConcatSubSample.', binConcatSubSample.', "Type", "Spearman");

%% Cluster Units Based on Correlations
threshold = 0.05;
A = corrMat > threshold;
A = A - diag(diag(A));

B = corrMat;
B = B - diag(diag(B));
B(B < threshold) = 0;
B(isnan(B)) = 0;

[communities, ~] = simple_louvain(B, threshold);
[counts, values] = groupcounts(communities);

% Identify clusters that have only one component
isolates = values((counts == 1));

% Restrict to communities with 5 or more units
validComms = values(counts>=5);
validCounts = counts(counts>=5);

[validCounts, idx] = sort(validCounts, "descend");
validComms = validComms(idx);

strengths = zeros(size(validCounts));

for i = 1:length(strengths)

    commID = validComms(i);
    toUse = communities == commID;
    subB = B(toUse, toUse);
    strengths(i) = mean(sum(subB,2));

end

[strengths, idx] = sort(strengths, "descend");
validCounts = validCounts(idx);
validComms = validComms(idx);

nComms = length(validComms);
cmap = colorcet('R1', 'N', nComms);

%% Display Ensembles
figure;
imagesc(corrMat);
colormap(colorcet('D01A')); clim([-1 1]);
for i = 1:length(validComms)
    xline(find(communities == validComms(i)), ...
        "Color", cmap(i,:));
end
title('Correlation Matrix with Ensembles');

%% Perform Acoustic Effects Analysis
numEnsembles = nComms;

activityWin = 10e-3; 
activityWinBins = round(activityWin / (binEdgesMotif(2) - binEdgesMotif(1)));

activityStep = 5e-3;
activityStepBins = round(activityStep / (binEdgesMotif(2) - binEdgesMotif(1)));

nTimes = size(specCorr,2);

% Build interleaved index mapping once
nSpecsOrig = size(specCorr, 3);

idx1 = 1:2:nSpecsOrig;
idx2 = fliplr(2:2:nSpecsOrig);
commonSize = min(length(idx1), length(idx2));
tmpIdx = [idx1(1:commonSize); idx2(1:commonSize)];
resultInterleave = tmpIdx(:).';
if length(idx1) > length(idx2)
    resultInterleave = [resultInterleave, idx1(end)];
end

numWins = floor((size(bins,1) - activityWinBins) / activityStepBins);

diffSpecMags = zeros(numEnsembles, numWins, nTimes);
diffSpecMagsPosShift = zeros(numEnsembles, numWins, nTimes); % Shift control
diffSpecMagsNegShift = zeros(numEnsembles, numWins, nTimes); % Shift control

% Convert spec once to single and reorder third dim
specSingle = single(specCorr(:,:,resultInterleave)); % freq x time x trials(interleaved)
nSpecs = size(specSingle,3);
splitPoint = floor(nSpecs/2);

pool = gcp('nocreate');
if isempty(pool) || pool.NumWorkers > 4
    delete(pool);
    parpool(4);
end

parfor u = 1:numEnsembles
    %units = find(ismember(communities, validComms(S(u,:))));
    units = find(communities == validComms(u));

    % Get activity over the units
    tActivity = squeeze(mean(bins(:,resultInterleave,units),3));
    
    % Do cumulative sum to 
    tActivityCum = [zeros(1, size(tActivity,2),'like',tActivity); cumsum(tActivity,1)];
    
    % Preallocate diffSpecs for this unit only (keeps memory low)
    diffSpecs = zeros(numWins, size(specCorr,1), size(specCorr,2), 'single');
    diffSpecsPosShift = zeros(numWins, size(specCorr,1), size(specCorr,2), 'single');
    diffSpecsNegShift = zeros(numWins, size(specCorr,1), size(specCorr,2), 'single');
    
    for i = 1:numWins
    
        winStart = (i-1) * activityStepBins + 1;
        winEnd = winStart + activityWinBins - 1;
        tActivityWin = squeeze(tActivityCum(winEnd+1,:) - tActivityCum(winStart,:));
        
        % sort activity and split
        [~, idxSorted] = sort(tActivityWin(:), 'ascend');
        
        % CONTROL
        for shift = [-1 0 1]

            highActMotifs = idxSorted(1:splitPoint);
            lowActMotifs  = idxSorted(splitPoint+1:end);

            highActMotifs = mod(highActMotifs-1+shift,nSpecs) + 1;
            lowActMotifs = mod(lowActMotifs-1+shift,nSpecs) + 1;
    
            if isempty(highActMotifs)
                highMean = zeros(size(specSingle,1), size(specSingle,2), 'single');
            else
                highMean = mean(specSingle(:,:,highActMotifs), 3);
            end
            if isempty(lowActMotifs)
                lowMean = zeros(size(specSingle,1), size(specSingle,2), 'single');
            else
                lowMean = mean(specSingle(:,:,lowActMotifs), 3);
            end
            
            switch shift
                case -1
                    diffSpecsNegShift(i,:,:) = highMean - lowMean; % freq x time
                case 0
                    diffSpecs(i,:,:) = highMean - lowMean;
                case 1
                    diffSpecsPosShift(i,:,:) = highMean - lowMean;
            end

        end

    end
    
    diffSpecMags(u, :, :) = squeeze(vecnorm(permute(diffSpecs, [1 3 2]), 2, 3)); % numWins x time
    diffSpecMagsPosShift(u, :, :) = squeeze(vecnorm(permute(diffSpecsPosShift, [1 3 2]), 2, 3)); % numWins x time
    diffSpecMagsNegShift(u, :, :) = squeeze(vecnorm(permute(diffSpecsNegShift, [1 3 2]), 2, 3)); % numWins x time
end

fprintf("Initial analysis done.\n");

%% Get Row Centres
rowCentres = zeros(numWins,1);

for i = 1:numWins
    winSamps = (1:activityWinBins) + (i-1) * activityStepBins;
    cent = mean(binCentres(winSamps));
    rowCentres(i) = cent;

end

%% Bootstrap for Normalization of Acoustic Differences
controlTrials = 1000;
nWins = size(diffSpecMags,1);
nFreqs = size(diffSpecMags,2);
nTimes = size(specCorr,2);

% Preallocate randDiffMags in single precision to save memory
randDiffMags = zeros(controlTrials, nTimes, 'single');

% Ensure a pool exists with desired number of workers (use local cluster threads/workers)
pool = gcp('nocreate');
if isempty(pool)
    parpool('local', 18);
elseif pool.NumWorkers ~= 18
    delete(pool);
    parpool('local', 18);
end

% Precompute constants
nSpecs = size(specCorr, 3);

% Convert specCorr to single to reduce memory (if not already)
% specCorr is freq x time x trials
specCorrSingle = single(specCorr);
% Permute to trials x freq x time for convenient indexing in parfor
specCorrPerm = permute(specCorrSingle, [3,1,2]); % trials x freq x time

% Work per window: avoid creating large temporaries; compute vector norms incrementally
controlTrialsLocal = controlTrials;

parfor t = 1:controlTrialsLocal
    rng(t + 1e4, 'combRecursive');

    indices = randperm(nSpecs, nSpecs);
    highIdx = indices(1:splitPoint);
    lowIdx = indices(splitPoint+1:end);

    highMean = squeeze(mean(specCorrPerm(highIdx,:,:),1)); % freq x time
    lowMean  = squeeze(mean(specCorrPerm(lowIdx,:,:),1));  % freq x time

    diffSpec = highMean - lowMean;                          % freq x time
    tmp = vecnorm(diffSpec, 2, 1).';                   % 1 x nFreqs
    randDiffMags(t,:) = tmp;
end

% Clear large temporaries to free memory
clear specCorrSingle specCorrPerm tmp diffSpec highMean lowMean

%% Z-Score Acoustic Differences
meanDiffs = mean(log(randDiffMags), 1);
stdDiffs = std(log(randDiffMags), 1);

meanDiffs = reshape(meanDiffs, [1, size(meanDiffs)]);
stdDiffs = reshape(stdDiffs, [1, size(stdDiffs)]);

diffSpecMagsZScored = (log(diffSpecMags) - meanDiffs) ./ stdDiffs;
diffSpecMagsPosShiftZScored = (log(diffSpecMagsPosShift) - meanDiffs) ./ stdDiffs;
diffSpecMagsNegShiftZScored = (log(diffSpecMagsNegShift) - meanDiffs) ./ stdDiffs;

fprintf("Normalization done.\n");

% Zero-out gaps
gapStarts = [1; tempOff(1:chopPoint)];
gapEnds = [tempOn(1:chopPoint); size(diffSpecMags, 3)];

for i = 1:length(gapStarts)

    diffSpecMagsZScored(:, :, gapStarts(i):gapEnds(i)) = NaN;
    diffSpecMagsPosShiftZScored(:, :, gapStarts(i):gapEnds(i)) = NaN;
    diffSpecMagsNegShiftZScored(:, :, gapStarts(i):gapEnds(i)) = NaN;

end

%% Plot Normalized Differences
% Per ensemble
dr = max(abs(diffSpecMagsZScored), [], "all");
for u = 1:size(diffSpecMags, 1)

    figure;
    imagesc(tspec*1e3, rowCentres*1e3, squeeze(diffSpecMagsZScored(u,:,:)));
    colormap(colorcet('D13'));
    hold on;
    plot([0 songLength*1e3], [0 songLength*1e3], 'k--');
    clim([-dr dr]);

end

% Mean across ensembles
totalDiffs = squeeze(mean(diffSpecMagsZScored,1));
totalDiffsPosShift = squeeze(mean(diffSpecMagsPosShiftZScored,1));
totalDiffsNegShift = squeeze(mean(diffSpecMagsNegShiftZScored,1));

figure;
imagesc(tspec*1e3, rowCentres*1e3, totalDiffs);
colormap(colorcet('D13'));
hold on;
plot([0 songLength*1e3], [0 songLength*1e3], 'k--');
clim([-dr dr]);

%% Shift Matrices
paddedMat = nan(size(totalDiffs,1), size(totalDiffs,2) * 3);
paddedMat(:,size(totalDiffs,2)+1:2*size(totalDiffs,2)) = totalDiffs;

paddedMatPosShift = nan(size(totalDiffs,1), size(totalDiffs,2) * 3);
paddedMatPosShift(:,size(totalDiffs,2)+1:2*size(totalDiffs,2)) = totalDiffsPosShift;

paddedMatNegShift = nan(size(totalDiffs,1), size(totalDiffs,2) * 3);
paddedMatNegShift(:,size(totalDiffs,2)+1:2*size(totalDiffs,2)) = totalDiffsNegShift;

colsPerRow = activityStep / (tspec(2)-tspec(1));

for i = 1:size(paddedMat,1)
    paddedMat(i,:) = circshift(paddedMat(i,:), -(i-1)*colsPerRow, 2);
    paddedMatPosShift(i,:) = circshift(paddedMatPosShift(i,:), -(i-1)*colsPerRow, 2);
    paddedMatNegShift(i,:) = circshift(paddedMatNegShift(i,:), -(i-1)*colsPerRow, 2);
end

figure;
imagesc(paddedMat);

%% Compute and Plot Difference Curves
rng(0);

sumEffect = mean(paddedMat, 1, "omitmissing");
sumEffect = sumEffect(1:2*length(tspec));
[ci, bootstat] = bootci(1e3, @bootFun, paddedMat);
ci = ci(:,1:2*length(tspec));
bootstat = bootstat(:,1:2*length(tspec));

sumEffectPosShift = mean(paddedMatPosShift, 1, "omitmissing");
sumEffectPosShift = sumEffectPosShift(1:2*length(tspec));
[ciPosShift, bootstatPosShift] = bootci(1e3, @bootFun, paddedMatPosShift);
ciPosShift = ciPosShift(:,1:2*length(tspec));
bootstatPosShift = bootstatPosShift(:,1:2*length(tspec));

sumEffectNegShift = mean(paddedMatNegShift, 1, "omitmissing");
sumEffectNegShift = sumEffectNegShift(1:2*length(tspec));
[ciNegShift, bootstatNegShift] = bootci(1e3, @bootFun, paddedMatNegShift);
ciNegShift = ciNegShift(:,1:2*length(tspec));
bootstatNegShift = bootstatNegShift(:,1:2*length(tspec));

tspecRes = tspec(2) - tspec(1);
halfAxis = tspec - rowCentres(1);
leftHalfAxis = fliplr(-tspecRes * (1:length(halfAxis))) + halfAxis(1);

tAxis = [leftHalfAxis, halfAxis];

figure;
hold on;

% True effect
xConf = [tAxis, fliplr(tAxis)];
yConf = [ci(2,:), fliplr(ci(1,:))];
fill(xConf, yConf, [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Pos Shift
xConf = [tAxis, fliplr(tAxis)];
yConf = [ciPosShift(2,:), fliplr(ciPosShift(1,:))];
fill(xConf, yConf, [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Neg shift
xConf = [tAxis, fliplr(tAxis)];
yConf = [ciNegShift(2,:), fliplr(ciNegShift(1,:))];
fill(xConf, yConf, [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

plot(tAxis, sumEffect);
plot(tAxis, sumEffectPosShift);
plot(tAxis, sumEffectNegShift);
xlim([-0.05 0.1]);

hold off;
xlabel("Time relative to LMAN activity (s)");
ylabel("Acoustic difference (a.u.)");

%% Save outputs
% General information
acousticEffectsStruct.bird = bird;
acousticEffectsStruct.tAxis = tAxis;

% Zero-shift effect
acousticEffectsStruct.effect = sumEffect;
acousticEffectsStruct.ci = ci;
acousticEffectsStruct.bootstat = bootstat;

% +1 motif shift effect
acousticEffectsStruct.effectPS = sumEffectPosShift;
acousticEffectsStruct.ciPS = ciPosShift;
acousticEffectsStruct.bootstatPS = bootstatPosShift;

% -1 motif shift effect
acousticEffectsStruct.effectNS = sumEffectNegShift;
acousticEffectsStruct.ciNS = ciPosShift;
acousticEffectsStruct.bootstatNS = bootstatNegShift;

% Generate filename
if splitMotifs
    outputName = sprintf("%s-V%i-acoustic-effects.mat", bird, motifVariant);
else
    outputName = sprintf("%s-acoustic-effects.mat", bird);
end

outputPath = fullfile("../presorted_data/acoustic-effect-curves/", ...
    outputName);

% Save result
if ~isfile(outputPath)
    save(outputPath, "acousticEffectsStruct");
    disp('File saved successfully');
else
    disp('File already exists. Save skipped');
end

%% Helper functions
function m = bootFun(x)

    m = mean(x, 1, "omitmissing");
    m(isnan(m)) = 0;

end