%% Run analysis on example bird
bird = "10977";
set(groot, 'defaultFigureVisible', 'off');
run("../preprocessing/acoustic_effects.m");
set(groot, 'defaultFigureVisible', 'on');

%% Correlation Matrix
u = 1;
units = find(communities == validComms(u));

corrMatDisplay = corrMat;
corrMatDisplay(logical(eye(size(corrMatDisplay)))) = 0;

figure;
imagesc(corrMatDisplay);
dr = max(abs(corrMatDisplay), [], "all");
colormap(colorcet('D06')); clim([-dr dr]);
axis square

for i = 1:length(validComms)
    xline(find(communities == validComms(i)), ...
        "Color", cmap(i,:));
end
colorbar();

%% Example Difference Spectrogram
contingentRow = 43;

% Get activity over the units
tActivity = squeeze(mean(bins(:,resultInterleave,units),3));

% Do cumulative sum to 
tActivityCum = [zeros(1, size(tActivity,2),'like',tActivity); cumsum(tActivity,1)];

i = contingentRow;
winStart = (i-1) * activityStepBins + 1;
winEnd = winStart + activityWinBins - 1;
tActivityWin = squeeze(tActivityCum(winEnd+1,:) - tActivityCum(winStart,:));

% sort activity and split
[~, idxSorted] = sort(tActivityWin(:), 'ascend');
highActMotifs = idxSorted(1:splitPoint);
lowActMotifs  = idxSorted(splitPoint+1:end);

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

diffSpec = highMean - lowMean;
diffSpec(isnan(diffSpec)) = 0;

% Display result
dr = max(abs(diffSpec), [], "all");
diffSpec = diffSpec / dr;

figure;
ax1 = subplot(2,1,1);
imagesc(tspec*1e3, fspec / 1e3, squeeze(mean(specCorr, 3)));
set(gca, "YDir", "normal");
xlabel('Time (ms)');
ylabel('Frequency (kHz)');
xlim([tspec(1)*1e3, tspec(tempOn(chopPoint+1))*1e3]);
colormap(ax1, colorcet('L01'));
colorbar(); % Just so they align


ax2 = subplot(2,1,2);
imagesc(tspec * 1e3, fspec / 1e3, diffSpec);
xlabel('Time (ms)');
ylabel('Frequency (kHz)');
set(gca, "YDir", "normal");
colormap(ax2, colorcet('D01A'));
clim([-1 1]);
xline(rowCentres(contingentRow)*1e3 + 5);
xline(rowCentres(contingentRow)*1e3 - 5);
xlim([tspec(1)*1e3, tspec(tempOn(chopPoint+1))*1e3]);
colorbar();

%% Full Difference Matrix
zScoredMat = squeeze(diffSpecMagsZScored(u,:,:));
zScoredMat(isnan(zScoredMat)) = 0;
dr = max(abs(zScoredMat), [], "all");

figure;
imagesc(tspec*1e3, rowCentres*1e3, zScoredMat);
xlabel('Time of acoustic difference (ms)');
ylabel('Time of LMAN activity');
xlim([tspec(1)*1e3, tspec(tempOn(chopPoint+1))*1e3]);
%ylim([tspec(1)*1e3, tspec(tempOff(chopPoint))*1e3]);
%colormap(colorcet('D13'));
colormap(colorcet('D01A'));
hold on;
plot([0 songLength*1e3], [0 songLength*1e3], 'k--');
plot([0 songLength*1e3] + 21.6, [0 songLength*1e3], 'g--');
clim([-dr dr]);
colorbar();

% Add line to mark location of previous
yline(rowCentres(contingentRow)*1e3);