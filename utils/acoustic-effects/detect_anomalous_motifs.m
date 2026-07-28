function loudDtw = detect_anomalous_motifs(upStack, tspec, noiseOn, noisePeriod)
%DETECT_ANOMALOUS_MOTIFS Summary of this function goes here
%   Detailed explanation goes here

loudStack = squeeze(sum(upStack, 1));

if noiseOn
    loudStack(noisePeriod, :) = 0.25 * max(loudStack(:));
end

loudStack = loudStack - repmat(mean(loudStack, 1), size(loudStack, 1), 1);
loudStack = loudStack./repmat(sqrt(sum(loudStack.^2, 1)), size(loudStack, 1), 1);
loudTemp = mean(loudStack, 2);
loudTemp = loudTemp-mean(loudTemp);
loudTemp = loudTemp/sqrt(sum(loudTemp.^2));

figure;
plot(tspec*1e3, loudTemp);
title('Loudness Template');
xlabel('Time (ms)');
ylabel('Loudness (au)');

% if noiseOn
%     xline(noiseStart*1e3, 'r--');
%     xline(noiseEnd*1e3, 'r--');
% end

loudCorr = zeros(2*numel(loudTemp)-1, size(loudStack, 2));
loudDtw = zeros(size(loudStack, 2), 1);
for i = 1:size(loudCorr, 2)
    loudCorr(:, i) = xcorr(loudStack(:, i), loudTemp);
    loudDtw(i) = dtw(loudStack(:, i), loudTemp);
end

figure;
histogram(loudDtw);
title('Distance from Template');
xlabel('DTW Distance');
ylabel('Counts');

end

