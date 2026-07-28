function specCorr = warp_audio(specWarp, onsets, offsets, ...
    tempOn, tempOff)
%WARPAUDIO Summary of this function goes here
%   Detailed explanation goes here

specCorr = zeros(size(specWarp));   % Initialize output

% For each motif, 
for i = 1:size(specCorr, 3)
  
    % FIRST SILENCE
    % If the first gap is longer than that in the template,
    % then crop off the beginning of the spectrogram
    if onsets(i, 1)>=tempOn(1)
        newStart = onsets(i, 1) - tempOn(1) + 1; % New start point
        specCorr(:, 1:tempOn(1)-1, i) = specWarp(:, ...
            newStart:onsets(i,1)-1, i);


    % Otherwise, pad the beginning with NaNs
    else
        padPoint = tempOn(1) - onsets(i,1);
        specCorr(:, 1:padPoint, i) = NaN;
        specCorr(:, padPoint+1:tempOn(1)-1, i) = specWarp(:, ...
            1:onsets(i, 1)-1, i);
    end

    
    % SYLLABLES
    % For each syllable (equivalently, onset)
    for j = 1:numel(tempOn)
        xq = linspace(onsets(i, j), offsets(i, j), tempOff(j)-tempOn(j)+1);
        specCorr(:, tempOn(j):tempOff(j), i) = interp1(specWarp(:, :, i)', ...
            xq)';
    end

    % GAPS
    % For each gap (equivalently, offset excluding the last)
    for j = 1:numel(tempOff)-1
        specCorr(:, tempOff(j)+1:tempOn(j+1)-1, i) = interp1(specWarp(:, :, i)', linspace(offsets(i, j)+1, onsets(i, j+1)-1, tempOn(j+1)-tempOff(j)-1))';
    end

    % Fill in the last silence
    if offsets(i, end)<=tempOff(end)
        specCorr(:, tempOff(end)+1:end, i) = specWarp(:, offsets(i, end)+1:offsets(i, end)+size(specWarp, 2)-tempOff(end), i);
    else
        specCorr(:, tempOff(end)+1:tempOff(end)+size(specWarp, 2)-offsets(i, end), i) = specWarp(:, offsets(i, end)+1:end, i);
        specCorr(:, tempOff(end)+size(specWarp, 2)-offsets(i, end)+1:end, i) = NaN;
    end

end
end

