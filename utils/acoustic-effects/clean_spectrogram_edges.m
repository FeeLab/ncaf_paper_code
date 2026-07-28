function upStackNew = clean_spectrogram_edges(upStack)
    % CLEAN_SPECTROGRAM_EDGES Function to clean the edges of a spectrogram
    %
    % Input Arguments:
    %     upStack - input spectrogram data
    %
    % Output Arguments:
    %     upStackNew - cleaned spectrogram data

    % Calculate the mean across the first dimension to get the loudness profile
    loudStack = squeeze(mean(upStack, 1));
    % Compute histogram counts and levels for the loudness profile
    [counts, levels] = histcounts(loudStack(:), ...
        linspace(0,max(loudStack(:)),1e5));
    % Find the index of the maximum count in the histogram
    [~, I] = max(counts);
    % Calculate the gap level based on the histogram
    gapL = 0.5 * levels(I) + 0.5 * levels(I+1);
    % Set the threshold for cleaning edges
    threshL = 5*gapL;
    
    % Loop through each column of the loudness profile
    for m = 1:size(loudStack, 2)
        % Clean the top edge if it exceeds the threshold
        if loudStack(1,m) > threshL
            firstEdge = find(loudStack(:,m) < threshL, 1, "first");
            loudStack(1:firstEdge, m) = gapL; % Set values to gapL
            upStack(:, 1:firstEdge, m) = gapL; % Update the original stack
        end
        % Clean the bottom edge if it exceeds the threshold
        if loudStack(end,m) > threshL
            lastEdge = find(loudStack(:,m) < threshL, 1, "last");
            loudStack(lastEdge:end, m) = gapL; % Set values to gapL
            upStack(:, lastEdge:end, m) = gapL; % Update the original stack
        end
    end
    
    % Return the cleaned spectrogram
    upStackNew = upStack;

end