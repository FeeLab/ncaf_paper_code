function MLBounds = find_lman_bounds_waveform(positions, posPeak, ...
    negPeak, waveWidth, defLMANBounds, defNotLMANBounds, C)
% COMMENT REQUIRED HERE!

    arguments
        positions (1,:) double        % Locations of units on shanks
        posPeak (1,:) double          % Positive peak amplitudes
        negPeak (1,:) double          % Negative peak amplitudes
        waveWidth (1,:) double        % Waveform spike widths
        defLMANBounds (1,2) double    % Lower & upper bounds of sure LMAN
        defNotLMANBounds (1,2) double % Lower & upper bounds of not LMAN
        C (2,2) double = [0 1; 1 0]   % Cost matrix
    end

%% Identify Training Units
defLMAN = positions <= defLMANBounds(2) & ...
    positions >= defLMANBounds(1);
defNot = positions > defNotLMANBounds(2) | ...
    positions < defNotLMANBounds(1);

%% Identify Data
XFull = [posPeak; negPeak; waveWidth].';
XFull = normalize(XFull, 1);
XTrain = [XFull(defLMAN, :); XFull(defNot, :)];
Y = [ones(sum(defLMAN), 1); zeros(sum(defNot), 1)];

%% Use SMOTE to Balance Data
rng(42);
[XTrain,Y,~,~] = smote(XTrain, [], 'Class', Y);

%% Train the SVM
Mdl = fitSVMPosterior(fitcsvm(XTrain,Y,"Cost",C));
[~, posteriorProbs] = predict(Mdl, XFull);

%% Find the optimal bounds by maximum likelihood

sortedPositions = sort(positions, 1, "ascend");

boundsToTry = (sortedPositions(1:end-1) + sortedPositions(2:end)) / 2;

searchSpace = -Inf(length(positions));

for i = 1:length(boundsToTry)
    for j = i:length(boundsToTry)

        lb = boundsToTry(i);
        ub = boundsToTry(j);

        searchSpace(i,j) = compute_log_likelihood([lb ub], ...
            positions, posteriorProbs);

    end
end

[~, idx] = max(searchSpace, [], "all");
[iMax, jMax] = ind2sub(size(searchSpace), idx);

MLBounds = [boundsToTry(iMax), boundsToTry(jMax)];

    function logLikelihood = compute_log_likelihood(bounds, ...
            positions, posteriorProbs)

        lowerBound = bounds(1);
        upperBound = bounds(2);

        % Find which units are within the LMAN bounds
        outsideBounds = positions > upperBound | ...
            positions < lowerBound;
        insideBounds = ~outsideBounds;

        % Compute log likelihood
        logLikelihood = sum(log(posteriorProbs(insideBounds,2))) + ...
            sum(log(posteriorProbs(outsideBounds,1)));

    end


end