function [rPair, CPair] = pairwiseCorrelations(fileName)
    
    load(fileName);
    unitSignal = full(unitSigSparse);
    
    
    tw = .020;
    Nmax = 100;
    
    filtSignal = smoothdata(unitSignal(:, 1:Nmax, ismember(unitNum, lmanNum(lmanSingleUnit))), 1, 'movmean', tw*np_fs);
    meanSignal = mean(filtSignal, 2);
    
    timeseries = zeros(size(filtSignal, 1)*Nmax, size(filtSignal, 3));
    for i = 1:Nmax
        timeseries((i-1)*size(filtSignal, 1)+1:i*size(filtSignal, 1), :) = filtSignal(:, i, :)-meanSignal;
    end
    
    N = size(timeseries, 2);
    corrVals = zeros(N, N);
    for i = 1:N
        for j = 1:N
            corrVals(i, j) = corr(timeseries(:, i), timeseries(:, j));
        end
    end

    unitLocs = unitCOM(ismember(unitNum, lmanNum(lmanSingleUnit)), :);
    
    distMatrix = zeros(size(corrVals));
    for i = 1:N
        for j = 1:N
            distMatrix(i, j) = sqrt(sum((unitLocs(i, :) - unitLocs(j, :)).^2));
        end
    end
    
    upperVals = logical(triu(ones(N), 1));
    
    rPair = distMatrix(upperVals);
    CPair = corrVals(upperVals);

end