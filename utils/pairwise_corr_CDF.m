function [corrCDF, xBin, yBin] = pairwise_corr_CDF(r, C)


    Nbin = 500;
    xEdge = linspace(0, 450, Nbin);
    xRes = mean(diff(xEdge));
    xBin = xEdge(1:end-1)+xRes;
    yEdge = linspace(min(C), max(C), Nbin);
    yRes = mean(diff(yEdge));
    yBin = yEdge(1:end-1)+yRes;
    [Ncounts, ~, ~] = histcounts2(C, r, yEdge, xEdge);
    smoothCounts = smoothdata2(Ncounts, 'gaussian', {Nbin/10 Nbin/10});
    corrCDF = cumsum(smoothCounts./repmat(sum(smoothCounts, 1), size(smoothCounts, 1), 1), 1);

end