function [a,b,c,d,e,mu,sig] = fit_recorded_phi_normed(unitDistTot, xdata)

    phiNorm = zeros(size(xdata, 1), 1);
    targetPhi = zeros(numel(unique(xdata(:, 5))), 1);
    for i = unique(xdata(:, 5))'
        targetPhi(i) = xdata(unitDistTot==0&xdata(:, 5)==i, 1);
        phiNorm(xdata(:, 5)==i) = xdata(xdata(:, 5)==i, 1)/targetPhi(i);
    end

    [f, ~, output] = fit(unitDistTot(unitDistTot>0), phiNorm(unitDistTot>0), 'exp1');
    
    a = f.a;
    b = f.b;
    
    
    fType = fittype("a*exp(b*x)+c");
    fResid = fit(unitDistTot(unitDistTot>0), output.residuals.^2, fType, 'Upper', [Inf 0 Inf], 'Lower', [0 -Inf 0], 'StartPoint', [.1 -.01 0.1]);

    c = fResid.a;
    d = fResid.b;
    e = fResid.c;

    mu = mean(targetPhi);
    sig = std(targetPhi);

    figure;
    hold on;
    plotDist = 0:400;
    plot(plotDist, a*exp(b*plotDist));

    plot(plotDist, a*exp(b*plotDist)+sqrt(c*exp(d*plotDist)+e), '--k');
    plot(plotDist, a*exp(b*plotDist)-sqrt(c*exp(d*plotDist)+e), '--k');
    
    
    scatter(unitDistTot, phiNorm, 'filled');
    xlabel('distance from target neuron');
    ylabel('noise correlation relative to target');
    set(gca, 'TickDir', 'out');
end

