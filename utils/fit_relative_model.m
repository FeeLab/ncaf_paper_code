function g = fit_relative_model(isiHist, binVals)
    
    fitStart = 0.003;
    fitEnd = 0.015;
    f = fit(binVals(binVals>=fitStart & binVals<=fitEnd)', isiHist(binVals>=fitStart & binVals<=fitEnd)', 'exp1');
    
    figure;
    subplot(1, 2, 1);
    plot(binVals, isiHist);
    hold on;
    plot(binVals, f(binVals));
    xlim([0 0.03]);
    ylim([1 Inf]);
    yscale('log');
    xlabel('time (s)')
    title('ISI distribution and exponential fit');
    
    plotMax = 0.01;
    plotVals = binVals(binVals<=plotMax)';
    w_t = isiHist(binVals<=plotMax)'./f(plotVals);
    
    ft = fittype(@(a, b, c, x) max(0, a*(1-exp(-b*(x-c)))));
    
    g = fit(plotVals, w_t, ft, 'Lower', [0 0 0], 'Upper', [2 Inf Inf], 'StartPoint', [1 3000 .001]);
    
    subplot(1, 2, 2);
    plot(plotVals, w_t);
    hold on;
    plot(plotVals, g(plotVals));
    xlim([0 0.01]);
    xlabel('time from previous spike (s)');
    ylabel('recovery function w_t');
    title("T_r = "+1000*g.c +" ms, 1/e = "+1/g.b*1000+" ms");

end