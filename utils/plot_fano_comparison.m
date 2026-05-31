function [binVals, Ncounts] = plot_fano_comparison(spikeCount, spikeCountSim, limVal, plotAx)

    jitterW = .01;

    [datax, datay] = fano_plot_count(spikeCount);
    datax = datax + (rand(numel(datax), 1)-0.5)*jitterW;
    datay = datay + (rand(numel(datax), 1)-0.5)*jitterW;


    [simx, simy] = fano_plot_count(spikeCountSim);
    simx = simx + (rand(numel(simx), 1)-0.5)*jitterW;
    simy = simy + (rand(numel(simy), 1)-0.5)*jitterW;

    if ~exist("limVal")
        limVal = ceil(max([datax; datay; simx; simy]));
    end
    
    
    binRes = 0.04;
    edgeVals = 0:binRes:limVal;
    binVals = edgeVals(1:end-1)+binRes/2;

    [Ncounts, ~, ~] = histcounts2(simy, simx, edgeVals, edgeVals, 'Normalization', 'probability');
    if ~isempty(plotAx)
        imagesc(plotAx, binVals, binVals, Ncounts);
        clim(plotAx, [0 prctile(Ncounts(:), 99.9)])
        colormap(plotAx, colorcet('L20'));
        colorbar(plotAx);
        set(plotAx, 'YDir', 'normal');
        hold on;
        scatter(plotAx, datax, datay, [], 'r', '.');
        
        plot(plotAx, [0, limVal], [0, limVal], '--w');
        xlim(plotAx, [0 limVal]);
        axis square;
        xlabel('mean spikes in contingent window');
        ylabel('variance of spikes in contingent window');
        
        for i = 0:limVal-1
            xval = (0:.01:1);
            yval = xval-xval.^2;
            fill(plotAx, xval+i, yval, 'w');
        end
    end
    
end