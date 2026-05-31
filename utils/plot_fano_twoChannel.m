function plot_fano_twoChannel(spikeCount, spikeCountSim1, spikeCountSim2)

    jitterW = .01;

    [datax, datay] = fano_plot_count(spikeCount);
    datax = datax + (rand(numel(datax), 1)-0.5)*jitterW;
    datay = datay + (rand(numel(datax), 1)-0.5)*jitterW;

    [simx1, simy1] = fano_plot_count(spikeCountSim1);
    simx1 = simx1 + (rand(numel(simx1), 1)-0.5)*jitterW;
    simy1 = simy1 + (rand(numel(simy1), 1)-0.5)*jitterW;

    [simx2, simy2] = fano_plot_count(spikeCountSim2);
    simx2 = simx2 + (rand(numel(simx2), 1)-0.5)*jitterW;
    simy2 = simy2 + (rand(numel(simy2), 1)-0.5)*jitterW;

    limVal = ceil(max([datax; datay; simx1; simy1; simx2; simy2]));
    
    binRes = 0.04;
    edgeVals = 0:binRes:limVal;
    binVals = edgeVals(1:end-1)+binRes/2;
    figure;

    gouldMap = colorcet('L20');
    blackVal(1, 1, :) = gouldMap(1, :);

    C = repmat(blackVal, numel(binVals), numel(binVals), 1);
    [Ncounts1, ~, ~] = histcounts2(simy1, simx1, edgeVals, edgeVals, 'Normalization', 'probability');
    C(:, :, 1) = C(:, :, 1) + (1-blackVal(1))*Ncounts1/prctile(Ncounts1(:), 99.9);
    [Ncounts2, ~, ~] = histcounts2(simy2, simx2, edgeVals, edgeVals, 'Normalization', 'probability');
    C(:, :, 2) = C(:, :, 2) + (1-blackVal(2))*Ncounts2/prctile(Ncounts2(:), 99.9)*.5;
    C(:, :, 3) = C(:, :, 3) + (1-blackVal(3))*Ncounts2/prctile(Ncounts2(:), 99.9);

    subplot(1, 3, 1);
    imshow(C, 'XData', [min(edgeVals) max(edgeVals)], 'YData', [min(edgeVals) max(edgeVals)]);

    set(gca, 'YDir', 'normal');
    hold on;
    scatter(datax, datay, [], 'g', '.');
    
    plot([0, limVal], [0, limVal], '--w');
    axis square;
    xlabel('mean spikes in contingent window');
    ylabel('variance of spikes in contingent window');
    xticks(0:limVal);
    yticks(0:limVal);
    set(gca, 'Visible', 'on');

    for i = 0:limVal-1
        xval = (0:.01:1);
        yval = xval-xval.^2;
        fill(xval+i, yval, [1 1 1]*.5);
    end

    ax2 = subplot(1, 3, 2);
    imagesc([1 1], [0 prctile(Ncounts1(:), 99.9)], (0:255)');
    colormap(ax2, [linspace(blackVal(1), 1, 256)' ones(256, 1)*blackVal(2) ones(256, 1)*blackVal(3)]);
    set(gca, 'YDir', 'normal');
    
    ax3 = subplot(1, 3, 3);
    imagesc([1 1], [0 prctile(Ncounts2(:), 99.9)], (0:255)');
    colormap(ax3, [ones(256, 1)*blackVal(1) linspace(blackVal(2), blackVal(2)+0.5*(1-blackVal(2)), 256)' linspace(blackVal(3), 1, 256)']);
    set(gca, 'YDir', 'normal');
end