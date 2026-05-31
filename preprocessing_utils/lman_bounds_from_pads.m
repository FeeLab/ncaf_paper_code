function lmanBounds = lman_bounds_from_pads(padSignal, trigDelay, file_path, np_fs, activeShank, plotAx)

    songMod = squeeze(sum(padSignal, 2));
    songMod = songMod(round([1:(trigDelay-.06)*np_fs (trigDelay+.04)*np_fs:size(songMod, 1)]), :);
    songMod = smoothdata(songMod, 'gaussian', 0.01*np_fs);
    padMod = std(songMod)./prctile(songMod,50);
    [~, removeI] = rmoutliers(padMod, "quartiles", 'ThresholdFactor', 4);
    padMod(removeI) = NaN;
    
    [h, edges] = histcounts(padMod, 100);
    thresh = min(edges)+(max(edges)-min(edges))*otsuthresh(h);

    file_name = name_from_path(file_path);
    npMeta = SGLX_readMeta.ReadMeta(strcat(file_name, '.imec0.ap.meta'), file_path);

    [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = geomMapToGeom(npMeta);
    
    xCoord = shankInd*shankPitch + xCoord;
    xCoord = xCoord-27; %offset from edge of NP2 probe
    
    padIndices = get_pad_indices(npMeta);

    shankI = unique(shankInd);
    shankMeans = zeros(numel(shankI), 1);
    for i = 1:numel(shankMeans)
        shankMeans(i) = mean(padMod(shankInd==shankI(i)), 'omitnan');
    end

    isActive = shankInd == activeShank;
    xActive = xCoord(isActive);
    yActive = yCoord(isActive);
    modActive =  padMod(isActive);
    activeCenter = (min(xActive)+max(xActive))/2;

    depthVals = unique(yActive);
    depthPitch = median(diff(depthVals));
    depthMod = zeros(size(depthVals));
    for i = 1:numel(depthMod)
        depthMod(i) = mean(modActive(yActive==depthVals(i)), 'omitnan');
    end

    depthMod = movmedian(depthMod, 3);
    depthMod = [0; depthMod; 0];
    depthVals = [min(depthVals)-depthPitch; depthVals; max(depthVals)+depthPitch];

    aboveThresh = single(depthMod-thresh>0);
    [~, locs, w] = findpeaks(aboveThresh, depthVals);
    [~, I] = sort(w);
    startDepth = locs(I(end))-depthPitch/2;
    endDepth = depthVals(find(depthVals>startDepth & aboveThresh<1, 1))-depthPitch/2;
    lmanBounds = [startDepth endDepth];

    padSize = 80;
    scatter3(plotAx, xCoord, yCoord, padIndices, padSize, padMod, "MarkerFaceColor", "flat", 'MarkerEdgeColor','w');
    view([0 90]);
    colormap(colorcet('L20'));
    clim([0 1]*prctile(modActive,99));
    axis equal;
    xlim([min(xCoord)-60 max(xCoord)+60]);
    ylim([min(yCoord)-20 max(yCoord)+20]);
    yline(lmanBounds, '--w');
    xlabel('AP distance (microns)');
    ylabel('DV distance (microns)');
    cb = colorbar;
    cb.Label.String = 'std of song-locked modulation relative to burst firing rate';
    isEven = mod(padIndices, 2)==0;
    textOffset = 20;
    text(xCoord(isEven)-textOffset, yCoord(isEven), num2str(padIndices(isEven)), 'HorizontalAlignment', 'right', 'Color', 'k');
    text(xCoord(~isEven)+textOffset, yCoord(~isEven), num2str(padIndices(~isEven)), 'HorizontalAlignment', 'left', 'Color', 'k');
    set(gca, 'Color', 'k');
end
