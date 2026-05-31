function isNoise = find_noise_motifs(upStack)

    maxL = max(squeeze(median(upStack, 1)));
    [N, edges] = histcounts(maxL, 100);
    thresh = min(edges)+(max(edges)-min(edges))*otsuthresh(N);

    isNoise = maxL>thresh;

end