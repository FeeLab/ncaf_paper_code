function signalDec = downsample_spikes(unitSignal, binW, pShift)

    signalMean = smoothdata(unitSignal, 1, 'gaussian', 5*binW/2.355);
    signalDec = downsample(signalMean, binW, pShift);

end