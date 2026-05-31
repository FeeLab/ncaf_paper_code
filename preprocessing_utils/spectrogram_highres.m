function [upStack, tspec] = spectrogram_highres(audio, songT, daq_fs, songLength)

    windowT = .01;
    wwindow = windowT * daq_fs;
    ovlp = .9;
    stepw = windowT*(1-ovlp);
    
    aSnip = zeros(round((songLength+windowT)*daq_fs), numel(songT));
    for i = 1:numel(songT)
        aSnip(:, i) = audio(round((songT(i)-windowT/2)*daq_fs)+(0:size(aSnip, 1)-1));
    end
    

    fvec = 500:10:8000;
    [s,f,tspec] = spectrogram(aSnip(:, 1), windowT*daq_fs, windowT*daq_fs*ovlp, fvec, round(daq_fs));
    upStack = zeros(numel(fvec), size(s, 2), numel(songT));
    parfor i = 1:numel(songT)
        [s,f,t] = spectrogram(aSnip(:, i), windowT*daq_fs, windowT*daq_fs*ovlp, fvec, round(daq_fs));
        s = flipud(abs(s));
        s = log(s - min(s(:)) + 1);
        upStack(:, :, i) = s;
    end

    aTemplate = mean(upStack, 3);
    figure;
    subplot(2, 1, 1);
    imagesc([0 size(aTemplate, 2)-1], [8000 500], aTemplate);
    set(gca, 'Ydir', 'normal');
    set(gca, 'TickDir', 'none');
    xlabel('time (ms)');
    ylabel('frequency (Hz)');
    subplot(2, 1, 2);
    imagesc(squeeze(sum(upStack, 1))');
    xlabel('time (ms)');
    ylabel('motif');
end