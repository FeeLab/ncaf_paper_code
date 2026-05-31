function [songStack, songLength, trigDelay] = define_song_edges(audio, trigT, daq_fs)

    w = 1;
    specL = 2*w*daq_fs+1;
    aSnip = zeros(specL, numel(trigT));
    for i = 1:numel(trigT)
        aSnip(:, i) = audio(round((trigT(i)-w)*daq_fs):round((trigT(i)+w)*daq_fs));
    end

    windowT = 10e-3;
    ovlp = .90;
    stepw = windowT*(1-ovlp);
    fvec = 500:10:8000;
    N = min(50, numel(trigT));
    upStack = zeros(numel(fvec), floor((specL-windowT*ovlp*daq_fs)/(stepw*daq_fs)), N);
    parfor i = 1:N
        [s,f,t] = spectrogram(aSnip(:, i), windowT*daq_fs, windowT*daq_fs*ovlp, fvec, round(daq_fs));
        s = flipud(abs(s));
        s = log(s - min(s(:)) + 1);
        upStack(:, :, i) = s;
    end

    imagesc(mean(upStack, 3));
    zoom xon;
    set(gcf, 'Visible', 'on');
    disp("zoom to desired location, then press any key to continue");
    pause;
    ax = gca;
    disp("select region for song template, then press any key to update");
    roi = drawrectangle(ax);
    pause;
    position = roi.Position;
    xmin = round(position(1));
    xmax = floor(position(1) + position(3));
    songL = xmax-xmin;
    songLength = songL*stepw;
    trigDelay = ((size(upStack, 2)-1)/2 - xmin)*stepw;

    aCrop = aSnip(((specL+1)/2 - round(trigDelay*daq_fs)) + (round(1-windowT/2*daq_fs):round((songLength+windowT/2-stepw)*daq_fs)), :);

    songStack = zeros(numel(fvec), songL, numel(trigT));
    [s,f,t] = spectrogram(aCrop(:, 1), windowT*daq_fs, windowT*daq_fs*ovlp, fvec, round(daq_fs));
    parfor i = 1:numel(trigT)
        [s,f,t] = spectrogram(aCrop(:, i), windowT*daq_fs, windowT*daq_fs*ovlp, fvec, round(daq_fs));
        s = flipud(abs(s));
        s = log(s - min(s(:)) + 1);
        songStack(:, :, i) = s;
    end

    imagesc(t*1000, fvec, mean(songStack, 3));
    xlabel('time (ms)');
    ylabel('frequency (Hz)');
    xline(trigDelay*1000, '--r');

end





