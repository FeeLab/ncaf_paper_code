function [upStack, fspec, tspec] = create_motif_spectrograms(audio, ...
    songT, daqFs, songLength)
    %CREATE_MOTIF_SPECTROGRAMS Create spectrograms for audio motifs
    %   Inputs:
    %       audio - audio signal data
    %       songT - time points for the motifs
    %       daqFs - data acquisition sampling frequency
    %       songLength - length of the audio snippet to analyze
    %   Outputs:
    %       upStack - 3D array of spectrograms for each motif
    %       fspec - frequency specification from the spectrogram
    %       tspec - time specification from the spectrogram

    % Preallocate a matrix to hold audio snippets for each motif
    aSnip = zeros(round(songLength*daqFs), numel(songT));
    % Extract audio snippets based on the specified time points
    for i = 1:numel(songT)
        aSnip(:, i) = audio(round(songT(i)*daqFs)+(0:size(aSnip, 1)-1));
    end
    % Compute the spectrogram for the first audio snippet
    [s,fspec,tspec] = my_spectrogram(aSnip(:, 1), daqFs);
    % Preallocate a 3D array to hold all spectrograms
    upStack = zeros(size(s,1), size(s,2), numel(songT));
    upStack(:,:,1) = s; % Store the first spectrogram
    % Compute spectrograms for the remaining audio snippets in parallel
    parfor i = 2:numel(songT)
        [s,~,~] = my_spectrogram(aSnip(:, i), daqFs);
        upStack(:, :, i) = s; % Store each spectrogram
    end
end

