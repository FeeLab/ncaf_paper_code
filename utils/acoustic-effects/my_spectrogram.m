function [s,f,t] = my_spectrogram(audio, Fs)
%MY_SPECTROGRAM Compute song spectrogram
%   Params:
%       audio   - Audio signal to compute spectrogram of
%       Fs      - Sampling frequency of audio (Hz)
%   Returns:
%       s   - Normalized spectrogram
%       f   - Frequency values (Hz)
%       t   - Time values (s)

% Constants
windowT = 10e-3;    % Time window
ovlp = .90;         % Overlap proportion
fvec = 500:10:8000; % Frequencies 

% Convert constants into samples
windowSamps = round(windowT * Fs);
ovlpSamps = round(windowT * ovlp * Fs);

% Compute the spectrogram
[s,f,t] = spectrogram(audio, windowSamps, ovlpSamps, fvec, Fs);

% Normalize and take the log to improve contrast
s = abs(s);
s = log(s-min(s(:))+1);

end

