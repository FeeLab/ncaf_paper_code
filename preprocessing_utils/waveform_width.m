function [widths, posPeak, negPeak] = waveform_width(templates, fs)
%WAVEFORM_WIDTH Summary of this function goes here
%   Detailed explanation goes here

    widths = zeros(1, size(templates, 1));
    posPeak = zeros(1, size(templates, 1));
    negPeak = zeros(1, size(templates, 1));

    for cl = 1:size(templates, 1)

        maxDev = squeeze(min(templates(cl, :, :), [], 2));
        
        [~, peakCh] = min(maxDev(:));

        waveform = -squeeze(templates(cl, :, peakCh));
        waveform = resample(waveform, 10, 1);

        [p, I1] = findpeaks(waveform, "SortStr", "descend");        
        waveform2 = -waveform(I1(1)+1:end);
        [t, I2] = findpeaks(waveform2, "SortStr", "descend");

        negPeak(cl) = p(1);
        posPeak(cl) = t(1);

        %[~, ~, w, ~] = findpeaks(waveform, ...
        %    "WidthReference", 'halfheight', ...
        %    "SortStr", "descend");
        %w = w(1);

        %widths(cl) = w / (fs * 10);
        if isempty(I2)
            widths(cl) = NaN;
        else
            widths(cl) = (I2(1)) / (fs * 10);
            %peakToTrough(cl) = p(1)/t(1);
            % plot(-waveform);
            % xline(I1(1));
            % xline(I1(1) + I2(1));
            % title(widths(cl) * 1e6);
            % pause;
        end

    end

end

