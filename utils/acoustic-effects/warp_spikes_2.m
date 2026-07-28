function spikeTWarp = warp_spikes_2(songSpikeTimes, onsetT, ...
    offsetT, tempOnSecs, tempOffSecs, premotorOffset)
%WARP_SPIKES Summary of this function goes here
%   Detailed explanation goes here

% Find those spikes that occur during songs

spikeTWarp = cell(size(songSpikeTimes));

% ons = numel(spikeTWarp);
% offs = numel(spikeTWarp);
% diffs = numel(spikeTWarp);

tempOnSecs = tempOnSecs - premotorOffset;
tempOffSecs = tempOffSecs - premotorOffset;
% 
% j = 1;

% For each motif, 
for m = 1:size(songSpikeTimes, 2)

    %onsetsSecs = tspec(squeeze(onsets(m, :)));
    %offsetsSecs = tspec(squeeze(offsets(m, :)));

    onsetsSecs = squeeze(onsetT(:,m)) - premotorOffset;
    offsetsSecs = squeeze(offsetT(:,m)) - premotorOffset;

    spikes = songSpikeTimes(:, m);
    spikesWarp = cell(size(spikes));

    % For each unit, 
    for n = 1:numel(spikes)
        
        spikesWarp{n} = zeros(size(spikes{n}));

        for i = 1:numel(spikes{n})

            s = spikes{n}(i); % The spike time

            if s < onsetsSecs(1)
                % Handle if it's in the first silence
                offset = s - onsetsSecs(1);
                sCorr = tempOnSecs(1) + offset;
                %if sCorr < 0 
                %    sCorr = NaN; 
                %end
            elseif s >= offsetsSecs(end)
                % Handle if it's in the last silence
                
                offset = s - offsetsSecs(end);
                sCorr = tempOffSecs(end) + offset;
                %if sCorr > songLength 
                %    sCorr = NaN; 
                %end

                %fprintf("Detected! Old time %f, new time %f\n", ...
                %    s, sCorr);

            else
                Ion = find(s >= onsetsSecs, 1, "last");
                Ioff = find(s >= offsetsSecs, 1, "last");
    
                % Syllable
                if any([Ion > Ioff, isempty(Ioff)])
                    offset = s - onsetsSecs(Ion);
    
                    l = offsetsSecs(Ion) - onsetsSecs(Ion);
                    lTemp = tempOffSecs(Ion) - tempOnSecs(Ion);
    
                    sCorr = tempOnSecs(Ion) + lTemp / l * (offset);
                    assert(sCorr >= tempOnSecs(Ion) & sCorr < tempOffSecs(Ion));
                % Gap
                else
                    offset = s - offsetsSecs(Ioff);
    
                    l = onsetsSecs(Ioff+1) - offsetsSecs(Ioff);
                    lTemp = tempOnSecs(Ioff+1) - tempOffSecs(Ioff);
    
                    sCorr = tempOffSecs(Ioff) + lTemp / l * (offset);
                    assert(sCorr >= tempOffSecs(Ioff) & sCorr < tempOnSecs(Ioff+1));
                end
    
            end
            spikesWarp{n}(i) = sCorr;

            % ons(j) = s;
            % offs(j) = sCorr;
            % diffs(j) = sCorr - s;
            % j = j+1;

        end
        
        outsideSong = isnan(spikesWarp{n});
        spikesWarp{n} = spikesWarp{n}(~outsideSong);
        % spikes{n} = spikes{n}(~outsideSong);
        %  if ~all(diff(spikesWarp{n})>0)
        %     figure;
        %     plot(diff(spikes{n}));
        %     hold on;
        %     plot(diff(spikesWarp{n}));
        % end
        
        assert(all(diff(spikesWarp{n})>=0));

    end
    
    % Spike times should still be monotonic
    %
    spikeTWarp(:, m) = spikesWarp;

end

% figure;
% subplot(2,1,1);
% scatter(ons, diffs);
% subplot(2,1,2);
% scatter(offs, diffs);

%spikeTWarp(:, end) = songSpikeTimes(:, end);

end

