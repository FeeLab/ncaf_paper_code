function [burstTimes, isiBurst, isiPoisson, fracBurst, ISICat] = find_bursts_in_signal(unitSignal, np_fs, trigDelay, cW, sThresh, isiThresh, blankW)
    
    tmin = round((trigDelay-cW-blankW)*np_fs);
    tmax = round((trigDelay+blankW)*np_fs);
    
    ISICat = [];
    sts = cell(size(unitSignal, 2), 1);
    spikeNum = cell(size(unitSignal, 2), 1);
    spikeMotif = [];
    spikeN = 1;
    for i = 1:size(unitSignal, 2)
        ISICat = [ISICat diff(find(unitSignal(1:tmin, i))/np_fs)'];
        ISICat = [ISICat diff(find(unitSignal(tmax:end, i))/np_fs)'];
    
        spikeTimes = find(unitSignal(:, i));
        spikeBefore = spikeTimes(spikeTimes<=tmin);
        spikeAfter = spikeTimes(spikeTimes>=tmax);
        sts{i} = [spikeBefore; spikeAfter];
    
        spikeMotif = [spikeMotif ones(1, max(0, numel(spikeBefore)-1))*i];
        spikeMotif = [spikeMotif ones(1, max(0, numel(spikeAfter)-1))*(i+0.5)];
    
        if numel(spikeBefore)>1
            spikeNum{i} = [spikeN:spikeN+numel(spikeBefore)-2 0];
            spikeN = spikeN + numel(spikeBefore)-1;
        elseif numel(spikeBefore)==1
            spikeNum{i} = 0;
        end
        if numel(spikeAfter)>1
            spikeNum{i} = [spikeNum{i} spikeN:spikeN+numel(spikeAfter)-2 0];
            spikeN = spikeN + numel(spikeAfter)-1;
        elseif numel(spikeAfter)==1
            spikeNum{i} = [spikeNum{i} 0];
        end
    
    end
    
    [archive_burst_RS,archive_burst_length,archive_burst_start]=rank_surprise_burst(ISICat, spikeMotif, isiThresh, sThresh);
    
    burstTimes = cell(size(unitSignal, 2), 1);
    for i = 1:numel(archive_burst_start)
        thisMotif = floor(spikeMotif(archive_burst_start(i)));
        burstI = sts{thisMotif}(spikeNum{thisMotif}==archive_burst_start(i));
        burstF = sts{thisMotif}(find(spikeNum{thisMotif}==archive_burst_start(i))+archive_burst_length(i)-1);
        burstTimes{thisMotif} = [burstTimes{thisMotif}; [burstI burstF]];
    end


    isiBurst = [];
    isiPoisson = ISICat;
    for i = 1:numel(archive_burst_start)
        isiBurst = [isiBurst, ISICat(archive_burst_start(i):(archive_burst_start(i)+archive_burst_length(i)-2))];
        isiPoisson(archive_burst_start(i):(archive_burst_start(i)+archive_burst_length(i)-2)) = NaN;
    end
    isiPoisson(isnan(isiPoisson))=[];

    fracBurst = sum(archive_burst_length)/sum(unitSignal([1:tmin tmax:end], :), 'all');
end