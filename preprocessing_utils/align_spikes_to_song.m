function [np_aligned, unitNum] = align_spikes_to_song(tempNum, songT, spikeT, daq_sync, np_sync, daq_fs, np_fs, songLength)
    
    unitNum = unique(tempNum);
    Ntemp = numel(unitNum);
    np_aligned = zeros(round(songLength*np_fs), numel(songT), Ntemp, 'logical');
    eventI = cell(Ntemp, 1);
    for i = 1:numel(eventI)
        eventI{i} = round(spikeT(tempNum==unitNum(i))*np_fs);
    end
    
    for j = 1:numel(songT)
        %{
        prevSync = find(diff(daq_sync>songT(j)));
        syncOffset = double(songT(j)-daq_sync(prevSync));
        [~, I] = min(abs(np_sync - daq_sync(prevSync)/daq_fs*np_fs));
        index = round(np_sync(I)*daq_fs/np_fs) + syncOffset;
        %}
   
        prevSync = find(diff(daq_sync>songT(j)));
        if isempty(prevSync)
            prevSync = 0; %may need to be 1/daq_fs
        end
        syncOffset = songT(j)-daq_sync(prevSync);
        [~, I] = min(abs(np_sync - daq_sync(prevSync)));
        index = round((np_sync(I)+syncOffset)*np_fs);
        for i = 1:Ntemp
            motifEvents = eventI{i}((eventI{i}-index)>=0 & (eventI{i}-index)<(songLength*np_fs))-index;
            np_aligned(motifEvents+1, j, i) = true;
        end
    
    end

end