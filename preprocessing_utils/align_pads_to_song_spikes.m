function [signal, padIndices] = align_pads_to_song_spikes(file_path, songT, daq_sync, np_sync, daq_fs, np_fs, songLength, ffactor, chanSubset)


    file_name = name_from_path(file_path);
    npMeta = SGLX_readMeta.ReadMeta(strcat(file_name, '.imec0.ap.meta'), file_path);

    nChan = str2double(npMeta.nSavedChans);
    if ~exist('chanSubset', 'var')
        chanSubset = 1:nChan-1;
    end
    nFileSamp = str2double(npMeta.fileSizeBytes)/(2*str2double(npMeta.nSavedChans));
    sizeA = [nChan, nFileSamp];
    m = memmapfile(fullfile(file_path, strcat(file_name, '.imec0.ap.bin')), 'Format', {'int16', sizeA, 'x'});


    nplength = round(songLength*np_fs);
    signal = zeros(nplength, numel(songT), numel(chanSubset), 'logical');


    Ni = numel(songT);
    parfor i = 1:Ni
        prevSync = find(diff(daq_sync>songT(i)));
        if isempty(prevSync)
            prevSync = 1;
        end
        syncOffset = songT(i)-daq_sync(prevSync);
        [~, I] = min(abs(np_sync - daq_sync(prevSync)));
        index = round((np_sync(I)+syncOffset)*np_fs);

        snippet = m.Data.x(chanSubset, index:index+nplength-1)';

        thresh = -ffactor * median(abs(snippet)/0.6745, 1);
        signal(:,i,:) = islocalmax(snippet < repmat(thresh, size(snippet, 1), 1), 'FlatSelection', 'center');
    end


    cellPad = extractBetween(npMeta.snsChanMap, "AP", ";");
    padIndices = zeros(numel(cellPad), 1);
    for i = 1:numel(cellPad)
        padIndices(i) = str2num(cellPad{i});
    end
end