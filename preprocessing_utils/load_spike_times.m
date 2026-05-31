function [tempNum, spikeT, isSingleUnit] = load_spike_times(file_path, daq_fs, np_fs)
    
    spikesamp = readNPY(fullfile(file_path, 'spike_times.npy'));
    
    tempNum = readNPY(fullfile(file_path, 'spike_clusters.npy'));
    spikeT = double(spikesamp)/np_fs;
    
    t = readtable(fullfile(file_path, 'cluster_group.tsv'), 'FileType', 'text');
    cId = t.cluster_id;
    
    %load single and multiunit
    %cGroup = strcmp(t.group, 'good') | strcmp(t.group, 'mua');

    %load single unit
    cGroup = strcmp(t.group, 'good');

    cGood = cId(cGroup);
    if ~isempty(cGood)
        spikeT = spikeT(ismember(tempNum, cGood));
        tempNum = tempNum(ismember(tempNum, cGood));
    end
    
end