function isSingleUnit = load_single_units(file_path, unitNum)

    t = readtable(fullfile(file_path, 'cluster_group.tsv'), 'FileType', 'text');
    cId = t.cluster_id;
    cGroup = strcmp(t.group, 'good');

    isSingleUnit = ismember(unitNum, cId(cGroup));
end