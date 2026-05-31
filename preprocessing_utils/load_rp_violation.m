function rp_violation = load_rp_violation(file_path, unitNum)

    t = readtable(fullfile(file_path, 'cluster_sliding_rp_violation.tsv'), 'FileType', 'text');
    cId = t.cluster_id;
    rp_violation_all = t.sliding_rp_violation;
    rp_violation = zeros(size(unitNum));
    for i = 1:numel(rp_violation)
        if ismember(unitNum(i), cId)
            rp_violation(i) = rp_violation_all(cId==unitNum(i));
        else
            rp_violation(i) = NaN;
        end
    end
end