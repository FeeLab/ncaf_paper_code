function padIndices = get_pad_indices(npMeta)
    cellPad = extractBetween(npMeta.snsChanMap, "AP", ";");
    padIndices = zeros(numel(cellPad), 1);
    for i = 1:numel(cellPad)
        padIndices(i) = str2num(cellPad{i});
    end
end