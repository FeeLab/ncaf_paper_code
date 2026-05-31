
function [isiVals, isiDist] = compute_ISI_distribution_all(unitSignal, np_fs)

    isiVals = cell(size(unitSignal, 2), size(unitSignal, 3));
    
    for i = 1:size(unitSignal, 3)
        for j = 1:size(unitSignal, 2)
            thisSignal = unitSignal(:, j, i);
            spikeT = find(thisSignal);
            isiVals{j, i} = diff(spikeT)/np_fs;
        end
    end

    isiDist = cell(size(isiVals, 2), 1);
    for i = 1:size(isiVals, 1)
        for j = 1:size(isiVals, 2)
            isiDist{j} = [isiDist{j}; isiVals{i, j}];
        end
    end

end
