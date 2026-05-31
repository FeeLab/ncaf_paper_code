function unitLoc = find_unit_locations(file_path, unitNum)

    %unit locations are calculated by the spikeinterface analyzer for
    %uncurated unit numberings. DO NOT USE for units that are the result of
    %merges in phy, as these do not correspond to original unit outputs

    unitLoc = readNPY(fullfile(file_path, 'analyzer', 'extensions', 'unit_locations', 'unit_locations.npy'));
    unitLoc = unitLoc(unitNum+1, :);
end