function waveforms = load_unit_waveforms(file_path, unitNum)

    waveforms = readNPY(fullfile(file_path, 'analyzer', 'extensions', 'templates', 'average.npy'));

    waveforms = waveforms(unitNum+1, :, :);
end