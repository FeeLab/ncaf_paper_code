function file_name = name_from_path(file_path)

    folders = split(file_path, filesep);
    run_name = folders{end-1};

    catgt_output = dir(fullfile(file_path, strcat(run_name, '_g*_t0.imec0.ap.bin')));
    if ~isempty(catgt_output)
        file_name = extractBefore(catgt_output(1).name, '.imec0.ap.bin');
    else
        catgt_output = dir(fullfile(file_path, strcat(run_name, '_g*_t0.nidq.bin')));
        file_name = extractBefore(catgt_output(1).name, '.nidq.bin');
    end

end