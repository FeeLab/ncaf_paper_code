function [np_sync, np_fs] = extract_np_sync_fast(file_path)

    file_name = name_from_path(file_path);
    npMeta = SGLX_readMeta.ReadMeta(strcat(file_name, '.imec0.ap.meta'), file_path);
    np_fs = str2num(npMeta.imSampRate);

    np_sync = readmatrix(fullfile(file_path, strcat(file_name, '.imec0.ap.xd_384_6_500.txt')));
end