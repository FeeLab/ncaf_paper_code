function [audio, daq_sync, trigT, daq_fs, isCatch] = load_nidaq(file_path, iTrig)

    file_name = name_from_path(file_path);

    daqMeta = SGLX_readMeta.ReadMeta(strcat(file_name, '.nidq.bin'), file_path);
    daq_fs = round(str2double(daqMeta.niSampRate));
    
    iAudio = 2;
    rawdata = SGLX_readMeta.ReadBin(0, Inf, daqMeta, strcat(file_name, '.nidq.bin'), file_path);
    corrdata = SGLX_readMeta.GainCorrectNI(rawdata, [iAudio], daqMeta);
    audio = corrdata(iAudio, :);
    iSync = 1;
    daq_sync = find(diff(rawdata(iSync,:)>max(rawdata(iSync,:))/2)>0)/daq_fs;
    
    if ~exist('iTrig')
        iTrig = 4;
    end

    trigBuff = 1;
    [~, trigT] = findpeaks(single(rawdata(iTrig,:)>max(rawdata(iTrig,:))/2), "WidthReference", "halfheight", "MaxPeakWidth", daq_fs*0.1);
    trigT = trigT/daq_fs;
    trigT = trigT(trigT>trigBuff & trigT<(numel(audio)/daq_fs-trigBuff));
    isCatch = rawdata(iTrig, :)>max(rawdata(iTrig, :)/2);
end