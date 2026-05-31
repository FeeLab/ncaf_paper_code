function shankID = assign_shanks_to_units(file_path, unitCOM)

    file_name = name_from_path(file_path);
    npMeta = SGLX_readMeta.ReadMeta(strcat(file_name, '.imec0.ap.meta'), file_path);

    [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = geomMapToGeom(npMeta);
  
    xCoord = shankInd*shankPitch + xCoord;
    xCoord = xCoord-27; %offset from edge of NP2 probe
    shankI = unique(shankInd);

    shankID = zeros(size(unitCOM, 1), 1);
    shankCenter = zeros(size(shankI));
    for i = 1:numel(shankI)
        shankCenter(i) = (min(xCoord(shankInd==shankI(i)))+max(xCoord(shankInd==shankI(i))))/2;
    end

    for i = 1:numel(shankID)
        [~,minI] = min(abs(unitCOM(i, 1)-shankCenter));
        shankID(i) = shankI(minI);
    end
end