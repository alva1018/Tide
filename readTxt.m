function [table_gps, gps_id, table_gal, gal_id, table_bds, bds_id] = readTxt(filePath)
    opts = detectImportOptions(filePath);
    opts.VariableNames = {'SBF_M','TOW','WN','SVID','SignalType','AntennaID','PR_m','L_cyc','Doppler_Hz','SNR','LockTime_s'};
    opts.VariableTypes = {'string','double','double','string','string','string','double','double','double','double','double'};
    table = readtable(filePath, opts);
    % DateTime conversion (from TOW and WN to datetime)
    table.TimeStamp = datetime(1980,1,6,0,0,0) + days(table.WN * 7) + seconds(table.TOW);

    %% GPS satellite filtering
    % find GPS satellites only
    isGPS = startsWith(table.SVID, 'G');
    isGPSL1 = contains(table.SignalType, 'GPS_L1CA');
    table_gps = table(isGPS & isGPSL1, :);
    % extract relevant columns
    table_gps.PRN = str2double(strip(strip(table_gps.SVID, "'"), 'left', 'G'));
    gps_id = unique(table_gps.PRN);

    %% Galileo satellite filtering
    % find Galileo satellites only
    isGal = startsWith(table.SVID, 'E');
    isGalE1 = contains(table.SignalType, 'GAL_E1BC');
    table_gal = table(isGal & isGalE1, :);
    % extract relevant columns
    table_gal.PRN = str2double(strip(strip(table_gal.SVID, "'"), 'left', 'E'));
    gal_id = unique(table_gal.PRN);

    %% BDS satellite filtering
    % find BDS satellites only
    % find Galileo satellites only
    isBds = startsWith(table.SVID, 'C');
    isBdsB1 = contains(table.SignalType, 'BDS_B1I');
    table_bds = table(isBds & isBdsB1, :);
    % extract relevant columns
    table_bds.PRN = str2double(strip(strip(table_bds.SVID, "'"), 'left', 'C'));
    bds_id = unique(table_bds.PRN);
end