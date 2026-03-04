close all; clc;
% clear;
%% data load
if exist('table_gps', 'var') && exist('gps_id', 'var') && ...
   exist('table_gal', 'var') && exist('gal_id', 'var') && ...
   exist('table_bds', 'var') && exist('bds_id', 'var')
    disp('Sat data already loaded.');
else
    disp('Loading Sat data...');
    [table_gps, gps_id, table_gal, gal_id, table_bds, bds_id] = readTxt("./log2.txt");
end

%% calculate Azimuth and Elevation angles
startTime = datetime(2026,1,22,3,30,0);
endTime = datetime(2026,1,22,10,30,0);
time = startTime : seconds(1) : endTime;
% calculate azimuth and elevation based on ephemeris data
site_pos_lla = [22.4302766120, 114.2089969609, 5.627475]; % PolyU site position
if exist('Az_gps', 'var') && exist('El_gps', 'var') ...
    && exist('Az_gal', 'var') && exist('El_gal', 'var') ...
    && exist('Az_bds', 'var') && exist('El_bds', 'var')
    disp('Azimuth and Elevation already calculated.');
else
    disp('Calculating Azimuth and Elevation...');
    [Az_gps, El_gps] = calculate_gps_az_el(time, gps_id, site_pos_lla);
    [Az_gal, El_gal] = calculate_gal_az_el(time, gal_id, site_pos_lla);
    [Az_bds, El_bds] = calculate_bds_az_el(time, bds_id, site_pos_lla);
end
% figure;
El_range = [5 25];
Az_range = [-10 50];
% plotSky(Az_gps, El_gps, gps_id, 'o', Az_range, El_range);
% figure;
% plotSky(Az_gal, El_gal, gal_id, 's', Az_range, El_range);
% figure;
% plotSky(Az_bds, El_bds, bds_id, '^', Az_range, El_range);

% Tide Polynomial Fitting
[tideTime, tideHeight, tideTime_real, tideHeight_real] = tideFitting();
% tideHeight is in mcd, convert to WGS84
tideHeight = tideHeight - 0.146 - 4.1; % mcd to HKPD coordinate
tideHeight_real = tideHeight_real - 0.146 - 4.1;
tideHeight = Hkpd2Wgs84(tideHeight, site_pos_lla);
tideHeight_real = Hkpd2Wgs84(tideHeight_real, site_pos_lla);

tide.tideTime = tideTime;
tide.tideHeight = tideHeight;
tide.tideTime_real = tideTime_real;
tide.tideHeight_real = tideHeight_real;

% results structure initialization
results = struct('PRN', {}, 'EstimatedTime', {}, 'EstimatedHeight', {}, ...
                 'SNR', {}, 'TimeData', {}, 'TideHeight', {}, 'ElData', {});
results_gps = repmat(results, length(gps_id), 1);
results_gal = repmat(results, length(gal_id), 1);
results_bds = repmat(results, length(bds_id), 1);
plotId = [];
%% Reflector height estimation
for i = 1:length(gps_id)
    if gps_id(i) ~= 1 && gps_id(i) ~= 2 && gps_id(i) ~= 27
        continue;
    end
    disp(['Processing GPS Satellite PRN ' num2str(gps_id(i)) ' ...']);
    satID = gps_id(i);
    % calculate estimated height using single antenna data
    idx = table_gps.PRN == satID;
    satData = table_gps(idx, :);
    timeData = satData.TimeStamp;
    snrData = satData.SNR;
    % extract EL and AZ data for the current satellite
    elData = El_gps(i, :);
    azData = Az_gps(i, :);
    % find the valid EL data points
    validIdx = nan(length(timeData), 1);
    for j = 1:length(timeData)
        t_idx = find(time == timeData(j));
        if ~isempty(t_idx)
            validIdx(j) = t_idx;
        end
    end
    elData = elData(validIdx);
    azData = azData(validIdx);
    % find the corresponding tide height
    validIdx_tide = nan(length(timeData), 1);
    for j = 1:length(timeData)
        t_idx = find(tideTime == timeData(j));
        if ~isempty(t_idx)
            validIdx_tide(j) = t_idx;
        end
    end
    current_tideHeight = tideHeight(validIdx_tide);
    % find the specific time points during the specific elevation angle range
    mask = (elData >= El_range(1)) & (elData <= El_range(2)) & ...
           ((azData >= mod(Az_range(1), 360)) | (azData <= Az_range(2)));
    if isempty(find(mask, 1))
        continue;
    end
    timeData = timeData(mask);
    t_last = seconds(timeData(end) - timeData(1));
    if t_last < 1800  % 0.5 hour
        % skip if the data duration is less than 0.5 hour
        continue;
    end
    snrData = snrData(mask);
    elData = elData(mask);
    elData = eleCorr(elData);
    current_tideHeight = current_tideHeight(mask);
    % convert dB-Hz to linear scale
    linSNR = 10.^(snrData / 20);

    % polyfit to remove the trend of elevation angle
    secondArray = seconds(timeData - timeData(1));
    p = polyfit(secondArray', linSNR, 2);
    trend = polyval(p, secondArray');
    % detrend the direct signal component (2nd order polynomial)
    detrendedLinSNR = linSNR - trend';
    % compute sine of elevation angle
    sinE = sind(elData);

    %% calculate the integer Ambiguity
    light_speed = 299792458; % m/s
    lambda = light_speed / 1575.42e6; % GPS L1 wavelength

    % %% Lomb-Scargle periodogram
    % % define the searching window of height
    % h_min = -20;
    % h_max = 0;
    % % convert height to frequency range
    % f_min = 2 * h_min / lambda;
    % f_max = 2 * h_max / lambda;
    % % frequency resolution
    % delta_h = 0.01;
    % h_vec = h_min : delta_h : h_max;
    % f_vec = 2 * h_vec / lambda;
    % 
    % % Set a window to calculate Lomb-Scargle periodogram
    % windowSize = 480; % number of data points in each window
    % stepSize = 50; % step size for moving the window
    % numPoints = length(detrendedLinSNR);
    % estHeight_all = [];
    % estTime_all = [];
    % est_PRN = [];
    % for startIdx = 1 : stepSize : (numPoints - windowSize + 1)
    %     endIdx = startIdx + windowSize - 1;
    %     meanIdx = round((startIdx + endIdx) / 2);
    %     currentTime = timeData(meanIdx);
    %     sinE_window = sinE(startIdx:endIdx);
    %     detrendedLinSNR_window = detrendedLinSNR(startIdx:endIdx);
    % 
    %     % Lomb-Scargle periodogram for the current window
    %     [Pxx_win, F_win] = lomb_scargle_custom(sinE_window, detrendedLinSNR_window, f_vec);
    % 
    %     % find the peak frequency in the current window
    %     [peakPower_win, peakIdx_win] = max(Pxx_win);
    %     peakFreq_win = F_win(peakIdx_win);
    % 
    %     estHeight_win = (peakFreq_win * lambda) / 2;
    %     estHeight_win = estHeight_win + site_pos_lla(3); % add antenna height offset
    %     estTime_all = [estTime_all; currentTime];
    %     estHeight_all = [estHeight_all; estHeight_win];
    % end
    % % save the result into the structure
    % results_gps(i).PRN = satID;
    % results_gps(i).EstimatedTime = estTime_all;
    % results_gps(i).EstimatedHeight = estHeight_all;
    % results_gps(i).SNR = snrData;
    % results_gps(i).TimeData = timeData;
    % results_gps(i).TideHeight = current_tideHeight; % transfer from HK mcd to wgs84
    % results_gps(i).ElData = elData;

    [N, frac_phase] = solve_ambiguity(current_tideHeight, sinE.', site_pos_lla(3), lambda);
    % Output estimated height based on integer ambiguity
    plotId = [plotId satID];
    H1(length(plotId)) = N2Height(timeData, detrendedLinSNR, sinE, N + frac_phase, lambda, site_pos_lla);
    
end

% plotSNR(results_gps);
plotHeight(H1, plotId, tide);

%% Reflector height estimation (for GAL satellites) 26 33
plotId = [];
for i = 1:length(gal_id)
    disp(['Processing GAL Satellite PRN ' num2str(gal_id(i)) ' ...']);
    satID = gal_id(i);
    % calculate estimated height using single antenna data
    idx = table_gal.PRN == satID;
    satData = table_gal(idx, :);
    timeData = satData.TimeStamp;
    snrData = satData.SNR;
    % extract EL and AZ data for the current satellite
    elData = El_gal(i, :);
    azData = Az_gal(i, :);
    % find the valid EL data points
    validIdx = nan(length(timeData), 1);
    for j = 1:length(timeData)
        t_idx = find(time == timeData(j));
        if ~isempty(t_idx)
            validIdx(j) = t_idx;
        end
    end
    elData = elData(validIdx);
    azData = azData(validIdx);
    % find the corresponding tide height
    validIdx_tide = nan(length(timeData), 1);
    for j = 1:length(timeData)
        t_idx = find(tideTime == timeData(j));
        if ~isempty(t_idx)
            validIdx_tide(j) = t_idx;
        end
    end
    current_tideHeight = tideHeight(validIdx_tide);
    % find the specific time points during the specific elevation angle range
    mask = (elData >= El_range(1)) & (elData <= El_range(2)) & ...
           ((azData >= mod(Az_range(1), 360)) | (azData <= Az_range(2)));
    if isempty(find(mask, 1))
        continue;
    end
    timeData = timeData(mask);
    t_last = seconds(timeData(end) - timeData(1));
    if t_last < 1800  % 0.5 hour
        % skip if the data duration is less than 0.5 hour
        continue;
    end
    snrData = snrData(mask);
    elData = elData(mask);
    current_tideHeight = current_tideHeight(mask);
    % convert dB-Hz to linear scale
    linSNR = 10.^(snrData / 20);

    % polyfit to remove the trend of elevation angle
    secondArray = seconds(timeData - timeData(1));
    p = polyfit(secondArray', linSNR, 2);
    trend = polyval(p, secondArray');
    % detrend the direct signal component (2nd order polynomial)
    detrendedLinSNR = linSNR - trend';
    % compute sine of elevation angle
    sinE = sind(elData);

    %% calculate the integer Ambiguity
    light_speed = 299792458; % m/s
    lambda = light_speed / 1575.42e6; % GAL E1 wavelength

    % %% Lomb-Scargle periodogram
    % % define the searching window of height
    % h_min = -20;
    % h_max = 0;
    % % convert height to frequency range
    % f_min = 2 * h_min / lambda;
    % f_max = 2 * h_max / lambda;
    % % frequency resolution
    % delta_h = 0.01;
    % h_vec = h_min : delta_h : h_max;
    % f_vec = 2 * h_vec / lambda;
    % 
    % % Set a window to calculate Lomb-Scargle periodogram
    % windowSize = 480; % number of data points in each window
    % stepSize = 50; % step size for moving the window
    % numPoints = length(detrendedLinSNR);
    % estHeight_all = [];
    % estTime_all = [];
    % est_PRN = [];
    % for startIdx = 1 : stepSize : (numPoints - windowSize + 1)
    %     endIdx = startIdx + windowSize - 1;
    %     meanIdx = round((startIdx + endIdx) / 2);
    %     currentTime = timeData(meanIdx);
    %     sinE_window = sinE(startIdx:endIdx);
    %     detrendedLinSNR_window = detrendedLinSNR(startIdx:endIdx);
    % 
    %     % Lomb-Scargle periodogram for the current window
    %     [Pxx_win, F_win] = lomb_scargle_custom(sinE_window, detrendedLinSNR_window, f_vec);
    % 
    %     % find the peak frequency in the current window
    %     [peakPower_win, peakIdx_win] = max(Pxx_win);
    %     peakFreq_win = F_win(peakIdx_win);
    % 
    %     estHeight_win = (peakFreq_win * lambda) / 2;
    %     estHeight_win = estHeight_win + site_pos_lla(3); % add antenna height offset
    %     estTime_all = [estTime_all; currentTime];
    %     estHeight_all = [estHeight_all; estHeight_win];
    % end
    % % save the result into the structure
    % results_gps(i).PRN = satID;
    % results_gps(i).EstimatedTime = estTime_all;
    % results_gps(i).EstimatedHeight = estHeight_all;
    % results_gps(i).SNR = snrData;
    % results_gps(i).TimeData = timeData;
    % results_gps(i).TideHeight = current_tideHeight; % transfer from HK mcd to wgs84
    % results_gps(i).ElData = elData;

    [N, frac_phase] = solve_ambiguity(current_tideHeight, sinE.', site_pos_lla(3), lambda);
    % Output estimated height based on integer ambiguity
    plotId = [plotId satID];
    H1(length(plotId)) = N2Height(timeData, detrendedLinSNR, sinE, N + frac_phase, lambda, site_pos_lla);
    
end

% plotSNR(results_gps);
plotHeight(H1, plotId, tide);

%% Reflector height estimation (for BDS satellites)
plotId = [];
for i = 1:length(bds_id)
    if bds_id(i) ~= 11 && bds_id(i) ~= 12 && bds_id(i) ~= 34 && bds_id(i) ~= 43
        continue;
    end
    disp(['Processing BDS Satellite PRN ' num2str(bds_id(i)) ' ...']);
    satID = bds_id(i);
    % calculate estimated height using single antenna data
    idx = table_bds.PRN == satID;
    satData = table_bds(idx, :);
    timeData = satData.TimeStamp;
    snrData = satData.SNR;
    % extract EL and AZ data for the current satellite
    elData = El_bds(i, :);
    azData = Az_bds(i, :);
    % find the valid EL data points
    validIdx = nan(length(timeData), 1);
    for j = 1:length(timeData)
        t_idx = find(time == timeData(j));
        if ~isempty(t_idx)
            validIdx(j) = t_idx;
        end
    end
    elData = elData(validIdx);
    azData = azData(validIdx);
    % find the corresponding tide height
    validIdx_tide = nan(length(timeData), 1);
    for j = 1:length(timeData)
        t_idx = find(tideTime == timeData(j));
        if ~isempty(t_idx)
            validIdx_tide(j) = t_idx;
        end
    end
    current_tideHeight = tideHeight(validIdx_tide);
    % find the specific time points during the specific elevation angle range
    mask = (elData >= El_range(1)) & (elData <= El_range(2)) & ...
           ((azData >= mod(Az_range(1), 360)) | (azData <= Az_range(2)));
    if isempty(find(mask, 1))
        continue;
    end
    timeData = timeData(mask);
    t_last = seconds(timeData(end) - timeData(1));
    if t_last < 1800  % 0.5 hour
        % skip if the data duration is less than 0.5 hour
        continue;
    end
    snrData = snrData(mask);
    elData = elData(mask);
    current_tideHeight = current_tideHeight(mask);
    % convert dB-Hz to linear scale
    linSNR = 10.^(snrData / 20);

    % polyfit to remove the trend of elevation angle
    secondArray = seconds(timeData - timeData(1));
    p = polyfit(secondArray', linSNR, 2);
    trend = polyval(p, secondArray');
    % detrend the direct signal component (2nd order polynomial)
    detrendedLinSNR = linSNR - trend';
    % compute sine of elevation angle
    sinE = sind(elData);

    %% calculate the integer Ambiguity
    light_speed = 299792458; % m/s
    lambda = light_speed / 1561.098e6; % BDS B1I wavelength

    % %% Lomb-Scargle periodogram
    % % define the searching window of height
    % h_min = -20;
    % h_max = 0;
    % % convert height to frequency range
    % f_min = 2 * h_min / lambda;
    % f_max = 2 * h_max / lambda;
    % % frequency resolution
    % delta_h = 0.01;
    % h_vec = h_min : delta_h : h_max;
    % f_vec = 2 * h_vec / lambda;
    % 
    % % Set a window to calculate Lomb-Scargle periodogram
    % windowSize = 480; % number of data points in each window
    % stepSize = 50; % step size for moving the window
    % numPoints = length(detrendedLinSNR);
    % estHeight_all = [];
    % estTime_all = [];
    % est_PRN = [];
    % for startIdx = 1 : stepSize : (numPoints - windowSize + 1)
    %     endIdx = startIdx + windowSize - 1;
    %     meanIdx = round((startIdx + endIdx) / 2);
    %     currentTime = timeData(meanIdx);
    %     sinE_window = sinE(startIdx:endIdx);
    %     detrendedLinSNR_window = detrendedLinSNR(startIdx:endIdx);
    % 
    %     % Lomb-Scargle periodogram for the current window
    %     [Pxx_win, F_win] = lomb_scargle_custom(sinE_window, detrendedLinSNR_window, f_vec);
    % 
    %     % find the peak frequency in the current window
    %     [peakPower_win, peakIdx_win] = max(Pxx_win);
    %     peakFreq_win = F_win(peakIdx_win);
    % 
    %     estHeight_win = (peakFreq_win * lambda) / 2;
    %     estHeight_win = estHeight_win + site_pos_lla(3); % add antenna height offset
    %     estTime_all = [estTime_all; currentTime];
    %     estHeight_all = [estHeight_all; estHeight_win];
    % end
    % % save the result into the structure
    % results_gps(i).PRN = satID;
    % results_gps(i).EstimatedTime = estTime_all;
    % results_gps(i).EstimatedHeight = estHeight_all;
    % results_gps(i).SNR = snrData;
    % results_gps(i).TimeData = timeData;
    % results_gps(i).TideHeight = current_tideHeight; % transfer from HK mcd to wgs84
    % results_gps(i).ElData = elData;

    [N, frac_phase, ~] = solve_ambiguity(current_tideHeight, sinE.', site_pos_lla(3), lambda);
    % Output estimated height based on integer ambiguity
    plotId = [plotId satID];
    H1(length(plotId)) = N2Height(timeData, detrendedLinSNR, sinE, N + frac_phase, lambda, site_pos_lla);
    
end

% plotSNR(results_gps);
plotHeight(H1, plotId, tide);