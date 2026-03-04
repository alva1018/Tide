function results = N2Height(time, snr_detrended, sinE, N, lambda, site_pos_lla)
% inputs:
%   time:           time sequence (datetime array)
%   snr_detrended:  snr after detrending
%   sinE:           sine of satellite elevation angle
%   N:              integer ambiguity + fractional phase (array)
%   lambda:        wavelength of the GNSS signal (m)
%  site_pos_lla:  [lat, lon, alt] of the receiver site (degrees, degrees, meters)
% outputs:
%   precise_heights: a Table containing [Time, Elevation, Type (Peak/Trough), Integer N, Precise Height]

    %% Data smoothing and Peak/Trough Detection
    results = struct('time', [], 'snr_detrended', [], ...
                     'snr_smooth', [], 'types_sorted', [], 'out_time', [],...
                     'out_h', [], 'out_N', [], 'events_sorted', [], ...
                     'N_origin', [],'N_origin_int', [], 'N_origin_int_time', [], 'N_integer', [], 'N_current', [], ...
                     'h_prior', [], 'sinE', []);
    % sampling frequency
    sampling_freq = 1; % Hz
    snr_smooth = sgolayfilt(snr_detrended, 2, 15);
    % zero-phase lowpass filter
    fc = 0.03;
    [b, a] = butter(2, fc / (sampling_freq / 2), 'low');
    snr_smooth = filtfilt(b, a, snr_smooth);
    
    min_dist_samples = round(22 / sampling_freq);
    amp_range = max(snr_smooth) - min(snr_smooth);
    min_prom = amp_range * 0.25;

    [~, locs_p] = findpeaks(snr_smooth, ...
        'MinPeakDistance', min_dist_samples, ...
        'MinPeakProminence', min_prom);
    [~, locs_t] = findpeaks(-snr_smooth, ...
        'MinPeakDistance', min_dist_samples, ...
        'MinPeakProminence', min_prom);

    events = [locs_p; locs_t];
    
    events(1:end-1)  = events(1:end-1) + 60;

    types = [ones(size(locs_p)); zeros(size(locs_t))];
    % sort peaks and troughs by time, 1 - peak, 0 - trough
    [events_sorted, sort_idx] = sort(events);
    types_sorted = types(sort_idx);

    % output result initialization
    results_count = length(events_sorted); % the number of detected events
    out_time = time(events_sorted); % time of events
    out_type = strings(results_count, 1); % "Peak" or "Trough"
    out_h = zeros(results_count, 1);
    out_N_value = zeros(results_count, 1);
    out_current_N = zeros(results_count, 1);
    out_sinE = zeros(results_count, 1);
    for i = 1 : results_count
        % current event point
        idx = events_sorted(i);
        curr_sinE = sinE(idx);
        current_N = N(idx);
        is_peak = types_sorted(i);

        if is_peak
            N_val = round(current_N);
            h_precise = (N_val * lambda) / (2 * curr_sinE);
            out_type(i) = "Peak (0)";
        else
            N_val = round(current_N - 0.5) + 0.5;
            h_precise = (N_val * lambda) / (2 * curr_sinE);
            out_type(i) = "Trough (π)";
        end
        out_current_N(i) = current_N;
        out_N_value(i) = N_val;
        out_sinE(i) = curr_sinE;
        % % limit exceedance check (half wavelength limit)
        % limit = lambda / (4 * curr_sinE);
        % % cycle slip warning
        % if abs(h_precise - out_h_prior(i)) > limit
        %     h_precise = NaN; 
        %     out_type(i) = out_type(i) + " [WARN]";
        % end
        out_h(i) = -h_precise + site_pos_lla(3); % convert to height above mean sea level (WGS84)
    end

    % time: original time array (1 Hz)
    results.time = time;
    results.snr_detrended = snr_detrended;
    results.snr_smooth = snr_smooth;
    % time: time of detected eventstime
    results.out_time = out_time;
    results.types_sorted = types_sorted;
    results.out_h = out_h;
    results.out_N = out_N_value;
    results.events_sorted = events_sorted;
    results.N_origin = N;
    results.N_current = out_current_N;
    results.N_integer = out_N_value;
    results.sinE = out_sinE;
    % calculate prior height
    h_prior = (N .* lambda) ./ (2 * sinE.');
    results.h_prior = -h_prior + site_pos_lla(3);

    % find the integer N / half-integer N in the original N array
    N_origin_int = nan(size(N));
    % datetime array
    N_origin_int_time = datetime(2024,1,22,0,0,0) + seconds(0:length(N)-1);
    rounded = round(N*2) / 2; % round to nearest 0.5
    rounded = unique(rounded); % get unique values to reduce computation
    for i = 1:length(rounded)
        target_num = rounded(i);
        [~, idx_rounded] = min(abs(N - target_num)); % find index of closest value in N
        N_origin_int(i) = N(idx_rounded);
        N_origin_int_time(i) = time(idx_rounded);
    end
    results.N_origin_int = N_origin_int;
    results.N_origin_int_time = N_origin_int_time;
end