function [time_fine, h_fitted, time_current_segment, height_current_segment] = tideFitting()
    % read raw tide data from text file
    filename = 'raw_forTidePred.txt'; 
    opts = detectImportOptions(filename);
    opts.VariableNames = {'time', 'height'};
    opts.VariableTypes = {'string', 'double'};
    T = readtable(filename, opts);
    T.TimeStamp = datetime(T.time, 'InputFormat', 'yyyy-MM-dd HH:mm');
    T.TimeStamp = T.TimeStamp - hours(8); % convert to UTC+0
    time = T.TimeStamp;
    h = T.height;

    % read the current tide data
    filename_current = "./raw_20260122_233836.txt";
    T_current = readtable(filename_current, opts);
    T_current.TimeStamp = datetime(T_current.time, 'InputFormat', 'yyyy-MM-dd HH:mm');
    T_current.TimeStamp = T_current.TimeStamp - hours(8); % convert to UTC+0
    time_current = T_current.TimeStamp;
    height_current = T_current.height;

    % main tidal constituents to consider
    M2 = 12.4206012; 
    S2 = 12.0000000; 
    K1 = 23.9344721; 
    O1 = 25.8193387; 
    N2 = 12.6583482;
    M4 = 6.2103006;
    Q1 = 26.868356;
    
    % transform periods to angular frequencies (rad/hour)
    constituents = struct();
    constituents.M2 = 2 * pi / M2;
    constituents.S2 = 2 * pi / S2;
    constituents.K1 = 2 * pi / K1;
    constituents.O1 = 2 * pi / O1;
    constituents.N2 = 2 * pi / N2;
    constituents.M4 = 2 * pi / M4;
    constituents.Q1 = 2 * pi / Q1;

    constituents.MS4 = 2 * pi / 6.1033393;
    constituents.M6 = 2 * pi / 4.1402004;
    constituents.MK3 = 2 * pi / 8.1771302;
    constituents.MN4 = 2 * pi / 6.2691739;
    
    names = fieldnames(constituents);
    num_cons = length(names);

    % time array in hours since start
    base_time = time(1);
    t_hours = hours(time - base_time);
    
    % construct design matrix A
    A = ones(length(t_hours), 1);
    for i = 1 : num_cons
        name = names{i};
        omega = constituents.(name);
        % add cos and sin pairs
        A = [A, cos(omega * t_hours), sin(omega * t_hours)];
    end

    % --- least squares solution ---
    % Check condition number to prevent unstable solutions due to short data duration
    if cond(A) > 1e4
        warning('Data duration too short, matrix condition number too large.');
    end

    % Ridge Regression for better stability
    % coeffs = A \ h;
    % coeffs = (A'A + lambda*I)^-1 * A'h
    lambda = 0.6; 
    XTX = A' * A;
    I = eye(size(XTX));
    I(1,1) = 0; % do not regularize the intercept term
    coeffs = (XTX + lambda * I) \ (A' * h);
    Z0 = coeffs(1);
    % 
    % coeffs = A \ h;
    % Z0 = coeffs(1);

    % --- output analysis results ---
    idx = 2;
    for i = 1:num_cons
        beta_c = coeffs(idx);
        beta_s = coeffs(idx + 1);
        amp = sqrt(beta_c^2 + beta_s^2);
        phase = atan2(beta_s, beta_c); % in radians
        phase_deg = rad2deg(phase);
        idx = idx + 2;
    end
    % reconstruct fitted tide over original time points
    endtime_fine = time(end) + hours(8);
    time_fine = time(1) : seconds(1) : endtime_fine;
    t_fine_hours = hours(time_fine - base_time);
    t_fine_hours = t_fine_hours.';
    % reconstruct fitted tide
    A_pred = ones(length(t_fine_hours), 1);
    for i = 1:num_cons
        name = names{i};
        omega = constituents.(name);
        A_pred = [A_pred, cos(omega * t_fine_hours), sin(omega * t_fine_hours)];
    end
    h_fitted = A_pred * coeffs;
    
    % extract the current tide fitting segment
    starttime_current = time(end) + seconds(1);
    endtime_current = time_current(end);
    % find indices in the current tide time range
    idx_start = find(time_current >= starttime_current, 1, 'first');
    idx_end = find(time_current <= endtime_current, 1, 'last');
    time_current_segment = time_current(idx_start:idx_end);
    height_current_segment = height_current(idx_start:idx_end);
    
    % calculate residual RMS
    idx = nan(length(time), 1);
    for j = 1:length(time)
        t_idx = find(time_fine == time(j));
        if ~isempty(t_idx)
            idx(j) = t_idx;
        end
    end
    rms_error = sqrt(mean((h - h_fitted(idx)).^2));
    figure;
    subplot(2,1,1);
    plot(time, h, 'b-', 'DisplayName', 'Raw Tide', 'LineWidth', 1.5); hold on;
    plot(time_fine, h_fitted, 'r--', 'DisplayName', 'Fitted', 'LineWidth', 1.5);
    xlim([time_fine(1), time_fine(end)]);
    ylabel('Height (m)');
    title(['Tidal Fit (RMS Error: ' num2str(rms_error, '%.3f') ' m)']);
    grid on;
    % plot current tide data segment
    plot(time_current_segment, height_current_segment, 'g-', 'DisplayName', 'Current Tide Data', 'LineWidth', 1.5);
    legend('show');

    subplot(2,1,2);
    plot(time, h - h_fitted(idx), 'k-', 'LineWidth', 1);
    xlim([time_fine(1), time_fine(end)]);
    yline(0, 'r--');
    ylabel('Residual (m)');
    xlabel('Time');
    title('Residuals (Raw - Fitted)');
    grid on;
end