function [az_result, el_result] = calculate_bds_az_el(time, bds_id, site_pos_lla)
    % file and site parameters
    rinex_file = '2026022_01D_MN.rnx';
    
    % Observation time in UTC (Year, Month, Day, Hour, Minute, Second)
    start_time = time(1);
    end_time = time(end);
    total_secs = seconds(end_time - start_time);
    % initialize result arrays
    num_sat = length(bds_id);
    az_result = nan(num_sat, total_secs + 1);
    el_result = nan(num_sat, total_secs + 1);
    % read rinex ephemeris data
    ephemeris = rinexread(rinex_file);
    % for BeiDou specific ephemeris
    eph = ephemeris.BeiDou;
    for i = 1 : num_sat
        prn = bds_id(i);
        
        current_eph = eph(eph.SatelliteID == prn,:);
        % iterate all seconds
        for sec = 0 : total_secs
            if mod(sec, 600) == 0
                disp(['Processing BDS C' num2str(prn) ', Time Elapsed: ' num2str(sec) ' seconds']);
            end
            current_time = start_time + seconds(sec);
            time_vec = current_eph.Time;
            % find the closest ephemeris for this time
            time_diffs = abs(time_vec - current_time);
            [min_diff, min_idx] = min(time_diffs);
            if min_diff > seconds(7200) % 2 hours threshold
                continue; % no valid ephemeris
            end
            selected_eph = current_eph(min_idx, :);
            
            % convert current time to GPS week and seconds
            [~, gps_seconds] = utc2gps(current_time);

            % BDT lags GPS time by 14 seconds (as of 2006)
            bdt_seconds = gps_seconds - 14; 
            
            % calculate satellite position in ECEF
            [sat_pos_ecef, ~] = sat_position(selected_eph, bdt_seconds, prn);
            
            % calculate azimuth and elevation
            [az, el, ~] = calculate_az_el(site_pos_lla, sat_pos_ecef);
            if el < 0
                continue; % satellite below horizon
            end
            az_result(i, sec + 1) = az;
            el_result(i, sec + 1) = el;
        end
    end
end

% calculate satellite position in ECEF coordinates
function [pos, clk_err] = sat_position(eph, t, prn)
    % GPS constants (Using Standard Pi)
    gpsPi          = 3.1415926535898; 
    
    %--- Constants for satellite position calculation (CGCS2000) --------------
    Omegae_dot     = 7.2921150e-5;     % Earth rotation rate for BDS, [rad/s]
    GM             = 3.986004418e14;   % GM for BeiDou (same as Galileo), [m^3/s^2]
    F              = -4.442807633e-10; % Constant, [sec/(meter)^(1/2)]
    
    % dt: time from ephemeris clock reference epoch
    Toc = eph.Time;
    [~, Toc_gps_seconds] = utc2gps(Toc);
    % Ensure Ephemeris reference time is also in BDT
    Toc_bdt_seconds = Toc_gps_seconds - 14; 
    
    dt = t - Toc_bdt_seconds;
    if dt > 302400, dt = dt - 604800; end
    if dt < -302400, dt = dt + 604800; end
    
    % Use TGD1 (for B1I) or TGD2 (for B2I). Usually TGD1 is default.
    % rinexread usually provides TGD1/TGD2. Assuming TGD1 exists or using TGD if mapped.
    if ismember('TGD1', eph.Properties.VariableNames)
        tgd_val = eph.TGD1;
    else
        tgd_val = eph.TGD; % Fallback
    end
    clk_err = eph.SVClockBias + eph.SVClockDrift * dt + eph.SVClockDriftRate * dt^2 - tgd_val;
    
    % corrected time
    t = t - clk_err;
    % calculate time from ephemeris reference epoch
    tk = t - eph.Toe;
    if tk > 302400, tk = tk - 604800; end
    if tk < -302400, tk = tk + 604800; end
    
    % A: semi-major axis
    A = eph.sqrtA ^ 2;
    % n0: computed mean motion
    n0 = sqrt(GM / A^3);
    % n: corrected mean motion
    n = n0 + eph.Delta_n;
    % Mean anomaly M
    M = eph.M0 + n * tk;
    M = rem(M + 2*gpsPi, 2*gpsPi);
    
    % E: Eccentric anomaly (Iterative solution)
    E = M;
    for k = 1:10
        E_old = E;
        E = M + eph.Eccentricity * sin(E);
        dE      = rem(E - E_old, 2*gpsPi);
        if abs(dE) < 1.e-12
            break;
        end
    end
    E = rem(E + 2*gpsPi, 2*gpsPi);
    
    % v: True anomaly
    v = atan2(sqrt(1 - eph.Eccentricity^2) * sin(E), cos(E)-eph.Eccentricity);
    
    % phi: Argument of latitude
    phi = v + eph.omega;
    phi = rem(phi, 2*gpsPi);
    % Corrections
    du = eph.Cus * sin(2*phi) + eph.Cuc * cos(2*phi);
    dr = eph.Crs * sin(2*phi) + eph.Crc * cos(2*phi);
    di = eph.Cis * sin(2*phi) + eph.Cic * cos(2*phi);
    
    u = phi + du;
    r = A * (1 - eph.Eccentricity * cos(E)) + dr;
    i = eph.i0 + eph.IDOT * tk + di;
    
    % the satellite position in orbital plane
    x_prime = r * cos(u);
    y_prime = r * sin(u);
    
    % Omega: Longitude of ascending node
    Omega = eph.OMEGA0 + (eph.OMEGA_DOT - Omegae_dot) * tk - Omegae_dot * eph.Toe;
    Omega = rem(Omega + 2*gpsPi, 2*gpsPi);

    % GEO Satellite Handling (PRN 1-5 and others defined as GEO)
    % BeiDou GEOs use a different coordinate transformation (Rotation around X axis)
    is_geo = (prn >= 1 && prn <= 5) || prn > 58; % Simple heuristic for GEOs
    
    if is_geo
        % --- GEO Algorithm (BDS-SIS-ICD) ---
        % 1. Calculate P_GEO in inertial frame
        Omega_geo = eph.OMEGA0 + eph.OMEGA_DOT * tk - Omegae_dot * eph.Toe;
        xg = x_prime * cos(Omega_geo) - y_prime * cos(i) * sin(Omega_geo);
        yg = x_prime * sin(Omega_geo) + y_prime * cos(i) * cos(Omega_geo);
        zg = y_prime * sin(i);
        
        % 2. Rotation Matrix Rx(-5 deg) and Rz(We*tk)
        % Coordinate rotation angle of the CGCS2000 coordinate system
        sin_5 = -0.0871557427476582; % sin(-5 deg)
        cos_5 = 0.9961946980917456;  % cos(-5 deg)
        wt = Omegae_dot * tk;
        
        % Apply rotation (Equation from BDS ICD 5.2.4)
        x = xg * cos(wt) + yg * sin(wt) * cos_5 + zg * sin(wt) * sin_5;
        y = -xg * sin(wt) + yg * cos(wt) * cos_5 + zg * cos(wt) * sin_5;
        z = -yg * sin_5 + zg * cos_5;
        
    else
        % --- MEO / IGSO Algorithm (Standard, Same as GPS) ---
        x = x_prime * cos(Omega) - y_prime * cos(i) * sin(Omega);
        y = x_prime * sin(Omega) + y_prime * cos(i) * cos(Omega);
        z = y_prime * sin(i);
    end
    
    pos = [x, y, z];
    
    % Relativistic correction
    dtr = F * eph.Eccentricity * eph.sqrtA * sin(E);
    
    % Satellite clock error
    % Re-apply clock calc with dtr and TGD1
    clk_err = eph.SVClockBias + eph.SVClockDrift * dt + eph.SVClockDriftRate * dt^2 + dtr - tgd_val;
end

% 4. Coordinate transformation (Kept same as original)
function [az, el, d] = calculate_az_el(site_lla, sat_ecef)
    % ... (No changes needed for coordinate frame conversion) ...
    % site_lla: [lat, lon, alt] (deg, deg, m)
    % sat_ecef: [x, y, z] (m)
    
    % 1. Station to ECEF
    phi = deg2rad(site_lla(1));
    lam = deg2rad(site_lla(2));
    h = site_lla(3);
    
    a = 6378137.0; % WGS84 semi-major axis (CGCS2000 is practically identical)
    f = 1/298.257223563;
    e2 = 2*f - f^2;
    
    N = a / sqrt(1 - e2 * sin(phi)^2);
    x_stn = (N + h) * cos(phi) * cos(lam);
    y_stn = (N + h) * cos(phi) * sin(lam);
    z_stn = (N * (1 - e2) + h) * sin(phi);
    stn_ecef = [x_stn, y_stn, z_stn];
    
    % Vector from station to satellite
    dx = sat_ecef - stn_ecef;
    d = norm(dx);
    
    % Rotation matrix from ECEF to ENU
    R = [-sin(lam),           cos(lam),          0;
         -sin(phi)*cos(lam), -sin(phi)*sin(lam), cos(phi);
          cos(phi)*cos(lam),  cos(phi)*sin(lam), sin(phi)];
      
    enu = R * dx';
    E = enu(1); N = enu(2); U = enu(3);
    
    % Calculate azimuth and elevation
    az = rad2deg(atan2(E, N));
    if az < 0, az = az + 360; end
    el = rad2deg(asin(U / d));
end