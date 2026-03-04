function hpol = plotSky(varargin)
%Function plots "sky view" from the receiver perspective.
%
%h = skyPlot(AZ, EL, PRN, line_style, mask_az, mask_el)
%
%   Inputs:
%       AZ              - contains satellite azimuth angles. It is a 2D
%                       matrix. One line contains data of one satellite.
%                       The columns are the calculated azimuth values.
%       EL              - contains satellite elevation angles. It is a 2D
%                       matrix. One line contains data of one satellite.
%                       The columns are the calculated elevations.
%       PRN             - a row vector containing PRN numbers of the
%                       satellites.
%       line_style      - line style of the plot. The same style will be
%                       used to plot all satellite positions (including
%                       color).
%       mask_az         - azimuth angles for masking certain areas (optional)
%       mask_el         - elevation angles for masking certain areas (optional)
%
%   Outputs:
%       h               - handle to the plot

%% Check arguments and sort them ==========================================
[hAxis, args, nargs] = axescheck(varargin{:});

if nargs < 3 || nargs > 6 || nargs == 5
    error('Requires 3 to 6 data arguments.')
elseif nargs == 3
    [az, el, prn]   = deal(args{1:3});
    line_style      = 'auto';
elseif nargs == 4
    [az, el, prn, line_style] = deal(args{1:4});
elseif nargs == 6
    [az, el, prn, line_style, mask_az, mask_el] = deal(args{1:6});
    %--- Apply masking ----------------------------------------------------
end

if ischar(az) || ischar(el) || ischar(prn)
    error('AZ and EL must be numeric.');
end

if ~isequal(size(az), size(el))
    error('AZ and EL must be same size.');
end

%% Prepare axis ===========================================================
hAxis = newplot(hAxis);

%--- Get x-axis text color so grid is in same color -----------------------
tc = get(hAxis, 'xcolor');

hold(hAxis, 'on');

%--- Plot white background ------------------------------------------------
rectangle('position', [-90, -90, 180, 180], ...
    'Curvature', [1 1], ...
    'facecolor', 'white', ...
    'edgecolor', 'black');

%% Mask areas (if specified) ===============================================
if nargs == 6 || nargs == 7
    %--- Convert mask elevation to spherical elevation ---------------------
    % note: mask_el and mask_az are the two-element vectors defining the
    % area to be masked
    % convert azimuth angles to XYZ coordinates angles
    mask_az = 90 - mask_az;
    % convert elevation angles to spherical elevation angles
    mask_el = abs(mask_el - 90);
    % create a closed sector to cover the masked area
    theta = linspace(mask_az(1)*pi/180, mask_az(2)*pi/180, 100);
    r_inner = mask_el(1);
    r_outer = mask_el(2);
    x_inner = r_inner * cos(theta);
    y_inner = r_inner * sin(theta);
    x_outer = r_outer * cos(fliplr(theta));
    y_outer = r_outer * sin(fliplr(theta));
    x_mask = [x_inner, x_outer];
    y_mask = [y_inner, y_outer];
    %--- Plot the masked area (filled with white color) ----------------------
    azel = fill(hAxis, x_mask, y_mask, 'blue', 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
    set(azel, 'FaceAlpha', 0.2); % set transparency

    % plot again (the same az, el ranges from 5 to 30 degrees)
    hold(hAxis, 'on');
    mask_el = [5 30];
    mask_el = abs(mask_el - 90);
    r_inner = mask_el(1);
    r_outer = mask_el(2);
    x_inner = r_inner * cos(theta);
    y_inner = r_inner * sin(theta);
    x_outer = r_outer * cos(fliplr(theta));
    y_outer = r_outer * sin(fliplr(theta));
    x_mask = [x_inner, x_outer];
    y_mask = [y_inner, y_outer];
    %--- Plot the masked area (filled with white color) ----------------------
    azel2 = fill(hAxis, x_mask, y_mask, 'blue', 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
    set(azel2, 'FaceAlpha', 0.2); % set transparency
end

%% Plot spokes ============================================================
%--- Find spoke angles ----------------------------------------------------
% Only 6 lines are needed to divide circle into 12 parts
th = (1:6) * 2*pi / 12;

%--- Convert spoke end point coordinate to Cartesian system ---------------
cst = cos(th); snt = sin(th);
cs = [cst; -cst];
sn = [snt; -snt];

%--- Plot the spoke lines -------------------------------------------------
line(90*sn, 90*cs, 'linestyle', '-', 'color', '#AFC7BF', 'linewidth', 0.3, ...
    'handlevisibility', 'on');

%% Annotate spokes in degrees =============================================
rt = 1.1 * 90;
for i = 1 : length(th)
    if i == length(th)
        % Write "North" at the top of the plot
        text(rt*snt(i), rt*cst(i), "S", ...
            'horizontalalignment', 'center', 'handlevisibility', 'off');
        % Write "South" at the bottom of the plot
        text(-rt*snt(i), -rt*cst(i), "N", ...
            'handlevisibility', 'off', 'horizontalalignment', 'center');
    elseif i == 3
        % Write "East" at the right side of the plot
        text(rt*snt(i), rt*cst(i), "E", ...
            'horizontalalignment', 'center', 'handlevisibility', 'off');
        % Write "West" at the left side of the plot
        text(-rt*snt(i), -rt*cst(i), "W", ...
            'handlevisibility', 'off', 'horizontalalignment', 'center');
    else
        %--- Write text in the first half of the plot -------------------------
        text(rt*snt(i), rt*cst(i), int2str(i*30), ...
            'horizontalalignment', 'center', 'handlevisibility', 'off');

        %--- Write text in the opposite half of the plot ----------------------
        text(-rt*snt(i), -rt*cst(i), int2str(180 + i*30), ...
            'handlevisibility', 'off', 'horizontalalignment', 'center');
    end
end

%% Plot elevation grid ====================================================
%--- Define a "unit" radius circle ----------------------------------------
th = 0 : pi/50 : 2*pi;
xunit = cos(th);
yunit = sin(th);

%--- Plot elevation grid lines and tick text ------------------------------
for elevation = 0 : 20 : 80
    %elevationSpherical = 90*cos((pi/180) * elevation);
    elevationSpherical = abs(elevation-90);  %sara p fix

    line(yunit * elevationSpherical, xunit * elevationSpherical, ...
        'linestyle', '-', 'color', '#AFC7BF', 'linewidth', 0.3, ...
        'handlevisibility', 'off');

    text(-elevationSpherical, 0, num2str(elevation) + "°", ...
        'BackgroundColor', 'none', 'horizontalalignment','center', ...
        'handlevisibility', 'off');
end

%--- Set view to 2-D ------------------------------------------------------
view(0, 90);

%--- Set axis limits ------------------------------------------------------
%save some space for the title
axis([-95 95 -90 101]);

%% Transform elevation angle to a distance to the center of the plot ------
%elSpherical = 90*cos(el * pi/180);
elSpherical = abs(el-90);  %sara p fix

%--- Transform data to Cartesian coordinates ------------------------------
yy = elSpherical .* cos(az * pi/180);
xx = elSpherical .* sin(az * pi/180);

%% Plot data on top of the grid ===========================================
if strcmp(line_style, 'auto')
    %--- Plot with "default" line style -----------------------------------
    hpol = plot(hAxis, xx', yy', '.-');
else
    %--- Plot with user specified line style ------------------------------
    % The same line style and color will be used for all satellites
    hpol = plot(hAxis, xx', yy', line_style);
end

%--- Place satellite PRN numbers at the latest position -------------------
for i = 1:length(prn)
    if(prn(i) ~= 0)
        % The empthy space is used to place the text a side of the last
        % point. This solution results in constant offset even if a zoom
        % is used.
        % find the last non-NaN point
        last_idx = find(~isnan(xx(i, :)), 1, 'last');
        text(xx(i, last_idx), yy(i, last_idx), ['  ', int2str(prn(i))], 'color', 'b');

        %--- Mark the last position of the satellite ------------------------------
        plot(hAxis, xx(i, last_idx), yy(i, last_idx), 'o', 'MarkerFaceColor', 'blue', ...
            'MarkerEdgeColor', 'k', 'MarkerSize', 8, ...
            'handlevisibility', 'off');
    end
end

%--- Make sure both axis have the same data aspect ratio ------------------
axis(hAxis, 'equal');

%--- Switch off the standard Cartesian axis -------------------------------
axis(hAxis, 'off');
