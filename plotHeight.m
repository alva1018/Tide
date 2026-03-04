function plotHeight(results, plotId, tide)
% plotHeight - Plot SNR peaks/troughs and height estimates

    hFig = figure('Color', 'w', 'Name', 'Height Figure');
    % create axes
    ax = axes('Position', [0.1, 0.2, 0.8, 0.7]); 
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    subplot(2,1,1);
    yyaxis left;
    xlabel('Time');
    ylabel('Detrended SNR');

    yyaxis right;
    ylabel('Height (m)');

    subplot(2,1,2);
    ylabel('N value');
    xlabel('Time');

    % create popup menu
    uicontrol('Parent', hFig, ...
        'Style', 'popupmenu', ...
        'String', string(plotId), ...       
        'Position', [50, 20, 100, 10], ... 
        'Callback', @updatePlot);    % call the nested function below
    
    % text label
    uicontrol('Parent', hFig, ...
        'Style', 'text', ...
        'String', 'Select Satellite:', ...
        'Position', [50, 35, 100, 20], ...
        'FontSize',10, ...
        'BackgroundColor', 'w', ...
        'HorizontalAlignment', 'left');
    
    % initial plot
    if ~isempty(results)
        firstSat = plotId(1);
        drawSatellite(firstSat);
    end

    %% inner functions
    function updatePlot(src, ~)
        % get the selected satellite ID
        val = src.Value;    
        str = src.String;
        selectedSat = str(val,:); 
        selectedID = str2double(selectedSat);
        drawSatellite(selectedID);
    end

    function drawSatellite(satName)
        % clear axes
        subplot(2,1,1);
        yyaxis('left');
        cla;
        yyaxis('right');
        cla;
        for i = 1:length(plotId)
            if plotId(i) == satName
                satData = results(i);
                break;
            end
            if i == length(results)
                satData = [];
            end
        end
        
        % redraw plot
        if ~isempty(satData)
            yyaxis('left');
            plot(satData.time, satData.snr_detrended, ...
                'Color', [0.8 0.8 0.8], 'DisplayName', 'Detrended SNR'); hold on;
            plot(satData.time, satData.snr_smooth, ...
                'k-', 'LineWidth', 1.5, 'DisplayName', 'Smoothed SNR');
            % peak plot
            peak_mask = satData.types_sorted == 1;
            if any(peak_mask)
                scatter(satData.out_time(peak_mask), satData.snr_smooth(satData.events_sorted(peak_mask)), ...
                    50, 'r^', 'filled', 'DisplayName', 'Peaks');
            end
            % trough plot
            trough_mask = satData.types_sorted == 0;
            if any(trough_mask)
                scatter(satData.out_time(trough_mask), satData.snr_smooth(satData.events_sorted(trough_mask)), ...
                    50, 'cv', 'filled', 'DisplayName', 'Troughs');
            end

            yyaxis('right');
            % % filter the height at trough points
            % idx = satData.types_sorted == 0;
            % plot_time = satData.out_time(idx);
            % plot_h = satData.out_h(idx);
            % plot(plot_time, plot_h, ...
            %     'bo-', 'LineWidth', 1.5, 'DisplayName', 'Estimated Height');
            plot(satData.out_time, satData.out_h, ...
                'bo-', 'LineWidth', 1.5, 'DisplayName', 'Estimated Height');
            hold on;
            % filter the tide data to match time range
            tide_mask = tide.tideTime >= min(satData.time) & tide.tideTime <= max(satData.time);
            plot(tide.tideTime(tide_mask), tide.tideHeight(tide_mask), ...
                'g-', 'LineWidth', 1.5, 'DisplayName', 'Tide Polynomial Height');
            tide_mask = tide.tideTime_real >= min(satData.time) & tide.tideTime_real <= max(satData.time);
            plot(tide.tideTime_real(tide_mask), tide.tideHeight_real(tide_mask), ...
                'm--', 'LineWidth', 1.5, 'DisplayName', 'Tide Real Height');
            % plot prior height
            plot(satData.time, satData.h_prior, ...
                'c-.', 'LineWidth', 1.5, 'DisplayName', 'Prior Height');
            legend('show', 'Location', 'southwest');
            title(['Satellite: ' string(satName)]);
        else
            title(['No Data for ' string(satName)]);
        end

        % subplot for N value
        subplot(2,1,2);
        yyaxis('left');
        cla;
        plot(satData.time, satData.N_origin, 'b-', 'LineWidth', 1.5, 'DisplayName', 'N Value Derived by Gauge'); hold on;
        plot(satData.out_time, satData.N_integer, 'ro-', 'LineWidth', 1.5, 'DisplayName', 'Integer N from Peaks/Troughs');
        plot(satData.out_time, satData.N_current, 'kx', 'LineWidth', 1.5, 'DisplayName', 'Corresponding N to Peaks/Troughs');
        % find the nan values in N_origin and plot them as magenta points
        nan_mask = ~isnan(satData.N_origin_int);
        N_origin_int_time = satData.N_origin_int_time(nan_mask);
        N_origin_int = satData.N_origin_int(nan_mask);
        plot(N_origin_int_time, N_origin_int, 'ms', 'MarkerFaceColor', 'm', 'DisplayName', 'Rounded N in Original N Array');
        xlabel('Time');
        ylabel('N Value');
        grid on;

        yyaxis('right');
        cla;
        plot(satData.out_time, satData.sinE, 'g', 'LineWidth', 1.5, 'DisplayName', 'sin(Elevation)');
        ylabel('sin(Elevation)');
        legend('show', 'Location', 'northeast');
    end
end