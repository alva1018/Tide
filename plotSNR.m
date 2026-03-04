function plotSNR (results)
% plotSNR - Plot data for multiple satellites with interactive GUI
    %% create GUI for interactive selection
    hFig = figure('Color', 'w', 'Name', 'Figure');
    % create axes
    ax = axes('Position', [0.1, 0.2, 0.85, 0.7]); 
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

    yyaxis left;
    xlabel(ax, 'Time');
    ylabel(ax, 'C/N0 (dB-Hz)');
    ylim(ax, [10, 55]);

    yyaxis right;
    ylabel(ax, 'Height (m)');

    % create popup menu
    uicontrol('Parent', hFig, ...
        'Style', 'popupmenu', ...
        'String', [results.PRN], ...       
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
        firstSat = [results.PRN];
        firstSat = firstSat(1);
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
        yyaxis(ax, 'left');
        cla(ax);
        yyaxis(ax, 'right');
        cla(ax);
        for i = 1:length(results)
            if results(i).PRN == satName
                satData = results(i);
                break;
            end
            if i == length(results)
                satData = [];
            end
        end
        
        % redraw plot
        if ~isempty(satData)
            yyaxis(ax, 'left');
            plot(satData.TimeData, satData.SNR, ...
                'DisplayName', "SNR", ...
                'LineWidth', 1.5, ...
                'Color', '#0072BD');
            yyaxis(ax, 'right');
            plot(satData.EstimatedTime, satData.EstimatedHeight, ...
                'DisplayName', "Estimated Height", ...
                'LineWidth', 1.5, ...
                'Color', '#D95319');
            hold(ax, 'on');
            plot(ax, satData.TimeData, satData.TideHeight, ...
                'DisplayName', "Tide Height", ...
                'LineWidth', 1.5, ...
                'Color', '#77AC30');
            legend(ax, 'show', 'Location', 'northeast');
            title(ax, ['Satellite: ' string(satName)]);
            % startTime = datetime(2026,1,22,3,30,0);
            % endTime = datetime(2026,1,22,10,30,0);
            % xlim(ax, [startTime, endTime]);
        else
            title(ax, ['No Data for ' string(satName)]);
        end
    end
end