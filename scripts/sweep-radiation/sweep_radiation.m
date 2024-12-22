close all; clear; clc;
script_directory = fileparts(matlab.desktop.editor.getActiveFilename) + "\";
addpath(genpath(script_directory))
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

%%
sweep_parameters = ["feed_length", "outlet_length", "polarizer_length", "chamfer_width", "chamfer_angle", "probe1_distance", "probe_length"];
for sweep_parameter = sweep_parameters
    data_directory = script_directory + strrep(sweep_parameter, "_", "-") + "\";
    sweep_values = importdata(data_directory + "sweep_values.txt")';
    files = dir(data_directory);
    files = {files(3:end).name};
    for i = 1:length(files)
        if strcmp(files{i}, "sweep_values.txt")
            files(i) = [];
        end
    end
    line_groups = cell(length(sweep_values), 1);

    figure("Name", strrep(sweep_parameter, "_", " ") + " sweep");
    tiles = tiledlayout(3, 2, "TileSpacing", "compact");
    for file = files
        [tit, sweep_count, sweep_value] = deal(0);
        frequency = [];
        data = {};
        fileID = fopen(string(data_directory + file), "r");
        while ~feof(fileID)
            line = fgetl(fileID);
            if contains(line, "Parameters")
                sweep_count = sweep_count + 1;
                sweep_value = regexp(lower(line), join(split(sweep_parameter, "_"), "") + "=(\d*\.?\d*)", "tokens");
                sweep_value = str2double(sweep_value{1}{1});
                sweep_index = find(sweep_values == sweep_value);
                continue
            end
            if contains(line, "Frequency / GHz")
                data{sweep_index} = [];
                if tit == 0
                    tit = split(line, sprintf("\t"));
                    end_of_title = regexp(tit{2}, " \(\d*\)", "start");
                    tit = tit{2}(2:end_of_title-1);
                end
                continue;
            end
            if ~startsWith(line, "#") && ~isempty(line)
                values = str2double(split(line, sprintf("\t")));
                if sweep_count == 1
                    frequency = [frequency; values(1)];
                end
                data{sweep_index} = [data{sweep_index}; values(2)];
            end
        end
        fclose(fileID);
        
        nexttile;
        hold on;
        for sweep = 1:length(sweep_values)
            handle = plot(frequency, data{sweep});
            if isempty(line_groups{sweep})
                line_groups{sweep} = handle;
            else
                line_groups{sweep} = [line_groups{sweep}, handle];
            end
        end
        hold off;
        title(tit);
        grid on;
        axis tight
    end

    xlabel(tiles, 'Frequency (GHz)');
    ylabel(tiles, 'Magnitude (-)');
    leg = legend(sweep_values + " mm", "Orientation", "vertical", 'ItemHitFcn', @(src, event) cb_legend(event));
    leg.UserData = line_groups;
    leg.Layout.Tile = "east";
end

%% Functions
function cb_legend(event)
    % disp('Legend item clicked!');
    clicked_line = event.Peer;
    line_groups = event.Source.UserData;

    for group_idx = 1:length(line_groups)
        group = line_groups{group_idx};
        if strcmp(group(end).DisplayName, clicked_line.DisplayName)
            % disp(['Group ' num2str(group_idx) ' toggled']);
            if strcmp(group(1).Visible, 'on')
                for line = group
                    line.Visible = 'off';
                end
            else
                for line = group
                    line.Visible = 'on';
                end
            end
            break;
        end
    end
end
