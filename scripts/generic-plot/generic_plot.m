close all; clear; clc;
script_directory = fileparts(matlab.desktop.editor.getActiveFilename) + "\";
addpath(genpath(script_directory))
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

%% Main
sweep_parameters = ["arc_angle_span"];
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

    for file = files
        file_index = find(strcmpi(files, file));
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
                data{sweep_index} = [];
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
        
        figure("Name", tit);
        hold on;
        for sweep = 1:length(sweep_values)
            handle = plot(frequency, data{sweep});
        end
        hold off;
        grid on;
        axis tight
        legend(sweep_values + " mm", "Orientation", "vertical");
        xlabel("Frequency [GHz]")
    end
end
