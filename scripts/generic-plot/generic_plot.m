close all; clear; clc;
addpath(genpath(fileparts(matlab.desktop.editor.getActiveFilename) + "\"))
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultAxesYAxisLocation', 'left');
disp('LaTeX interpreter and left-side Y-axis set as default.');

%% Polarizer - circular (cross-section)
sweep_values = importdata("polarizer_circular_sweep_values.txt")';
[freq, polarizer_length] = get_sweep_plot_data("polarizer_circular_length_for_90deg.txt", "chamferWidth", sweep_values);
[~, amplitude_ratio] = get_sweep_plot_data("polarizer_circular_amplitude_ratio.txt", "chamferWidth", sweep_values);

figure("Name", "Polarizer - circular (cross-section)");
tiledlayout(2, 1, "TileSpacing", "compact");

nexttile;
hold on;
for sweep = 1:length(sweep_values)
    plot(freq, polarizer_length{sweep});
end
hold off
grid on;
box on;
xlim([min(freq), max(freq)]);
ylabel("$L_\perp\ [\mathrm{mm}]$")

nexttile;
hold on
for sweep = 1:length(sweep_values)
    plot(freq, amplitude_ratio{sweep});
end
hold off;
grid on;
box on;
xlim([min(freq), max(freq)]);
xlabel("$f\ [\mathrm{GHz}]$")
ylabel("$E_2/E_1\ [-]$")

leg = legend("$" + sweep_values + "\ \mathrm{mm}$", "Orientation", "horizontal");
leg.Layout.Tile = "north";

%% Polarizer - square (cross-section)
sweep_values = importdata("polarizer_square_sweep_values.txt")';
[freq, polarizer_length] = get_sweep_plot_data("polarizer_square_length_for_90deg.txt", "chamferWidth", sweep_values);
[~, amplitude_ratio] = get_sweep_plot_data("polarizer_square_amplitude_ratio.txt", "chamferWidth", sweep_values);

figure("Name", "Polarizer - square (cross-section)");
tiledlayout(2, 1, "TileSpacing", "compact");

nexttile;
hold on;
for sweep = 1:length(sweep_values)
    plot(freq, polarizer_length{sweep});
end
hold off
grid on;
box on;
xlim([min(freq), max(freq)]);
ylabel("$L_\perp\ [\mathrm{mm}]$")

nexttile;
hold on
for sweep = 1:length(sweep_values)
    plot(freq, amplitude_ratio{sweep});
end
hold off;
grid on;
box on;
xlim([min(freq), max(freq)]);
xlabel("$f\ [\mathrm{GHz}]$")
ylabel("$E_2/E_1\ [-]$")

leg = legend("$" + sweep_values + "\ \mathrm{mm}$", "Orientation", "horizontal");
leg.Layout.Tile = "north";

%% Polarizers (cutoff frequencies)
[chamfer_width, cutoff_frequency] = get_single_plot_data("polarizer_circular_cutoff_frequency.txt");

figure("Name", "Polarizers (cutoff frequencies)");
tiledlayout(2, 1, "TileSpacing", "compact");

nexttile;
hold on
for mode = 1:2
    plot(chamfer_width, cutoff_frequency{mode});
end
hold off;
grid on;
box on;
xlim([min(chamfer_width), max(chamfer_width)]);
ylabel("$f_{\mathrm{cutoff}}\ [\mathrm{GHz}]$")
legend(["Mode 1", "Mode 2"], "Location", "east")
title("Circular waveguide")

[chamfer_width, cutoff_frequency] = get_single_plot_data("polarizer_square_cutoff_frequency.txt");

nexttile;
hold on
for mode = 1:2
    plot(chamfer_width, cutoff_frequency{mode});
end
hold off;
grid on;
box on;
xlim([min(chamfer_width), max(chamfer_width)]);
xlabel("$w\ [\mathrm{mm}]$")
ylabel("$f_{\mathrm{cutoff}}\ [\mathrm{GHz}]$")
legend(["Mode 1", "Mode 2"], "Location", "east")
title("Square waveguide")

%% Functions - get_sweep_plot_data
function [xdata, ydata] = get_sweep_plot_data(file, sweep_parameter, sweep_values)
    [sweep_count, sweep_value] = deal(0);
    xdata = [];
    ydata = {};
    fileID = fopen(string(file), "r");
    while ~feof(fileID)
        line = fgetl(fileID);
        if contains(line, "Parameters")
            sweep_count = sweep_count + 1;
            sweep_value = regexp(lower(line), lower(sweep_parameter) + "=(\d*\.?\d*)", "tokens");
            sweep_value = str2double(sweep_value{1}{1});
            sweep_index = find(sweep_values == sweep_value);
            ydata{sweep_index} = [];
            continue;
        end
        if ~startsWith(line, "#") && ~isempty(line)
            values = str2double(split(line, sprintf("\t")));
            if sweep_count == 1
                xdata = [xdata; values(1)];
            end
            ydata{sweep_index} = [ydata{sweep_index}; values(2)];
        end
    end
    fclose(fileID);
end

%% Functions - get_single_plot_data
function [xdata, ydata] = get_single_plot_data(file)
    curve = 0;
    xdata = [];
    ydata = {};
    fileID = fopen(string(file), "r");
    while ~feof(fileID)
        line = fgetl(fileID);
        if startsWith(line, "#---")
            curve = curve + 1;
            ydata{curve} = [];
            continue
        end
        if ~startsWith(line, "#") && ~isempty(line)
            values = str2double(split(line, sprintf("\t")));
            if curve == 1
                xdata = [xdata; values(1)];
            end
            ydata{curve} = [ydata{curve}; values(2)];
        end
    end
    fclose(fileID);
end
