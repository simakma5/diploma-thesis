% Default figure position: [680, 458, 560, 420]
close all; clear; clc;
addpath(genpath(fileparts(matlab.desktop.editor.getActiveFilename) + "\"))
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultAxesYAxisLocation', 'left');
disp('LaTeX interpreter and left-side Y-axis set as default.');

%% Eigenmode analysis - square (cross-section)
sweep_values = importdata("square_sweep_values.txt")';
[freq, polarizer_length] = get_sweep_plot_data("square_length_for_90deg.txt", "chamferWidth", sweep_values);
[~, amplitude_ratio] = get_sweep_plot_data("square_amplitude_ratio.txt", "chamferWidth", sweep_values);

figure("Name", "Eigenmode analysis - square (cross-section)");
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
ylabel("$E_2/E_1\ [\mathrm{dB}]$")

leg = legend("$" + sweep_values + "\ \mathrm{mm}$", "Orientation", "horizontal");
leg.Layout.Tile = "north";

saveas(gcf, fullfile(pwd, '\latex\src\square_polarizer_cross_section.svg'), 'svg')

%% Eigenmode analysis - circular (cross-section)
sweep_values = importdata("circular_sweep_values.txt")';
[freq, polarizer_length] = get_sweep_plot_data("circular_length_for_90deg.txt", "chamferWidth", sweep_values);
[~, amplitude_ratio] = get_sweep_plot_data("circular_amplitude_ratio.txt", "chamferWidth", sweep_values);

figure("Name", "Eigenmode analysis - circular (cross-section)");
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
ylabel("$E_2/E_1\ [\mathrm{dB}]$")

leg = legend("$" + sweep_values + "\ \mathrm{mm}$", "Orientation", "horizontal");
leg.Layout.Tile = "north";

saveas(gcf, fullfile(pwd, '\latex\src\circular_polarizer_cross_section.svg'), 'svg')

%% Eigenmode analysis (cutoff frequencies)
[chamfer_width, cutoff_frequency] = get_single_plot_data("square_cutoff_frequency.txt");

figure("Name", "Eigenmode analysis (cutoff frequencies)");
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
title("Square waveguide")

[chamfer_width, cutoff_frequency] = get_single_plot_data("circular_cutoff_frequency.txt");

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
title("Circular waveguide")

saveas(gcf, fullfile(pwd, '\latex\src\polarizers_cutoff_frequencies.svg'), 'svg')

%% Polarizer - farfield (test)
% farfield1 = process_farfield_data("radiation1.txt");
% farfield2 = process_farfield_data("radiation2.txt");
% 
% figure("Name", "Polarizer - farfield (test)");
% cst_polarplot(farfield1.Theta, [farfield1.Directivity, farfield2.Directivity], [-20, 15])
% legend(["Mode 1", "Mode2"])

%% Polarizer - radiation patterns
[theta, rhcp1] = get_single_plot_data("rhcp1.txt");
[~, lhcp1] = get_single_plot_data("lhcp1.txt");
[~, rhcp2] = get_single_plot_data("rhcp2.txt");
[~, lhcp2] = get_single_plot_data("lhcp2.txt");
[freq, axial_ratio1] = get_single_plot_data("axial_ratio1.txt");
% [~, axial_ratio2] = get_single_plot_data("axial_ratio2.txt");

figure("Name", "Polarizer - radiation patterns");
polartiles = tiledlayout(3, 2, "TileSpacing", "compact");
rlim_all = [-50, 10];

pax = polaraxes(polartiles);
pax.Layout.Tile = 1;
cst_polarplot(theta, rhcp1, rlim_all, pax)
pax.Title.String = "RHCP (mode 1)";

pax = polaraxes(polartiles);
pax.Layout.Tile = 2;
cst_polarplot(theta, lhcp1, rlim_all, pax)
pax.Title.String = "LHCP (mode 1)";

pax = polaraxes(polartiles);
pax.Layout.Tile = 3;
cst_polarplot(theta, rhcp2, rlim_all, pax)
pax.Title.String = "RHCP (mode 2)";

pax = polaraxes(polartiles);
pax.Layout.Tile = 4;
cst_polarplot(theta, lhcp2, rlim_all, pax)
pax.Title.String = "LHCP (mode 2)";

nexttile([1, 2]);
hold on
plot(freq, axial_ratio1)
% plot(freq, axial_ratio2)
hold off
grid on;
box on;
xlim([min(freq), max(freq)]);
xlabel("$f\ [\mathrm{GHz}]$")
ylabel("$\mathrm{AR}\ [\mathrm{dB}]$")
title("Axial ratio (identical for both modes)")

set(gcf, "Position", [680, 160, 560, 700])
saveas(gcf, fullfile(pwd, '\latex\src\polarizer_radiation.svg'), 'svg')

%% Feed - single (reflection)
[freq, reflection] = get_single_plot_data("single_feed_reflection.txt");

figure("Name", "Feed - single (reflection)");
plot(freq, reflection)
grid on;
box on;
xlim([min(freq), max(freq)]);
xlabel("$f\ [\mathrm{GHz}]$")
ylabel("$S_{11}\ [\mathrm{dB}]$")
set(gcf, "Position", [680, 458, 560, 260])
saveas(gcf, fullfile(pwd, '\latex\src\single_feed_reflection.svg'), 'svg')

%% Feed - grating distance sweep (reflection)
sweep_values = importdata("grating_distance_sweep_values.txt")';
[freq, reflection] = get_sweep_plot_data("grating_distance_sweep_data.txt", "gratingDistance", sweep_values);

figure("Name", "Feed - grating distance sweep (reflection)");
hold on;
for sweep = 1:length(sweep_values)
    plot(freq, reflection{sweep});
end
hold off;
grid on;
box on;
xlim([min(freq), max(freq)]);
xlabel("$f\ [\mathrm{GHz}]$")
ylabel("$S_{11}\ [\mathrm{dB}]$")
sweep_values_legend = ["$\lambda_{\mathrm g}/4$", "$3\lambda_{\mathrm g}/8$", "$\lambda_{\mathrm g}/2$", "$5\lambda_{\mathrm g}/8$", "$3\lambda_{\mathrm g}/4$"];
legend(sweep_values_legend, "Orientation", "horizontal", "Location", "northoutside");

set(gcf, "Position", [680, 458, 560, 420])  % [680, 458, 560, 420]
saveas(gcf, fullfile(pwd, '\latex\src\grating_distance_sweep.svg'), 'svg')

%% Feed - single with grating (reflection)
[freq, reflection] = get_single_plot_data("single_feed_with_grating_reflection.txt");

figure("Name", "Feed - single with grating (reflection)");
plot(freq, reflection)
grid on;
box on;
xlim([min(freq), max(freq)]);
xlabel("$f\ [\mathrm{GHz}]$")
ylabel("$S_{11}\ [\mathrm{dB}]$")
set(gcf, "Position", [680, 458, 560, 260])
saveas(gcf, fullfile(pwd, '\latex\src\single_feed_with_grating_reflection.svg'), 'svg')

%% Feed - dual with guideline values (S-parameters)
sdata = sparameters("dual_feed.s2p");
freq = sdata.Frequencies*1e-9;
s11 = 20*log10(abs(squeeze(sdata.Parameters(1, 1, :))));
s21 = 20*log10(abs(squeeze(sdata.Parameters(2, 1, :))));
s12 = 20*log10(abs(squeeze(sdata.Parameters(1, 2, :))));
s22 = 20*log10(abs(squeeze(sdata.Parameters(2, 2, :))));

figure("Name", "Feed - dual with guideline values (S-parameters)");
tiledlayout(2, 1, "TileSpacing", "compact");

nexttile;
hold on;
plot(freq, s11)
plot(freq, s22)
hold off;
grid on;
box on;
xlim([min(freq), max(freq)]);
ylabel("$[\mathrm{dB}]$")

nexttile;
hold on;
plot(0,0)
plot(0,0)
plot(freq, s21)
plot(freq, s12)
hold off;
grid on;
box on;
xlim([min(freq), max(freq)]);
xlabel("$f\ [\mathrm{GHz}]$")
ylabel("$[\mathrm{dB}]$")

leg = legend(["$S_{11}$", "$S_{22}$", "$S_{21}$", "$S_{12}$"], "Orientation", "horizontal");
leg.Layout.Tile = "north";

saveas(gcf, fullfile(pwd, '\latex\src\dual_feed_sparameters.svg'), 'svg')

%% Functions - get_sweep_plot_data
function [xdata, ydata] = get_sweep_plot_data(file, sweep_parameter, sweep_values)
    sweep_count = 0;
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
    if length(ydata) == 1
        ydata = ydata{1};
    end
end

%% Function - process_farfield_data
function T = process_farfield_data(filename)
    fid = fopen(filename, 'r');
    headerLine = fgetl(fid);
    fclose(fid);
    headerLine = strtrim(regexprep(headerLine, '\s*\[.*?\]', ''));
    headerLine = regexprep(headerLine, 'Abs\(Dir.*?\)', 'Directivity');
    headerLine = regexprep(headerLine, 'Abs\(Left.*?\)', 'LHCP');
    headerLine = regexprep(headerLine, 'Abs\(Right.*?\)', 'RHCP');
    headerLine = regexprep(headerLine, 'Ax\.Ratio.*?', 'AR');
    headerLine = regexprep(headerLine, 'Phase\(.*?\)', 'REMOVE');
    variableNames = strsplit(headerLine);
    T = readtable(filename, 'FileType', 'text', 'Delimiter', ' ', ...
                  'MultipleDelimsAsOne', true, 'HeaderLines', 2, ...
                  'ReadVariableNames', false, 'TreatAsEmpty', {'NA'});
    if any(isnan(T{:, 1}))
        T(:, 1) = [];
    end
    for i = 1:length(variableNames)
        if i > length(variableNames)
            break
        end
        if variableNames(i) == "REMOVE"
            variableNames(i) = [];
            T(:, i) = [];
        end
    end
    T.Properties.VariableNames = variableNames;
    T.Theta(182:end) = (181:359)';
    T(361, :) = T(1, :);
end

%% Function - cst_polarplot
function cst_polarplot(theta, rdata, rlim, pax)
    if nargin < 4
        pax = polaraxes;
    end
    hold on
    for i = 1:size(rdata, 2)
        polarplot(deg2rad(theta), rdata(:, i))
    end
    hold off
    pax.ThetaZeroLocation = 'top';
    pax.RLim = rlim;
    pax.RAxisLocation = 0;
end