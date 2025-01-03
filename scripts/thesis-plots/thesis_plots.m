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

%% Eigenmode analysis: square (cross-section)
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

%% Eigenmode analysis: circular (cross-section)
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

%% Eigenmode analysis: cutoff frequencies
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

%% Polarizer: radiation patterns
[theta, rhcp1] = get_single_plot_data("rhcp1.txt");
[~, lhcp1] = get_single_plot_data("lhcp1.txt");
[~, rhcp2] = get_single_plot_data("rhcp2.txt");
[~, lhcp2] = get_single_plot_data("lhcp2.txt");
[freq, axial_ratio1] = get_single_plot_data("axial_ratio1.txt");
% [~, axial_ratio2] = get_single_plot_data("axial_ratio2.txt");

figure("Name", "Polarizer - radiation patterns");
tiles = tiledlayout(3, 2, "TileSpacing", "compact");
rlim_all = [-50, 10];

pax = polaraxes(tiles);
pax.Layout.Tile = 1;
cst_polarplot(theta, rhcp1, rlim_all, pax)
pax.Title.String = "RHCP (mode 1)";

pax = polaraxes(tiles);
pax.Layout.Tile = 2;
cst_polarplot(theta, lhcp1, rlim_all, pax)
pax.Title.String = "LHCP (mode 1)";

pax = polaraxes(tiles);
pax.Layout.Tile = 3;
cst_polarplot(theta, rhcp2, rlim_all, pax)
pax.Title.String = "RHCP (mode 2)";

pax = polaraxes(tiles);
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

set(gcf, "Position", [680, 178, 560, 700])
saveas(gcf, fullfile(pwd, '\latex\src\polarizer_radiation.svg'), 'svg')

%% Feed: single (reflection)
[freq, reflection] = get_single_plot_data("single_feed_reflection.txt");

figure("Name", "Feed - single (reflection)");
plot(freq, reflection)
grid on;
box on;
xlim([min(freq), max(freq)]);
xlabel("$f\ [\mathrm{GHz}]$")
ylabel("$S_{11}\ [\mathrm{dB}]$")
set(gcf, "Position", [680, 618, 560, 260])
saveas(gcf, fullfile(pwd, '\latex\src\single_feed_reflection.svg'), 'svg')

%% Feed: grating distance sweep (reflection)
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

saveas(gcf, fullfile(pwd, '\latex\src\grating_distance_sweep.svg'), 'svg')

%% Feed: single with grating (reflection)
[freq, reflection] = get_single_plot_data("single_feed_with_grating_reflection.txt");

figure("Name", "Feed - single with grating (reflection)");
plot(freq, reflection)
grid on;
box on;
xlim([min(freq), max(freq)]);
xlabel("$f\ [\mathrm{GHz}]$")
ylabel("$S_{11}\ [\mathrm{dB}]$")
set(gcf, "Position", [680, 618, 560, 260])
saveas(gcf, fullfile(pwd, '\latex\src\single_feed_with_grating_reflection.svg'), 'svg')

%% Feed: dual with guideline values (S-parameters)
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

%% Antenna: conical horns comparison
[theta, magus_gain_pattern] = get_single_plot_data("magus_gain_pattern.txt");
[~, final_gain_pattern] = get_single_plot_data("final_gain_pattern.txt");
[freq_gain, magus_gain] = get_single_plot_data("magus_gain.txt");
[~, final_gain] = get_single_plot_data("final_gain.txt");
[freq_reflection, magus_reflection] = get_single_plot_data("magus_reflection.txt");
[~, final_reflection] = get_single_plot_data("final_reflection.txt");

figure("Name", "Antenna - conical horns comparison")
tiles = tiledlayout(3, 2, "TileSpacing", "compact");
rlim_all = [-40, 20];

pax = polaraxes(tiles);
pax.Layout.Tile = 1;
pax.Layout.TileSpan = [1,2];
hold on
cst_polarplot(theta, magus_gain_pattern, rlim_all, pax)
cst_polarplot(theta, final_gain_pattern, rlim_all, pax)
hold off
pax.Title.String = "Gain ($f = 5.2\ \mathrm{GHz}$)";

nexttile([1, 2]);
hold on
plot(freq_gain, magus_gain)
plot(freq_gain, final_gain)
hold off
grid on;
box on;
xlim([min(freq_gain), max(freq_gain)]);
xlabel("$f\ [\mathrm{GHz}]$")
ylabel("$G\ [\mathrm{dB}]$")
title("Gain (boresight)")

nexttile([1, 2]);
hold on
plot(freq_reflection, magus_reflection)
plot(freq_reflection, final_reflection)
hold off
grid on;
box on;
xlim([min(freq_reflection), max(freq_reflection)]);
xlabel("$f\ [\mathrm{GHz}]$")
ylabel("$S_{11}\ [\mathrm{dB}]$")
title("Reflection")

leg = legend(["Antenna Magus", "Finalized"], "Orientation", "horizontal");
leg.Layout.Tile = "north";

set(gcf, "Position", [680, 78, 560, 800])
saveas(gcf, fullfile(pwd, '\latex\src\conical_horns_comparison.svg'), 'svg')

%% Final structure: grating vs MVP comparison (S-parameters)
grating = sparameters("final_grating.s2p");
freq = grating.Frequencies*1e-9;
grating_s11 = 20*log10(abs(squeeze(grating.Parameters(1, 1, :))));
grating_s21 = 20*log10(abs(squeeze(grating.Parameters(2, 1, :))));
grating_s12 = 20*log10(abs(squeeze(grating.Parameters(1, 2, :))));
grating_s22 = 20*log10(abs(squeeze(grating.Parameters(2, 2, :))));
mvp = sparameters("final.s2p");
mvp_s11 = 20*log10(abs(squeeze(mvp.Parameters(1, 1, :))));
mvp_s21 = 20*log10(abs(squeeze(mvp.Parameters(2, 1, :))));
mvp_s12 = 20*log10(abs(squeeze(mvp.Parameters(1, 2, :))));
mvp_s22 = 20*log10(abs(squeeze(mvp.Parameters(2, 2, :))));

figure("Name", "Final structure: grating vs MVP comparison (S-parameters)")
tiledlayout(2, 2, "TileSpacing", "compact")

nexttile;
hold on;
plot(freq, grating_s11)
plot(freq, mvp_s11)
hold off;
grid on;
box on;
xlim([min(freq), max(freq)]);
xlabel("$f\ [\mathrm{GHz}]$")
ylabel("$S_{11}\ [\mathrm{dB}]$")

nexttile;
hold on;
plot(freq, grating_s12)
plot(freq, mvp_s12)
hold off;
grid on;
box on;
xlim([min(freq), max(freq)]);
xlabel("$f\ [\mathrm{GHz}]$")
ylabel("$S_{12}\ [\mathrm{dB}]$")

nexttile;
hold on;
plot(freq, grating_s21)
plot(freq, mvp_s21)
hold off;
grid on;
box on;
xlim([min(freq), max(freq)]);
xlabel("$f\ [\mathrm{GHz}]$")
ylabel("$S_{21}\ [\mathrm{dB}]$")

nexttile;
hold on;
plot(freq, grating_s22)
plot(freq, mvp_s22)
hold off;
grid on;
box on;
xlim([min(freq), max(freq)]);
xlabel("$f\ [\mathrm{GHz}]$")
ylabel("$S_{22}\ [\mathrm{dB}]$")

leg = legend(["Grating", "No grating"], "Orientation", "horizontal");
leg.Layout.Tile = "north";
saveas(gcf, fullfile(pwd, '\latex\src\grating_vs_mvp_sparameters.svg'), 'svg')

%% Final - grating vs MVP comparison (gain)
% [theta, grating_gain] = get_single_plot_data("final_grating_gain_5G2.txt");
% [~, mvp_gain] = get_single_plot_data("final_gain_5G2.txt");
% 
% figure("Name", "Final - grating vs MVP comparison (gain)")
% tiles = tiledlayout(1, 2, "TileSpacing", "compact");
% rlim_all = [-40, 20];
% 
% pax = polaraxes(tiles);
% pax.Layout.Tile = 1;
% hold on
% cst_polarplot(theta, grating_gain{1}, rlim_all, pax)
% cst_polarplot(theta, mvp_gain{1}, rlim_all, pax)
% hold off
% pax.Title.String = "Port 1";
% 
% pax = polaraxes(tiles);
% pax.Layout.Tile = 2;
% hold on
% cst_polarplot(theta, grating_gain{2}, rlim_all, pax)
% cst_polarplot(theta, mvp_gain{2}, rlim_all, pax)
% hold off
% pax.Title.String = "Port 2";
% 
% leg = legend(["Grating", "No grating"], "Orientation", "horizontal");
% leg.Layout.Tile = "north";
% saveas(gcf, fullfile(pwd, '\latex\src\grating_vs_mvp_gain.svg'), 'svg')

%% Final - grating vs MVP comparison (polarizations)
% [freq, grating_ar] = get_single_plot_data("final_grating_boresight_ar.txt");
% [freq_mvp, mvp_ar] = get_single_plot_data("final_boresight_ar.txt");
% [~, grating_pol] = get_single_plot_data("final_grating_boresight_polarizations.txt");
% [~, mvp_pol] = get_single_plot_data("final_boresight_polarizations.txt");
% 
% figure("Name", "Feed - grating vs MVP comparison (polarizations)");
% tiles = tiledlayout(3, 1, "TileSpacing", "compact");
% 
% nexttile;
% hold on;
% plot(freq, grating_ar{1})
% plot(freq, grating_ar{2})
% plot(freq_mvp, mvp_ar{2}, "LineStyle", "--", "Color", [.4 .4 .4])
% hold off;
% grid on;
% box on;
% xlim([min(freq), max(freq)]);
% ylim([0, 7]);
% ylabel("$\mathrm{AR}\ [\mathrm{dB}]$")
% title("Axial ratio")
% 
% nexttile;
% hold on;
% plot(freq, grating_pol{3})
% plot(freq, grating_pol{2})
% plot(freq_mvp, mvp_pol{2}, "LineStyle", "--", "Color", [.4 .4 .4])
% hold off;
% grid on;
% box on;
% xlim([min(freq), max(freq)]);
% ylabel("$G\ [\mathrm{dBi}]$")
% title("Co-polarization")
% 
% nexttile;
% hold on;
% plot(freq, grating_pol{1})
% plot(freq, grating_pol{4})
% plot(freq_mvp, mvp_pol{1}, "LineStyle", "--", "Color", [.4 .4 .4])
% hold off;
% grid on;
% box on;
% xlim([min(freq), max(freq)]);
% xlabel("$f\ [\mathrm{GHz}]$")
% ylabel("$G\ [\mathrm{dBi}]$")
% title("Cross-polarization")
% 
% leg = legend(["Grating (port 1)", "Grating (port 2)", "No grating"], "Orientation", "horizontal");
% leg.Layout.Tile = "north";
% set(gcf, "Position", [680, 178, 560, 800])
% saveas(gcf, fullfile(pwd, '\latex\src\grating_vs_mvp_polarization.svg'), 'svg')

%% Measurement: S-parameters
sim_data = sparameters("final.s2p");
sim_freq = sim_data.Frequencies(:)*1e-9;
sim_s11 = 20*log10(abs(squeeze(sim_data.Parameters(1, 1, :))));
sim_s21 = 20*log10(abs(squeeze(sim_data.Parameters(2, 1, :))));
sim_s12 = 20*log10(abs(squeeze(sim_data.Parameters(1, 2, :))));
sim_s22 = 20*log10(abs(squeeze(sim_data.Parameters(2, 2, :))));
for meas = ["meas1", "meas2"]
    meas_data = sparameters(meas + "_sparameters.s2p");
    meas_freq = meas_data.Frequencies(:)*1e-9;
    meas_s11 = 20*log10(abs(squeeze(meas_data.Parameters(1, 1, :))));
    meas_s21 = 20*log10(abs(squeeze(meas_data.Parameters(2, 1, :))));
    meas_s12 = 20*log10(abs(squeeze(meas_data.Parameters(1, 2, :))));
    meas_s22 = 20*log10(abs(squeeze(meas_data.Parameters(2, 2, :))));
    
    figure("Name", "Measurement " + meas{1}(end) + ": S-parameters")
    tiledlayout(2, 2, "TileSpacing", "compact")
    
    nexttile;
    hold on
    plot(meas_freq, meas_s11)
    % plot(sim_freq, sim_s11)
    hold off
    grid on;
    box on;
    xlim([min(meas_freq), max(meas_freq)]);
    xlabel("$f\ [\mathrm{GHz}]$")
    ylabel("$S_{11}\ [\mathrm{dB}]$")
    
    nexttile;
    hold on
    plot(meas_freq, meas_s21)
    % plot(sim_freq, sim_s21)
    hold off
    grid on;
    box on;
    xlim([min(meas_freq), max(meas_freq)]);
    xlabel("$f\ [\mathrm{GHz}]$")
    ylabel("$S_{12}\ [\mathrm{dB}]$")
    
    nexttile;
    hold on
    plot(meas_freq, meas_s21)
    % plot(sim_freq, sim_s21)
    hold off
    grid on;
    box on;
    xlim([min(meas_freq), max(meas_freq)]);
    xlabel("$f\ [\mathrm{GHz}]$")
    ylabel("$S_{21}\ [\mathrm{dB}]$")
    
    nexttile;
    hold on
    plot(meas_freq, meas_s22)
    % plot(sim_freq, sim_s22)
    hold off
    grid on;
    box on;
    xlim([min(meas_freq), max(meas_freq)]);
    xlabel("$f\ [\mathrm{GHz}]$")
    ylabel("$S_{22}\ [\mathrm{dB}]$")
    
    % leg = legend(["Measurement", "Simulation"], "Orientation", "horizontal");
    % leg.Layout.Tile = "north";
    saveas(gcf, fullfile(pwd, '\latex\src\' + meas + '_sparameters.svg'), 'svg')
end

%% Measurement: boresight radiation
[freq, sim_gain] = get_single_plot_data("final_boresight_gain.txt");
[~, sim_ar] = get_single_plot_data("final_boresight_ar.txt");
for meas = ["meas1", "meas2"]
    [~, gain1] = get_umt_data(meas + "_port1_gain.dat", "Gain");
    [~, gain2] = get_umt_data(meas + "_port2_gain.dat", "Gain");
    [~, ar1] = get_umt_data(meas + "_port1_axial_ratio.dat", "Gain");
    [~, ar2] = get_umt_data(meas + "_port2_axial_ratio.dat", "Gain");
    
    figure("Name", "Measurement " + meas{1}(end) + ": boresight radiation");
    tiledlayout(2, 2, "TileSpacing", "compact")

    nexttile;
    hold on
    plot(freq, gain1)
    % plot(freq, sim_gain)
    hold off
    grid on;
    box on;
    xlim([min(freq), max(freq)]);
    ylim([10, 25])
    ylabel("$G\ [\mathrm{dBi}]$")
    title("Port 1 boresight gain")

    nexttile;
    hold on
    plot(freq, gain2)
    % plot(freq, sim_gain)
    hold off
    grid on;
    box on;
    xlim([min(freq), max(freq)]);
    ylim([10, 25])
    title("Port 2 boresight gain")

    nexttile;
    hold on
    plot(freq, ar1)
    % plot(freq, sim_ar)
    hold off
    grid on;
    box on;
    xlim([min(freq), max(freq)]);
    ylim([0, 6])
    xlabel("$f\ [\mathrm{GHz}]$")
    ylabel("$\mathrm{AR}\ [\mathrm{dB}]$")
    title("Port 1 boresight axial ratio")

    nexttile;
    hold on
    plot(freq, ar2)
    % plot(freq, sim_ar)
    hold off
    grid on;
    box on;
    xlim([min(freq), max(freq)]);
    ylim([0, 6])
    xlabel("$f\ [\mathrm{GHz}]$")
    title("Port 2 boresight axial ratio")

    % leg = legend(["Measurement", "Simulation"], "Orientation", "horizontal");
    % leg.Layout.Tile = "north";
    saveas(gcf, fullfile(pwd, '\latex\src\' + meas + '_boresight_radiation.svg'), 'svg')
end

%% Measurement - radiation
% [sim_theta, sim_port1_copol] = get_single_plot_data("final_port1_copol_5G0Hz.txt");
for meas = ["meas1", "meas2"]
    for angle_var = ["elevation", "azimuth"]
        files = [
            "_port1_" + angle_var + "_5G0Hz.dat", "_port2_" + angle_var + "_5G0Hz.dat", ...
            "_port1_" + angle_var + "_5G3Hz.dat", "_port2_" + angle_var + "_5G3Hz.dat", ...
            "_port1_" + angle_var + "_5G5Hz.dat", "_port2_" + angle_var + "_5G5Hz.dat"
        ];
        titles =  [
            "Port 1 ($5\ \mathrm{GHz}$)", "Port 2 ($5\ \mathrm{GHz}$)", ...
            "Port 1 ($5.3\ \mathrm{GHz}$)", "Port 2 ($5.3\ \mathrm{GHz}$)", ...
            "Port 1 ($5.5\ \mathrm{GHz}$)", "Port 2 ($5.5\ \mathrm{GHz}$)"
        ];
        
        figure("Name", "Measurement " + meas{1}(end) + ": " + angle_var + " " + "pattern")
        tiles = tiledlayout(3, 2, "TileSpacing", "compact");
        rlim = [-40, 20];
        for i = 1:6
            pax = polaraxes(tiles);
            [angle, gains] = get_umt_data(meas + files(i), ["Co-pol.", "Cross-pol."]);
            polarplot(deg2rad(angle), gains{1}, deg2rad(angle), gains{2})
            pax.ThetaZeroLocation = 'top';
            pax.ThetaLim = [-90, 90];
            pax.RLim = rlim;
            pax.Layout.Tile = i;
            pax.Title.String = titles(i);
        end
        leg = legend(["Co-polarization", "Cross-polarization"], "Orientation", "Horizontal");
        leg.Layout.Tile = "north";
        set(gcf, "Position", [680, 318, 560, 700])
        saveas(gcf, fullfile(pwd, '\latex\src\' + meas + '_' + angle_var + '_radiation.svg'), 'svg')
    end
end

%% Functions
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
    pax.ThetaLim = [-180, 180];
    pax.RLim = rlim;
    pax.RAxisLocation = 0;
end

function [xdata, ydata] = get_umt_data(filename, labels)
    xdata = [];
    ydata = cell(1, length(labels));
    for label = labels
        ydata{labels==label} = [];
    end
    fileID = fopen(string(filename), "r");
    while ~feof(fileID)
        line = fgetl(fileID);
        values = split(line, sprintf("\t"));
        if ~strcmp(values{3}, " ")
            for label = labels
                if strcmp(values{1}, label) || strcmp(values{1}(2:end), label)
                    xdata = [xdata; double(string(values{3}))];
                    ydata{labels==label} = [ydata{labels==label}; double(string(values{4}))];
                end
            end
        end
    end
    fclose(fileID);
    xdata = unique(xdata);
    if length(ydata) == 1
        ydata = ydata{1};
    end
end
