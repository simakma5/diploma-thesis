close all; clear; clc
directory = fullfile(pwd, '/scripts/postprocessing/');
addpath(genpath(directory))
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
c = 299792458;
epsilon0 = 8.85418782e-12;
mu0 = 1/(c^2*epsilon0);

%% TODO add legends
filename = "feedlength_sweep";
directory = directory + strrep(filename, "_", "-") + "\";
extension = '.s2p';

if isfile(directory + filename + extension)
    data = sparameters(filename + extension);
    figure;
    rfplot(data)
else
    if isfile(directory + filename + "_1" + extension)
        freq = sparameters(filename + "_1" + extension).Frequencies;
    else
        throw(MException("Files not found!"))
    end
    s11 = zeros(1, length(freq));
    s21 = zeros(1, length(freq));
    s12 = zeros(1, length(freq));
    s22 = zeros(1, length(freq));
    result_count = 1;
    while isfile(directory + filename + "_" + result_count + extension)
        data = sparameters(filename + "_" + result_count + extension);
        s11(result_count, :) = squeeze(data.Parameters(1, 1, :));
        s21(result_count, :) = squeeze(data.Parameters(2, 1, :));
        s12(result_count, :) = squeeze(data.Parameters(1, 2, :));
        s22(result_count, :) = squeeze(data.Parameters(2, 2, :));
        result_count = result_count + 1;
    end

    figure;
    tiles = tiledlayout(2, 2, "TileSpacing", "compact");

    nexttile;
    hold on
    for i = 1:result_count-1
        plot(freq*1e-9, 10*log10(abs(s11(i, :))))
    end
    hold off
    ylabel("Magnitude [dB]")
    title("|S11|")
    
    nexttile;
    hold on
    for i = 1:result_count-1
        plot(freq*1e-9, 10*log10(abs(s12(i, :))))
    end
    hold off
    title("|S12|")
    
    nexttile;
    hold on
    for i = 1:result_count-1
        plot(freq*1e-9, 10*log10(abs(s21(i, :))))
    end
    hold off
    xlabel("Frequency [GHz]")
    ylabel("Magnitude [dB]")
    title("|S21|")
    
    nexttile;
    hold on
    for i = 1:result_count-1
        plot(freq*1e-9, 10*log10(abs(s22(i, :))))
    end
    hold off
    xlabel("Frequency [GHz]")
    title("|S22|")
end
