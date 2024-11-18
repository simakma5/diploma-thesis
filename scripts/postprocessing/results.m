close all; clear; clc; addpath(genpath(fullfile(pwd, '/Uni/diploma_thesis/scripts/')))
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
c = 299792458;
epsilon0 = 8.85418782e-12;
mu0 = 1/(c^2*epsilon0);

filePath = matlab.desktop.editor.getActiveFilename;
mfilePath = mfilename("fullpath");
%%
opts = detectImportOptions("probe12_phase.txt", "VariableNamesLine", 2);
T = readtable("probe12_phase.txt", opts);

figure("Name", "E-field(0 0 PolarizerLength [1(2)] probe")
x = T.x__Frequency_GHz_;
y1 = T.E_Field_00100__X__1_2___2__PhaseInDegrees_*pi/180;
y2 = unwrap(y1);
plot(x, y1)