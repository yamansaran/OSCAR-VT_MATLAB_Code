%% plot_temp_envelope.m
% Efficiently plots temperature min/max envelope from a large CSV.
% Bins data into N windows, extracts extremes per bin.

clear; clc; close all;

%% --- User Config ---
filename = 'OSCAR_SEET_Temperature.csv';  % <-- your CSV path
nBins    = 2000;                % number of time bins (tune for resolution vs speed)

% Date range filter (set both to empty [] to plot everything)
tStart = datetime(2030, 12, 20, 17, 0, 0);
tEnd   = datetime(2032, 4, 17, 19, 0, 0);

%% --- Read Data ---
opts = detectImportOptions(filename);
opts.VariableNames = {'Time_UTCG', 'Temperature_degC', 'SolarFlux_Wm2', 'SolarIntensity_pct'};
opts = setvartype(opts, 'Time_UTCG', 'char');

fprintf('Reading CSV...\n');
T = readtable(filename, opts);
fprintf('Done. %d rows loaded.\n', height(T));

% Parse timestamps
t = datetime(T.Time_UTCG, 'InputFormat', 'dd MMM yyyy HH:mm:ss.SSS');
temp = T.Temperature_degC;

% Apply date range filter
rangeMask = true(size(t));
if ~isempty(tStart), rangeMask = rangeMask & (t >= tStart); end
if ~isempty(tEnd),   rangeMask = rangeMask & (t <= tEnd);   end
t    = t(rangeMask);
temp = temp(rangeMask);
fprintf('After filtering: %d rows in range.\n', numel(t));

%% --- Bin into time windows, extract min/max ---
nPts   = numel(t);
binIdx = min(floor((0:nPts-1)' / nPts * nBins) + 1, nBins);

tBinCenter = NaT(nBins, 1);
tempMin    = nan(nBins, 1);
tempMax    = nan(nBins, 1);

for k = 1:nBins
    mask = (binIdx == k);
    tBinCenter(k) = t(find(mask, 1, 'first')) + (t(find(mask, 1, 'last')) - t(find(mask, 1, 'first'))) / 2;
    tempMin(k)    = min(temp(mask));
    tempMax(k)    = max(temp(mask));
end

%% --- Plot ---
figure('Position', [100 100 1200 500], 'Color', 'w');

% Filled envelope
fill([tBinCenter; flipud(tBinCenter)], ...
     [tempMin; flipud(tempMax)], ...
     [0.2 0.5 0.9], ...
     'FaceAlpha', 0.35, 'EdgeColor', 'none','HandleVisibility', 'off');
hold on;

% Min/Max lines
plot(tBinCenter, tempMax, '-', 'Color', [0.85 0.2 0.2], 'LineWidth', 1.2, 'DisplayName', 'T_{max}');
plot(tBinCenter, tempMin, '-', 'Color', [0.2 0.3 0.85], 'LineWidth', 1.2, 'DisplayName', 'T_{min}');

hold off;
grid on;
xlabel('Time (UTC)');
ylabel('Temperature (\circC)');
title('Spacecraft Temperature Envelope');
legend('Location', 'best');
set(gca, 'FontSize', 11);

fprintf('Plot complete.\n');



