%% Load satellite data with preserved column names
opts = detectImportOptions('satellites.xlsx');
opts.VariableNamingRule = 'preserve';  % Keep original names
data = readtable('satellites.xlsx', opts);

% Trim whitespace from all column names
data.Properties.VariableNames = strtrim(data.Properties.VariableNames);

fprintf('Loaded %d rows from Excel\n', height(data));

%% Create satellite structure array
satellites = struct([]);

for i = 1:height(data)
    satellites(i).Name = char(data.("Name of Satellite, Alternate Names")(i));
    satellites(i).OfficialName = char(data.("Current Official Name of Satellite")(i));
    satellites(i).Country = char(data.("Country/Org of UN Registry")(i));
    satellites(i).CountryOfOperator = char(data.("Country of Operator/Owner")(i));
    satellites(i).Operator = char(data.("Operator/Owner")(i));
    satellites(i).Users = char(data.("Users")(i));
    satellites(i).Purpose = char(data.("Purpose")(i));
    satellites(i).DetailedPurpose = char(data.("Detailed Purpose")(i));
    satellites(i).ClassOfOrbit = char(data.("Class of Orbit")(i));
    satellites(i).TypeOfOrbit = char(data.("Type of Orbit")(i));
    
    % Numeric fields - use direct assignment, not str2double
    % If the value is already numeric, just use it. If it's text/missing, handle it
    val = data.("Longitude of GEO (degrees)")(i);
    satellites(i).LongitudeGEO_deg = convertToNumeric(val);
    
    val = data.("Perigee (km)")(i);
    satellites(i).Perigee_km = convertToNumeric(val);
    
    val = data.("Apogee (km)")(i);
    satellites(i).Apogee_km = convertToNumeric(val);
    
    val = data.("Eccentricity")(i);
    satellites(i).Eccentricity = convertToNumeric(val);
    
    val = data.("Inclination (degrees)")(i);
    satellites(i).Inclination_deg = convertToNumeric(val);
    
    val = data.("Period (minutes)")(i);
    satellites(i).Period_min = convertToNumeric(val);
    
    val = data.("Launch Mass (kg.)")(i);
    satellites(i).LaunchMass_kg = convertToNumeric(val);
    
    val = data.("Dry Mass (kg.)")(i);
    satellites(i).DryMass_kg = convertToNumeric(val);
    
    val = data.("Power (watts)")(i);
    satellites(i).Power_watts = convertToNumeric(val);
    
    % Calculate mean altitude
    if ~isnan(satellites(i).Perigee_km) && ~isnan(satellites(i).Apogee_km)
        satellites(i).MeanAltitude_km = (satellites(i).Perigee_km + satellites(i).Apogee_km) / 2;
    else
        satellites(i).MeanAltitude_km = NaN;
    end
    
    % Additional fields
    satellites(i).DateOfLaunch = char(data.("Date of Launch")(i));
    
    val = data.("Expected Lifetime (yrs.)")(i);
    satellites(i).ExpectedLifetime_yrs = convertToNumeric(val);
    
    satellites(i).Contractor = char(data.("Contractor")(i));
    satellites(i).CountryOfContractor = char(data.("Country of Contractor")(i));
    satellites(i).LaunchSite = char(data.("Launch Site")(i));
    satellites(i).LaunchVehicle = char(data.("Launch Vehicle")(i));
    satellites(i).COSPAR = char(data.("COSPAR Number")(i));
    
    val = data.("NORAD Number")(i);
    satellites(i).NORAD_Number = convertToNumeric(val);
    
    satellites(i).Comments = char(data.("Comments")(i));
    satellites(i).SourceUsedForOrbitalData = char(data.("Source Used for Orbital Data")(i));
    
    % Price field - initially empty
    satellites(i).Price = NaN;
end

fprintf('Created structure array with %d satellites\n', length(satellites));

%% Debug: Check some altitudes
fprintf('\nFirst 10 satellites - checking altitudes:\n');
for i = 1:min(10, length(satellites))
    fprintf('Satellite %d: %s\n', i, satellites(i).Name);
    fprintf('  Perigee: %.2f km, Apogee: %.2f km, Mean: %.2f km\n', ...
        satellites(i).Perigee_km, satellites(i).Apogee_km, satellites(i).MeanAltitude_km);
end

% Count how many have valid altitudes
valid_altitudes = sum(~isnan([satellites.MeanAltitude_km]));
fprintf('\nSatellites with valid altitudes: %d out of %d\n', valid_altitudes, length(satellites));

% Show altitude range
altitudes = [satellites.MeanAltitude_km];
altitudes = altitudes(~isnan(altitudes));
if ~isempty(altitudes)
    fprintf('Altitude range: %.2f to %.2f km\n', min(altitudes), max(altitudes));
end

%% Define orbital bins
% LEO: 160-2000 km in 20 km bands
LEO_edges = 160:20:2000;
LEO_labels = cell(length(LEO_edges)-1, 1);
for i = 1:length(LEO_edges)-1
    LEO_labels{i} = sprintf('LEO_%d_%d', LEO_edges(i), LEO_edges(i+1));
end

% MEO: 2000-35586 km in 200 km bands
MEO_edges = 2000:200:35800;
MEO_labels = cell(length(MEO_edges)-1, 1);
for i = 1:length(MEO_edges)-1
    MEO_labels{i} = sprintf('MEO_%d_%d', MEO_edges(i), MEO_edges(i+1));
end

% GEO: 35586-35986 km (single bin)
GEO_edges = [35586, 35986];
GEO_labels = {'GEO_35586_35986'};

% HEO: >35986 km in 1000 km bands
HEO_edges = 35986:1000:100000;
HEO_labels = cell(length(HEO_edges)-1, 1);
for i = 1:length(HEO_edges)-1
    HEO_labels{i} = sprintf('HEO_%d_%d', HEO_edges(i), HEO_edges(i+1));
end

% Combine all bins
all_edges = [LEO_edges, MEO_edges(2:end), GEO_edges(2:end), HEO_edges(2:end)];
all_labels = [LEO_labels; MEO_labels; GEO_labels; HEO_labels];

%% Create binned structure
orbital_bins = struct();
for i = 1:length(all_labels)
    orbital_bins.(all_labels{i}) = [];
end

%% Sort satellites into bins
binned_count = 0;
skipped_count = 0;

for i = 1:length(satellites)
    alt = satellites(i).MeanAltitude_km;
    
    if isnan(alt)
        skipped_count = skipped_count + 1;
        continue;
    end
    
    bin_idx = discretize(alt, all_edges);
    
    if ~isnan(bin_idx) && bin_idx <= length(all_labels)
        bin_name = all_labels{bin_idx};
        orbital_bins.(bin_name) = [orbital_bins.(bin_name), satellites(i)];
        binned_count = binned_count + 1;
    end
end

fprintf('\nBinning results:\n');
fprintf('  Successfully binned: %d\n', binned_count);
fprintf('  Skipped (NaN altitude): %d\n', skipped_count);

%% Create summary statistics for each bin
bin_summary = struct();
bin_names = fieldnames(orbital_bins);

for i = 1:length(bin_names)
    bin_name = bin_names{i};
    bin_sats = orbital_bins.(bin_name);
    
    bin_summary.(bin_name).count = length(bin_sats);
    
    if length(bin_sats) > 0
        masses = [bin_sats.LaunchMass_kg];
        masses = masses(~isnan(masses));
        
        bin_summary.(bin_name).total_mass = sum(masses);
        bin_summary.(bin_name).mean_mass = mean(masses);
        bin_summary.(bin_name).median_mass = median(masses);
    else
        bin_summary.(bin_name).total_mass = 0;
        bin_summary.(bin_name).mean_mass = 0;
        bin_summary.(bin_name).median_mass = 0;
    end
end

%% Save results
save('satellite_bins.mat', 'orbital_bins', 'bin_summary', 'satellites');

%% Export summary to Excel
num_bins = length(bin_names);
bin_name_array = cell(num_bins, 1);
count_array = zeros(num_bins, 1);
total_mass_array = zeros(num_bins, 1);
mean_mass_array = zeros(num_bins, 1);
median_mass_array = zeros(num_bins, 1);

for i = 1:num_bins
    bin_name = bin_names{i};
    bin_name_array{i} = bin_name;
    count_array(i) = bin_summary.(bin_name).count;
    total_mass_array(i) = bin_summary.(bin_name).total_mass;
    mean_mass_array(i) = bin_summary.(bin_name).mean_mass;
    median_mass_array(i) = bin_summary.(bin_name).median_mass;
end

summary_table = table(bin_name_array, count_array, total_mass_array, ...
    mean_mass_array, median_mass_array, ...
    'VariableNames', {'BinName', 'Count', 'TotalMass_kg', 'MeanMass_kg', 'MedianMass_kg'});

writetable(summary_table, 'orbital_bins_summary.xlsx');

%% Create line plot
counts = summary_table.Count;
bin_labels = summary_table.BinName;

non_empty_idx = counts > 0;
counts_filtered = counts(non_empty_idx);
bin_labels_filtered = bin_labels(non_empty_idx);

if isempty(counts_filtered)
    fprintf('No satellites found in any bins!\n');
    return;
end

altitude_centers = zeros(length(bin_labels_filtered), 1);
for i = 1:length(bin_labels_filtered)
    parts = split(bin_labels_filtered{i}, '_');
    if length(parts) >= 3
        altitude_centers(i) = (str2double(parts{2}) + str2double(parts{3})) / 2;
    end
end

figure('Position', [100, 100, 1400, 600]);
plot(altitude_centers, counts_filtered, 'b-', 'LineWidth', 2);
hold on;
plot(altitude_centers, counts_filtered, 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
title('Number of Satellites by Altitude (LEO to GEO)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Altitude (km)', 'FontSize', 12);
ylabel('Number of Satellites', 'FontSize', 12);
grid on;

% Set x-axis to only show up to end of GEO
xlim([0 36000]);  % This line limits the plot view

ymax = max(counts_filtered);
if ymax > 0
    ymax = ymax * 1.1;
    ylim([0 ymax]);
    plot([2000 2000], [0 ymax], 'r--', 'LineWidth', 2);
    plot([35586 35586], [0 ymax], 'g--', 'LineWidth', 2);
    text(1000, ymax*0.95, 'LEO', 'FontSize', 12, 'FontWeight', 'bold');
    text(18000, ymax*0.95, 'MEO', 'FontSize', 12, 'FontWeight', 'bold');
    text(35700, ymax*0.95, 'GEO', 'FontSize', 12, 'FontWeight', 'bold');
end
hold off;

fprintf('\nChart created! Check the figure window.\n');
fprintf('Summary saved to orbital_bins_summary.xlsx\n');
fprintf('Data saved to satellite_bins.mat\n');

%% Price
fprintf('\nPopulating prices...\n');

% Heuristic 1: Starlink satellites
starlink_count = 0;
for i = 1:length(satellites)
    % Check if name contains "Starlink" (case-insensitive)
    if contains(satellites(i).Name, 'Starlink', 'IgnoreCase', true)
        mass = satellites(i).LaunchMass_kg;
        
        % Check mass and assign price
        if ~isnan(mass)
            if ismember(mass, [227, 260, 290, 300])
                satellites(i).Price = 800000;  % $800k
                starlink_count = starlink_count + 1;
            elseif ismember(mass, [750, 800])
                satellites(i).Price = 2000000;  % $2M
                starlink_count = starlink_count + 1;
                %fprintf('\nExpensive Starlink');
            end
        end
    end
end

fprintf('  Priced %d Starlink satellites\n', starlink_count);

% Heuristic 1b: OneWeb satellites (all altitudes)
oneweb_count = 0;
for i = 1:length(satellites)
    % Check if name contains "OneWeb" (case-insensitive)
    if contains(satellites(i).Name, 'OneWeb', 'IgnoreCase', true)
        satellites(i).Price = 2400000;  % $2.4M
        oneweb_count = oneweb_count + 1;
    end
end

fprintf('  Priced %d OneWeb satellites\n', oneweb_count);



% LEO
% Power law based on:
%
% $90,000 for 2kg for launch (SATCatalogue)
%
% $275,000 for 50kg for launch (SATCatalogue)
%
% $13,000,000 for 500 kg (IDA)
%
% $50,000,000 for 1000 kg (McKinsey)
%
% $125,000,000 for 2500 kg (McKinsey)
%
% 
%
%
% GEO
%
% $500,000,000 for 3000 kg (ARSAT-1)
%
% $700,000,000 for 4000 kg (MIT)
%
% $800,000,000 for 6000 kg (ARABSAT-6A)
%
% $ Science sats cost way, way more
%
%
%
%
%
% Heuristic 2: LEO satellites with piecewise mass function
leo_count = 0;
for i = 1:length(satellites)
    % Only apply to LEO satellites that don't already have a price
    % Exclude both Starlink and OneWeb from this heuristic
    if strcmp(satellites(i).ClassOfOrbit, 'LEO') && isnan(satellites(i).Price) && ...
       ~contains(satellites(i).Name, 'Starlink', 'IgnoreCase', true) && ...
       ~contains(satellites(i).Name, 'OneWeb', 'IgnoreCase', true)
        
        mass = satellites(i).LaunchMass_kg;
        
        if ~isnan(mass) && mass >= 1
            if mass < 50
                satellites(i).Price = 76900 * mass^0.243;
                leo_count = leo_count + 1;
            elseif mass < 500
                satellites(i).Price = 535 * mass^1.674;
                leo_count = leo_count + 1;
            elseif mass < 1000
                % FIXED: More reasonable pricing for mid-size LEO satellites
                satellites(i).Price = mass * 80000;  % $80k per kg
                leo_count = leo_count + 1;
            elseif mass < 2500
                satellites(i).Price = 50000 * mass;  % $50k per kg
                leo_count = leo_count + 1;
            end
        end
    end
end

fprintf('  Priced %d LEO satellites using mass heuristic\n', leo_count);
%% Populate Price field for GEO satellites
fprintf('\nPricing GEO satellites...\n');

geo_count = 0;
for i = 1:length(satellites)
    % Only price satellites in GEO that don't already have a price
    if strcmp(satellites(i).ClassOfOrbit, 'GEO') && isnan(satellites(i).Price)
        satellites(i).Price = GEOsatPricer(satellites(i));
        if ~isnan(satellites(i).Price)
            geo_count = geo_count + 1;
        end
    end
end

fprintf('  Priced %d GEO satellites\n', geo_count);

% Count total satellites with prices
priced_count = sum(~isnan([satellites.Price]));
fprintf('Total satellites with prices: %d out of %d\n', priced_count, length(satellites));

% Count total satellites with prices
priced_count = sum(~isnan([satellites.Price]));
fprintf('Total satellites with prices: %d out of %d\n', priced_count, length(satellites));

%% Plot satellite prices by altitude (color-coded by orbit class)

altitudes_for_price_plot = [satellites.MeanAltitude_km];
prices_for_price_plot = [satellites.Price];

valid_idx = ~isnan(altitudes_for_price_plot) & ~isnan(prices_for_price_plot);
altitudes_valid = altitudes_for_price_plot(valid_idx);
prices_valid = prices_for_price_plot(valid_idx);

% Get orbit classes for coloring
orbit_classes = {satellites(valid_idx).ClassOfOrbit};

if ~isempty(altitudes_valid)
    figure('Position', [100, 100, 1400, 600]);
    
    % Plot each orbit class with different color
    hold on;
    leo_idx = strcmp(orbit_classes, 'LEO');
    meo_idx = strcmp(orbit_classes, 'MEO');
    geo_idx = strcmp(orbit_classes, 'GEO');
    
    if any(leo_idx)
        scatter(altitudes_valid(leo_idx), prices_valid(leo_idx), 20, 'b', 'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'LEO');
    end
    if any(meo_idx)
        scatter(altitudes_valid(meo_idx), prices_valid(meo_idx), 20, 'r', 'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'MEO');
    end
    if any(geo_idx)
        scatter(altitudes_valid(geo_idx), prices_valid(geo_idx), 20, 'g', 'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'GEO');
    end
    
    title('Satellite Prices by Altitude', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Altitude (km)', 'FontSize', 12);
    ylabel('Price ($) - log scale', 'FontSize', 12);
    set(gca, 'YScale', 'log');
    grid on;
    xlim([0 1500]);
    legend('Location', 'best');
    hold off;
    
    fprintf('\nPrice vs Altitude plot created!\n');
    fprintf('Plotted %d satellites with both altitude and price data\n', length(altitudes_valid));
end



%% Plot average mass by altitude with log x-axis

% Use bins
avg_mass_by_bin = zeros(length(bin_labels_filtered), 1);
altitude_centers_for_mass = zeros(length(bin_labels_filtered), 1);

for i = 1:length(bin_labels_filtered)
    bin_name = bin_labels_filtered{i};
    bin_sats = orbital_bins.(bin_name);
    
    if length(bin_sats) > 0
        masses = [bin_sats.LaunchMass_kg];
        masses = masses(~isnan(masses));
        
        if ~isempty(masses)
            avg_mass_by_bin(i) = mean(masses);
        else
            avg_mass_by_bin(i) = NaN;
        end
    else
        avg_mass_by_bin(i) = NaN;
    end
    
    % Get altitude center
    parts = split(bin_labels_filtered{i}, '_');
    if length(parts) >= 3
        altitude_centers_for_mass(i) = (str2double(parts{2}) + str2double(parts{3})) / 2;
    end
end

% Remove NaN values
valid_mass_idx = ~isnan(avg_mass_by_bin);
altitude_centers_mass = altitude_centers_for_mass(valid_mass_idx);
avg_mass_filtered = avg_mass_by_bin(valid_mass_idx);

if ~isempty(avg_mass_filtered)
    figure('Position', [100, 100, 1400, 600]);
    plot(altitude_centers_mass, avg_mass_filtered, 'r-', 'LineWidth', 2);
    hold on;
    plot(altitude_centers_mass, avg_mass_filtered, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
    
    title('Average Satellite Mass by Altitude', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Altitude (km) - log scale', 'FontSize', 12);
    ylabel('Average Launch Mass (kg)', 'FontSize', 12);
    grid on;
    
    % Set x-axis to log scale
    set(gca, 'XScale', 'log');
    xlim([200 36000]);  % Start at 200 to avoid log(0) issues
    
    % Add orbit boundary lines
    yrange = ylim;
    plot([2000 2000], yrange, 'k--', 'LineWidth', 2);
    plot([35586 35586], yrange, 'k--', 'LineWidth', 2);
    
    % Adjust text positions for log scale
    text(500, yrange(2)*0.95, 'LEO', 'FontSize', 12, 'FontWeight', 'bold');
    text(10000, yrange(2)*0.95, 'MEO', 'FontSize', 12, 'FontWeight', 'bold');
    text(36000, yrange(2)*0.95, 'GEO', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    hold off;
    
    fprintf('\nAverage mass by altitude plot created!\n');
else
    fprintf('\nNo mass data available to plot\n');
end

%% Plot average prices for specific altitude bands

% Define the bands we want to analyze
band_500_600_starlink = [];
band_500_600_non_starlink = [];
band_oneweb_low = [];  % 500-560 km
band_oneweb_mid = [];  % 590-733 km
band_oneweb_high = []; % 977-1215 km
band_1200_1400 = [];
band_GEO = [];

% Collect satellites in each band
for i = 1:length(satellites)
    alt = satellites(i).MeanAltitude_km;
    price = satellites(i).Price;
    
    % Skip if no valid altitude or price
    if isnan(alt) || isnan(price)
        continue;
    end
    
    % OneWeb 500-560 km band
    if alt >= 500 && alt <= 560 && contains(satellites(i).Name, 'OneWeb', 'IgnoreCase', true)
        band_oneweb_low = [band_oneweb_low, price];
    end
    
    % OneWeb 590-733 km band
    if alt >= 590 && alt <= 733 && contains(satellites(i).Name, 'OneWeb', 'IgnoreCase', true)
        band_oneweb_mid = [band_oneweb_mid, price];
    end
    
    % OneWeb 977-1215 km band
    if alt >= 977 && alt <= 1215 && contains(satellites(i).Name, 'OneWeb', 'IgnoreCase', true)
        band_oneweb_high = [band_oneweb_high, price];
    end
    
    % 500-600 km band (Starlink vs non-Starlink, excluding OneWeb counted above)
    if alt >= 500 && alt <= 600
        if contains(satellites(i).Name, 'Starlink', 'IgnoreCase', true)
            band_500_600_starlink = [band_500_600_starlink, price];
        elseif ~contains(satellites(i).Name, 'OneWeb', 'IgnoreCase', true)
            band_500_600_non_starlink = [band_500_600_non_starlink, price];
        end
    end
    
    % 1200-1400 km band (excluding OneWeb)
    if alt >= 1200 && alt <= 1400 && ~contains(satellites(i).Name, 'OneWeb', 'IgnoreCase', true)
        band_1200_1400 = [band_1200_1400, price];
    end
    
    % GEO band (35586-35986 km)
    if alt >= 35586 && alt <= 35986
        band_GEO = [band_GEO, price];
    end
end

% Calculate averages
avg_500_600_starlink = mean(band_500_600_starlink);
avg_500_600_non_starlink = mean(band_500_600_non_starlink);
avg_oneweb_low = mean(band_oneweb_low);
avg_oneweb_mid = mean(band_oneweb_mid);
avg_oneweb_high = mean(band_oneweb_high);
avg_1200_1400 = mean(band_1200_1400);
avg_GEO = mean(band_GEO);

% Prepare data for plotting
categories = {'500-600 km\n(Starlink)', '500-600 km\n(Non-Starlink)', ...
              '500-560 km\n(OneWeb)', '590-733 km\n(OneWeb)', '977-1215 km\n(OneWeb)', ...
              '1200-1400 km', 'GEO'};
averages = [avg_500_600_starlink, avg_500_600_non_starlink, ...
            avg_oneweb_low, avg_oneweb_mid, avg_oneweb_high, ...
            avg_1200_1400, avg_GEO];
counts = [length(band_500_600_starlink), length(band_500_600_non_starlink), ...
          length(band_oneweb_low), length(band_oneweb_mid), length(band_oneweb_high), ...
          length(band_1200_1400), length(band_GEO)];

% Remove NaN values (bands with no data)
valid_idx = ~isnan(averages);
categories = categories(valid_idx);
averages = averages(valid_idx);
counts = counts(valid_idx);

if ~isempty(averages)
    figure('Position', [100, 100, 1400, 600]);
    
    % Create bar chart
    bar_handle = bar(averages, 'FaceColor', 'flat');
    
    % Color code the bars
    colors = [0.2 0.4 0.8;   % Blue for Starlink (500-600)
              0.8 0.4 0.2;   % Orange for Non-Starlink (500-600)
              0.3 0.7 0.3;   % Green for OneWeb (500-560)
              0.2 0.6 0.2;   % Medium Green for OneWeb (590-733)
              0.1 0.5 0.1;   % Dark Green for OneWeb (977-1215)
              0.6 0.4 0.7;   % Purple for 1200-1400
              0.9 0.3 0.3];  % Red for GEO
    bar_handle.CData = colors(valid_idx, :);
    
    % Labels and formatting
    title('Average Satellite Price by Altitude Band (Log Scale)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Average Price ($)', 'FontSize', 12);
    set(gca, 'XTickLabel', categories);
    set(gca, 'YScale', 'log');  % Set logarithmic scale
    xtickangle(20);
    grid on;
    
    % Add value labels on top of bars
    hold on;
    for i = 1:length(averages)
        if averages(i) >= 1e6
            label_text = sprintf('$%.2fM\n(n=%d)', averages(i)/1e6, counts(i));
        else
            label_text = sprintf('$%.0fk\n(n=%d)', averages(i)/1e3, counts(i));
        end
        text(i, averages(i), label_text, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
    end
    hold off;
    
    fprintf('\nAverage prices by altitude band:\n');
    for i = 1:length(categories)
        fprintf('  %s: $%.2f M (%d satellites)\n', ...
            strrep(categories{i}, sprintf('\n'), ' '), averages(i)/1e6, counts(i));
    end
else
    fprintf('\nNo price data available for the specified bands\n');
end





%% Debug: Check what's in the 1200-1400 band
fprintf('\nDebugging 1200-1400 km band...\n');

debug_1200_1400 = [];
for i = 1:length(satellites)
    alt = satellites(i).MeanAltitude_km;
    price = satellites(i).Price;
    
    if ~isnan(alt) && ~isnan(price) && alt >= 1200 && alt <= 1400 && ...
       ~contains(satellites(i).Name, 'OneWeb', 'IgnoreCase', true)
        debug_1200_1400 = [debug_1200_1400; i];
        fprintf('  %s: Alt=%.1f km, Mass=%.0f kg, Price=$%.2fM, Purpose=%s, Operator=%s\n', ...
            satellites(i).Name, alt, satellites(i).LaunchMass_kg, price/1e6, ...
            satellites(i).Purpose, satellites(i).Operator);
    end
end

fprintf('\nTotal in 1200-1400 band (non-OneWeb): %d\n', length(debug_1200_1400));
fprintf('Average price: $%.2fM\n', mean([satellites(debug_1200_1400).Price])/1e6);
fprintf('Median price: $%.2fM\n', median([satellites(debug_1200_1400).Price])/1e6);

%% Helper function to convert values to numeric
function num = convertToNumeric(val)
    if isnumeric(val)
        num = double(val);
    elseif ischar(val) || isstring(val)
        num = str2double(val);
    else
        num = NaN;
    end
end