%% Constellation Servicing Delta-V Analysis
% Calculates the total delta-V cost of visiting every satellite in a constellation
% Uses optimal phasing maneuvers within planes and plane changes between planes
%
% Author: Generated for OSCAR@VT analysis
% Date: January 2025

clear; clc; close all;

% Start logging all command window output to file
diary('constellation_output.txt');
diary on;

fprintf('=======================================================================\n');
fprintf('CONSTELLATION SERVICING DELTA-V ANALYSIS\n');
fprintf('Run date: %s\n', datestr(now));
fprintf('=======================================================================\n\n');

%% ===================== USER CONFIGURATION =====================
% Time constraint for visiting entire constellation (seconds)
% Set to Inf for unconstrained (minimum delta-V regardless of time)
TOTAL_TIME_CONSTRAINT = 30 * 24 * 3600;  % 30 days in seconds

% If total time is too constrained, fall back to per-phase constraints
TIME_PER_PLANE_SHIFT = 2 * 24 * 3600;    % 2 days for plane changes
TIME_PER_SATELLITE = 0.5 * 24 * 3600;    % 12 hours per satellite phasing

% RAAN tolerance for grouping satellites into planes (degrees)
% GPS planes are nominally 60 deg apart, but drift occurs
% Use ~15 deg to capture satellites in the same nominal plane
RAAN_TOLERANCE = 15.0;

% Earth parameters
MU_EARTH = 3.986004418e14;  % Earth gravitational parameter (m^3/s^2)
R_EARTH = 6378.137e3;       % Earth equatorial radius (m)

% Select constellation to analyze
% Options: 'gpsLegacy', 'gpsDecomissioned', 'iridiumCurrent', ...
                  %'iridiumDecomAlive', 'globalstar', 'goes', 'wgs', ...
                  %'tdrss', 'jpss', 'poesItos'
CONSTELLATION = 'gpsLegacy';

%% ===================== LOAD SATELLITE DATA =====================
csvFile = [CONSTELLATION, '.csv'];
fprintf('Loading data from %s...\n', csvFile);

satellites = loadJCATFormat(csvFile);
fprintf('Loaded %d satellites from %s\n\n', length(satellites), CONSTELLATION);

% Display satellite summary (will show detailed list after plane grouping)
fprintf('Processing satellites...\n\n');

%% ===================== GROUP INTO ORBITAL PLANES =====================
planes = groupIntoPlanes(satellites, RAAN_TOLERANCE);
fprintf('Identified %d orbital planes:\n', length(planes));
fprintf('%s\n', repmat('=', 1, 85));

for i = 1:length(planes)
    fprintf('\nPLANE %d: RAAN = %.2f deg, Inc = %.2f deg (%d satellites)\n', ...
        i, planes(i).raan, planes(i).inc, length(planes(i).satellites));
    fprintf('%s\n', repmat('-', 1, 85));
    fprintf('  %-35s %12s %12s %12s\n', 'Satellite Name', 'Alt (km)', 'Inc (deg)', 'Ecc');
    fprintf('  %s\n', repmat('-', 1, 80));
    
    for j = 1:length(planes(i).satellites)
        sat = planes(i).satellites(j);
        alt = (sat.a / 1e3) - 6378.137;
        fprintf('  %-35s %12.1f %12.2f %12.6f\n', ...
            sat.name, alt, sat.i, sat.e);
    end
end
fprintf('\n%s\n\n', repmat('=', 1, 85));

%% ===================== OPTIMIZE PLANE VISIT ORDER =====================
% Use nearest-neighbor heuristic for plane ordering (minimizes plane change dV)
planeOrder = optimizePlaneOrder(planes);
fprintf('Optimal plane visit order: ');
fprintf('%d ', planeOrder);
fprintf('\n\n');

%% ===================== CALCULATE DELTA-V =====================
% Distribute time budget
nPlanes = length(planes);
nSatellites = length(satellites);
nPlaneChanges = nPlanes - 1;

if ~isinf(TOTAL_TIME_CONSTRAINT)
    % Allocate time proportionally
    timeForPhasing = TOTAL_TIME_CONSTRAINT * 0.7;  % 70% for phasing
    timeForPlaneChanges = TOTAL_TIME_CONSTRAINT * 0.3;  % 30% for plane changes
    
    timePerSatPhasing = timeForPhasing / max(nSatellites - nPlanes, 1);
    timePerPlaneChange = timeForPlaneChanges / max(nPlaneChanges, 1);
else
    timePerSatPhasing = TIME_PER_SATELLITE;
    timePerPlaneChange = TIME_PER_PLANE_SHIFT;
end

fprintf('Time budget allocation:\n');
fprintf('  Time per satellite phasing: %.2f hours\n', timePerSatPhasing/3600);
fprintf('  Time per plane change: %.2f hours\n', timePerPlaneChange/3600);
fprintf('\n');

% Calculate delta-V for each segment
results = struct();
results.phasingDV = zeros(nPlanes, 1);
results.planeChangeDV = zeros(nPlanes - 1, 1);
results.planeDetails = cell(nPlanes, 1);

totalPhasingDV = 0;
totalPlaneChangeDV = 0;

for i = 1:nPlanes
    planeIdx = planeOrder(i);
    plane = planes(planeIdx);
    
    % Calculate phasing delta-V within this plane
    [phasingDV, visitOrder, phasingDetails] = calculatePhasingDV(...
        plane.satellites, timePerSatPhasing, MU_EARTH);
    
    results.phasingDV(i) = phasingDV;
    results.planeDetails{i} = phasingDetails;
    totalPhasingDV = totalPhasingDV + phasingDV;
    
    fprintf('Plane %d (RAAN = %.2f deg, Inc = %.2f deg):\n', planeIdx, plane.raan, plane.inc);
    fprintf('  Satellites visited: %d\n', length(plane.satellites));
    fprintf('  Phasing delta-V: %.2f m/s\n', phasingDV);
    
    % Calculate plane change delta-V to next plane
    if i < nPlanes
        nextPlaneIdx = planeOrder(i + 1);
        nextPlane = planes(nextPlaneIdx);
        
        % Get representative orbital elements from current and next plane
        currentSat = plane.satellites(end);  % Last satellite visited
        nextSat = nextPlane.satellites(1);   % First satellite in next plane
        
        planeChangeDV = calculatePlaneChangeDV(currentSat, nextSat, MU_EARTH);
        results.planeChangeDV(i) = planeChangeDV;
        totalPlaneChangeDV = totalPlaneChangeDV + planeChangeDV;
        
        fprintf('  Plane change to next: %.2f m/s\n', planeChangeDV);
    end
    fprintf('\n');
end

%% ===================== SUMMARY RESULTS =====================
totalDV = totalPhasingDV + totalPlaneChangeDV;

fprintf('========================================\n');
fprintf('DELTA-V SUMMARY FOR %s\n', upper(CONSTELLATION));
fprintf('========================================\n');
fprintf('Total satellites:           %8d\n', nSatellites);
fprintf('Total orbital planes:       %8d\n', nPlanes);
fprintf('Total phasing delta-V:      %8.2f m/s\n', totalPhasingDV);
fprintf('Total plane change delta-V: %8.2f m/s\n', totalPlaneChangeDV);
fprintf('----------------------------------------\n');
fprintf('TOTAL DELTA-V:              %8.2f m/s (%.2f km/s)\n', totalDV, totalDV/1000);
fprintf('========================================\n');

% Store results
results.constellation = CONSTELLATION;
results.totalPhasingDV = totalPhasingDV;
results.totalPlaneChangeDV = totalPlaneChangeDV;
results.totalDV = totalDV;
results.nSatellites = nSatellites;
results.nPlanes = nPlanes;
results.planeOrder = planeOrder;

%% ===================== SAVE OUTPUTS =====================
% Save figure
figFileName = sprintf('%s_analysis.png', CONSTELLATION);
figure('Position', [100, 100, 800, 600]);
bar([totalPhasingDV, totalPlaneChangeDV, totalDV] / 1000);
set(gca, 'XTickLabel', {'Phasing', 'Plane Change', 'Total'});
ylabel('Delta-V (km/s)');
title(sprintf('%s Constellation Delta-V Breakdown', upper(CONSTELLATION)));
grid on;
saveas(gcf, figFileName);
fprintf('\nFigure saved as %s\n', figFileName);

% Stop logging
fprintf('\n=======================================================================\n');
fprintf('Analysis complete. Output saved to constellation_output.txt\n');
fprintf('=======================================================================\n');
diary off;

%% ===================== FUNCTION DEFINITIONS =====================

function satellites = loadJCATFormat(filename)
    % Load satellite data from JCAT format CSV file
    % Fetches accurate RAAN, AoP, and MA from CelesTrak TLEs using NORAD ID
    %
    % The file uses DispPeri/DispApo (km altitude) and DispInc (deg) for orbit params
    % RAAN and other elements are fetched from TLEs for accuracy
    
    % Read the CSV file
    opts = detectImportOptions(filename);
    opts.VariableNamingRule = 'preserve';
    
    % Suppress datetime ambiguity warnings by treating date columns as text
    for i = 1:length(opts.VariableTypes)
        if strcmp(opts.VariableTypes{i}, 'datetime')
            opts.VariableTypes{i} = 'char';
        end
    end
    
    data = readtable(filename, opts);
    
    % Get column names (handle potential variations)
    colNames = data.Properties.VariableNames;
    
    % Find relevant columns
    nameCol = findColumn(colNames, {'Name', 'name'});
    jcatCol = findColumn(colNames, {'#JCAT', 'JCAT', 'NORAD', 'CatalogNumber'});
    
    % Try DispPeri/DispApo first, fall back to UNPerigee/UNApogee
    periCol = findColumn(colNames, {'DispPeri', 'UNPerigee'});
    apoCol = findColumn(colNames, {'DispApo', 'UNApogee'});
    incCol = findColumn(colNames, {'DispInc', 'UNInc'});
    
    % Check for RAAN column (may exist if user added it)
    raanCol = findColumn(colNames, {'RAAN', 'raan', 'DispRAAN', 'UNRaan'});
    
    nSats = height(data);
    satellites = struct('name', cell(1, nSats), 'a', cell(1, nSats), ...
                        'e', cell(1, nSats), 'i', cell(1, nSats), ...
                        'raan', cell(1, nSats), 'aop', cell(1, nSats), ...
                        'ma', cell(1, nSats), 'norad', cell(1, nSats));
    
    % Collect NORAD IDs for batch TLE fetch
    noradIDs = [];
    validRows = [];
    
    fprintf('Extracting NORAD IDs from catalog...\n');
    
    for k = 1:nSats
        % Get perigee and apogee altitudes (km)
        peri = getNumericValue(data, k, periCol);
        apo = getNumericValue(data, k, apoCol);
        inc = getNumericValue(data, k, incCol);
        
        % Skip invalid entries
        if isnan(peri) || isnan(apo) || isnan(inc) || peri <= 0 || apo <= 0
            continue;
        end
        
        % Get NORAD ID from JCAT column (format: S##### or just #####)
        if ~isempty(jcatCol)
            jcatVal = data{k, jcatCol};
            if iscell(jcatVal)
                jcatStr = jcatVal{1};
            else
                jcatStr = char(jcatVal);
            end
            noradID = extractNoradID(jcatStr);
        else
            noradID = NaN;
        end
        
        if ~isnan(noradID)
            noradIDs(end+1) = noradID;
            validRows(end+1) = k;
        end
    end
    
    % Fetch TLEs from CelesTrak
    fprintf('Fetching TLEs from CelesTrak for %d satellites...\n', length(noradIDs));
    tleData = fetchTLEsFromCelesTrak(noradIDs);
    
    % Now build satellite structures with TLE data
    validCount = 0;
    
    for idx = 1:length(validRows)
        k = validRows(idx);
        noradID = noradIDs(idx);
        
        validCount = validCount + 1;
        
        % Get satellite name
        if ~isempty(nameCol)
            nameVal = data{k, nameCol};
            if iscell(nameVal)
                satellites(validCount).name = nameVal{1};
            else
                satellites(validCount).name = char(nameVal);
            end
        else
            satellites(validCount).name = sprintf('SAT-%d', noradID);
        end
        
        % Clean up name (remove brackets)
        satellites(validCount).name = strrep(satellites(validCount).name, '[', '');
        satellites(validCount).name = strrep(satellites(validCount).name, ']', '');
        
        % Store NORAD ID
        satellites(validCount).norad = noradID;
        
        % Check if we have TLE data for this satellite
        fieldName = sprintf('n%d', noradID);
        if isfield(tleData, fieldName)
            tle = tleData.(fieldName);
            
            % Use TLE orbital elements (more accurate)
            satellites(validCount).a = tle.a;
            satellites(validCount).e = tle.e;
            satellites(validCount).i = tle.i;
            satellites(validCount).raan = tle.raan;
            satellites(validCount).aop = tle.aop;
            satellites(validCount).ma = tle.ma;
            
            % Debug: show first few satellites
            if validCount <= 5
                fprintf('  Applied TLE for NORAD %d: RAAN=%.2f\n', noradID, tle.raan);
            end
        else
            fprintf('  WARNING: No TLE found for NORAD %d (%s)\n', noradID, satellites(validCount).name);
            % Fall back to CSV data
            peri = getNumericValue(data, k, periCol);
            apo = getNumericValue(data, k, apoCol);
            inc = getNumericValue(data, k, incCol);
            
            R_EARTH_KM = 6378.137;
            r_peri = (peri + R_EARTH_KM) * 1e3;
            r_apo = (apo + R_EARTH_KM) * 1e3;
            
            satellites(validCount).a = (r_peri + r_apo) / 2;
            satellites(validCount).e = (r_apo - r_peri) / (r_apo + r_peri);
            satellites(validCount).i = inc;
            
            % Check for RAAN in CSV
            if ~isempty(raanCol)
                raan = getNumericValue(data, k, raanCol);
                satellites(validCount).raan = raan;
            else
                satellites(validCount).raan = 0;
                fprintf('  WARNING: No TLE found for NORAD %d, RAAN set to 0\n', noradID);
            end
            
            satellites(validCount).aop = 0;
            satellites(validCount).ma = mod(validCount * 45, 360);
        end
    end
    
    % Trim to valid satellites
    satellites = satellites(1:validCount);
    
    fprintf('Successfully loaded %d satellites with orbital elements\n', validCount);
end

function noradID = extractNoradID(jcatStr)
    % Extract NORAD catalog number from JCAT string (e.g., "S36585" -> 36585)
    jcatStr = strtrim(char(jcatStr));
    
    % Remove leading 'S' if present
    if ~isempty(jcatStr) && upper(jcatStr(1)) == 'S'
        jcatStr = jcatStr(2:end);
    end
    
    % Extract numeric part
    numStr = regexp(jcatStr, '\d+', 'match', 'once');
    if ~isempty(numStr)
        noradID = str2double(numStr);
    else
        noradID = NaN;
    end
end

function tleData = fetchTLEsFromCelesTrak(noradIDs)
    % Fetch TLEs from CelesTrak for given NORAD IDs
    % Uses local cache file to avoid hitting CelesTrak rate limits
    % Cache is updated once per day
    
    cacheFile = 'tle_cache.txt';
    tleData = struct();
    nSats = length(noradIDs);
    
    % Check if cache exists and is from today
    cacheValid = false;
    cachedTLEs = struct();
    
    if isfile(cacheFile)
        fid = fopen(cacheFile, 'r');
        if fid ~= -1
            firstLine = fgetl(fid);
            fclose(fid);
            
            % Check date in first line (format: "CACHE_DATE: YYYY-MM-DD")
            if startsWith(firstLine, 'CACHE_DATE:')
                cacheDateStr = strtrim(firstLine(12:end));
                todayStr = datestr(now, 'yyyy-mm-dd');
                
                if strcmp(cacheDateStr, todayStr)
                    cacheValid = true;
                    fprintf('TLE cache is current (dated %s). Loading from cache...\n', cacheDateStr);
                    cachedTLEs = loadTLECache(cacheFile);
                    fprintf('Loaded %d TLEs from cache.\n', length(fieldnames(cachedTLEs)));
                else
                    fprintf('TLE cache is outdated (cache: %s, today: %s). Will fetch new TLEs.\n', ...
                        cacheDateStr, todayStr);
                end
            end
        end
    else
        fprintf('No TLE cache file found. Will fetch TLEs from CelesTrak.\n');
    end
    
    % Determine which satellites need to be fetched
    if cacheValid
        missingIDs = [];
        for i = 1:nSats
            fieldName = sprintf('n%d', noradIDs(i));
            if ~isfield(cachedTLEs, fieldName)
                missingIDs(end+1) = noradIDs(i);
            end
        end
        tleData = cachedTLEs;
        fprintf('%d of %d satellites found in cache. Need to fetch %d.\n', ...
            nSats - length(missingIDs), nSats, length(missingIDs));
    else
        missingIDs = noradIDs;
    end
    
    % Fetch missing TLEs from CelesTrak
    if ~isempty(missingIDs)
        fprintf('Fetching %d TLEs from CelesTrak...\n', length(missingIDs));
        
        options = weboptions('Timeout', 30);
        fetchedCount = 0;
        failedCount = 0;
        
        for i = 1:length(missingIDs)
            noradID = missingIDs(i);
            url = sprintf('https://celestrak.org/NORAD/elements/gp.php?CATNR=%d&FORMAT=TLE', noradID);
            
            try
                tleText = webread(url, options);
                lines = strsplit(tleText, '\n');
                
                if length(lines) >= 3
                    line0 = strtrim(lines{1});
                    line1 = strtrim(lines{2});
                    line2 = strtrim(lines{3});
                    
                    if length(line1) >= 69 && line1(1) == '1' && ...
                       length(line2) >= 69 && line2(1) == '2'
                        
                        inc = str2double(line2(9:16));
                        raan = str2double(line2(18:25));
                        eccStr = ['0.' line2(27:33)];
                        ecc = str2double(eccStr);
                        aop = str2double(line2(35:42));
                        ma = str2double(line2(44:51));
                        n = str2double(line2(53:63));
                        
                        mu = 3.986004418e14;
                        n_rad_s = n * 2 * pi / 86400;
                        a = (mu / n_rad_s^2)^(1/3);
                        
                        fieldName = sprintf('n%d', noradID);
                        tleData.(fieldName).a = a;
                        tleData.(fieldName).e = ecc;
                        tleData.(fieldName).i = inc;
                        tleData.(fieldName).raan = raan;
                        tleData.(fieldName).aop = aop;
                        tleData.(fieldName).ma = ma;
                        tleData.(fieldName).name = line0;
                        tleData.(fieldName).line1 = line1;
                        tleData.(fieldName).line2 = line2;
                        
                        fetchedCount = fetchedCount + 1;
                        
                        if mod(fetchedCount, 10) == 0
                            fprintf('  Fetched %d/%d TLEs...\n', fetchedCount, length(missingIDs));
                        end
                    else
                        failedCount = failedCount + 1;
                        fprintf('  Invalid TLE format for NORAD %d\n', noradID);
                    end
                else
                    failedCount = failedCount + 1;
                    fprintf('  No TLE data returned for NORAD %d\n', noradID);
                end
                
            catch ME
                failedCount = failedCount + 1;
                fprintf('  Failed to fetch NORAD %d: %s\n', noradID, ME.message);
            end
            
            if i < length(missingIDs)
                pause(0.2);
            end
        end
        
        fprintf('TLE fetch complete: %d successful, %d failed\n', fetchedCount, failedCount);
        
        % Save updated cache
        saveTLECache(cacheFile, tleData);
    end
end

function cachedTLEs = loadTLECache(cacheFile)
    % Load TLEs from cache file
    cachedTLEs = struct();
    
    fid = fopen(cacheFile, 'r');
    if fid == -1
        return;
    end
    
    % Skip first line (date)
    fgetl(fid);
    
    while ~feof(fid)
        line0 = fgetl(fid);  % Name line or separator
        if ~ischar(line0) || isempty(strtrim(line0)) || startsWith(line0, '---')
            continue;
        end
        
        line1 = fgetl(fid);  % TLE line 1
        line2 = fgetl(fid);  % TLE line 2
        
        if ~ischar(line1) || ~ischar(line2)
            break;
        end
        
        line0 = strtrim(line0);
        line1 = strtrim(line1);
        line2 = strtrim(line2);
        
        % Validate and parse TLE
        if length(line1) >= 69 && line1(1) == '1' && ...
           length(line2) >= 69 && line2(1) == '2'
            
            noradID = str2double(line1(3:7));
            
            inc = str2double(line2(9:16));
            raan = str2double(line2(18:25));
            eccStr = ['0.' line2(27:33)];
            ecc = str2double(eccStr);
            aop = str2double(line2(35:42));
            ma = str2double(line2(44:51));
            n = str2double(line2(53:63));
            
            mu = 3.986004418e14;
            n_rad_s = n * 2 * pi / 86400;
            a = (mu / n_rad_s^2)^(1/3);
            
            fieldName = sprintf('n%d', noradID);
            cachedTLEs.(fieldName).a = a;
            cachedTLEs.(fieldName).e = ecc;
            cachedTLEs.(fieldName).i = inc;
            cachedTLEs.(fieldName).raan = raan;
            cachedTLEs.(fieldName).aop = aop;
            cachedTLEs.(fieldName).ma = ma;
            cachedTLEs.(fieldName).name = line0;
            cachedTLEs.(fieldName).line1 = line1;
            cachedTLEs.(fieldName).line2 = line2;
        end
    end
    
    fclose(fid);
end

function saveTLECache(cacheFile, tleData)
    % Save TLEs to cache file
    fid = fopen(cacheFile, 'w');
    if fid == -1
        fprintf('WARNING: Could not save TLE cache file.\n');
        return;
    end
    
    % Write date header
    fprintf(fid, 'CACHE_DATE: %s\n', datestr(now, 'yyyy-mm-dd'));
    fprintf(fid, '-----------------------------------------------------------\n');
    
    % Write each TLE
    fields = fieldnames(tleData);
    for i = 1:length(fields)
        tle = tleData.(fields{i});
        if isfield(tle, 'line1') && isfield(tle, 'line2')
            fprintf(fid, '%s\n', tle.name);
            fprintf(fid, '%s\n', tle.line1);
            fprintf(fid, '%s\n', tle.line2);
        end
    end
    
    fclose(fid);
    fprintf('TLE cache saved to %s (%d TLEs)\n', cacheFile, length(fields));
end

function col = findColumn(colNames, possibleNames)
    % Find column index from list of possible names
    col = [];
    for i = 1:length(possibleNames)
        idx = find(strcmpi(colNames, possibleNames{i}), 1);
        if ~isempty(idx)
            col = idx;
            return;
        end
    end
end

function val = getNumericValue(data, row, col)
    % Safely extract numeric value from table
    if isempty(col)
        val = NaN;
        return;
    end
    
    rawVal = data{row, col};
    
    if iscell(rawVal)
        rawVal = rawVal{1};
    end
    
    if isnumeric(rawVal)
        val = rawVal;
    elseif ischar(rawVal) || isstring(rawVal)
        % Remove commas and whitespace, then convert
        rawVal = strtrim(char(rawVal));
        rawVal = strrep(rawVal, ',', '');
        rawVal = strrep(rawVal, ' ', '');
        if isempty(rawVal) || strcmp(rawVal, '-') || strcmp(rawVal, '*')
            val = NaN;
        else
            val = str2double(rawVal);
        end
    else
        val = NaN;
    end
end

function planes = groupIntoPlanes(satellites, raanTolerance)
    % Group satellites into orbital planes based on RAAN and inclination clustering
    
    nSats = length(satellites);
    raans = [satellites.raan];
    incs = [satellites.i];
    
    % Sort by RAAN
    [~, sortIdx] = sort(raans);
    
    planes = struct('raan', {}, 'inc', {}, 'satellites', {});
    assigned = false(1, nSats);
    
    for i = 1:nSats
        idx = sortIdx(i);
        if assigned(idx)
            continue;
        end
        
        % Start new plane
        planeIdx = length(planes) + 1;
        planes(planeIdx).raan = satellites(idx).raan;
        planes(planeIdx).inc = satellites(idx).i;
        planes(planeIdx).satellites = satellites(idx);
        assigned(idx) = true;
        
        % Find all satellites within tolerance (both RAAN and inclination)
        for j = i+1:nSats
            jdx = sortIdx(j);
            if assigned(jdx)
                continue;
            end
            
            raanDiff = abs(satellites(jdx).raan - planes(planeIdx).raan);
            raanDiff = min(raanDiff, 360 - raanDiff);  % Handle wraparound
            
            incDiff = abs(satellites(jdx).i - planes(planeIdx).inc);
            
            if raanDiff <= raanTolerance && incDiff <= 2.0  % 2 deg inc tolerance
                planes(planeIdx).satellites(end+1) = satellites(jdx);
                assigned(jdx) = true;
            end
        end
        
        % Update plane RAAN/inc to mean of included satellites
        planeRaans = [planes(planeIdx).satellites.raan];
        planeIncs = [planes(planeIdx).satellites.i];
        planes(planeIdx).raan = mean(planeRaans);
        planes(planeIdx).inc = mean(planeIncs);
    end
end

function order = optimizePlaneOrder(planes)
    % Optimize plane visit order using nearest-neighbor heuristic
    % Minimizes total plane change delta-V (considers both RAAN and inclination)
    
    nPlanes = length(planes);
    if nPlanes <= 1
        order = 1:nPlanes;
        return;
    end
    
    visited = false(1, nPlanes);
    order = zeros(1, nPlanes);
    
    % Start with first plane
    current = 1;
    order(1) = current;
    visited(current) = true;
    
    for i = 2:nPlanes
        minCost = Inf;
        nextPlane = -1;
        
        for j = 1:nPlanes
            if visited(j)
                continue;
            end
            
            % Calculate approximate plane change cost
            cost = estimatePlaneChangeCost(planes(current), planes(j));
            
            if cost < minCost
                minCost = cost;
                nextPlane = j;
            end
        end
        
        order(i) = nextPlane;
        visited(nextPlane) = true;
        current = nextPlane;
    end
end

function cost = estimatePlaneChangeCost(plane1, plane2)
    % Estimate relative cost of plane change (for ordering optimization)
    
    % RAAN difference
    raanDiff = abs(plane2.raan - plane1.raan);
    raanDiff = min(raanDiff, 360 - raanDiff);
    
    % Inclination difference
    incDiff = abs(plane2.inc - plane1.inc);
    
    % Combined cost (weighted)
    % RAAN changes at same inclination are generally cheaper than inc changes
    cost = raanDiff + 3 * incDiff;
end

function [totalDV, visitOrder, details] = calculatePhasingDV(satellites, timePerManeuver, mu)
    % Calculate delta-V for visiting all satellites within a "plane"
    % Accounts for actual RAAN, inclination, altitude, and phase differences
    % between each pair of satellites
    
    nSats = length(satellites);
    
    if nSats <= 1
        totalDV = 0;
        visitOrder = 1;
        details = struct('dv', 0, 'from', '', 'to', '');
        return;
    end
    
    % Sort satellites by mean anomaly for initial visiting order
    % (This is a heuristic - optimal ordering is a TSP problem)
    mas = [satellites.ma];
    [~, visitOrder] = sort(mas);
    
    totalDV = 0;
    details = struct('dv', {}, 'from', {}, 'to', {}, 'phaseAngle', {}, ...
                     'raanChange', {}, 'incChange', {}, 'altChange', {});
    
    for i = 1:nSats-1
        currentIdx = visitOrder(i);
        nextIdx = visitOrder(i+1);
        
        currentSat = satellites(currentIdx);
        nextSat = satellites(nextIdx);
        
        % Calculate full transfer delta-V accounting for all orbital differences
        [dv, dvBreakdown] = calculateTransferDV(currentSat, nextSat, timePerManeuver, mu);
        
        totalDV = totalDV + dv;
        
        details(i).dv = dv;
        details(i).from = currentSat.name;
        details(i).to = nextSat.name;
        details(i).phaseAngle = dvBreakdown.phaseAngle;
        details(i).raanChange = dvBreakdown.raanChange;
        details(i).incChange = dvBreakdown.incChange;
        details(i).altChange = dvBreakdown.altChange;
    end
end

function [totalDV, breakdown] = calculateTransferDV(sat1, sat2, transferTime, mu)
    % Calculate delta-V to transfer from sat1 to sat2
    % Accounts for: phase difference, RAAN difference, inclination difference, altitude difference
    
    % Initialize breakdown
    breakdown = struct('phaseAngle', 0, 'raanChange', 0, 'incChange', 0, 'altChange', 0);
    
    % 1. Calculate plane change component (RAAN and inclination)
    i1 = deg2rad(sat1.i);
    i2 = deg2rad(sat2.i);
    raan1 = deg2rad(sat1.raan);
    raan2 = deg2rad(sat2.raan);
    
    % Angle between orbital planes using spherical trig
    dRAAN = raan2 - raan1;
    cosTheta = cos(i1)*cos(i2) + sin(i1)*sin(i2)*cos(dRAAN);
    theta = acos(max(-1, min(1, cosTheta)));  % Total plane change angle
    
    breakdown.raanChange = rad2deg(abs(dRAAN));
    if breakdown.raanChange > 180
        breakdown.raanChange = 360 - breakdown.raanChange;
    end
    breakdown.incChange = abs(sat2.i - sat1.i);
    
    % Velocity for plane change (at average altitude)
    a_avg = (sat1.a + sat2.a) / 2;
    v_avg = sqrt(mu / a_avg);
    
    % Plane change delta-V (single impulse approximation)
    if theta > 0.001  % More than ~0.06 degrees
        dv_plane = 2 * v_avg * sin(theta / 2);
    else
        dv_plane = 0;
    end
    
    % 2. Calculate altitude change component (Hohmann transfer)
    breakdown.altChange = (sat2.a - sat1.a) / 1000;  % km
    
    if abs(sat1.a - sat2.a) > 1000  % More than 1 km difference
        dv_altitude = calculateHohmannDV(sat1.a, sat2.a, mu);
    else
        dv_altitude = 0;
    end
    
    % 3. Calculate phasing component
    % Phase angle difference (accounting for different RAANs affecting relative position)
    phaseAngle = sat2.ma - sat1.ma;
    if phaseAngle < 0
        phaseAngle = phaseAngle + 360;
    end
    % Take the shorter path
    if phaseAngle > 180
        phaseAngle = 360 - phaseAngle;
    end
    breakdown.phaseAngle = phaseAngle;
    
    phaseAngleRad = deg2rad(phaseAngle);
    
    if phaseAngleRad > 0.01  % More than ~0.6 degrees
        dv_phase = calculatePhasingManeuverDV(a_avg, phaseAngleRad, transferTime, mu);
    else
        dv_phase = 5;  % Minimum rendezvous delta-V
    end
    
    % 4. Combine delta-V components
    % For small plane changes, can combine with phasing maneuver
    % For large plane changes, they add more directly
    if theta < deg2rad(5)  % Less than 5 degree plane change
        % Small plane change - can combine efficiently
        totalDV = sqrt(dv_phase^2 + dv_plane^2) + dv_altitude;
    else
        % Large plane change - maneuvers are more sequential
        totalDV = dv_phase + dv_plane + dv_altitude;
    end
    
    % Minimum delta-V for proximity operations
    totalDV = max(totalDV, 5);
end

function dv = calculatePhasingManeuverDV(a, phaseAngle, transferTime, mu)
    % Calculate delta-V for a phasing maneuver
    % Uses a Hohmann-like transfer orbit to change phase
    
    % Current orbital period
    T = 2 * pi * sqrt(a^3 / mu);
    
    % Current velocity
    v_circ = sqrt(mu / a);
    
    % Number of orbits during transfer time
    nOrbits = transferTime / T;
    
    if nOrbits < 0.5
        nOrbits = 0.5;  % Minimum half orbit
    end
    
    % For phasing, we need to complete different number of orbits than target
    if phaseAngle <= pi
        % Need to speed up (lower orbit) to catch up
        targetOrbits = nOrbits + phaseAngle / (2*pi);
        T_transfer = transferTime / targetOrbits;
    else
        % Go the short way around
        phaseAngle = 2*pi - phaseAngle;
        targetOrbits = nOrbits + phaseAngle / (2*pi);
        T_transfer = transferTime / targetOrbits;
    end
    
    % Semi-major axis of transfer orbit
    a_transfer = (mu * (T_transfer / (2*pi))^2)^(1/3);
    
    % Ensure transfer orbit is valid (perigee above Earth surface)
    r_min = 2 * a_transfer - a;
    if r_min < 6478e3  % ~100 km altitude minimum
        % Need to use multiple orbits
        nOrbits = nOrbits * 2;
        targetOrbits = nOrbits + phaseAngle / (2*pi);
        T_transfer = (transferTime * 2) / targetOrbits;
        a_transfer = (mu * (T_transfer / (2*pi))^2)^(1/3);
    end
    
    % Delta-V for transfer
    if a_transfer < a
        % Lower orbit - current point is apogee of transfer
        r_a = a;
        v_transfer_apo = sqrt(mu * (2/r_a - 1/a_transfer));
        dv1 = abs(v_circ - v_transfer_apo);
    else
        % Higher orbit - current point is perigee of transfer
        r_p = a;
        v_transfer_peri = sqrt(mu * (2/r_p - 1/a_transfer));
        dv1 = abs(v_transfer_peri - v_circ);
    end
    
    dv = 2 * dv1;  % Symmetric for return
    
    % Add minimum delta-V for proximity operations
    dv = dv + 5;  % 5 m/s minimum for rendezvous
end

function dv = calculatePlaneChangeDV(sat1, sat2, mu)
    % Calculate delta-V for plane change between two satellites
    % Accounts for inclination and RAAN changes
    
    % Convert to radians
    i1 = deg2rad(sat1.i);
    i2 = deg2rad(sat2.i);
    raan1 = deg2rad(sat1.raan);
    raan2 = deg2rad(sat2.raan);
    
    % Calculate the angle between orbital planes using spherical trig
    dRAAN = raan2 - raan1;
    
    cosTheta = cos(i1)*cos(i2) + sin(i1)*sin(i2)*cos(dRAAN);
    theta = acos(max(-1, min(1, cosTheta)));  % Clamp for numerical stability
    
    % Velocity at the maneuvering point (use average altitude)
    a_avg = (sat1.a + sat2.a) / 2;
    v = sqrt(mu / a_avg);
    
    % Delta-V for plane change
    dv = 2 * v * sin(theta / 2);
    
    % If there's also an altitude change, combine
    if abs(sat1.a - sat2.a) > 1e3
        dv_altitude = calculateHohmannDV(sat1.a, sat2.a, mu);
        dv = sqrt(dv^2 + dv_altitude^2);
    end
end

function dv = calculateHohmannDV(a1, a2, mu)
    % Calculate Hohmann transfer delta-V between two circular orbits
    
    r1 = a1;
    r2 = a2;
    a_t = (r1 + r2) / 2;
    
    v1 = sqrt(mu / r1);
    v2 = sqrt(mu / r2);
    
    v_t1 = sqrt(mu * (2/r1 - 1/a_t));
    v_t2 = sqrt(mu * (2/r2 - 1/a_t));
    
    dv1 = abs(v_t1 - v1);
    dv2 = abs(v2 - v_t2);
    
    dv = dv1 + dv2;
end
