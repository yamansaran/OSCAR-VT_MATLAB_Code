%% tle_covariance_estimation.m
% Estimates position uncertainty covariance from TLE-to-TLE residuals (SGP4).
% Requires: MATLAB Aerospace Toolbox + Satellite Communications Toolbox
%
% Pipeline:
%   1. Read satellite list from Selected.csv (JCAT column -> NORAD IDs)
%   2. Load historical TLEs from local cache (populated by fetch_tles.py)
%   3. For each satellite, propagate consecutive TLE pairs via SGP4
%   4. Compute ECI residuals, rotate to RTN, build covariance
%   5. Compute Mahalanobis sigma risk for a given relative position
%

clear; clc; close all;

%% ======================== USER INPUTS ========================
csvFile = 'Selected.csv';

% Cache directory (populated by fetch_tles.py)
cacheDir = 'tle_covariance_cache';

% TLE pair gap window [hours]
minGapHrs = 1;       % skip near-duplicate epochs
maxGapHrs = 48;      % don't extrapolate too far

% Full 3x3 covariance or diagonal only
useFullCovariance = true;

%% ======================== 1. PARSE CSV ========================
sats = parseSelectedCSV(csvFile);
nSats = numel(sats);
fprintf('Read %d satellites from %s:\n', nSats, csvFile);
for s = 1:nSats
    fprintf('  NORAD %05d  %s\n', sats(s).noradId, sats(s).name);
end

%% ======================== 2. LOAD CACHED TLEs ========================
if ~isfolder(cacheDir)
    error(['Cache directory "%s" not found.\n' ...
           'Run fetch_tles.py first to download TLEs from Space-Track.'], cacheDir);
end

for s = 1:nSats
    cacheFile = fullfile(cacheDir, sprintf('%05d.tle', sats(s).noradId));
    if ~isfile(cacheFile)
        warning('No cached TLE file for NORAD %05d. Run fetch_tles.py first.', ...
                sats(s).noradId);
    end
    sats(s).tleFile = cacheFile;
end

%% ======================== 3-5. PER-SATELLITE COVARIANCE ========================
results = struct('noradId',{}, 'name',{}, 'nTLE',{}, 'nSamples',{}, ...
                 'mu_rtn',{}, 'sigma_rtn',{}, 'Sigma_target',{}, ...
                 'residuals_rtn',{}, 'gapHours',{});

for s = 1:nSats
    fprintf('\n========== NORAD %05d  %s ==========\n', ...
            sats(s).noradId, sats(s).name);

    tles = parseTLEFile(sats(s).tleFile);
    nTLE = numel(tles);
    fprintf('  Parsed %d TLE sets\n', nTLE);

    if nTLE < 2
        warning('Not enough TLEs for NORAD %d, skipping.', sats(s).noradId);
        continue;
    end

    % Sort by epoch
    [~, idx] = sort([tles.epochDatenum]);
    tles = tles(idx);

    fprintf('  Epoch range: %s to %s\n', ...
        datestr(tles(1).epochDatenum,   'yyyy-mm-dd HH:MM'), ...
        datestr(tles(end).epochDatenum, 'yyyy-mm-dd HH:MM'));

    % --- Residual loop ---
    res_rtn  = [];
    res_eci  = [];
    gapH     = [];
    nFail    = 0;

    for i = 1:nTLE-1
        j = i + 1;

        dtHrs = (tles(j).epochDatenum - tles(i).epochDatenum) * 24;
        if dtHrs < minGapHrs || dtHrs > maxGapHrs
            continue;
        end

        targetEpoch = tles(j).epochDatetime;

        [r_prop, ~]     = sgp4StateAtEpoch(tles(i), targetEpoch);
        [r_ref,  v_ref] = sgp4StateAtEpoch(tles(j), targetEpoch);

        if isempty(r_prop) || isempty(r_ref)
            nFail = nFail + 1;
            if nFail == 1
                warning('SGP4 propagation failed. Check that Satellite Communications Toolbox is installed.');
            end
            if nFail >= 5 && isempty(res_rtn)
                warning('First %d SGP4 calls all failed — skipping satellite.', nFail);
                break;
            end
            continue;
        end

        dr_eci = r_prop - r_ref;
        dr_rtn = eci2rtn(dr_eci, r_ref, v_ref);

        res_eci = [res_eci; dr_eci(:)'];   %#ok<AGROW>
        res_rtn = [res_rtn; dr_rtn(:)'];   %#ok<AGROW>
        gapH    = [gapH; dtHrs];            %#ok<AGROW>
    end

    if nFail > 0
        fprintf('  (%d pairs failed SGP4 propagation)\n', nFail);
    end

    nSamp = size(res_rtn, 1);
    fprintf('  %d residual samples\n', nSamp);

    if nSamp < 3
        warning('Too few samples for NORAD %d, skipping.', sats(s).noradId);
        continue;
    end

    mu_val  = mean(res_rtn, 1);
    sig_val = std(res_rtn, 0, 1);

    if useFullCovariance
        Sig = cov(res_rtn);
    else
        Sig = diag(sig_val.^2);
    end

    fprintf('  Mean   R: %+8.4f   T: %+8.4f   N: %+8.4f  km\n', mu_val);
    fprintf('  Sigma  R: %8.4f    T: %8.4f    N: %8.4f   km\n', sig_val);

    % Store
    r = numel(results) + 1;
    results(r).noradId       = sats(s).noradId;
    results(r).name          = sats(s).name;
    results(r).nTLE          = nTLE;
    results(r).nSamples      = nSamp;
    results(r).mu_rtn        = mu_val;
    results(r).sigma_rtn     = sig_val;
    results(r).Sigma_target  = Sig;
    results(r).residuals_rtn = res_rtn;
    results(r).gapHours      = gapH;
end

nResults = numel(results);
if nResults == 0
    error('No valid covariance results. Check TLE files and toolbox installation.');
end

%% ======================== 6. SUMMARY TABLE ========================
fprintf('\n\n============ COVARIANCE SUMMARY (km) ============\n');
fprintf('%-20s  %8s  %8s  %8s  %6s\n', 'Satellite', 'sig_R', 'sig_T', 'sig_N', 'N');
fprintf('%s\n', repmat('-', 1, 60));
for r = 1:nResults
    fprintf('%-20s  %8.4f  %8.4f  %8.4f  %6d\n', ...
        results(r).name, results(r).sigma_rtn, results(r).nSamples);
end

%% ======================== 7. PLOTS ========================
labels = {'Radial','Along-Track','Cross-Track'};
colors = lines(nResults);

for r = 1:nResults
    figure('Name', results(r).name, 'Position', [80+40*r 80+40*r 1100 400]);
    for k = 1:3
        subplot(1,3,k);
        histogram(results(r).residuals_rtn(:,k), 'Normalization','pdf', ...
                  'FaceColor', colors(r,:), 'FaceAlpha', 0.7);
        hold on;
        xv = linspace(min(results(r).residuals_rtn(:,k)), ...
                       max(results(r).residuals_rtn(:,k)), 200);
        plot(xv, normpdf(xv, results(r).mu_rtn(k), results(r).sigma_rtn(k)), ...
             'r-', 'LineWidth', 1.5);
        xlabel([labels{k} ' (km)']);
        ylabel('PDF');
        title(sprintf('%s: \\sigma=%.3f km', labels{k}, results(r).sigma_rtn(k)));
        grid on;
    end
    sgtitle(sprintf('%s (NORAD %d) — RTN Residuals', ...
            results(r).name, results(r).noradId));
end

% Combined error-growth plot
figure('Name','Error Growth','Position',[100 550 900 400]);
hold on;
for r = 1:nResults
    plot(results(r).gapHours, vecnorm(results(r).residuals_rtn, 2, 2), ...
         'o', 'MarkerSize', 4, 'DisplayName', results(r).name, ...
         'Color', colors(r,:));
end
xlabel('Propagation gap (hours)');
ylabel('||residual||_{RTN} (km)');
title('Position Error Growth — All Satellites');
legend('Location','best');
grid on;

%% ======================== 8. SAVE RESULTS ========================
saveFile = 'tle_covariance_results.mat';
save(saveFile, 'results');
fprintf('\nResults saved to %s\n', saveFile);
fprintf('Run sigma_risk_analysis.m to compute Mahalanobis distances.\n');

%% =====================================================================
%%  LOCAL FUNCTIONS
%% =====================================================================

function sats = parseSelectedCSV(filename)
%PARSESELECTEDCSV  Read Selected.csv and extract NORAD IDs + names.

    fid = fopen(filename, 'r');
    headerLine = fgetl(fid);
    fclose(fid);

    headers = strsplit(headerLine, ',');
    headers = strtrim(headers);

    jcatCol = find(contains(headers, 'JCAT', 'IgnoreCase', true), 1);
    nameCol = find(strcmpi(headers, 'Name'), 1);

    if isempty(jcatCol)
        error('Could not find JCAT column in %s', filename);
    end
    if isempty(nameCol)
        nameCol = jcatCol;
    end

    fid = fopen(filename, 'r');
    fgetl(fid);   % skip header
    sats = struct('noradId',{}, 'name',{});
    idx = 0;
    while ~feof(fid)
        line = fgetl(fid);
        if isempty(strtrim(line)) || line(1) == -1
            continue;
        end
        fields = strsplit(line, ',');
        if numel(fields) < max(jcatCol, nameCol)
            continue;
        end

        jcatStr  = strtrim(fields{jcatCol});
        noradStr = regexprep(jcatStr, '^[A-Za-z]+', '');
        noradId  = str2double(noradStr);

        if isnan(noradId)
            continue;
        end

        idx = idx + 1;
        sats(idx).noradId = noradId;
        sats(idx).name    = strtrim(fields{nameCol});
    end
    fclose(fid);
end


function tles = parseTLEFile(filename)
%PARSETLEFILE  Read a TLE file of arbitrary length and return struct array.
%   Handles 3-line (name+L1+L2) and 2-line (L1+L2) formats.
%   No limit on number of TLEs — reads entire file.

    lines = readlines(filename);
    lines = strtrim(lines);
    lines(lines == "") = [];

    nLines = numel(lines);

    % Pre-allocate generously, trim at end
    maxTLEs = ceil(nLines / 2);
    tles = struct('name', cell(1, maxTLEs), ...
                  'line1', cell(1, maxTLEs), ...
                  'line2', cell(1, maxTLEs), ...
                  'epochDatenum', cell(1, maxTLEs), ...
                  'epochDatetime', cell(1, maxTLEs));

    k = 1;
    idx = 0;
    while k <= nLines - 1
        if k <= nLines - 2 && ...
           startsWith(lines(k+1), '1 ') && startsWith(lines(k+2), '2 ')
            nameLine = lines(k);
            l1 = char(lines(k+1));
            l2 = char(lines(k+2));
            k = k + 3;
        elseif startsWith(lines(k), '1 ') && startsWith(lines(k+1), '2 ')
            nameLine = "UNKNOWN";
            l1 = char(lines(k));
            l2 = char(lines(k+1));
            k = k + 2;
        else
            k = k + 1;
            continue;
        end

        epochYr  = str2double(l1(19:20));
        epochDay = str2double(l1(21:32));

        if epochYr < 57
            fullYr = 2000 + epochYr;
        else
            fullYr = 1900 + epochYr;
        end
        dn = datenum(fullYr, 1, 0) + epochDay;

        idx = idx + 1;
        tles(idx).name          = char(nameLine);
        tles(idx).line1         = l1;
        tles(idx).line2         = l2;
        tles(idx).epochDatenum  = dn;
        tles(idx).epochDatetime = datetime(dn, 'ConvertFrom','datenum');
    end

    % Trim unused pre-allocation
    tles = tles(1:idx);
end


function [r_eci_km, v_eci_kms] = sgp4StateAtEpoch(tleSt, targetDatetime)
%SGP4STATEATEPOCH  Propagate a single TLE to targetDatetime via SGP4.
%   Uses Satellite Communications Toolbox satelliteScenario.

    r_eci_km  = [];
    v_eci_kms = [];

    tmpFile = [tempname, '.tle'];
    fid = fopen(tmpFile, 'w');
    fprintf(fid, '%s\n', tleSt.name);
    fprintf(fid, '%s\n', tleSt.line1);
    fprintf(fid, '%s',   tleSt.line2);
    fclose(fid);

    cleanup = onCleanup(@() safeDelete(tmpFile));

    try
        t0 = targetDatetime;
        t1 = t0 + seconds(1);
        sc  = satelliteScenario(t0, t1, 1);
        sat = satellite(sc, tmpFile);

        [r, v] = states(sat, 'CoordinateFrame', 'inertial');

        r_eci_km  = r(:,1) / 1e3;
        v_eci_kms = v(:,1) / 1e3;

        delete(sc);
    catch ME
        warning('sgp4StateAtEpoch: %s', ME.message);
    end
end


function dr_rtn = eci2rtn(dr_eci, r_eci, v_eci)
%ECI2RTN  Rotate an ECI vector into RTN (Hill) frame.

    r_eci  = r_eci(:);
    v_eci  = v_eci(:);
    dr_eci = dr_eci(:);

    R_hat = r_eci / norm(r_eci);
    N_hat = cross(r_eci, v_eci);
    N_hat = N_hat / norm(N_hat);
    T_hat = cross(N_hat, R_hat);

    C = [R_hat'; T_hat'; N_hat'];
    dr_rtn = C * dr_eci;
end


function safeDelete(f)
    if isfile(f)
        delete(f);
    end
end