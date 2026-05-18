%% sigma_risk_analysis.m
% Loads TLE covariance results and computes Mahalanobis sigma risk
% for a variety of relative position vectors and inspector covariances.
%
% Run tle_covariance_estimation.m first to generate tle_covariance_results.mat
%

clear; clc; close all;

%% ======================== LOAD RESULTS ========================
matFile = 'tle_covariance_results.mat';
if ~isfile(matFile)
    error('Results file "%s" not found. Run tle_covariance_estimation.m first.', matFile);
end
load(matFile, 'results');
nTargets = numel(results);

fprintf('Loaded covariance results for %d satellites:\n', nTargets);
for t = 1:nTargets
    fprintf('  [%d] %-20s  sig_R=%.4f  sig_T=%.4f  sig_N=%.4f km  (%d samples)\n', ...
        t, results(t).name, results(t).sigma_rtn, results(t).nSamples);
end



%% ======================== USER INPUTS ========================

% --- Inspector (own) covariance options to sweep [km^2] ---
% Each row: [sig_R^2, sig_T^2, sig_N^2] (diagonal elements)
% Add as many rows as you want to test
Sigma_own_options = [
    0.01,  0.04,  0.005;    % tight nav (100m R, 200m T, 70m N)
    %0.05,  0.10,  0.02;     % moderate nav
    %0.25,  1.00,  0.10;     % loose nav
    %0.00,  0.00,  0.00;     % no inspector uncertainty (target-only)
];
nSigmaOwn = size(Sigma_own_options, 1);

% Labels for each Sigma_own row (for table/plot readability)
Sigma_own_labels = {
    'Tight nav'
    'Moderate nav'
    'Loose nav'
    'Target only'
};

% --- Relative position vectors to sweep [km in RTN] ---
% Each row: [R, T, N]
r_vectors = [
    0.0,  0.3,  0.0;    % 500m along-track
    %0.0,  1.0,  0.0;    % 1 km along-track
    %0.0,  2.0,  0.0;    % 2 km along-track
    %0.0,  5.0,  0.0;    % 5 km along-track
    %0.5,  0.0,  0.0;    % 500m radial
    %1.0,  0.0,  0.0;    % 1 km radial
    %0.0,  0.0,  0.5;    % 500m cross-track
    %0.0,  0.0,  1.0;    % 1 km cross-track
    %0.5,  0.5,  0.5;    % 500m all axes
    %1.0,  1.0,  1.0;    % 1 km all axes
];
nR = size(r_vectors, 1);

% Keep-out zone center [km in RTN]
r_safe = [0; 0; 0];

% --- SENSOR DEFINITIONS ---
% Toggle sensors on/off, set noise parameters.
% These are used by both the inverse solver and the sensor fusion section.
% NOTE: Don't enable rangefinder + LIDAR simultaneously (double-counts range).

sensors = struct();

% Sensor 1: Rangefinder (measures range only)
sensors(1).name    = 'Rangefinder';
sensors(1).enabled = true;
sensors(1).type    = 'range';          % 1 measurement: range
sensors(1).sigma   = 0.003;            % 1-sigma range noise [km] (= 10 m)

% Sensor 2: LIDAR (measures range + azimuth + elevation)
sensors(2).name    = 'LIDAR';
sensors(2).enabled = false;
sensors(2).type    = 'lidar';          % 3 measurements: range, azimuth, elevation
sensors(2).sigma   = [0.005, ...       % range noise [km] (= 1 m)
                      deg2rad(0.01), ... % azimuth noise [rad] (= 0.1 deg)
                      deg2rad(0.01)];    % elevation noise [rad] (= 0.1 deg)

% Sensor 3: Star tracker relative bearing (angles only, no range)
sensors(3).name    = 'Star Tracker (bearing)';
sensors(3).enabled = false;
sensors(3).type    = 'angles_only';    % 2 measurements: azimuth, elevation
sensors(3).sigma   = [deg2rad(0.05), ... % azimuth noise [rad]
                      deg2rad(0.05)];    % elevation noise [rad]

%% ======================== COMPUTE SIGMA RISK ========================
% sigma_risk(i, j, t) = Mahalanobis distance for
%   r_vector(i,:), Sigma_own_option(j,:), target satellite t

% Pre-compute dr vectors: nR x 3
dr_all = r_vectors - r_safe';   % nR x 3

% Vectorized computation per target and Sigma_own
sigma_risk = zeros(nR, nSigmaOwn, nTargets);


for t = 1:nTargets
    Sig_target = results(t).Sigma_target;   % 3x3

    for j = 1:nSigmaOwn
        Sig_own = diag(Sigma_own_options(j,:));
        Sig_combined = Sig_own + Sig_target;        % 3x3
        % Assume 2-band matching
        Sig_combined = Sig_own + Sig_own;         % 3x3
        Sig_inv = inv(Sig_combined);                  %#ok<MINV>

        % Vectorized over all r: dr_all is nR x 3
        % quad(i) = dr_all(i,:) * Sig_inv * dr_all(i,:)'
        quad = sum((dr_all * Sig_inv) .* dr_all, 2);  % nR x 1
        sigma_risk(:, j, t) = sqrt(quad);
    end
end

%% ======================== DISPLAY TABLES ========================
fprintf('\n');
for t = 1:nTargets
    fprintf('============================================================\n');
    fprintf('  %s (NORAD %d)\n', results(t).name, results(t).noradId);
    fprintf('  Target covariance diag: [%.4f, %.4f, %.4f] km^2\n', ...
        results(t).Sigma_target(1,1), results(t).Sigma_target(2,2), results(t).Sigma_target(3,3));
    fprintf('============================================================\n');

    % Header
    fprintf('%-24s', 'r_RTN [km]');
    for j = 1:nSigmaOwn
        fprintf('  %14s', Sigma_own_labels{j});
    end
    fprintf('\n');
    fprintf('%s\n', repmat('-', 1, 24 + nSigmaOwn*16));

    % Rows
    for i = 1:nR
        rStr = sprintf('[%5.2f %5.2f %5.2f]', r_vectors(i,:));
        fprintf('%-24s', rStr);
        for j = 1:nSigmaOwn
            val = sigma_risk(i, j, t);
            % Color code: mark dangerous values
            if val < 1
                marker = ' !!';
            elseif val < 3
                marker = '  *';
            else
                marker = '   ';
            end
            fprintf('  %10.4f%s', val, marker);
        end
        fprintf('\n');
    end
    fprintf('\n');
end

fprintf('Legend: !! = inside 1-sigma (high risk),  * = inside 3-sigma (caution)\n\n');

%% ======================== HEATMAP PLOTS ========================
% One heatmap per target satellite, axes = r_vector vs Sigma_own

rLabels = cell(nR, 1);
for i = 1:nR
    rLabels{i} = sprintf('[%.1f %.1f %.1f]', r_vectors(i,:));
end

figPos = [418 182 700 500];   % centered on 1536x864

for t = 1:nTargets
    figure('Name', sprintf('Sigma Risk — %s', results(t).name), ...
           'Position', figPos);

    data = sigma_risk(:, :, t);   % nR x nSigmaOwn

    imagesc(data);
    colormap(flipud(hot));   % red = low sigma (dangerous), yellow/white = safe
    colorbar;
    caxis([0, max(data(:))*1.1]);

    set(gca, 'XTick', 1:nSigmaOwn, 'XTickLabel', Sigma_own_labels, ...
             'YTick', 1:nR, 'YTickLabel', rLabels);
    xlabel('Inspector Covariance');
    ylabel('Relative Position r_{RTN} [km]');
    title(sprintf('Mahalanobis \\sigma-risk — %s (NORAD %d)', ...
          results(t).name, results(t).noradId));

    % Annotate cells with values
    for i = 1:nR
        for j = 1:nSigmaOwn
            val = data(i,j);
            if val < 3
                clr = 'w';
            else
                clr = 'k';
            end
            text(j, i, sprintf('%.2f', val), ...
                 'HorizontalAlignment','center', 'Color', clr, 'FontSize', 8);
        end
    end
end

%% ======================== SWEEP: ALONG-TRACK DISTANCE ========================
% Continuous sigma risk vs distance along one axis (e.g. along-track)

dist_km = linspace(0.1, 10, 200)';   % 100m to 10km

figPos = [318 182 900 500];   % centered on 1536x864

figure('Name','Sigma Risk vs Along-Track Distance', ...
       'Position', figPos);
hold on;
plotStyles = {'-', '--', ':', '-.'};
colors = lines(nTargets);

legendEntries = {};
for t = 1:nTargets
    Sig_target = results(t).Sigma_target;

    for j = 1:nSigmaOwn
        Sig_own = diag(Sigma_own_options(j,:));
        Sig_combined = Sig_own + Sig_target;
        Sig_inv = inv(Sig_combined);                  %#ok<MINV>

        % dr = [0, d, 0] for along-track sweep
        dr_sweep = [zeros(numel(dist_km),1), dist_km, zeros(numel(dist_km),1)];
        quad = sum((dr_sweep * Sig_inv) .* dr_sweep, 2);
        sr = sqrt(quad);

        styleIdx = mod(j-1, numel(plotStyles)) + 1;
        plot(dist_km, sr, plotStyles{styleIdx}, ...
             'Color', colors(t,:), 'LineWidth', 1.3);

        legendEntries{end+1} = sprintf('%s / %s', results(t).name, Sigma_own_labels{j}); %#ok<AGROW>
    end
end

yline(1, 'r--', '1\sigma', 'LineWidth', 1, 'LabelHorizontalAlignment','left');
yline(3, 'b--', '3\sigma', 'LineWidth', 1, 'LabelHorizontalAlignment','left');

xlabel('Along-Track Distance (km)');
ylabel('\sigma_{risk} (Mahalanobis Distance)');
title('Sigma Risk vs Along-Track Separation');
legend(legendEntries, 'Location', 'eastoutside', 'FontSize', 7);
grid on;
hold off;

%% ======================== SWEEP: RADIAL DISTANCE ========================
figure('Name','Sigma Risk vs Radial Distance', ...
       'Position', figPos);
hold on;
legendEntries = {};
for t = 1:nTargets
    Sig_target = results(t).Sigma_target;
    for j = 1:nSigmaOwn
        Sig_own = diag(Sigma_own_options(j,:));
        Sig_combined = Sig_own + Sig_target;
        Sig_inv = inv(Sig_combined);                  %#ok<MINV>

        dr_sweep = [dist_km, zeros(numel(dist_km),1), zeros(numel(dist_km),1)];
        quad = sum((dr_sweep * Sig_inv) .* dr_sweep, 2);
        sr = sqrt(quad);

        styleIdx = mod(j-1, numel(plotStyles)) + 1;
        plot(dist_km, sr, plotStyles{styleIdx}, ...
             'Color', colors(t,:), 'LineWidth', 1.3);
        legendEntries{end+1} = sprintf('%s / %s', results(t).name, Sigma_own_labels{j}); %#ok<AGROW>
    end
end

yline(1, 'r--', '1\sigma', 'LineWidth', 1, 'LabelHorizontalAlignment','left');
yline(3, 'b--', '3\sigma', 'LineWidth', 1, 'LabelHorizontalAlignment','left');

xlabel('Radial Distance (km)');
ylabel('\sigma_{risk} (Mahalanobis Distance)');
title('Sigma Risk vs Radial Separation');
legend(legendEntries, 'Location', 'eastoutside', 'FontSize', 7);
grid on;
hold off;

fprintf('Done. Figures show sigma risk heatmaps and distance sweeps.\n');

%% ======================== INVERSE: SIGMA -> DISTANCE ========================
% Given desired sigma risk thresholds, compute the required separation
% distance along each pure axis (R-only, T-only, N-only).
%
% WITHOUT sensors: d = sigma_desired / sqrt(Sig_inv(k,k))  (closed-form)
% WITH sensors:    solve numerically via fzero, because H depends on d

% --- Desired sigma thresholds to solve for ---
sigma_thresholds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
nThresh = numel(sigma_thresholds);

axisLabels = {'R (Radial)', 'T (Along-Track)', 'N (Cross-Track)'};

% Unit vectors for each pure-axis direction
axisUnitVecs = {[1;0;0], [0;1;0], [0;0;1]};

fprintf('\n\n============================================================\n');
fprintf('  INVERSE: Required separation distance for desired sigma\n');
fprintf('  (comparing: without sensors vs with sensors)\n');
fprintf('============================================================\n');

for t = 1:nTargets
    Sig_target = results(t).Sigma_target;

    fprintf('\n--- %s (NORAD %d) ---\n', results(t).name, results(t).noradId);

    for j = 1:nSigmaOwn
        Sig_own = diag(Sigma_own_options(j,:));
        Sig_combined = Sig_own + Sig_target;
        Sig_inv = inv(Sig_combined);                  %#ok<MINV>

        fprintf('\n  Inspector: %s\n', Sigma_own_labels{j});

        % Header
        fprintf('  %-20s', 'Axis');
        for s = 1:nThresh
            fprintf('  %18s', sprintf('%d-sigma', sigma_thresholds(s)));
        end
        fprintf('\n');
        fprintf('  %-20s', '');
        for s = 1:nThresh
            fprintf('   %7s / %7s', 'no sns', 'w/ sns');
        end
        fprintf('\n');
        fprintf('  %s\n', repmat('-', 1, 20 + nThresh*20));

        for k = 1:3
            fprintf('  %-20s', axisLabels{k});
            uv = axisUnitVecs{k};

            for s = 1:nThresh
                % --- Without sensors (closed-form) ---
                d_no = sigma_thresholds(s) / sqrt(Sig_inv(k,k));

                % --- With sensors (numerical solve) ---
                % sigma_risk(d) = sqrt( (d*uv)' * Sig_new(d)^-1 * (d*uv) )
                % Solve: sigma_risk(d) - sigma_target = 0
                sigmaFun = @(d) sigmaRiskAtDist(d, uv, Sig_combined, sensors) ...
                               - sigma_thresholds(s);

                % Bracket: start from closed-form solution, search outward
                try
                    d_with = fzero(sigmaFun, [1e-6, d_no * 2]);
                catch
                    % If bracket fails, try wider range
                    try
                        d_with = fzero(sigmaFun, [1e-6, d_no * 10]);
                    catch
                        d_with = NaN;  % no solution found
                    end
                end

                % Format output
                fprintf('   %s / %s', fmtDist(d_no), fmtDist(d_with));
            end
            fprintf('\n');
        end
    end
end

fprintf('\nInterpretation: distance needed along each axis (with others = 0)\n');
fprintf('to achieve the given sigma risk level.\n');
fprintf('"no sns" = nav+TLE only,  "w/ sns" = with onboard sensors\n');

%% ======================== INVERSE: PLOTS ========================

figPos_wide = [218 182 1100 500];   % centered on 1536x864
figPos_std  = [268 182 1000 500];   % centered on 1536x864

% --- Bar chart: required distance per axis, with/without sensors ---
for t = 1:nTargets
    Sig_target = results(t).Sigma_target;

    figure('Name', sprintf('Required Distance — %s', results(t).name), ...
           'Position', figPos_wide);

    for j = 1:nSigmaOwn
        Sig_own = diag(Sigma_own_options(j,:));
        Sig_combined = Sig_own + Sig_target;
        Sig_inv = inv(Sig_combined);                  %#ok<MINV>

        % Without sensors (closed-form)
        d_no = zeros(3, nThresh);
        for k = 1:3
            d_no(k, :) = sigma_thresholds / sqrt(Sig_inv(k,k));
        end

        % With sensors (numerical)
        d_with = zeros(3, nThresh);
        for k = 1:3
            uv = axisUnitVecs{k};
            for s = 1:nThresh
                sigmaFun = @(d) sigmaRiskAtDist(d, uv, Sig_combined, sensors) ...
                               - sigma_thresholds(s);
                try
                    d_with(k,s) = fzero(sigmaFun, [1e-6, d_no(k,s)*2]);
                catch
                    try
                        d_with(k,s) = fzero(sigmaFun, [1e-6, d_no(k,s)*10]);
                    catch
                        d_with(k,s) = NaN;
                    end
                end
            end
        end

        subplot(1, nSigmaOwn, j);
        % Group: for each sigma threshold, show 6 bars (3 axes × 2 modes)
        % Reshape into grouped bar format
        barData = zeros(nThresh, 6);
        for s = 1:nThresh
            barData(s,:) = [d_no(1,s), d_with(1,s), ...   % R no/with
                            d_no(2,s), d_with(2,s), ...   % T no/with
                            d_no(3,s), d_with(3,s)] * 1000; % N no/with, in meters
        end

        b = bar(categorical(compose('%d\\sigma', sigma_thresholds)), barData);
        % Solid = no sensors, hatched/lighter = with sensors
        b(1).FaceColor = [0.85 0.33 0.10];   b(2).FaceColor = [0.85 0.33 0.10]*0.5 + 0.5;
        b(3).FaceColor = [0.00 0.45 0.74];   b(4).FaceColor = [0.00 0.45 0.74]*0.5 + 0.5;
        b(5).FaceColor = [0.47 0.67 0.19];   b(6).FaceColor = [0.47 0.67 0.19]*0.5 + 0.5;

        ylabel('Distance (m)');
        title(Sigma_own_labels{j}, 'FontSize', 10);
        legend({'R','R+sns','T','T+sns','N','N+sns'}, ...
               'Location','northwest', 'FontSize', 6);
        grid on;
    end
    sgtitle(sprintf('Required Separation — %s (NORAD %d)\n(dark=no sensors, light=with sensors)', ...
            results(t).name, results(t).noradId));
end

% --- Continuous: distance vs sigma (0.5 to 6 sigma) per axis ---
sigma_sweep = linspace(0.5, 6, 200);
axisColors = [0.85 0.33 0.10;   % R
              0.00 0.45 0.74;   % T
              0.47 0.67 0.19];  % N

for t = 1:nTargets
    Sig_target = results(t).Sigma_target;

    figure('Name', sprintf('Distance vs Sigma — %s', results(t).name), ...
           'Position', figPos_std);

    legendEntries = {};
    hold on;

    for j = 1:nSigmaOwn
        Sig_own = diag(Sigma_own_options(j,:));
        Sig_combined = Sig_own + Sig_target;
        Sig_inv = inv(Sig_combined);                  %#ok<MINV>

        styleIdx = mod(j-1, numel(plotStyles)) + 1;

        for k = 1:3
            uv = axisUnitVecs{k};

            % Without sensors (closed-form)
            d_sweep_no = sigma_sweep / sqrt(Sig_inv(k,k));

            % With sensors (numerical for each sigma)
            d_sweep_with = zeros(size(sigma_sweep));
            for s = 1:numel(sigma_sweep)
                sigmaFun = @(d) sigmaRiskAtDist(d, uv, Sig_combined, sensors) ...
                               - sigma_sweep(s);
                try
                    d_sweep_with(s) = fzero(sigmaFun, [1e-6, d_sweep_no(s)*2]);
                catch
                    try
                        d_sweep_with(s) = fzero(sigmaFun, [1e-6, d_sweep_no(s)*10]);
                    catch
                        d_sweep_with(s) = NaN;
                    end
                end
            end

            plot(sigma_sweep, d_sweep_no, plotStyles{styleIdx}, ...
                 'Color', axisColors(k,:), 'LineWidth', 1.5);
            plot(sigma_sweep, d_sweep_with, plotStyles{styleIdx}, ...
                 'Color', axisColors(k,:)*0.5 + 0.5, 'LineWidth', 1.5);

            legendEntries{end+1} = sprintf('%s / %s', axisLabels{k}, Sigma_own_labels{j}); %#ok<AGROW>
            legendEntries{end+1} = sprintf('%s / %s +sns', axisLabels{k}, Sigma_own_labels{j}); %#ok<AGROW>
        end
    end

    for s = [1, 3]
        xline(s, 'k:', sprintf('%d\\sigma', s), 'LineWidth', 0.8, ...
              'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
    end

    xlabel('\sigma_{risk} Threshold');
    ylabel('Required Separation Distance (km)');
    title(sprintf('Min Separation vs \\sigma — %s (NORAD %d)\n(dark=no sensors, light=with sensors)', ...
          results(t).name, results(t).noradId));
    legend(legendEntries, 'Location', 'eastoutside', 'FontSize', 6);
    grid on;
    hold off;
end

% --- Polar / radar-style: 3-axis footprint at each sigma level ---
figPos_env = [443 157 650 550];   % centered on 1536x864

for t = 1:nTargets
    Sig_target = results(t).Sigma_target;

    figure('Name', sprintf('Safety Envelope — %s', results(t).name), ...
           'Position', figPos_env);

    % Use first Sigma_own option for the envelope plot
    j_plot = 1;
    Sig_own = diag(Sigma_own_options(j_plot,:));
    Sig_combined = Sig_own + Sig_target;
    Sig_inv = inv(Sig_combined);                      %#ok<MINV>

    % For each sigma threshold, draw an ellipse in R-T plane and R-N plane
    theta = linspace(0, 2*pi, 360);
    envelopeColors = parula(nThresh + 1);

    % --- R-T plane (subplot 1) ---
    subplot(1,2,1);
    hold on;
    for s = 1:nThresh
        sig = sigma_thresholds(s);
        % Ellipse in R-T: extract 2x2 submatrix [R,T]
        Sig2 = Sig_combined([1,2], [1,2]);
        [V, D] = eig(Sig2);
        radii = sig * sqrt(diag(D));
        ellipse_pts = V * diag(radii) * [cos(theta); sin(theta)];

        plot(ellipse_pts(2,:), ellipse_pts(1,:), '-', ...
             'Color', envelopeColors(s,:), 'LineWidth', 1.5);
    end
    plot(0, 0, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
    xlabel('Along-Track T (km)');
    ylabel('Radial R (km)');
    title('R-T Plane');
    axis equal; grid on;
    legend([compose('%d\\sigma', sigma_thresholds), {'Target'}], ...
           'Location', 'best', 'FontSize', 7);
    hold off;

    % --- R-N plane (subplot 2) ---
    subplot(1,2,2);
    hold on;
    for s = 1:nThresh
        sig = sigma_thresholds(s);
        Sig2 = Sig_combined([1,3], [1,3]);
        [V, D] = eig(Sig2);
        radii = sig * sqrt(diag(D));
        ellipse_pts = V * diag(radii) * [cos(theta); sin(theta)];

        plot(ellipse_pts(2,:), ellipse_pts(1,:), '-', ...
             'Color', envelopeColors(s,:), 'LineWidth', 1.5);
    end
    plot(0, 0, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
    xlabel('Cross-Track N (km)');
    ylabel('Radial R (km)');
    title('R-N Plane');
    axis equal; grid on;
    legend([compose('%d\\sigma', sigma_thresholds), {'Target'}], ...
           'Location', 'best', 'FontSize', 7);
    hold off;

    sgtitle(sprintf('Uncertainty Envelope — %s / %s\n(NORAD %d)', ...
            results(t).name, Sigma_own_labels{j_plot}, results(t).noradId));
end

%% ======================== SAFETY ENVELOPE: SENSOR-FUSED ONLY ========================
figPos_env = [443 157 650 550];

for t = 1:nTargets
    Sig_target = results(t).Sigma_target;

    figure('Name', sprintf('Safety Envelope (Sensors) — %s', results(t).name), ...
           'Position', figPos_env);

    j_plot = 1;
    Sig_own = diag(Sigma_own_options(j_plot,:));
    Sig_combined = Sig_own + Sig_target;

    % Representative position for sensor Jacobian
    dr_repr = [0; 1; 0];
    for i = 1:nR
        if norm(r_vectors(i,:)) > 0
            dr_repr = dr_all(i,:)';
            break;
        end
    end

    Sig_fused = applySensorFusion(Sig_combined, dr_repr, sensors);

    theta = linspace(0, 2*pi, 360);
    envelopeColors = parula(nThresh + 1);

    % --- R-T plane ---
    subplot(1,2,1); hold on;
    for s = 1:nThresh
        sig = sigma_thresholds(s);
        Sig2 = Sig_fused([1,2],[1,2]);
        [V, D] = eig(Sig2);
        radii = sig * sqrt(diag(D));
        pts = V * diag(radii) * [cos(theta); sin(theta)];
        plot(pts(2,:), pts(1,:), '-', ...
             'Color', envelopeColors(s,:), 'LineWidth', 1.5);
    end
    plot(0, 0, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
    xlabel('Along-Track T (km)');
    ylabel('Radial R (km)');
    title('R-T Plane');
    axis equal; grid on;
    legend([compose('%d\\sigma', sigma_thresholds), {'Target'}], ...
           'Location', 'best', 'FontSize', 7);
    hold off;

    % --- R-N plane ---
    subplot(1,2,2); hold on;
    for s = 1:nThresh
        sig = sigma_thresholds(s);
        Sig2 = Sig_fused([1,3],[1,3]);
        [V, D] = eig(Sig2);
        radii = sig * sqrt(diag(D));
        pts = V * diag(radii) * [cos(theta); sin(theta)];
        plot(pts(2,:), pts(1,:), '-', ...
             'Color', envelopeColors(s,:), 'LineWidth', 1.5);
    end
    plot(0, 0, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
    xlabel('Cross-Track N (km)');
    ylabel('Radial R (km)');
    title('R-N Plane');
    axis equal; grid on;
    legend([compose('%d\\sigma', sigma_thresholds), {'Target'}], ...
           'Location', 'best', 'FontSize', 7);
    hold off;

    sgtitle(sprintf('Sensor-Fused Uncertainty Envelope — %s / %s\n(NORAD %d)  at r=[%.1f, %.1f, %.1f] km', ...
            results(t).name, Sigma_own_labels{j_plot}, results(t).noradId, dr_repr));
end

%% ========================================================================
%%  SENSOR FUSION: Improved covariance with proximity sensors
%% ========================================================================
% Uses sensor definitions from USER INPUTS section above.

fprintf('\n\n============================================================\n');
fprintf('  SENSOR FUSION: Covariance improvement with onboard sensors\n');
fprintf('============================================================\n');

enabledNames = {};
for si = 1:numel(sensors)
    if sensors(si).enabled
        enabledNames{end+1} = sensors(si).name; %#ok<AGROW>
    end
end
fprintf('Enabled sensors: %s\n', strjoin(enabledNames, ', '));

% --- Compute improved sigma risk for each (r_vector, Sigma_own, target) ---
sigma_risk_improved = zeros(nR, nSigmaOwn, nTargets);

for t = 1:nTargets
    Sig_target = results(t).Sigma_target;

    for j = 1:nSigmaOwn
        Sig_own = diag(Sigma_own_options(j,:));
        Sig_old = Sig_own + Sig_target;

        for i = 1:nR
            dr = dr_all(i,:)';   % [R; T; N]

            % Apply sensor fusion
            Sig_new = applySensorFusion(Sig_old, dr, sensors);

            % Mahalanobis distance with improved covariance
            sigma_risk_improved(i, j, t) = sqrt(dr' / Sig_new * dr);
        end
    end
end

% --- Display comparison table ---
fprintf('\n');
for t = 1:nTargets
    fprintf('============================================================\n');
    fprintf('  %s (NORAD %d) — With sensor fusion\n', ...
        results(t).name, results(t).noradId);
    fprintf('============================================================\n');

    for j = 1:nSigmaOwn
        fprintf('\n  Inspector: %s\n', Sigma_own_labels{j});

        fprintf('  %-24s  %12s  %12s  %10s\n', ...
            'r_RTN [km]', 'No sensors', 'W/ sensors', 'Improvement');
        fprintf('  %s\n', repmat('-', 1, 62));

        for i = 1:nR
            rStr = sprintf('[%5.2f %5.2f %5.2f]', r_vectors(i,:));
            old_val = sigma_risk(i, j, t);
            new_val = sigma_risk_improved(i, j, t);
            improv  = (new_val / old_val - 1) * 100;

            fprintf('  %-24s  %10.4f    %10.4f    %+8.1f%%\n', ...
                rStr, old_val, new_val, improv);
        end
    end
    fprintf('\n');
end

% --- Covariance ellipse comparison: before vs after sensors ---
figPos_fusion = [318 207 900 450];   % centered on 1536x864

for t = 1:nTargets
    Sig_target = results(t).Sigma_target;
    j_plot = 1;   % which Sigma_own to use for ellipse plot
    Sig_own = diag(Sigma_own_options(j_plot,:));
    Sig_old = Sig_own + Sig_target;

    % Pick a representative position for the sensor Jacobian
    % Use 1km along-track as default
    dr_repr = [0; 1; 0];
    % If there's a nonzero r_vector, use the first one
    for i = 1:nR
        if norm(r_vectors(i,:)) > 0
            dr_repr = dr_all(i,:)';
            break;
        end
    end

    Sig_new = applySensorFusion(Sig_old, dr_repr, sensors);

    figure('Name', sprintf('Sensor Fusion Envelope — %s', results(t).name), ...
           'Position', figPos_fusion);
    theta = linspace(0, 2*pi, 360);
    sig_level = 3;   % draw 3-sigma ellipses

    % --- R-T plane ---
    subplot(1,2,1); hold on;

    % Before
    Sig2_old = Sig_old([1,2],[1,2]);
    [V,D] = eig(Sig2_old);
    radii = sig_level * sqrt(diag(D));
    pts = V * diag(radii) * [cos(theta); sin(theta)];
    plot(pts(2,:), pts(1,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Nav + TLE only');

    % After
    Sig2_new = Sig_new([1,2],[1,2]);
    [V,D] = eig(Sig2_new);
    radii = sig_level * sqrt(diag(D));
    pts = V * diag(radii) * [cos(theta); sin(theta)];
    plot(pts(2,:), pts(1,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'With sensors');

    plot(dr_repr(2), dr_repr(1), 'k^', 'MarkerSize', 10, ...
         'MarkerFaceColor','k', 'DisplayName', 'Inspector pos');
    plot(0, 0, 'rx', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Target');

    xlabel('Along-Track T (km)'); ylabel('Radial R (km)');
    title(sprintf('R-T Plane (%d\\sigma)', sig_level));
    legend('Location','best','FontSize',7);
    axis equal; grid on; hold off;

    % --- R-N plane ---
    subplot(1,2,2); hold on;

    Sig2_old = Sig_old([1,3],[1,3]);
    [V,D] = eig(Sig2_old);
    radii = sig_level * sqrt(diag(D));
    pts = V * diag(radii) * [cos(theta); sin(theta)];
    plot(pts(2,:), pts(1,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Nav + TLE only');

    Sig2_new = Sig_new([1,3],[1,3]);
    [V,D] = eig(Sig2_new);
    radii = sig_level * sqrt(diag(D));
    pts = V * diag(radii) * [cos(theta); sin(theta)];
    plot(pts(2,:), pts(1,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'With sensors');

    plot(dr_repr(3), dr_repr(1), 'k^', 'MarkerSize', 10, ...
         'MarkerFaceColor','k', 'DisplayName', 'Inspector pos');
    plot(0, 0, 'rx', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Target');

    xlabel('Cross-Track N (km)'); ylabel('Radial R (km)');
    title(sprintf('R-N Plane (%d\\sigma)', sig_level));
    legend('Location','best','FontSize',7);
    axis equal; grid on; hold off;

    sgtitle(sprintf('Sensor Fusion: Uncertainty Reduction — %s / %s\n(at r=[%.1f, %.1f, %.1f] km)', ...
            results(t).name, Sigma_own_labels{j_plot}, dr_repr));
end

% --- Sweep: sigma risk vs along-track distance, with & without sensors ---
figure('Name', 'Sensor Fusion: Along-Track Sweep', ...
       'Position', figPos_std);
hold on;
legendEntries = {};
dist_km_sweep = linspace(0.1, 10, 200)';

for t = 1:nTargets
    Sig_target = results(t).Sigma_target;
    j_sweep = 1;   % use first Sigma_own for sweep
    Sig_own = diag(Sigma_own_options(j_sweep,:));
    Sig_old = Sig_own + Sig_target;
    Sig_inv_old = inv(Sig_old);                       %#ok<MINV>

    % Without sensors (vectorized)
    dr_sweep = [zeros(numel(dist_km_sweep),1), dist_km_sweep, zeros(numel(dist_km_sweep),1)];
    quad_old = sum((dr_sweep * Sig_inv_old) .* dr_sweep, 2);
    sr_old = sqrt(quad_old);

    % With sensors (loop — H changes at each position)
    sr_new = zeros(size(dist_km_sweep));
    for i = 1:numel(dist_km_sweep)
        dr_i = [0; dist_km_sweep(i); 0];
        Sig_new = applySensorFusion(Sig_old, dr_i, sensors);
        sr_new(i) = sqrt(dr_i' / Sig_new * dr_i);
    end

    plot(dist_km_sweep, sr_old, '-',  'Color', colors(t,:), 'LineWidth', 1.5);
    plot(dist_km_sweep, sr_new, '--', 'Color', colors(t,:), 'LineWidth', 1.5);

    legendEntries{end+1} = sprintf('%s (no sensors)', results(t).name);   %#ok<AGROW>
    legendEntries{end+1} = sprintf('%s (with sensors)', results(t).name); %#ok<AGROW>
end

yline(1, 'r:', '1\sigma', 'LineWidth', 0.8, 'LabelHorizontalAlignment','left');
yline(3, 'b:', '3\sigma', 'LineWidth', 0.8, 'LabelHorizontalAlignment','left');

xlabel('Along-Track Distance (km)');
ylabel('\sigma_{risk} (Mahalanobis Distance)');
title(sprintf('Sensor Fusion Effect on Sigma Risk (%s)', Sigma_own_labels{j_sweep}));
legend(legendEntries, 'Location', 'eastoutside', 'FontSize', 7);
grid on; hold off;

% --- Sweep: radial distance ---
figure('Name', 'Sensor Fusion: Radial Sweep', ...
       'Position', figPos_std);
hold on;
legendEntries = {};

for t = 1:nTargets
    Sig_target = results(t).Sigma_target;
    j_sweep = 1;
    Sig_own = diag(Sigma_own_options(j_sweep,:));
    Sig_old = Sig_own + Sig_target;
    Sig_inv_old = inv(Sig_old);                       %#ok<MINV>

    dr_sweep = [dist_km_sweep, zeros(numel(dist_km_sweep),1), zeros(numel(dist_km_sweep),1)];
    quad_old = sum((dr_sweep * Sig_inv_old) .* dr_sweep, 2);
    sr_old = sqrt(quad_old);

    sr_new = zeros(size(dist_km_sweep));
    for i = 1:numel(dist_km_sweep)
        dr_i = [dist_km_sweep(i); 0; 0];
        Sig_new = applySensorFusion(Sig_old, dr_i, sensors);
        sr_new(i) = sqrt(dr_i' / Sig_new * dr_i);
    end

    plot(dist_km_sweep, sr_old, '-',  'Color', colors(t,:), 'LineWidth', 1.5);
    plot(dist_km_sweep, sr_new, '--', 'Color', colors(t,:), 'LineWidth', 1.5);

    legendEntries{end+1} = sprintf('%s (no sensors)', results(t).name);   %#ok<AGROW>
    legendEntries{end+1} = sprintf('%s (with sensors)', results(t).name); %#ok<AGROW>
end

yline(1, 'r:', '1\sigma', 'LineWidth', 0.8, 'LabelHorizontalAlignment','left');
yline(3, 'b:', '3\sigma', 'LineWidth', 0.8, 'LabelHorizontalAlignment','left');

xlabel('Radial Distance (km)');
ylabel('\sigma_{risk} (Mahalanobis Distance)');
title(sprintf('Sensor Fusion Effect on Sigma Risk (%s)', Sigma_own_labels{j_sweep}));
legend(legendEntries, 'Location', 'eastoutside', 'FontSize', 7);
grid on; hold off;

fprintf('\nSensor fusion analysis complete.\n');

%% ======================== SWEEP: RSS DISTANCE (SENSOR-FUSED ONLY) ========================
figure('Name', 'Sigma Risk vs RSS Distance (Sensor-Fused)', ...
       'Position', figPos_std);
hold on;
legendEntries = {};

% Directions to sweep (unit vectors in RTN)
sweepDirs = {[1;0;0], [0;1;0], [0;0;1], [1;1;1]/sqrt(3)};
sweepLabels = {'Radial', 'Along-Track', 'Cross-Track', 'Equal RTN'};
sweepColors = [0.85 0.33 0.10;
               0.00 0.45 0.74;
               0.47 0.67 0.19;
               0.49 0.18 0.56];

dist_rss = linspace(0.1, 10, 200)';

for t = 1:nTargets
    Sig_target = results(t).Sigma_target;

    for j = 1:nSigmaOwn
        Sig_own = diag(Sigma_own_options(j,:));
        Sig_old = Sig_own + Sig_target;

        styleIdx = mod(j-1, numel(plotStyles)) + 1;

        for d = 1:numel(sweepDirs)
            uv = sweepDirs{d};
            sr_fused = zeros(size(dist_rss));

            for i = 1:numel(dist_rss)
                dr_i = dist_rss(i) * uv;
                Sig_new = applySensorFusion(Sig_old, dr_i, sensors);
                sr_fused(i) = sqrt(dr_i' / Sig_new * dr_i);
            end

            plot(dist_rss, sr_fused, plotStyles{styleIdx}, ...
                 'Color', sweepColors(d,:), 'LineWidth', 1.5);
            legendEntries{end+1} = sprintf('%s / %s / %s', ...
                results(t).name, Sigma_own_labels{j}, sweepLabels{d}); %#ok<AGROW>
        end
    end
end

yline(1, 'r:', '1\sigma', 'LineWidth', 0.8, 'LabelHorizontalAlignment','left');
yline(3, 'b:', '3\sigma', 'LineWidth', 0.8, 'LabelHorizontalAlignment','left');

xlabel('RSS Distance (km)');
ylabel('\sigma_{risk} (Mahalanobis Distance)');
title('Sigma Risk vs RSS Distance (Sensor-Fused)');
legend(legendEntries, 'Location', 'eastoutside', 'FontSize', 7);
grid on; hold off;

%% =====================================================================
%%  LOCAL FUNCTIONS
%% =====================================================================

function sr = sigmaRiskAtDist(d, unitVec, Sig_old, sensors)
%SIGMARISKATDIST  Compute sigma risk at distance d along a unit vector.
%   Used by fzero for the inverse problem (find d given desired sigma).
%
%   Inputs:
%     d       - scalar distance [km]
%     unitVec - 3x1 unit direction in RTN
%     Sig_old - 3x3 prior covariance [km^2]
%     sensors - sensor struct array
%
%   Output:
%     sr      - scalar Mahalanobis sigma risk

    dr = d * unitVec(:);
    Sig_new = applySensorFusion(Sig_old, dr, sensors);
    sr = sqrt(dr' / Sig_new * dr);
end


function s = fmtDist(d_km)
%FMTDIST  Format a distance as meters or km for table output.
    if isnan(d_km)
        s = '    N/A';
    elseif d_km < 1
        s = sprintf('%5.0f m', d_km * 1000);
    else
        s = sprintf('%4.2fkm', d_km);
    end
end


function Sig_new = applySensorFusion(Sig_old, r_rtn, sensors)
%APPLYSENSORFUSION  Information-filter covariance update with sensors.
%
%   Sig_new = ( Sig_old^-1 + sum_i( H_i' * R_i^-1 * H_i ) )^-1
%
%   Inputs:
%     Sig_old  - 3x3 prior covariance [km^2] in RTN
%     r_rtn    - 3x1 relative position [km] in RTN [R; T; N]
%     sensors  - struct array with fields: enabled, type, sigma
%
%   Output:
%     Sig_new  - 3x3 updated covariance [km^2] in RTN

    R_val = r_rtn(1);
    T_val = r_rtn(2);
    N_val = r_rtn(3);
    r = norm(r_rtn);

    % Accumulate information from all enabled sensors
    info_gain = zeros(3);   % sum of H' * R^-1 * H

    for si = 1:numel(sensors)
        if ~sensors(si).enabled
            continue;
        end

        switch sensors(si).type
            case 'range'
                % Rangefinder: measures range rho = ||r||
                % H = [R/r, T/r, N/r]  (1x3)
                % R_sensor = sigma_range^2  (scalar)

                if r < 1e-10, continue; end

                H = [R_val/r, T_val/r, N_val/r];  % 1x3
                R_sensor = sensors(si).sigma(1)^2;  % scalar
                info_gain = info_gain + H' * (1/R_sensor) * H;

            case 'lidar'
                % LIDAR: measures (range, azimuth, elevation)
                %   rho   = sqrt(R^2 + T^2 + N^2)
                %   theta = atan2(T, R)         (azimuth in R-T plane)
                %   phi   = atan2(N, sqrt(R^2+T^2))  (elevation)
                %
                % H is 3x3 Jacobian

                if r < 1e-10, continue; end

                rho_RT = sqrt(R_val^2 + T_val^2);

                if rho_RT < 1e-10
                    % Pure cross-track: azimuth/elevation singular
                    % Fall back to range-only contribution
                    H = [R_val/r, T_val/r, N_val/r];
                    R_sensor = sensors(si).sigma(1)^2;
                    info_gain = info_gain + H' * (1/R_sensor) * H;
                    continue;
                end

                H = zeros(3);

                % Row 1: d(rho)/d[R,T,N]
                H(1,:) = [R_val/r, T_val/r, N_val/r];

                % Row 2: d(theta)/d[R,T,N]
                H(2,:) = [-T_val/(rho_RT^2), R_val/(rho_RT^2), 0];

                % Row 3: d(phi)/d[R,T,N]
                H(3,:) = [-R_val*N_val/(r^2 * rho_RT), ...
                          -T_val*N_val/(r^2 * rho_RT), ...
                           rho_RT / r^2];

                R_sensor = diag(sensors(si).sigma.^2);  % 3x3
                info_gain = info_gain + H' / R_sensor * H;

            case 'angles_only'
                % Angles-only sensor: measures (azimuth, elevation)
                % Same as LIDAR rows 2-3, no range measurement.

                if r < 1e-10, continue; end

                rho_RT = sqrt(R_val^2 + T_val^2);

                if rho_RT < 1e-10
                    continue;   % singular at pure cross-track
                end

                H = zeros(2, 3);

                % Row 1: d(theta)/d[R,T,N]
                H(1,:) = [-T_val/(rho_RT^2), R_val/(rho_RT^2), 0];

                % Row 2: d(phi)/d[R,T,N]
                H(2,:) = [-R_val*N_val/(r^2 * rho_RT), ...
                          -T_val*N_val/(r^2 * rho_RT), ...
                           rho_RT / r^2];

                R_sensor = diag(sensors(si).sigma.^2);  % 2x2
                info_gain = info_gain + H' / R_sensor * H;

            otherwise
                warning('Unknown sensor type: %s', sensors(si).type);
        end
    end

    % Information filter update
    Sig_new = inv(inv(Sig_old) + info_gain);           %#ok<MINV>
end

