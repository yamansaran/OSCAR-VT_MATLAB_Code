%% J2 Relative RAAN Maneuver for Satellite Rendezvous
% Find altitude change to match target satellite's RAAN using differential J2 precession

clear; clc;

%% User Inputs - Chaser Satellite (YOU)
h_chaser = 854;              % Your current altitude [km]
inc_chaser = 98.55;          % Your inclination [deg]
ecc_chaser = 0.0;            % Your eccentricity [-]
RAAN_chaser_initial = 108.64; % Your current RAAN [deg]

%% User Inputs - Target Satellite
h_target_sat = 854;          % Target satellite altitude [km] (assume same for now)
inc_target = 98.55;          % Target inclination [deg] (assume same)
ecc_target = 0.0;            % Target eccentricity [-]
RAAN_target = 53.0;          % Target satellite RAAN [deg]

%% Maneuver Constraints
max_time = 365;              % Maximum time for maneuver [days]
thrust_accel = 1e-6;         % Thrust acceleration [km/s^2]
h_maneuver_min = 300;        % Minimum altitude during maneuver [km]
h_maneuver_max = 2000;       % Maximum altitude during maneuver [km]

%% Constants
mu = 398600.4418;            % Earth gravitational parameter [km^3/s^2]
J2 = 1.08263e-3;             % J2 coefficient [-]
R_E = 6378.137;              % Earth equatorial radius [km]

%% Conversions
i_chaser_rad = deg2rad(inc_chaser);
i_target_rad = deg2rad(inc_target);
a_chaser = R_E + h_chaser;
a_target_sat = R_E + h_target_sat;

%% Calculate Relative RAAN
% Shortest angular distance between the two planes
RAAN_diff = RAAN_target - RAAN_chaser_initial;

% Normalize to [-180, 180]
if RAAN_diff > 180
    RAAN_diff = RAAN_diff - 360;
elseif RAAN_diff < -180
    RAAN_diff = RAAN_diff + 360;
end

fprintf('======================================\n');
fprintf('RELATIVE RAAN CHANGE PROBLEM\n');
fprintf('======================================\n');
fprintf('Chaser (you):  RAAN = %.2f°, h = %.0f km\n', RAAN_chaser_initial, h_chaser);
fprintf('Target:        RAAN = %.2f°, h = %.0f km\n', RAAN_target, h_target_sat);
fprintf('Relative RAAN: %.2f°\n\n', RAAN_diff);

%% J2 Precession Rate Function
J2_precession_rate = @(a, i_rad, ecc) -1.5 * J2 * R_E^2 * sqrt(mu) * cos(i_rad) / ...
                                       (a^(7/2) * (1 - ecc^2)^2);

%% Natural drift rates
Omega_dot_chaser = J2_precession_rate(a_chaser, i_chaser_rad, ecc_chaser);
Omega_dot_target = J2_precession_rate(a_target_sat, i_target_rad, ecc_target);

% Relative drift rate (how fast planes are separating/converging)
Omega_dot_relative = Omega_dot_chaser - Omega_dot_target;

fprintf('NATURAL PRECESSION RATES:\n');
fprintf('Chaser rate:   %.6f deg/day\n', rad2deg(Omega_dot_chaser)*86400);
fprintf('Target rate:   %.6f deg/day\n', rad2deg(Omega_dot_target)*86400);
fprintf('Relative rate: %.6f deg/day\n\n', rad2deg(Omega_dot_relative)*86400);

% Check if planes are naturally converging or diverging
if abs(Omega_dot_relative) < 1e-10
    fprintf('⚠️  Planes have identical precession rates - no natural convergence!\n');
    fprintf('    Maneuver is required to create differential precession.\n\n');
elseif sign(Omega_dot_relative) == sign(RAAN_diff)
    fprintf('⚠️  Planes are naturally DIVERGING!\n');
    fprintf('    Without maneuver, youre moving AWAY from target.\n\n');
else
    t_natural_convergence = abs(deg2rad(RAAN_diff) / Omega_dot_relative / 86400);
    fprintf('✓ Planes are naturally CONVERGING!\n');
    fprintf('  Time to natural rendezvous: %.1f days\n\n', t_natural_convergence);
    
    if t_natural_convergence <= max_time
        fprintf('✓ No maneuver needed! Just wait %.1f days.\n', t_natural_convergence);
        return;
    else
        fprintf('⚠️  Natural convergence takes too long (%.1f days > %.0f days max)\n', ...
                t_natural_convergence, max_time);
        fprintf('    Maneuver needed to accelerate convergence.\n\n');
    end
end

%% Differential Precession Strategy
% We need to create a differential precession rate to close the RAAN gap
% Strategy: temporarily change altitude to create different precession rate

fprintf('COMPUTING OPTIMAL MANEUVER...\n\n');

% Function to calculate time to close RAAN gap at given maneuver altitude
calc_rendezvous_time = @(h_maneuver) calculate_rendezvous_time(h_maneuver, ...
    h_chaser, RAAN_diff, thrust_accel, i_chaser_rad, ecc_chaser, ...
    i_target_rad, a_target_sat, R_E, J2, mu);

% Search for altitude that minimizes total time
h_search = linspace(h_maneuver_min, h_maneuver_max, 200);
rendezvous_times = zeros(size(h_search));

for k = 1:length(h_search)
    rendezvous_times(k) = calc_rendezvous_time(h_search(k));
end

% Filter out invalid solutions (NaN or negative times)
valid_idx = isfinite(rendezvous_times) & rendezvous_times > 0;
if ~any(valid_idx)
    fprintf('❌ No valid solution found in altitude range [%.0f, %.0f] km\n', ...
            h_maneuver_min, h_maneuver_max);
    fprintf('Try: 1) Increasing max_time, 2) Expanding altitude range, or 3) Higher thrust\n');
    return;
end

[min_time, idx] = min(rendezvous_times(valid_idx));
h_search_valid = h_search(valid_idx);
h_optimal = h_search_valid(idx);

% Refine solution
options = optimset('Display', 'off');
h_optimal = fminbnd(calc_rendezvous_time, max(h_maneuver_min, h_optimal-50), ...
                    min(h_maneuver_max, h_optimal+50), options);

%% Calculate Detailed Results
[total_time, t_spiral_out, t_coast, t_spiral_back, RAAN_chaser_final, ...
 RAAN_during_spiral_out, RAAN_during_coast, RAAN_during_spiral_back] = ...
    calculate_rendezvous_details(h_optimal, h_chaser, RAAN_diff, thrust_accel, ...
    i_chaser_rad, ecc_chaser, i_target_rad, a_target_sat, R_E, J2, mu);

delta_h = h_optimal - h_chaser;
a_maneuver = R_E + h_optimal;

Omega_dot_chaser_natural = J2_precession_rate(a_chaser, i_chaser_rad, ecc_chaser);
Omega_dot_chaser_maneuver = J2_precession_rate(a_maneuver, i_chaser_rad, ecc_chaser);

% Calculate what target RAAN will be at rendezvous
RAAN_target_at_rendezvous = RAAN_target + rad2deg(Omega_dot_target * total_time);
RAAN_target_at_rendezvous = mod(RAAN_target_at_rendezvous, 360);

%% Display Results
fprintf('======================================\n');
fprintf('OPTIMAL RENDEZVOUS SOLUTION\n');
fprintf('======================================\n\n');

fprintf('MANEUVER ALTITUDE:\n');
fprintf('  Current altitude:     %.2f km\n', h_chaser);
fprintf('  Maneuver altitude:    %.2f km\n', h_optimal);
fprintf('  Altitude change:      %+.2f km\n', delta_h);
fprintf('  Total time:           %.2f days\n\n', total_time/86400);

fprintf('TIME BREAKDOWN:\n');
fprintf('  Spiral out:           %.2f days (%.1f%%)\n', ...
        t_spiral_out/86400, 100*t_spiral_out/total_time);
fprintf('  Coast (differential): %.2f days (%.1f%%)\n', ...
        t_coast/86400, 100*t_coast/total_time);
fprintf('  Spiral back:          %.2f days (%.1f%%)\n\n', ...
        t_spiral_back/86400, 100*t_spiral_back/total_time);

fprintf('RAAN EVOLUTION (Chaser):\n');
fprintf('  Initial:              %.2f°\n', RAAN_chaser_initial);
fprintf('  After spiral out:     %.2f°\n', RAAN_chaser_initial + rad2deg(RAAN_during_spiral_out));
fprintf('  After coast:          %.2f°\n', RAAN_chaser_initial + rad2deg(RAAN_during_spiral_out + RAAN_during_coast));
fprintf('  Final (at rendezvous):%.2f°\n', mod(RAAN_chaser_final, 360));
fprintf('  Total change:         %.2f°\n\n', rad2deg(RAAN_during_spiral_out + RAAN_during_coast + RAAN_during_spiral_back));

fprintf('TARGET RAAN:\n');
fprintf('  Initial:              %.2f°\n', RAAN_target);
fprintf('  At rendezvous:        %.2f°\n', RAAN_target_at_rendezvous);
fprintf('  Target drift:         %.2f°\n\n', RAAN_target_at_rendezvous - RAAN_target);

fprintf('VERIFICATION:\n');
fprintf('  Chaser RAAN at end:   %.2f°\n', mod(RAAN_chaser_final, 360));
fprintf('  Target RAAN at end:   %.2f°\n', RAAN_target_at_rendezvous);
fprintf('  Residual error:       %.4f°\n\n', mod(RAAN_chaser_final - RAAN_target_at_rendezvous + 180, 360) - 180);

fprintf('PRECESSION RATES:\n');
fprintf('  Chaser at %.0f km:    %.6f deg/day\n', h_chaser, rad2deg(Omega_dot_chaser_natural)*86400);
fprintf('  Chaser at %.0f km:    %.6f deg/day\n', h_optimal, rad2deg(Omega_dot_chaser_maneuver)*86400);
fprintf('  Target (constant):    %.6f deg/day\n', rad2deg(Omega_dot_target)*86400);
fprintf('  Differential rate:    %.6f deg/day\n\n', rad2deg(Omega_dot_chaser_maneuver - Omega_dot_target)*86400);

fprintf('MANEUVER SEQUENCE:\n');
fprintf('  Day 0:   Start at h=%.0f km, RAAN=%.2f°\n', h_chaser, RAAN_chaser_initial);
fprintf('  Day %.1f: Reach h=%.0f km via spiral\n', t_spiral_out/86400, h_optimal);
fprintf('  Day %.1f: Complete coast at h=%.0f km\n', (t_spiral_out+t_coast)/86400, h_optimal);
fprintf('  Day %.1f: Return to h=%.0f km, RAAN≈%.2f° (matches target!)\n', ...
          total_time/86400, h_chaser, mod(RAAN_chaser_final, 360));
fprintf('======================================\n');

%% Visualization
figure('Position', [100 100, 1400, 900]);

% Plot 1: Rendezvous time vs maneuver altitude
subplot(2,2,1);
plot(h_search(valid_idx), rendezvous_times(valid_idx)/86400, 'b-', 'LineWidth', 2);
hold on;
plot(h_optimal, min_time/86400, 'ro', 'MarkerSize', 12, 'LineWidth', 2);
grid on;
xlabel('Maneuver Altitude [km]');
ylabel('Total Rendezvous Time [days]');
title('Rendezvous Time vs Maneuver Altitude');
legend('Total Time', 'Optimal Solution', 'Location', 'best');

% Plot 2: Altitude profile
subplot(2,2,2);
time_profile = [0, t_spiral_out, t_spiral_out+t_coast, total_time]/86400;
alt_profile = [h_chaser, h_optimal, h_optimal, h_chaser];
plot(time_profile, alt_profile, 'b-', 'LineWidth', 2);
grid on;
xlabel('Time [days]');
ylabel('Altitude [km]');
title('Chaser Altitude Profile');
xlim([0 total_time/86400]);

% Plot 3: RAAN evolution for both satellites
subplot(2,2,3);
% Chaser RAAN
RAAN_chaser_profile = [RAAN_chaser_initial, ...
                       RAAN_chaser_initial + rad2deg(RAAN_during_spiral_out), ...
                       RAAN_chaser_initial + rad2deg(RAAN_during_spiral_out + RAAN_during_coast), ...
                       RAAN_chaser_final];

% Target RAAN (linear drift)
RAAN_target_profile = RAAN_target + rad2deg(Omega_dot_target) * [0, t_spiral_out, t_spiral_out+t_coast, total_time];

plot(time_profile, RAAN_chaser_profile, 'b-', 'LineWidth', 2);
hold on;
plot(time_profile, RAAN_target_profile, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time [days]');
ylabel('RAAN [deg]');
title('RAAN Evolution: Chaser vs Target');
legend('Chaser (you)', 'Target Satellite', 'Location', 'best');
xlim([0 total_time/86400]);

% Plot 4: Differential precession rate vs altitude
subplot(2,2,4);
h_range = linspace(h_maneuver_min, h_maneuver_max, 500);
a_range = R_E + h_range;
Omega_dot_chaser_range = arrayfun(@(a) J2_precession_rate(a, i_chaser_rad, ecc_chaser), a_range);
diff_rate = Omega_dot_chaser_range - Omega_dot_target;

plot(h_range, rad2deg(diff_rate)*86400, 'b-', 'LineWidth', 2);
hold on;
yline(0, 'k--', 'No Differential', 'LineWidth', 1);
plot(h_chaser, rad2deg(Omega_dot_chaser_natural - Omega_dot_target)*86400, ...
     'gs', 'MarkerSize', 12, 'LineWidth', 2);
plot(h_optimal, rad2deg(Omega_dot_chaser_maneuver - Omega_dot_target)*86400, ...
     'ro', 'MarkerSize', 12, 'LineWidth', 2);
grid on;
xlabel('Chaser Altitude [km]');
ylabel('Differential RAAN Rate [deg/day]');
title('Differential Precession Rate (Chaser - Target)');
legend('Diff Rate', 'Zero Line', sprintf('Natural (%.0f km)', h_chaser), ...
       sprintf('Maneuver (%.0f km)', h_optimal), 'Location', 'best');

%% Helper Functions

function time = calculate_rendezvous_time(h_maneuver, h_initial, RAAN_diff_deg, ...
    thrust_accel, i_chaser_rad, ecc_chaser, i_target_rad, a_target, R_E, J2, mu)
    
    [time, ~, ~, ~, ~, ~, ~, ~] = calculate_rendezvous_details(h_maneuver, ...
        h_initial, RAAN_diff_deg, thrust_accel, i_chaser_rad, ecc_chaser, ...
        i_target_rad, a_target, R_E, J2, mu);
end

function [total_time, t_out, t_coast, t_back, RAAN_chaser_final, ...
          RAAN_out, RAAN_coast, RAAN_back] = ...
    calculate_rendezvous_details(h_maneuver, h_initial, RAAN_diff_deg, ...
    thrust_accel, i_chaser_rad, ecc_chaser, i_target_rad, a_target, R_E, J2, mu)
    
    a_initial = R_E + h_initial;
    a_maneuver = R_E + h_maneuver;
    delta_h = h_maneuver - h_initial;
    
    % J2 precession rate
    Omega_dot = @(a, i_rad, ecc) -1.5 * J2 * R_E^2 * sqrt(mu) * cos(i_rad) / ...
                                  (a^(7/2) * (1 - ecc^2)^2);
    
    % Target precession rate (constant)
    Omega_dot_target = Omega_dot(a_target, i_target_rad, 0);
    
    % Transfer times
    delta_v = abs(sqrt(mu) * delta_h / (a_initial^(3/2)));
    t_transfer = delta_v / thrust_accel;
    
    t_out = t_transfer;
    t_back = t_transfer;
    
    % Integrate RAAN during spiral out
    n_steps = 100;
    a_spiral_out = linspace(a_initial, a_maneuver, n_steps);
    dt_out = t_out / (n_steps - 1);
    RAAN_out = 0;
    RAAN_target_during_out = 0;
    for k = 1:n_steps-1
        a_avg = (a_spiral_out(k) + a_spiral_out(k+1)) / 2;
        RAAN_out = RAAN_out + Omega_dot(a_avg, i_chaser_rad, ecc_chaser) * dt_out;
        RAAN_target_during_out = RAAN_target_during_out + Omega_dot_target * dt_out;
    end
    
    % Integrate RAAN during spiral back
    a_spiral_back = linspace(a_maneuver, a_initial, n_steps);
    dt_back = t_back / (n_steps - 1);
    RAAN_back = 0;
    RAAN_target_during_back = 0;
    for k = 1:n_steps-1
        a_avg = (a_spiral_back(k) + a_spiral_back(k+1)) / 2;
        RAAN_back = RAAN_back + Omega_dot(a_avg, i_chaser_rad, ecc_chaser) * dt_back;
        RAAN_target_during_back = RAAN_target_during_back + Omega_dot_target * dt_back;
    end
    
    % Calculate required coast time to close RAAN gap
    Omega_dot_maneuver = Omega_dot(a_maneuver, i_chaser_rad, ecc_chaser);
    differential_rate = Omega_dot_maneuver - Omega_dot_target;
    
    % Remaining RAAN to close after spirals
    RAAN_diff_rad = deg2rad(RAAN_diff_deg);
    RAAN_closed_during_spirals = (RAAN_out - RAAN_target_during_out) + ...
                                  (RAAN_back - RAAN_target_during_back);
    RAAN_remaining = RAAN_diff_rad - RAAN_closed_during_spirals;
    
    % Coast time needed
    if abs(differential_rate) < 1e-12
        t_coast = inf; % No differential precession at this altitude
    else
        t_coast = RAAN_remaining / differential_rate;
    end
    
    if t_coast < 0
        t_coast = NaN; % Wrong direction
    end
    
    % RAAN changes during coast
    RAAN_coast = differential_rate * t_coast;
    
    % Total time
    total_time = t_out + t_coast + t_back;
    
    % Final chaser RAAN
    RAAN_chaser_final = RAAN_out + RAAN_coast + RAAN_back;
end