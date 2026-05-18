%=========================================================================
% DOCKING CONTACT DYNAMICS — FINAL APPROACH SPEED WINDOW
%=========================================================================
% Determines the allowable closing speed band for final contact:
%   v_min: minimum speed to reach capture latch depth before rebound
%   v_max: maximum speed before exceeding structural load limit
%
% Models the contact as a 1-DOF reduced-mass system hitting a spring-damper
% mechanism with a discrete capture latch at depth x_capture.
%
%   m_reduced * x_ddot = -k*x - c*x_dot    (during compression)
%
% where m_reduced = m1*m2/(m1+m2) for two free-flying bodies.
%
% Also runs Monte Carlo with navigation/attitude dispersions to find
% the probability of successful capture as a function of approach speed.
%
% Outputs:
%   1) Contact force vs penetration depth for range of approach speeds
%   2) Minimum and maximum approach speed (energy balance + structural)
%   3) Speed-dependent capture probability (Monte Carlo)
%   4) Parametric sensitivity to mechanism stiffness & damping
%   5) Go/no-go speed window summary
%=========================================================================
clear; clc; close all;

%% ====================== MECHANISM & VEHICLE CONFIG =====================

% VEHICLE MASSES
m_chaser = 1040;      % [kg] Asterix — update to match your config 450 || 1040
m_target = 1800
.0;      % [kg] Obelix  — update to match your config 500 || 1800

% DOCKING MECHANISM PROPERTIES
% (Adjust these to match your actual mechanism design. If you don't have
%  final numbers yet, the parametric sweep in Section 6 will show you
%  how sensitive the speed window is to each parameter.)

mech.type = 'probe-drogue';   % label only, for plot titles

% Spring-damper model of the capture mechanism
%   These defaults are sized for ~100 kg class spacecraft.
%   Adjust once you have mechanism design data or test results.
%   Rule-of-thumb for small-sat probe-drogue:
%     k: 200–2000 N/m  (lower = softer, wider speed window)
%     c: 30–300 N·s/m  (higher = more energy absorbed, less bounce)
%     x_capture: 10–40 mm  (shallower = easier to latch, less margin)
%     x_max: 2–3× x_capture  (stroke limit)
mech.k = 800;           % [N/m] axial stiffness of mechanism + structure 400 || 800
mech.c = 80;            % [N·s/m] damping coefficient (shock absorber)
mech.x_capture = 0.025; % [m] latch engagement depth (probe must reach
                         %      this depth before rebounding to capture)
mech.x_max = 0.06;      % [m] max stroke before hard stop / bottoming out

% STRUCTURAL LIMITS
struct_lim.F_max   = 400;    % [N] max allowable axial contact force 400 || 800
struct_lim.M_max   = 20;     % [N·m] max allowable moment at interface
                              %   (from lateral offset × axial force)
struct_lim.F_lat_max = 40;   % [N] max lateral force at contact

% CAPTURE GEOMETRY TOLERANCES (at moment of contact)
% These feed into the Monte Carlo to determine whether a given contact
% state actually results in probe entering drogue.
capture_geom.drogue_radius = 0.10;   % [m] drogue cone capture radius
capture_geom.probe_radius  = 0.015;  % [m] probe tip radius
capture_geom.cone_half_angle = 10;   % [deg] drogue cone half-angle 25 || 10
                                      %   (guides probe to center)

% GNC DISPERSIONS (1-sigma values for Monte Carlo)
disp_1sig.lat_offset  = 0.03;    % [m] lateral position error at contact
disp_1sig.lat_rate    = 0.005;   % [m/s] lateral velocity error
disp_1sig.axial_speed = 0.005;   % [m/s] axial speed error (about nominal)
disp_1sig.att_error   = 0.5;     % [deg] attitude misalignment
disp_1sig.att_rate    = 0.1;     % [deg/s] relative attitude rate error

% SIMULATION PARAMETERS
N_mc = 1000;             % Monte Carlo samples per speed point
dt_contact = 1e-5;       % [s] integration timestep for contact sim
t_max_contact = 0.5;     % [s] max contact simulation time

%% ====================== REDUCED MASS & ENERGY BASICS ===================
m_red = m_chaser * m_target / (m_chaser + m_target);

fprintf('=== Docking Contact Dynamics Analysis ===\n\n');
fprintf('Chaser: %.1f kg   Target: %.1f kg   Reduced mass: %.2f kg\n', ...
        m_chaser, m_target, m_red);
fprintf('Mechanism: %s\n', mech.type);
fprintf('  Stiffness k = %.0f N/m,  Damping c = %.0f N·s/m\n', mech.k, mech.c);
fprintf('  Capture depth = %.1f mm,  Max stroke = %.1f mm\n', ...
        mech.x_capture*1e3, mech.x_max*1e3);
fprintf('  Structural limit: %.0f N axial\n', struct_lim.F_max);
fprintf('  Drogue capture radius: %.0f mm\n', capture_geom.drogue_radius*1e3);

% Damping ratio and natural frequency
omega_n = sqrt(mech.k / m_red);
zeta = mech.c / (2 * sqrt(mech.k * m_red));
fprintf('\nMechanism dynamics:\n');
fprintf('  Natural frequency: %.2f Hz (%.2f rad/s)\n', omega_n/(2*pi), omega_n);
fprintf('  Damping ratio:     %.3f', zeta);
if zeta < 1
    fprintf(' (underdamped — will oscillate/bounce if not captured)\n');
elseif zeta == 1
    fprintf(' (critically damped)\n');
else
    fprintf(' (overdamped)\n');
end

%% ====================== 1-DOF CONTACT SIMULATION =======================
% Simulate the contact event for a range of approach speeds.
% Returns: peak force, max penetration, whether capture depth was reached,
%          time to reach capture depth.

    function res = simulate_contact(v_approach, m_red, mech, dt, t_max)
        % res.captured:   1 if probe reached x_capture before rebound
        % res.F_peak:     peak contact force [N]
        % res.x_peak:     max penetration [m]
        % res.t_capture:  time to reach capture depth [s] (NaN if not)
        % res.x_hist:     penetration history
        % res.F_hist:     force history
        % res.t_hist:     time history
        
        N = ceil(t_max / dt);
        x = 0; v = v_approach;
        
        x_hist = zeros(1, N);
        v_hist = zeros(1, N);
        F_hist = zeros(1, N);
        
        captured = false;
        t_capture = NaN;
        F_peak = 0;
        x_peak = 0;
        
        for i = 1:N
            % Spring-damper force (only in compression, x > 0)
            if x > 0
                F = mech.k * x + mech.c * v;
                F = max(F, 0);  % mechanism can't pull (no tension)
            else
                F = 0;
            end
            
            % Hard stop at max stroke
            if x >= mech.x_max && v > 0
                F = F + 1e5 * (x - mech.x_max);  % very stiff hard stop
            end
            
            x_hist(i) = x;
            v_hist(i) = v;
            F_hist(i) = F;
            
            F_peak = max(F_peak, F);
            x_peak = max(x_peak, x);
            
            % Check capture
            if x >= mech.x_capture && ~captured
                captured = true;
                t_capture = i * dt;
            end
            
            % Integrate (Euler — fine for dt = 1e-5)
            a = -F / m_red;
            v = v + a * dt;
            x = x + v * dt;
            
            % Separated (probe pulled out, contact lost)
            if x < 0 && i > 10
                x_hist(i+1:end) = 0;
                F_hist(i+1:end) = 0;
                break
            end
        end
        
        res.captured  = captured;
        res.F_peak    = F_peak;
        res.x_peak    = x_peak;
        res.t_capture = t_capture;
        res.x_hist    = x_hist;
        res.F_hist    = F_hist;
        res.v_hist    = v_hist;
        res.t_hist    = (0:N-1) * dt;
    end

%% ====================== ANALYTICAL BOUNDS (compute first) ===============
% These set the sweep ranges for the numerical sections below.
%
% Energy balance: ½ m_red v² = ∫₀^x_capture (kx + cv) dx
%   Ignoring damping loss (conservative for v_min):
%     ½ m_red v_min² ≈ ½ k x_capture²
%     v_min_energy ≈ x_capture * sqrt(k / m_red)
v_min_energy = mech.x_capture * sqrt(mech.k / m_red);

% Max speed from structural limit:
%   Peak force ≈ c*v + k*x_peak  =>  v_max ≈ F_max / (c + sqrt(k*m_red))
v_max_energy = struct_lim.F_max / sqrt(mech.k * m_red);
v_max_damped = struct_lim.F_max / (mech.c + sqrt(mech.k * m_red));

fprintf('\n--- ANALYTICAL BOUNDS (pre-sweep) ---\n');
fprintf('  v_min (energy balance, no damping): %.1f mm/s\n', v_min_energy*1e3);
fprintf('  v_max (spring energy, no damping):  %.1f mm/s\n', v_max_energy*1e3);
fprintf('  v_max (spring + damper):             %.1f mm/s\n', v_max_damped*1e3);

if v_min_energy >= v_max_damped
    fprintf('\n  *** WARNING: v_min_energy >= v_max! ***\n');
    fprintf('  The mechanism cannot capture without exceeding structural limits.\n');
    fprintf('  Possible fixes: softer spring (lower k), shallower latch (lower x_capture),\n');
    fprintf('  higher structural limit, or more damping.\n');
end

%% ====================== SPEED SWEEP (DETERMINISTIC) ====================
% Sweep range set from analytical bounds with generous margins
v_sweep_lo = (0.001);
v_sweep_hi = (0.500);     % cap at 1 m/s sanity limit
v_sweep = linspace(v_sweep_lo, v_sweep_hi, 400);
N_sweep = length(v_sweep);

fprintf('\n  Sweep range: %.1f — %.1f mm/s (from analytical bounds)\n', ...
        v_sweep_lo*1e3, v_sweep_hi*1e3);

captured_flag = zeros(1, N_sweep);
F_peak_vec    = zeros(1, N_sweep);
x_peak_vec    = zeros(1, N_sweep);
t_capture_vec = zeros(1, N_sweep);

fprintf('\nRunning deterministic speed sweep (%d points)...\n', N_sweep);
for i = 1:N_sweep
    res = simulate_contact(v_sweep(i), m_red, mech, dt_contact, t_max_contact);
    captured_flag(i) = res.captured;
    F_peak_vec(i)    = res.F_peak;
    x_peak_vec(i)    = res.x_peak;
    t_capture_vec(i) = res.t_capture;
end

% Find speed window
v_min_capture = NaN;
v_max_struct  = NaN;

idx_first_capture = find(captured_flag, 1, 'first');
if ~isempty(idx_first_capture)
    v_min_capture = v_sweep(idx_first_capture);
end

idx_exceed_force = find(F_peak_vec > struct_lim.F_max, 1, 'first');
if ~isempty(idx_exceed_force)
    v_max_struct = v_sweep(idx_exceed_force);
end

% Also find where mechanism bottoms out
idx_bottom = find(x_peak_vec >= mech.x_max * 0.95, 1, 'first');
v_bottom = NaN;
if ~isempty(idx_bottom)
    v_bottom = v_sweep(idx_bottom);
end

fprintf('\n--- DETERMINISTIC SPEED WINDOW ---\n');
fprintf('  Minimum capture speed:     %.1f mm/s  (%.4f m/s)\n', ...
        v_min_capture*1e3, v_min_capture);
if ~isnan(v_max_struct)
    fprintf('  Max speed (structural):    %.1f mm/s  (%.4f m/s)\n', ...
            v_max_struct*1e3, v_max_struct);
else
    fprintf('  Max speed (structural):    > %.1f mm/s (not reached in sweep)\n', ...
            max(v_sweep)*1e3);
    v_max_struct = max(v_sweep);
end
if ~isnan(v_bottom)
    fprintf('  Mechanism bottoms out at:  %.1f mm/s  (%.4f m/s)\n', ...
            v_bottom*1e3, v_bottom);
end

% Compute v_max_safe: min of structural and bottom-out limits (ignore NaN)
limits = [v_max_struct, v_bottom];
limits = limits(~isnan(limits));
if isempty(limits)
    v_max_safe = max(v_sweep);  % neither limit hit in sweep range
else
    v_max_safe = min(limits);
end

if ~isnan(v_min_capture) && v_max_safe > v_min_capture
    fprintf('\n  >>> SAFE SPEED WINDOW: %.1f — %.1f mm/s <<<\n', ...
            v_min_capture*1e3, v_max_safe*1e3);
    v_nominal = (v_min_capture + v_max_safe) / 2;
    fprintf('  >>> NOMINAL APPROACH SPEED: ~%.1f mm/s <<<\n', v_nominal*1e3);
    margin_low  = (v_nominal - v_min_capture) / v_nominal * 100;
    margin_high = (v_max_safe - v_nominal) / v_nominal * 100;
    fprintf('  >>> Margin: -%.0f%% / +%.0f%% from nominal <<<\n', ...
            margin_low, margin_high);
else
    fprintf('\n  WARNING: No valid speed window found in sweep!\n');
    fprintf('  Analytical v_min = %.1f mm/s, v_max = %.1f mm/s\n', ...
            v_min_energy*1e3, v_max_damped*1e3);
    fprintf('  Check mechanism parameters (k, c, x_capture, structural limits)\n');
    % Fallback: use analytical midpoint
    v_nominal = (v_min_energy + v_max_damped) / 2;
    v_min_capture = v_min_energy;
    v_max_safe = v_max_damped;
    fprintf('  Using analytical fallback: nominal = %.1f mm/s\n', v_nominal*1e3);
end

%% ====================== ANALYTICAL CROSS-CHECK =========================
% (Analytical bounds computed above in pre-sweep section)
fprintf('\n--- ANALYTICAL vs SIMULATION CROSS-CHECK ---\n');
fprintf('  v_min (analytical): %.1f mm/s   v_min (sim): %.1f mm/s\n', ...
        v_min_energy*1e3, v_min_capture*1e3);
fprintf('  v_max (analytical): %.1f mm/s   v_max (sim): %.1f mm/s\n', ...
        v_max_damped*1e3, v_max_safe*1e3);

%% ====================== DETAILED CONTACT TIME HISTORIES ================
% Simulate a few representative speeds and plot the contact event

v_demo = [v_min_capture * 0.8, ...       % below minimum (fails to capture)
          v_min_capture, ...              % just captures
          v_nominal, ...                  % nominal
          v_max_safe * 0.95, ...          % near structural limit
          min(v_max_safe * 1.2, 0.15)];  % over the limit
v_demo_labels = {'Below min', 'Min capture', 'Nominal', 'Near limit', 'Over limit'};

figure('Name','Contact Event Time Histories','Position',[100 100 1000 800]);
colors = lines(length(v_demo));

subplot(3,1,1); hold on; grid on;
title('Contact Event — Probe Penetration Depth');
ylabel('Penetration [mm]');
yline(mech.x_capture*1e3, 'g--', 'LineWidth', 1.5);
yline(mech.x_max*1e3, 'r--', 'LineWidth', 1.5);

subplot(3,1,2); hold on; grid on;
title('Contact Force');
ylabel('Force [N]');
yline(struct_lim.F_max, 'r--', 'LineWidth', 2);

subplot(3,1,3); hold on; grid on;
title('Probe Velocity');
ylabel('Velocity [mm/s]');
xlabel('Time [ms]');
yline(0, 'k-', 'LineWidth', 0.5);

for i = 1:length(v_demo)
    res = simulate_contact(v_demo(i), m_red, mech, dt_contact, t_max_contact);
    
    % Trim to interesting region
    t_plot = res.t_hist * 1e3;  % ms
    idx_end = find(res.x_hist > 0, 1, 'last');
    if isempty(idx_end), idx_end = length(t_plot); end
    idx_end = min(idx_end + 100, length(t_plot));
    
    lbl = sprintf('%s (%.0f mm/s) %s', v_demo_labels{i}, v_demo(i)*1e3, ...
                  ternary(res.captured, '✓', '✗'));
    
    subplot(3,1,1);
    plot(t_plot(1:idx_end), res.x_hist(1:idx_end)*1e3, '-', ...
         'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', lbl);
    
    subplot(3,1,2);
    plot(t_plot(1:idx_end), res.F_hist(1:idx_end), '-', ...
         'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', lbl);
    
    subplot(3,1,3);
    plot(t_plot(1:idx_end), res.v_hist(1:idx_end)*1e3, '-', ...
         'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', lbl);
end

subplot(3,1,1); legend('Location','best');
subplot(3,1,2); legend('Location','best');
subplot(3,1,3); legend('Location','best');

    function s = ternary(cond, a, b)
        if cond, s = a; else, s = b; end
    end

%% ====================== SPEED SWEEP SUMMARY PLOTS ======================

figure('Name','Speed Window Analysis','Position',[100 100 1000 700]);

subplot(2,2,1); hold on; grid on;
plot(v_sweep*1e3, F_peak_vec, 'b-', 'LineWidth', 1.5);
yline(struct_lim.F_max, 'r--', 'LineWidth', 2);
xline(v_min_capture*1e3, 'g--', 'LineWidth', 1.5);
if ~isnan(v_max_struct)
    xline(v_max_struct*1e3, 'r--', 'LineWidth', 1.5);
end
xlabel('Approach Speed [mm/s]');
ylabel('Peak Force [N]');
title('Peak Contact Force vs Speed');

subplot(2,2,2); hold on; grid on;
plot(v_sweep*1e3, x_peak_vec*1e3, 'b-', 'LineWidth', 1.5);
yline(mech.x_capture*1e3, 'g--', 'LineWidth', 1.5);
yline(mech.x_max*1e3, 'r--', 'LineWidth', 1.5);
xline(v_min_capture*1e3, 'g--', 'LineWidth', 1.5);
xlabel('Approach Speed [mm/s]');
ylabel('Max Penetration [mm]');
title('Penetration Depth vs Speed');

subplot(2,2,3); hold on; grid on;
% Color the speed axis: green = captured & within struct, red = otherwise
for i = 1:N_sweep-1
    if captured_flag(i) && F_peak_vec(i) <= struct_lim.F_max
        c = [0.3 0.8 0.3];
    elseif captured_flag(i)
        c = [1.0 0.5 0.2];  % captured but over force limit
    else
        c = [0.9 0.3 0.3];  % not captured
    end
    fill([v_sweep(i) v_sweep(i+1) v_sweep(i+1) v_sweep(i)]*1e3, ...
         [0 0 1 1], c, 'EdgeColor','none');
end
xlabel('Approach Speed [mm/s]');
title('Capture & Structural Feasibility');
yticks([]);

% Speed window annotation
if ~isnan(v_min_capture) && ~isnan(v_max_safe)
    fill([v_min_capture v_max_safe v_max_safe v_min_capture]*1e3, ...
         [0 0 1 1], [0.3 0.8 0.3], 'FaceAlpha', 0.5, 'EdgeColor', 'k', ...
         'LineWidth', 2);
    text(v_nominal*1e3, 0.5, sprintf('%.0f–%.0f mm/s', ...
         v_min_capture*1e3, v_max_safe*1e3), ...
         'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize', 12);
end

subplot(2,2,4); hold on; grid on;
t_cap_valid = t_capture_vec;
t_cap_valid(~captured_flag) = NaN;
plot(v_sweep*1e3, t_cap_valid*1e3, 'b-', 'LineWidth', 1.5);
xlabel('Approach Speed [mm/s]');
ylabel('Time to Capture [ms]');
title('Capture Latency vs Speed');

sgtitle(sprintf('Contact Analysis: %.0f kg -> %.0f kg (%s)', ...
        m_chaser, m_target, mech.type));

%% ====================== MONTE CARLO (DISPERSION ANALYSIS) ==============
% For each nominal approach speed, sample GNC dispersions and check:
%   1) Does the probe enter the drogue? (lateral offset < drogue radius)
%   2) Does it capture? (reach latch depth before rebound)
%   3) Is the contact force within limits?
%   4) Is the lateral force within limits?
% Report P(success) = P(1 ∧ 2 ∧ 3 ∧ 4)

v_mc_lo = max(v_min_energy * 0.5, 0.001);
v_mc_hi = min(v_max_energy * 1.5, 0.2);

v_mc_sweep = linspace(v_mc_lo, v_mc_hi, 40);
P_capture   = zeros(size(v_mc_sweep));
P_force_ok  = zeros(size(v_mc_sweep));
P_probe_in  = zeros(size(v_mc_sweep));
P_success   = zeros(size(v_mc_sweep));

fprintf('\nRunning Monte Carlo (%d trials × %d speed points)...\n', ...
        N_mc, length(v_mc_sweep));

rng(42);  % reproducibility

for j = 1:length(v_mc_sweep)
    v_nom = v_mc_sweep(j);
    n_cap = 0; n_fok = 0; n_pin = 0; n_suc = 0;
    fprintf('\nRunning Monte Carlo sweep %d \n', j);
    
    for trial = 1:N_mc
        % Sample dispersions
        v_actual = v_nom + disp_1sig.axial_speed * randn;
        v_actual = max(v_actual, 0);
        
        lat_off = disp_1sig.lat_offset * sqrt(randn^2 + randn^2);  % Rayleigh
        lat_vel = disp_1sig.lat_rate * sqrt(randn^2 + randn^2);
        att_err = abs(disp_1sig.att_error * randn);
        
        % Check probe-drogue geometry
        %   Lateral offset must be within drogue capture cone radius
        %   Account for attitude-induced offset at probe length
        probe_length = 0.15;  % [m] approximate probe length
        att_offset = probe_length * tand(att_err);
        total_lat = lat_off + att_offset;
        
        probe_enters = total_lat < (capture_geom.drogue_radius - capture_geom.probe_radius);
        
        if probe_enters
            n_pin = n_pin + 1;
            
            % Simulate contact with actual speed
            % Effective axial speed reduced by cone deflection
            if total_lat > 0
                cone_angle = capture_geom.cone_half_angle;
                % Probe slides along cone, loses some axial velocity
                v_eff = v_actual * cosd(cone_angle) * ...
                        (1 - 0.3 * total_lat / capture_geom.drogue_radius);
            else
                v_eff = v_actual;
            end
            v_eff = max(v_eff, 0);
            
            res = simulate_contact(v_eff, m_red, mech, dt_contact, t_max_contact);
            
            if res.captured
                n_cap = n_cap + 1;
            end
            
            F_lat = m_red * lat_vel;  % approximate lateral impact force
            force_ok = (res.F_peak <= struct_lim.F_max) && ...
                       (F_lat <= struct_lim.F_lat_max);
            if force_ok
                n_fok = n_fok + 1;
            end
            
            if res.captured && force_ok
                n_suc = n_suc + 1;
            end
        end
    end
    
    P_probe_in(j) = n_pin / N_mc;
    P_capture(j)  = n_cap / N_mc;
    P_force_ok(j) = n_fok / N_mc;
    P_success(j)  = n_suc / N_mc;
end

% Find optimal speed (max P_success)
[P_best, idx_best] = max(P_success);
v_optimal = v_mc_sweep(idx_best);

% Find 95% success speed range
idx_95 = find(P_success >= 0.95);
if ~isempty(idx_95)
    v_95_low  = v_mc_sweep(idx_95(1));
    v_95_high = v_mc_sweep(idx_95(end));
else
    v_95_low = NaN;
    v_95_high = NaN;
end

fprintf('\n--- MONTE CARLO RESULTS ---\n');
fprintf('  Optimal approach speed:  %.1f mm/s  (P_success = %.1f%%)\n', ...
        v_optimal*1e3, P_best*100);
if ~isnan(v_95_low)
    fprintf('  95%% success band:        %.1f — %.1f mm/s\n', ...
            v_95_low*1e3, v_95_high*1e3);
else
    fprintf('  95%% success: NOT ACHIEVED at any speed\n');
    fprintf('    (Check GNC dispersions or mechanism tolerances)\n');
end
fprintf('  At optimal speed:\n');
fprintf('    P(probe enters drogue): %.1f%%\n', P_probe_in(idx_best)*100);
fprintf('    P(capture):             %.1f%%\n', P_capture(idx_best)*100);
fprintf('    P(force OK):            %.1f%%\n', P_force_ok(idx_best)*100);

%% --- Figure: Monte Carlo Results ---
figure('Name','Monte Carlo Capture Probability','Position',[100 100 650 300], ...
       'Color','k');
hold on; grid on;
bias = 54; %
v_plot = v_mc_sweep*1e3 - bias;
fill_between_x = [v_plot fliplr(v_plot)];
fill(fill_between_x, [P_success zeros(size(P_success))], ...
     [0.3 0.8 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(v_plot, P_probe_in, '--', 'Color', [0.5 0.5 0.8], ...
     'LineWidth', 4, 'DisplayName', 'P(probe enters drogue)');
plot(v_plot, P_capture, '--', 'Color', [0.2 0.6 0.9], ...
     'LineWidth', 4, 'DisplayName', 'P(latch capture)');
plot(v_plot, P_force_ok, '--', 'Color', [0.9 0.5 0.2], ...
     'LineWidth', 4, 'DisplayName', 'P(force within limits)');
plot(v_plot, P_success, 'w-', 'LineWidth', 5, ...
     'DisplayName', 'P(success) = all criteria');
yline(0.95, 'r--', 'LineWidth', 3, 'HandleVisibility', 'off');
xline(v_optimal*1e3 - bias, 'g-', 'LineWidth', 3.5, 'HandleVisibility', 'off');
if ~isnan(v_95_low)
    fill([v_95_low v_95_high v_95_high v_95_low]*1e3 - bias, ...
         [0 0 1 1], [0.3 0.9 0.3], 'FaceAlpha', 0.1, 'EdgeColor', 'g', ...
         'LineWidth', 3, 'HandleVisibility', 'off');
end
set(gca, 'FontSize', 16, 'LineWidth', 2, ...
    'Color', 'k', 'XColor', 'w', 'YColor', 'w', ...
    'GridColor', [0.4 0.4 0.4], 'GridAlpha', 0.7);
xlabel('Approach Speed [mm/s]', 'Color', 'w');
ylabel('Probability', 'Color', 'w');
title(sprintf('Monte Carlo Capture Probability (%d trials/point)', N_mc), 'Color', 'w');
lgd = legend('Location','best');
lgd.TextColor = 'w';
lgd.Color = [0.15 0.15 0.15];
lgd.EdgeColor = 'w';


%xticks_current = get(gca, 'XTick');
%set(gca, 'XTickLabel', xticks_current - 20);



%% ====================== PARAMETRIC SENSITIVITY =========================
% How does the speed window change with mechanism stiffness and damping?

k_sweep = linspace(500, 10000, 30);
c_sweep = linspace(20, 500, 30);
[KK, CC] = meshgrid(k_sweep, c_sweep);

v_min_map = zeros(size(KK));
v_max_map = zeros(size(KK));
window_map = zeros(size(KK));

fprintf('\nRunning parametric sensitivity (k × c sweep)...\n');
for ik = 1:length(k_sweep)
    for ic = 1:length(c_sweep)
        mech_test = mech;
        mech_test.k = k_sweep(ik);
        mech_test.c = c_sweep(ic);
        
        m_r = m_red;
        
        % Quick binary search for v_min
        vlo = 0.001; vhi = 0.15;
        for iter = 1:20
            vm = (vlo + vhi) / 2;
            r = simulate_contact(vm, m_r, mech_test, dt_contact, t_max_contact);
            if r.captured
                vhi = vm;
            else
                vlo = vm;
            end
        end
        v_min_map(ic, ik) = vhi;
        
        % Quick binary search for v_max (structural)
        vlo = vhi; vhi = 0.20;
        found_max = false;
        for iter = 1:20
            vm = (vlo + vhi) / 2;
            r = simulate_contact(vm, m_r, mech_test, dt_contact, t_max_contact);
            if r.F_peak > struct_lim.F_max
                vhi = vm;
                found_max = true;
            else
                vlo = vm;
            end
        end
        if found_max
            v_max_map(ic, ik) = vhi;
        else
            v_max_map(ic, ik) = 0.20;
        end
        
        window_map(ic, ik) = v_max_map(ic, ik) - v_min_map(ic, ik);
        window_map(ic, ik) = max(window_map(ic, ik), 0);
    end
end

figure('Name','Parametric Sensitivity','Position',[100 100 1000 500]);

subplot(1,2,1);
contourf(KK, CC, v_min_map * 1e3, 20);
colorbar;
hold on;
plot(mech.k, mech.c, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
xlabel('Stiffness k [N/m]');
ylabel('Damping c [N·s/m]');
title('Minimum Capture Speed [mm/s]');

subplot(1,2,2);
contourf(KK, CC, window_map * 1e3, 20);
colorbar;
hold on;
plot(mech.k, mech.c, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
xlabel('Stiffness k [N/m]');
ylabel('Damping c [N·s/m]');
title('Speed Window Width [mm/s]');

sgtitle('Sensitivity to Mechanism Properties');

%% ====================== ENERGY PARTITION ANALYSIS ======================
% Where does the kinetic energy go during contact?

fprintf('\n--- ENERGY PARTITION AT NOMINAL SPEED (%.1f mm/s) ---\n', v_nominal*1e3);

res_nom = simulate_contact(v_nominal, m_red, mech, dt_contact, t_max_contact);
KE_initial = 0.5 * m_red * v_nominal^2;

% Energy stored in spring at max compression
PE_spring = 0.5 * mech.k * res_nom.x_peak^2;

% Energy dissipated by damper (integrate c*v^2 dt)
E_damped = 0;
for i = 1:length(res_nom.v_hist)-1
    if res_nom.x_hist(i) > 0
        E_damped = E_damped + mech.c * res_nom.v_hist(i)^2 * dt_contact;
    end
end

fprintf('  Initial KE:           %.4f J  (½ × %.2f kg × (%.4f m/s)²)\n', ...
        KE_initial, m_red, v_nominal);
fprintf('  Peak spring PE:       %.4f J  (%.1f%% of KE)\n', ...
        PE_spring, PE_spring/KE_initial*100);
fprintf('  Damper dissipation:   %.4f J  (%.1f%% of KE)\n', ...
        E_damped, E_damped/KE_initial*100);
fprintf('  Rebound KE:           %.4f J  (%.1f%% of KE)\n', ...
        KE_initial - PE_spring - E_damped, ...
        (KE_initial - PE_spring - E_damped)/KE_initial*100);
fprintf('  Peak contact force:   %.1f N  (limit: %.0f N)\n', ...
        res_nom.F_peak, struct_lim.F_max);
fprintf('  Capture latency:      %.1f ms\n', res_nom.t_capture*1e3);

%% ====================== FINAL SUMMARY ==================================
fprintf('\n================================================================\n');
fprintf('  DOCKING APPROACH SPEED RECOMMENDATION\n');
fprintf('================================================================\n');
fprintf('  Vehicle pair:   %.0f kg (chaser) -> %.0f kg (target)\n', m_chaser, m_target);
fprintf('  Reduced mass:   %.1f kg\n', m_red);
fprintf('  Mechanism:      %s (k=%.0f N/m, c=%.0f N·s/m)\n', ...
        mech.type, mech.k, mech.c);
fprintf('  ----------------------------------------------------------------\n');
fprintf('  DETERMINISTIC WINDOW:\n');
fprintf('    Min speed (capture):     %6.1f mm/s\n', v_min_capture*1e3);
fprintf('    Max speed (structural):  %6.1f mm/s\n', v_max_safe*1e3);
fprintf('    Nominal:                 %6.1f mm/s\n', v_nominal*1e3);
fprintf('  ----------------------------------------------------------------\n');
fprintf('  MONTE CARLO (with GNC dispersions):\n');
fprintf('    Optimal speed:           %6.1f mm/s  (P=%.1f%%)\n', ...
        v_optimal*1e3, P_best*100);
if ~isnan(v_95_low)
    fprintf('    95%% success band:        %6.1f — %.1f mm/s\n', ...
            v_95_low*1e3, v_95_high*1e3);
end
fprintf('  ----------------------------------------------------------------\n');
fprintf('  ANALYTICAL ESTIMATES:\n');
fprintf('    v_min (energy):          %6.1f mm/s\n', v_min_energy*1e3);
fprintf('    v_max (F_max/impedance): %6.1f mm/s\n', v_max_damped*1e3);
fprintf('  ----------------------------------------------------------------\n');
fprintf('  RECOMMENDED APPROACH SPEED:  ~%.0f mm/s  (%.3f m/s)\n', ...
        v_optimal*1e3, v_optimal);
fprintf('  with GNC accuracy sufficient for ±%.0f mm lateral, ±%.1f deg att\n', ...
        disp_1sig.lat_offset*1e3, disp_1sig.att_error);
fprintf('================================================================\n');
fprintf('Done.\n');