%% ========================================================================
%  REACTION WHEEL SIZING ANALYSIS
%  Disturbance Torque Integration for Momentum Storage Sizing
%  ========================================================================
%  Computes disturbance torques on a satellite modeled as two solid cubes
%  connected by a thin hollow cylinder, with two flat-plate solar panels
%  extending from opposite faces (±Y) of one cube.
%
%  APPROACH: The satellite is assumed to maintain a fixed attitude relative
%  to LVLH (nadir-pointing) at all times — i.e., reaction wheels perfectly
%  absorb all disturbance torques. This gives the true momentum storage
%  requirement for reaction wheel sizing between desaturation events.
%
%  The code computes:
%    - Disturbance torques in body frame at each timestep
%    - Cumulative angular impulse (integral of torque) = RW momentum
%    - Peak-to-peak momentum swing per orbit (cyclic storage)
%    - Secular momentum accumulation rate (desaturation budget)
%    - Per-source breakdown for design trades
%
%  Body frame:  X = along the axis connecting the two cubes
%               Y = lateral (panel span direction)
%               Z = completes right-hand set
%
%  Layout along X-axis:
%    Cube A: [0, a]
%    Cylinder: [a, a + L_cyl]
%    Cube B: [a + L_cyl, a + L_cyl + b]
%    Solar panels extend from Cube B in ±Y directions
%
%  Perturbations modeled:
%    1. Gravity-gradient torque
%    2. Solar radiation pressure (SRP) torque
%    3. Aerodynamic drag torque (free-molecular flow)
%    4. Residual magnetic dipole torque
%
%  Author : OSCAR@VT - Yaman Saran
%  ========================================================================
clear; clc; close all;

%% =====================  USER-DEFINED PARAMETERS  ========================

% ----- Cube A (no panels) -----
cubeA.side = 1.10;          % [m] side length
cubeA.mass = 447;           % [kg]

% ----- Cube B (panels attach here) -----
cubeB.side = 1.10;          % [m] side length
cubeB.mass = 580;           % [kg]

% ----- Connecting cylinder (thin hollow tube) -----
tube.length    = 0.08;      % [m] length along X-axis
tube.R_outer   = 0.5;       % [m] outer radius
tube.R_inner   = 0.475;     % [m] inner radius
tube.mass      = 7;        % [kg]

% ----- Solar panels (two identical, extending ±Y from Cube B) -----
panel.span     = 4.00;      % [m] extent in ±Y away from cube face
panel.length   = 1.00;      % [m] extent along X (chord on panel)
panel.mass     = 10;        % [kg] mass of each panel
panel.angle    = 0;         % [deg] dihedral: tilt about X-axis (0 = flat)
panel.angle_y  = 90;         % [deg] cant: tilt about Y-axis (0 = flat)
panel.x_center = 2.28;      % [m] absolute X-coordinate of panel hinge center
                             %     leave [] to default to center of Cube B

% ----- Surface optical properties -----
cube_rho_s   = 0.10;   cube_rho_d   = 0.10;   cube_alpha = 0.80; %#ok<NASGU>
panel_rho_s  = 0.05;   panel_rho_d  = 0.05;   panel_alpha = 0.90; %#ok<NASGU>
tube_rho_s   = 0.10;   tube_rho_d   = 0.10;   tube_alpha = 0.80; %#ok<NASGU>

% ----- Orbit parameters -----
orbit.altitude = 850;        % [km] circular orbit altitude
orbit.incl     = 98.8;       % [deg] SSO

% ----- Desired body attitude relative to LVLH -----
%  Euler angles (3-2-1) from LVLH to body frame.
%  [0; 0; 0] = nadir-pointing, X-body along velocity, Z-body toward Earth
att.roll  = 0;     % [deg] rotation about X_LVLH
att.pitch = 0;     % [deg] rotation about Y_LVLH
att.yaw   = 0;     % [deg] rotation about Z_LVLH

% ----- Magnetic residual dipole (A·m²) in body frame -----
mag.dipole_body = [0.1; 0.05; 0.02];   % [A·m²]

% ----- Simulation settings -----
sim.n_orbits  = 5;          % number of orbits to simulate
sim.pts_per_orbit = 10000;     % time points per orbit

% ----- Desaturation cycle -----
desat.period_orbits = 99;     % desaturation assumed every N orbits

% ----- Environment constants -----
const.mu    = 3.986004418e14;  % [m³/s²]
const.Re    = 6371e3;          % [m]
const.P_sun = 4.56e-6;        % [N/m²] SRP at 1 AU
const.c     = 2.998e8;        % [m/s]
const.B0    = 7.94e15;        % [T·m³] Earth dipole moment

%% ====================  DERIVED QUANTITIES  ==============================
a  = cubeA.side;
b  = cubeB.side;
Lc = tube.length;
Ro = tube.R_outer;
Ri = tube.R_inner;

% X-coordinates of each component
x_cubeA_start = 0;
x_cubeA_end   = a;              %#ok<NASGU>
x_tube_start  = a;              %#ok<NASGU>
x_tube_end    = a + Lc;         %#ok<NASGU>
x_cubeB_start = a + Lc;
x_cubeB_end   = a + Lc + b;    %#ok<NASGU>

% Centroids in body frame
r_cubeA = [a/2;           0; 0];
r_tube  = [a + Lc/2;      0; 0];
r_cubeB = [a + Lc + b/2;  0; 0];

% Panel position and orientation
if isempty(panel.x_center)
    panel_x_center = x_cubeB_start + b/2;
else
    panel_x_center = panel.x_center;
end

alpha_p = deg2rad(panel.angle);
beta_p  = deg2rad(panel.angle_y);
Sp = panel.span;

% Rotation matrices for panel dihedral + cant
Rx_a  = [1, 0, 0; 0, cos(alpha_p), -sin(alpha_p); 0, sin(alpha_p), cos(alpha_p)];
Rx_na = [1, 0, 0; 0, cos(alpha_p),  sin(alpha_p); 0,-sin(alpha_p), cos(alpha_p)];
Ry_b  = [cos(beta_p), 0, sin(beta_p); 0, 1, 0; -sin(beta_p), 0, cos(beta_p)];

R_panelP = Ry_b * Rx_a;
R_panelN = Ry_b * Rx_na;

span_dirP = R_panelP * [0; 1; 0];
span_dirN = R_panelN * [0;-1; 0];

r_panelP = [panel_x_center; b/2; 0] + (Sp/2) * span_dirP;
r_panelN = [panel_x_center;-b/2; 0] + (Sp/2) * span_dirN;

% Orbit
orbit.R      = const.Re + orbit.altitude * 1e3;
orbit.n      = sqrt(const.mu / orbit.R^3);
orbit.V      = sqrt(const.mu / orbit.R);
orbit.period = 2*pi / orbit.n;

[rho_atm, ~] = atm_density_exp(orbit.altitude);

% Fixed body-to-LVLH DCM from desired Euler angles
C_bl_fixed = euler321_to_dcm(deg2rad(att.roll), deg2rad(att.pitch), deg2rad(att.yaw));

%% ==============  MASS PROPERTIES  =======================================

I_cubeA = (cubeA.mass * a^2 / 6) * eye(3);
I_cubeB = (cubeB.mass * b^2 / 6) * eye(3);

I_tube_xx = tube.mass * (Ro^2 + Ri^2) / 2;
I_tube_yy = tube.mass * (3*(Ro^2 + Ri^2) + Lc^2) / 12;
I_tube = diag([I_tube_xx, I_tube_yy, I_tube_yy]);

Lp = panel.length;
I_panel_local = panel.mass / 12 * diag([Sp^2, Lp^2, Lp^2 + Sp^2]);
I_panelP_body = R_panelP * I_panel_local * R_panelP';
I_panelN_body = R_panelN * I_panel_local * R_panelN';

M_total = cubeA.mass + tube.mass + cubeB.mass + 2*panel.mass;
r_cm = (cubeA.mass * r_cubeA + tube.mass * r_tube + cubeB.mass * r_cubeB ...
        + panel.mass * r_panelP + panel.mass * r_panelN) / M_total;

bodies = {I_cubeA,       cubeA.mass, r_cubeA;
          I_tube,        tube.mass,  r_tube;
          I_cubeB,       cubeB.mass, r_cubeB;
          I_panelP_body, panel.mass, r_panelP;
          I_panelN_body, panel.mass, r_panelN};

I_total = zeros(3);
for k = 1:size(bodies,1)
    I_k = bodies{k,1};
    m_k = bodies{k,2};
    d_k = bodies{k,3} - r_cm;
    I_total = I_total + I_k + m_k * (dot(d_k,d_k)*eye(3) - d_k*d_k');
end

fprintf('====== COMPOSITE MASS PROPERTIES ======\n');
fprintf('Total mass        : %.2f kg\n', M_total);
fprintf('CM location (body): [%.4f, %.4f, %.4f] m\n', r_cm);
fprintf('Inertia tensor (kg·m²):\n');
disp(I_total);

%% ==============  SURFACE MODEL  =========================================
surfaces = build_surface_model_dual(cubeA, cubeB, tube, panel, ...
    a, b, Lc, Ro, r_cm, alpha_p, beta_p, Sp, ...
    x_cubeA_start, x_cubeB_start, panel_x_center, ...
    cube_rho_s, cube_rho_d, panel_rho_s, panel_rho_d, tube_rho_s, tube_rho_d);
surfaces.I_total = I_total;

%% ==============  3D GEOMETRY VISUALIZATION  =============================
plot_satellite_geometry_dual(a, b, Lc, Ro, Ri, panel, alpha_p, beta_p, Sp, r_cm, I_total, ...
    x_cubeA_start, x_cubeB_start, panel_x_center);

%% ==============  DISTURBANCE TORQUE COMPUTATION (FIXED ATTITUDE)  =======
%  The satellite holds a fixed orientation relative to LVLH. At each
%  timestep we compute the body-to-inertial DCM from the orbit geometry
%  and the fixed body-to-LVLH rotation, then evaluate all disturbance
%  torques in the body frame. No attitude dynamics are propagated.

T_total_sim = orbit.period * sim.n_orbits;
N = sim.n_orbits * sim.pts_per_orbit;
t = linspace(0, T_total_sim, N)';
dt = t(2) - t(1);

fprintf('\nComputing disturbance torques for %.1f orbits (%.1f hours)...\n', ...
    sim.n_orbits, T_total_sim/3600);
fprintf('Attitude held at LVLH + [roll=%.1f°, pitch=%.1f°, yaw=%.1f°]\n', ...
    att.roll, att.pitch, att.yaw);

T_gg   = zeros(N,3);
T_srp  = zeros(N,3);
T_aero = zeros(N,3);
T_mag  = zeros(N,3);
T_tot  = zeros(N,3);

for k = 1:N
    theta = orbit.n * t(k);   % true anomaly / orbit angle

    % --- LVLH-to-Inertial DCM ---
    %  LVLH: X = velocity, Y = -orbit-normal (for inclined), Z = nadir
    %  For a circular orbit in the equatorial plane simplification:
    %    r_hat_i = [cos(theta); sin(theta); 0]
    %    v_hat_i = [-sin(theta); cos(theta); 0]
    %    h_hat_i = [0; 0; 1]  (for i=0)
    %
    %  For inclined orbit, we rotate about the inertial Z by RAAN and
    %  then account for inclination. For SSO we use a simplified model
    %  where the orbit plane is fixed during the short simulation.

    % Inertial directions
    r_hat_i = [cos(theta); sin(theta); 0];
    v_hat_i = [-sin(theta); cos(theta); 0];

    % LVLH frame (X=velocity, Y=orbit-normal×velocity, Z=nadir)
    %  Here orbit normal = +Z_inertial for i=0 simplification
    %  (The inclination affects the magnetic field model but the LVLH
    %   construction is consistent for any equatorial-plane orbit.)
    z_lvlh_i = -r_hat_i;       % nadir
    x_lvlh_i = v_hat_i;        % velocity direction
    y_lvlh_i = cross(z_lvlh_i, x_lvlh_i);
    y_lvlh_i = y_lvlh_i / norm(y_lvlh_i);

    C_li = [x_lvlh_i, y_lvlh_i, z_lvlh_i]';   % LVLH-to-Inertial = rows are LVLH axes in inertial

    % Body-to-Inertial: C_bi = C_bl * C_li
    C_bi = C_bl_fixed * C_li;

    % --- Compute torques using this fixed-attitude DCM ---
    r_hat_b = C_bi * r_hat_i;
    sun_hat_i = [1; 0; 0];      % simplified fixed sun direction
    sun_hat_b = C_bi * sun_hat_i;
    v_hat_b = C_bi * v_hat_i;

    % 1. Gravity-Gradient
    Tgg = 3 * orbit.n^2 * cross(r_hat_b, I_total * r_hat_b);

    % 2 & 3. SRP and Aero
    Tsrp  = [0;0;0];
    Taero = [0;0;0];
    for s = 1:length(surfaces.list)
        surf = surfaces.list(s);
        n_hat = surf.normal;
        A     = surf.area;
        r_cp  = surf.r_cp;
        rho_s = surf.rho_s;
        rho_d = surf.rho_d;

        cos_sun = dot(sun_hat_b, n_hat);
        if cos_sun > 0
            f_srp = -const.P_sun * A * ( ...
                     (1 - rho_s)*sun_hat_b + ...
                     2*(rho_s*cos_sun + rho_d/3)*n_hat ) * cos_sun;
            Tsrp = Tsrp + cross(r_cp, f_srp);
        end

        cos_v = dot(v_hat_b, n_hat);
        if cos_v > 0
            Cd = 2.2;
            f_aero = -0.5 * rho_atm * orbit.V^2 * Cd * A * cos_v * v_hat_b;
            Taero = Taero + cross(r_cp, f_aero);
        end
    end

    % 4. Magnetic torque
    B_i = earth_mag_field(r_hat_i, orbit.R, const);
    B_b = C_bi * B_i;
    Tmag = cross(mag.dipole_body, B_b);

    % Store
    T_gg(k,:)   = Tgg';
    T_srp(k,:)  = Tsrp';
    T_aero(k,:) = Taero';
    T_mag(k,:)  = Tmag';
    T_tot(k,:)  = (Tgg + Tsrp + Taero + Tmag)';
end

%% ==============  CUMULATIVE ANGULAR IMPULSE  ============================
%  H(t) = integral_0^t T(tau) dtau   [N·m·s]
%  This is the momentum that reaction wheels must store.

H_gg   = cumtrapz(t, T_gg);
H_srp  = cumtrapz(t, T_srp);
H_aero = cumtrapz(t, T_aero);
H_mag  = cumtrapz(t, T_mag);
H_tot  = cumtrapz(t, T_tot);

%% ==============  PER-ORBIT ANALYSIS  ====================================
%  Break the simulation into individual orbits and compute:
%    - Per-orbit peak-to-peak momentum swing (cyclic component)
%    - Per-orbit net secular momentum accumulation
%    - Worst-case single-orbit storage requirement

orbit_idx = zeros(sim.n_orbits, 2);   % start/end indices per orbit
H_pp_orbit   = zeros(sim.n_orbits, 3);  % peak-to-peak per orbit
H_net_orbit  = zeros(sim.n_orbits, 3);  % net accumulation per orbit
T_rms_orbit  = zeros(sim.n_orbits, 3);  % RMS torque per orbit
T_peak_orbit = zeros(sim.n_orbits, 3);  % peak absolute torque per orbit (per axis)

for orb = 1:sim.n_orbits
    t_start = (orb - 1) * orbit.period;
    t_end   = orb * orbit.period;
    idx = find(t >= t_start & t <= t_end);
    orbit_idx(orb,:) = [idx(1), idx(end)];

    % Momentum within this orbit (re-zeroed at orbit start)
    H_orb = cumtrapz(t(idx) - t(idx(1)), T_tot(idx,:));

    H_pp_orbit(orb,:)  = max(H_orb) - min(H_orb);
    H_net_orbit(orb,:) = H_orb(end,:);
    T_rms_orbit(orb,:) = sqrt(mean(T_tot(idx,:).^2, 1));
    T_peak_orbit(orb,:) = max(abs(T_tot(idx,:)), [], 1);
end

% Worst-case orbit for sizing
[~, worst_orb_x] = max(H_pp_orbit(:,1));
[~, worst_orb_y] = max(H_pp_orbit(:,2));
[~, worst_orb_z] = max(H_pp_orbit(:,3));

% Desaturation interval analysis
desat_period_s = desat.period_orbits * orbit.period;
n_desat_intervals = floor(T_total_sim / desat_period_s);
H_desat_max = zeros(1,3);

for di = 1:n_desat_intervals
    t_start = (di - 1) * desat_period_s;
    t_end_d = di * desat_period_s;
    idx = find(t >= t_start & t <= t_end_d);
    H_interval = cumtrapz(t(idx) - t(idx(1)), T_tot(idx,:));
    pp = max(H_interval) - min(H_interval);
    H_desat_max = max(H_desat_max, pp);
end



%% ==============  PRINT RESULTS  =========================================
fprintf('\n====================================================================\n');
fprintf('         REACTION WHEEL SIZING SUMMARY\n');
fprintf('====================================================================\n');
fprintf('Orbit:  %.0f km altitude, %.1f° inclination (period = %.1f min)\n', ...
    orbit.altitude, orbit.incl, orbit.period/60);
fprintf('Attitude: LVLH-fixed [roll=%.1f°, pitch=%.1f°, yaw=%.1f°]\n', ...
    att.roll, att.pitch, att.yaw);
fprintf('Simulation: %d orbits (%.1f hours)\n', sim.n_orbits, T_total_sim/3600);
fprintf('Desaturation interval: %d orbit(s)\n', desat.period_orbits);
fprintf('\n--- Average Torques (over full simulation) ---\n');
T_mean = mean(T_tot, 1);
fprintf('  Mean T_x = %+.4e  N·m\n', T_mean(1));
fprintf('  Mean T_y = %+.4e  N·m\n', T_mean(2));
fprintf('  Mean T_z = %+.4e  N·m\n', T_mean(3));

fprintf('\n--- RMS Torques (over full simulation) ---\n');
T_rms = sqrt(mean(T_tot.^2, 1));
fprintf('  RMS  T_x = %.4e  N·m\n', T_rms(1));
fprintf('  RMS  T_y = %.4e  N·m\n', T_rms(2));
fprintf('  RMS  T_z = %.4e  N·m\n', T_rms(3));

fprintf('\n--- Peak Torque per Orbit (worst-case instantaneous) ---\n');
T_peak_worst = max(T_peak_orbit, [], 1);
[~, peak_orb_x] = max(T_peak_orbit(:,1));
[~, peak_orb_y] = max(T_peak_orbit(:,2));
[~, peak_orb_z] = max(T_peak_orbit(:,3));
fprintf('  Peak T_x = %.4e  N·m  (orbit %d)\n', T_peak_worst(1), peak_orb_x);
fprintf('  Peak T_y = %.4e  N·m  (orbit %d)\n', T_peak_worst(2), peak_orb_y);
fprintf('  Peak T_z = %.4e  N·m  (orbit %d)\n', T_peak_worst(3), peak_orb_z);
fprintf('  Peak |T|  = %.4e  N·m  (max instantaneous magnitude)\n', ...
    max(vecnorm(T_tot, 2, 2)));

fprintf('\n--- Per-Source Mean Torque Magnitudes ---\n');
fprintf('  Gravity-Gradient: |T| = %.4e  N·m (mean)\n', mean(vecnorm(T_gg,2,2)));
fprintf('  SRP:              |T| = %.4e  N·m (mean)\n', mean(vecnorm(T_srp,2,2)));
fprintf('  Aerodynamic:      |T| = %.4e  N·m (mean)\n', mean(vecnorm(T_aero,2,2)));
fprintf('  Magnetic dipole:  |T| = %.4e  N·m (mean)\n', mean(vecnorm(T_mag,2,2)));

fprintf('\n--- Worst-Case Per-Orbit Peak-to-Peak Momentum Storage ---\n');
H_pp_worst = max(H_pp_orbit, [], 1);
fprintf('  H_x = %.4e  N·m·s  (orbit %d)\n', H_pp_worst(1), worst_orb_x);
fprintf('  H_y = %.4e  N·m·s  (orbit %d)\n', H_pp_worst(2), worst_orb_y);
fprintf('  H_z = %.4e  N·m·s  (orbit %d)\n', H_pp_worst(3), worst_orb_z);
fprintf('  |H|  = %.4e  N·m·s  (RSS of worst-case axes)\n', norm(H_pp_worst));

fprintf('\n--- Momentum Storage per Desaturation Interval (%d orbit(s)) ---\n', ...
    desat.period_orbits);
fprintf('  H_x = %.4e  N·m·s\n', H_desat_max(1));
fprintf('  H_y = %.4e  N·m·s\n', H_desat_max(2));
fprintf('  H_z = %.4e  N·m·s\n', H_desat_max(3));
fprintf('  |H|  = %.4e  N·m·s  (RSS)\n', norm(H_desat_max));

fprintf('\n--- Secular Momentum Accumulation Rate ---\n');
H_secular_rate = H_tot(end,:) / T_total_sim;
fprintf('  dH/dt_x = %+.4e  N·m  (avg secular torque)\n', H_secular_rate(1));
fprintf('  dH/dt_y = %+.4e  N·m\n', H_secular_rate(2));
fprintf('  dH/dt_z = %+.4e  N·m\n', H_secular_rate(3));
fprintf('  Momentum per orbit (secular): [%.4e, %.4e, %.4e] N·m·s\n', ...
    H_secular_rate * orbit.period);

fprintf('\n--- RECOMMENDED MINIMUM RW CAPACITY (with 50%% margin) ---\n');
H_required = H_desat_max * 1.5;
fprintf('  Momentum storage:\n');
fprintf('    H_x >= %.4f  N·m·s\n', H_required(1));
fprintf('    H_y >= %.4f  N·m·s\n', H_required(2));
fprintf('    H_z >= %.4f  N·m·s\n', H_required(3));
fprintf('    |H| >= %.4f  N·m·s  (per-axis RSS with margin)\n', norm(H_required));
T_required = T_peak_worst * 1.5;
fprintf('  Torque authority:\n');
fprintf('    T_x >= %.4e  N·m\n', T_required(1));
fprintf('    T_y >= %.4e  N·m\n', T_required(2));
fprintf('    T_z >= %.4e  N·m\n', T_required(3));
fprintf('    |T| >= %.4e  N·m  (peak RSS with margin)\n', ...
    max(vecnorm(T_tot, 2, 2)) * 1.5);

% Secular momentum per desaturation interval
H_secular_per_interval = H_secular_rate * desat.period_orbits * orbit.period;  % [N·m·s]

% Time allocated for desaturation maneuver
t_desat = 0.25 * orbit.period;  % e.g., quarter orbit

% Required magnetorquer torque
tau_mag_req = abs(H_secular_per_interval') / t_desat;  % [N·m] per axis
tau_mag_req(2) = max(H_pp_orbit(:,2)) / (0.25 * orbit.period); % Y axis peak as proxy

fprintf("tau_mag_req = %.5f", tau_mag_req);

%% ==============  PLOTS  =================================================
t_orb = t / orbit.period;   % time in orbit units
t_hr  = t / 3600;

% =========================================================================
% FIGURE 1: Disturbance Torque Profiles (2 orbits shown)
% =========================================================================
n_show = min(2, sim.n_orbits);
mask = t_orb <= n_show;

figure('Name','Disturbance Torques','Position',[100 100 1200 800]);

subplot(2,2,1);
plot(t_orb(mask), T_gg(mask,:)*1e6, 'LineWidth', 1.2);
xlabel('Time [orbits]'); ylabel('Torque [\muN\cdotm]');
title('Gravity-Gradient'); legend('T_x','T_y','T_z'); grid on;

subplot(2,2,2);
plot(t_orb(mask), T_srp(mask,:)*1e6, 'LineWidth', 1.2);
xlabel('Time [orbits]'); ylabel('Torque [\muN\cdotm]');
title('Solar Radiation Pressure'); legend('T_x','T_y','T_z'); grid on;

subplot(2,2,3);
plot(t_orb(mask), T_aero(mask,:)*1e6, 'LineWidth', 1.2);
xlabel('Time [orbits]'); ylabel('Torque [\muN\cdotm]');
title('Aerodynamic Drag'); legend('T_x','T_y','T_z'); grid on;

subplot(2,2,4);
plot(t_orb(mask), T_mag(mask,:)*1e6, 'LineWidth', 1.2);
xlabel('Time [orbits]'); ylabel('Torque [\muN\cdotm]');
title('Residual Magnetic Dipole'); legend('T_x','T_y','T_z'); grid on;

sgtitle(sprintf('Disturbance Torques at Fixed LVLH Attitude  [%.0f km, %.1f° SSO]', ...
    orbit.altitude, orbit.incl));

% =========================================================================
% FIGURE 2: Total Torque Magnitude + Per-Component
% =========================================================================
figure('Name','Total Torque','Position',[100 100 1000 500]);
tiledlayout(1,2);

nexttile;
plot(t_orb(mask), T_tot(mask,:)*1e6, 'LineWidth', 1.2);
xlabel('Time [orbits]'); ylabel('Torque [\muN\cdotm]');
title('Total Disturbance Torque (components)');
legend('T_x','T_y','T_z'); grid on;

nexttile;
T_mag_sources = [vecnorm(T_gg(mask,:),2,2), vecnorm(T_srp(mask,:),2,2), ...
                 vecnorm(T_aero(mask,:),2,2), vecnorm(T_mag(mask,:),2,2)] * 1e6;
area(t_orb(mask), T_mag_sources, 'LineWidth', 0.5);
xlabel('Time [orbits]'); ylabel('|T| [\muN\cdotm]');
title('Torque Magnitude by Source (stacked)');
legend('Gravity-Gradient','SRP','Aerodynamic','Magnetic'); grid on;
alpha(0.7);

sgtitle('Total Disturbance Torque Breakdown');

% =========================================================================
% FIGURE 3: Cumulative Angular Impulse (Momentum Envelope)
% =========================================================================
figure('Name','Momentum Envelope','Position',[100 100 1200 700]);
tiledlayout(2,2);

nexttile([1,2]);
plot(t_orb, H_tot, 'LineWidth', 1.5);
xlabel('Time [orbits]'); ylabel('H [N\cdotm\cdots]');
title('Cumulative Angular Impulse — Total (RW Stored Momentum if No Desaturation)');
legend('H_x','H_y','H_z','Location','best'); grid on;

nexttile;
plot(t_orb, H_gg, '-', 'LineWidth', 1.0);
hold on;
plot(t_orb, H_srp, '--', 'LineWidth', 1.0);
xlabel('Time [orbits]'); ylabel('H [N\cdotm\cdots]');
title('GG (solid) + SRP (dashed)');
legend('GG_x','GG_y','GG_z','SRP_x','SRP_y','SRP_z','Location','best'); grid on;

nexttile;
plot(t_orb, H_aero, '-', 'LineWidth', 1.0);
hold on;
plot(t_orb, H_mag, '--', 'LineWidth', 1.0);
xlabel('Time [orbits]'); ylabel('H [N\cdotm\cdots]');
title('Aero (solid) + Mag (dashed)');
legend('Aero_x','Aero_y','Aero_z','Mag_x','Mag_y','Mag_z','Location','best'); grid on;

sgtitle('Momentum Accumulation by Disturbance Source');

% =========================================================================
% FIGURE 4: Per-Orbit Momentum Budget (Bar Chart)
% =========================================================================
figure('Name','Per-Orbit Momentum & Torque','Position',[100 100 1000 850]);
tiledlayout(3,1);

nexttile;
bar(1:sim.n_orbits, T_peak_orbit * 1e6, 'grouped');
xlabel('Orbit Number'); ylabel('Peak |T| [\muN\cdotm]');
title('Peak Instantaneous Torque per Orbit (RW Torque Authority Required)');
legend('T_x','T_y','T_z','Location','best'); grid on;

nexttile;
bar(1:sim.n_orbits, H_pp_orbit, 'grouped');
xlabel('Orbit Number'); ylabel('Peak-to-Peak H [N\cdotm\cdots]');
title('Cyclic Momentum Swing per Orbit (RW Storage Required)');
legend('H_x','H_y','H_z','Location','best'); grid on;

nexttile;
bar(1:sim.n_orbits, H_net_orbit, 'grouped');
xlabel('Orbit Number'); ylabel('Net \DeltaH [N\cdotm\cdots]');
title('Net Secular Momentum Accumulation per Orbit (Desaturation Budget)');
legend('\DeltaH_x','\DeltaH_y','\DeltaH_z','Location','best'); grid on;

sgtitle('Per-Orbit Reaction Wheel Analysis');

% =========================================================================
% FIGURE 5: Single Worst-Case Orbit Detail
% =========================================================================
% Show the orbit with largest total peak-to-peak momentum
[~, worst_orb] = max(vecnorm(H_pp_orbit, 2, 2));
idx_w = orbit_idx(worst_orb,1):orbit_idx(worst_orb,2);
t_orb_local = (t(idx_w) - t(idx_w(1))) / orbit.period;

H_worst = cumtrapz(t(idx_w) - t(idx_w(1)), T_tot(idx_w,:));

figure('Name','Worst-Case Orbit','Position',[100 100 1000 600]);
tiledlayout(2,1);

nexttile;
plot(t_orb_local, T_tot(idx_w,:)*1e6, 'LineWidth', 1.3);
hold on;
colors = get(gca, 'ColorOrder');
ax_labels = {'x','y','z'};
for ax = 1:3
    pk = T_peak_orbit(worst_orb, ax) * 1e6;
    yline(pk, '--', 'Color', colors(ax,:), 'LineWidth', 1);
    yline(-pk, '--', 'Color', colors(ax,:), 'LineWidth', 1);
end
xlim([t_orb_local(1), t_orb_local(round(end/2))]);
xlabel('Orbit Fraction'); ylabel('Torque [\muN\cdotm]');
title(sprintf('Disturbance Torque — Orbit %d (Worst-Case)  |  Peak lines shown dashed', worst_orb));
legend('T_x','T_y','T_z'); grid on;

nexttile;
plot(t_orb_local, H_worst, 'LineWidth', 1.5);
hold on;
for ax = 1:3
    yline(max(H_worst(:,ax)), '--', 'Color', get(gca,'ColorOrder')); %#ok<NOPTS>
    yline(min(H_worst(:,ax)), '--', 'Color', get(gca,'ColorOrder'));
end
xlabel('Orbit Fraction'); ylabel('H [N\cdotm\cdots]');
title('Accumulated Momentum in Worst-Case Orbit (RW Storage)');
legend('H_x','H_y','H_z','Location','best'); grid on;

sgtitle(sprintf('Worst-Case Single Orbit Detail — Orbit %d', worst_orb));

% =========================================================================
% FIGURE 6: Momentum Budget Breakdown (End-of-Sim Bar Chart)
% =========================================================================
figure('Name','Momentum Budget Breakdown','Position',[100 100 700 450]);
H_final_abs = [abs(H_gg(end,:)); abs(H_srp(end,:)); abs(H_aero(end,:)); abs(H_mag(end,:))];
bar(H_final_abs);
set(gca, 'XTickLabel', {'Gravity-Gradient','SRP','Aerodynamic','Magnetic'});
ylabel('|H| at End of Sim [N\cdotm\cdots]');
title(sprintf('Net Angular Impulse by Source (%d orbits)', sim.n_orbits));
legend('X','Y','Z'); grid on;

% =========================================================================
% FIGURE 7: Momentum Magnitude with Desaturation Resets
% =========================================================================
figure('Name','Desaturation Cycle View','Position',[100 100 1000 450]);

% Simulate desaturation resets: reset momentum integral at each desat event
H_desat_sim = zeros(N, 3);
last_desat_idx = 1;
desat_times = [];

for k = 1:N
    if t(k) - t(last_desat_idx) >= desat_period_s && k > 1
        last_desat_idx = k;
        desat_times = [desat_times; t(k)]; %#ok<AGROW>
    end
    if k == last_desat_idx
        H_desat_sim(k,:) = [0, 0, 0];
    else
        H_desat_sim(k,:) = H_desat_sim(k-1,:) + T_tot(k,:) * dt;
    end
end

plot(t_orb, vecnorm(H_desat_sim, 2, 2), 'LineWidth', 1.5, 'Color', [0.8 0.15 0.15]);
hold on;
for dti = 1:length(desat_times)
    xline(desat_times(dti)/orbit.period, '--', 'Color', [0.3 0.3 0.8], 'LineWidth', 1);
end
yline(norm(H_desat_max), '-.k', sprintf('Max = %.4e N·m·s', norm(H_desat_max)), ...
    'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
yline(norm(H_required), '--r', sprintf('Required (1.5×) = %.4e N·m·s', norm(H_required)), ...
    'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
xlabel('Time [orbits]'); ylabel('|H_{stored}| [N\cdotm\cdots]');
title(sprintf('RW Momentum Magnitude with Desaturation Every %d Orbit(s)', desat.period_orbits));
grid on;

fprintf('\nDone. %d figures generated.\n', 7);


% Y axis allocation bug
tau_mag_req(2) = max(H_pp_orbit(:,2)) / (0.25 * orbit.period);

%% ========================================================================
%                        LOCAL FUNCTIONS
%% ========================================================================

function C = euler321_to_dcm(phi, theta, psi)
    % 3-2-1 Euler angles (roll, pitch, yaw) to DCM
    % C rotates from the reference frame to the body frame
    cphi = cos(phi); sphi = sin(phi);
    cth  = cos(theta); sth = sin(theta);
    cpsi = cos(psi); spsi = sin(psi);

    C = [cth*cpsi,                cth*spsi,                -sth;
         sphi*sth*cpsi-cphi*spsi, sphi*sth*spsi+cphi*cpsi, sphi*cth;
         cphi*sth*cpsi+sphi*spsi, cphi*sth*spsi-sphi*cpsi, cphi*cth];
end

function surfaces = build_surface_model_dual(cubeA, cubeB, tube, panel, ...
    a, b, Lc, Ro, r_cm, alpha_p, beta_p, Sp, ...
    x_A0, x_B0, panel_xc, ...
    c_rs, c_rd, p_rs, p_rd, t_rs, t_rd)

    list = [];

    % ---- Cube A: 6 faces ----
    face_normals = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
    face_centers = [x_A0+a, 0, 0;   x_A0, 0, 0;
                    x_A0+a/2, a/2, 0;  x_A0+a/2,-a/2, 0;
                    x_A0+a/2, 0, a/2;  x_A0+a/2, 0,-a/2];
    A_face_A = a^2;
    for i = 1:6
        s.normal = face_normals(i,:)';
        s.area   = A_face_A;
        s.r_cp   = face_centers(i,:)' - r_cm;
        s.rho_s  = c_rs;
        s.rho_d  = c_rd;
        list = [list, s]; %#ok<AGROW>
    end

    % ---- Cube B: 6 faces ----
    face_centers_B = [x_B0+b, 0, 0;  x_B0, 0, 0;
                      x_B0+b/2, b/2, 0;  x_B0+b/2,-b/2, 0;
                      x_B0+b/2, 0, b/2;  x_B0+b/2, 0,-b/2];
    A_face_B = b^2;
    for i = 1:6
        s.normal = face_normals(i,:)';
        s.area   = A_face_B;
        s.r_cp   = face_centers_B(i,:)' - r_cm;
        s.rho_s  = c_rs;
        s.rho_d  = c_rd;
        list = [list, s]; %#ok<AGROW>
    end

    % ---- Connecting cylinder: discretize azimuthally ----
    n_az = 12;
    dth = 2*pi / n_az;
    x_tube_mid = a + Lc/2;
    for i = 1:n_az
        th = (i - 0.5) * dth;
        s.normal = [0; cos(th); sin(th)];
        s.area   = Ro * dth * Lc;
        s.r_cp   = [x_tube_mid; Ro*cos(th); Ro*sin(th)] - r_cm;
        s.rho_s  = t_rs;
        s.rho_d  = t_rd;
        list = [list, s]; %#ok<AGROW>
    end

    % ---- Solar panels ----
    Lp = panel.length;
    A_panel = Lp * Sp;

    ca = cos(alpha_p);  sa = sin(alpha_p);
    cb = cos(beta_p);   sb = sin(beta_p);
    Ry_b  = [cb, 0, sb; 0, 1, 0; -sb, 0, cb];
    Rx_a  = [1, 0, 0; 0, ca, -sa; 0, sa, ca];
    Rx_na = [1, 0, 0; 0, ca,  sa; 0,-sa, ca];
    R_pY = Ry_b * Rx_a;
    R_nY = Ry_b * Rx_na;

    span_pY   = R_pY * [0;1;0];
    normal_pY = R_pY * [0;0;1];
    r_pp = [panel_xc; b/2; 0] + (Sp/2) * span_pY - r_cm;

    s.normal = normal_pY;   s.area = A_panel;  s.r_cp = r_pp;
    s.rho_s = p_rs;  s.rho_d = p_rd;
    list = [list, s];
    s.normal = -normal_pY;  s.area = A_panel;  s.r_cp = r_pp;
    list = [list, s];

    span_nY   = R_nY * [0;-1;0];
    normal_nY = R_nY * [0;0;1];
    r_pn = [panel_xc;-b/2; 0] + (Sp/2) * span_nY - r_cm;

    s.normal = normal_nY;   s.area = A_panel;  s.r_cp = r_pn;
    s.rho_s = p_rs;  s.rho_d = p_rd;
    list = [list, s];
    s.normal = -normal_nY;  s.area = A_panel;  s.r_cp = r_pn;
    list = [list, s];

    surfaces.list = list;
    surfaces.I_total = [];
end

function plot_satellite_geometry_dual(a, b, Lc, Ro, Ri, panel, alpha_p, beta_p, Sp, r_cm, I_total, ...
    x_A0, x_B0, panel_xc) %#ok<INUSL>

    figure('Name','Satellite Geometry','Position',[80 80 1200 800]);
    hold on;

    col_cubeA  = [0.55 0.65 0.75];
    col_cubeB  = [0.65 0.55 0.70];
    col_tube   = [0.50 0.50 0.50];
    col_panel  = [0.15 0.20 0.55];
    col_back   = [0.75 0.75 0.70];
    col_cm     = [1.0  0.15 0.15];
    col_axes   = [0.9 0.1 0.1; 0.1 0.7 0.1; 0.1 0.2 0.9];

    draw_cube(x_A0, -a/2, -a/2, a, a, a, col_cubeA);
    draw_cube(x_B0, -b/2, -b/2, b, b, b, col_cubeB);

    n_circ = 40;
    theta_c = linspace(0, 2*pi, n_circ);
    x_ends = [a, a + Lc];
    [Theta, Xc] = meshgrid(theta_c, x_ends);
    Yc = Ro * cos(Theta);
    Zc = Ro * sin(Theta);
    surf(Xc, Yc, Zc, 'FaceColor', col_tube, 'EdgeColor', 'none', ...
        'FaceAlpha', 0.7, 'FaceLighting', 'gouraud');

    Lp = panel.length;
    ca = cos(alpha_p);  sa = sin(alpha_p);
    cb = cos(beta_p);   sb = sin(beta_p);
    Ry_b  = [cb, 0, sb; 0, 1, 0; -sb, 0, cb];
    Rx_a  = [1, 0, 0; 0, ca, -sa; 0, sa, ca];
    Rx_na = [1, 0, 0; 0, ca,  sa; 0,-sa, ca];

    for sign = [1, -1]
        hinge = [panel_xc; sign*b/2; 0];
        if sign > 0
            R_p = Ry_b * Rx_a;
        else
            R_p = Ry_b * Rx_na;
        end
        corners_local = [
            -Lp/2, 0, 0;
             Lp/2, 0, 0;
             Lp/2, sign*Sp, 0;
            -Lp/2, sign*Sp, 0;
        ]';
        corners_body = R_p * corners_local + hinge;

        fill3(corners_body(1,:), corners_body(2,:), corners_body(3,:), col_panel, ...
            'FaceAlpha', 0.9, 'EdgeColor', [0.1 0.1 0.4], 'LineWidth', 1.5);
        n_face = R_p * [0;0;1];
        offset = 0.005 * n_face;
        fill3(corners_body(1,:) - offset(1), ...
              corners_body(2,:) - offset(2), ...
              corners_body(3,:) - offset(3), col_back, ...
            'FaceAlpha', 0.7, 'EdgeColor', 'none');
    end

    plot3(r_cm(1), r_cm(2), r_cm(3), 'o', 'Color', col_cm, ...
        'MarkerSize', 10, 'MarkerFaceColor', col_cm, 'LineWidth', 2);
    text(r_cm(1), r_cm(2), r_cm(3) - 0.35, '  CM', ...
        'Color', col_cm, 'FontSize', 11, 'FontWeight', 'bold');

    total_len = a + Lc + b;
    ax_len = total_len * 0.35;
    labels = {'X_{body}', 'Y_{body}', 'Z_{body}'};
    dirs = eye(3);
    for i = 1:3
        quiver3(r_cm(1), r_cm(2), r_cm(3), ...
                dirs(1,i)*ax_len, dirs(2,i)*ax_len, dirs(3,i)*ax_len, 0, ...
                'Color', col_axes(i,:), 'LineWidth', 2.2, 'MaxHeadSize', 0.4);
        text(r_cm(1) + dirs(1,i)*ax_len*1.12, ...
             r_cm(2) + dirs(2,i)*ax_len*1.12, ...
             r_cm(3) + dirs(3,i)*ax_len*1.12, ...
             labels{i}, 'Color', col_axes(i,:), 'FontSize', 11, 'FontWeight', 'bold');
    end

    [V_p, D_p] = eig(I_total);
    I_p = diag(D_p);
    [I_p, idx] = sort(I_p);
    V_p = V_p(:, idx);
    pa_len = ax_len * 0.7;
    pa_col = [1.0 0.5 0.0; 0.6 0.0 0.8; 0.0 0.7 0.7];
    h_leg = gobjects(3,1);
    for i = 1:3
        v = V_p(:,i);
        quiver3(r_cm(1), r_cm(2), r_cm(3), ...
                v(1)*pa_len, v(2)*pa_len, v(3)*pa_len, 0, ...
                'Color', pa_col(i,:), 'LineWidth', 1.5, 'LineStyle', '--', 'MaxHeadSize', 0.3);
        h_leg(i) = plot3(NaN,NaN,NaN,'--','Color',pa_col(i,:),'LineWidth',1.5);
    end

    z_dim = -max(a,b)/2 - 0.3;
    plot3([x_A0, x_A0+a], [0,0], [z_dim, z_dim], 'k-', 'LineWidth', 1);
    text((x_A0+a)/2, 0, z_dim - 0.15, sprintf('%.2f m', a), ...
        'HorizontalAlignment','center','FontSize',8,'Color',[0.3 0.3 0.3]);
    plot3([a, a+Lc], [0,0], [z_dim, z_dim], 'k-', 'LineWidth', 1);
    text(a+Lc/2, 0, z_dim - 0.15, sprintf('%.2f m', Lc), ...
        'HorizontalAlignment','center','FontSize',8,'Color',[0.3 0.3 0.3]);
    plot3([x_B0, x_B0+b], [0,0], [z_dim, z_dim], 'k-', 'LineWidth', 1);
    text(x_B0+b/2, 0, z_dim - 0.15, sprintf('%.2f m', b), ...
        'HorizontalAlignment','center','FontSize',8,'Color',[0.3 0.3 0.3]);

    light('Position',[2*(a+Lc+b), -3*(a+Lc+b), 2*(a+Lc+b)],'Style','infinite');
    light('Position',[-(a+Lc+b), 2*(a+Lc+b), (a+Lc+b)],'Style','infinite','Color',[0.3 0.3 0.3]);
    material([0.6 0.7 0.3 10]);

    axis equal; grid on;
    xlabel('X_{body} [m]'); ylabel('Y_{body} [m]'); zlabel('Z_{body} [m]');
    title(sprintf(['Satellite Geometry  |  CubeA: %.2f m  |  Tube: %.2f m \\times \\varnothing%.2f m' ...
           '  |  CubeB: %.2f m  |  Panels: %.1f \\times %.1f m'], ...
           a, Lc, 2*Ro, b, panel.span, panel.length), 'FontSize', 11);
    view(135, 25);

    info_str = {sprintf('I_1 = %.2f  kg{\\cdot}m^2  (min)', I_p(1)), ...
                sprintf('I_2 = %.2f  kg{\\cdot}m^2  (mid)', I_p(2)), ...
                sprintf('I_3 = %.2f  kg{\\cdot}m^2  (max)', I_p(3))};
    legend(h_leg, info_str, 'Location', 'northeast', 'FontSize', 9);
    set(gcf, 'Color', 'w');
    rotate3d on;
end

function draw_cube(x0, y0, z0, dx, dy, dz, col)
    verts = [x0,    y0,    z0;
             x0+dx, y0,    z0;
             x0+dx, y0+dy, z0;
             x0,    y0+dy, z0;
             x0,    y0,    z0+dz;
             x0+dx, y0,    z0+dz;
             x0+dx, y0+dy, z0+dz;
             x0,    y0+dy, z0+dz];
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 3 4 8 7; 1 4 8 5; 2 3 7 6];
    patch('Vertices', verts, 'Faces', faces, 'FaceColor', col, ...
          'EdgeColor', [0.3 0.3 0.3], 'FaceAlpha', 0.85, ...
          'FaceLighting', 'gouraud', 'LineWidth', 0.8);
end

function [rho, H] = atm_density_exp(alt_km)
    alt_ref = [200, 300, 400, 500, 600, 700, 800, 900, 1000];
    rho_ref = [2.53e-10, 6.24e-12, 3.73e-13, 4.89e-14, ...
               9.15e-15, 2.39e-15, 8.19e-16, 3.58e-16, 1.84e-16];
    H_ref   = [37.5, 53.6, 58.5, 63.8, 71.8, 88.7, 124.6, 181.1, 268.0];
    if alt_km < alt_ref(1)
        rho = rho_ref(1); H = H_ref(1);
    elseif alt_km >= alt_ref(end)
        rho = rho_ref(end); H = H_ref(end);
    else
        idx = find(alt_ref <= alt_km, 1, 'last');
        H = H_ref(idx);
        rho = rho_ref(idx) * exp(-(alt_km - alt_ref(idx)) / H);
    end
end

function B = earth_mag_field(r_hat, R_orbit, const)
    m_hat = [0; 0; 1];
    B0_r = const.B0 / R_orbit^3;
    B = B0_r * (3 * dot(m_hat, r_hat) * r_hat - m_hat);
end