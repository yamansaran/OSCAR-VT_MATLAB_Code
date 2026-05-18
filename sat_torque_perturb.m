%% ========================================================================
%  SATELLITE ROTATIONAL TENDENCY ANALYSIS
%  Dual-Cube + Connecting Cylinder + Two Solar Panels
%  ========================================================================
%  Computes and propagates attitude disturbance torques on a satellite
%  modeled as two solid cubes connected by a thin hollow cylinder, with
%  two flat-plate solar panels extending from opposite faces (±Y) of
%  one cube.
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
%  Outputs include cumulative angular impulse (integral of torque) for
%  reaction wheel momentum storage sizing.
%
%  Author : OSCAR@VT
%  Yaman Saran
%  ========================================================================
clear; clc; close all;

%% =====================  USER-DEFINED PARAMETERS  ========================

% ----- Cube A (no panels) -----
cubeA.side = 1.10;          % [m] side length
cubeA.mass = 450;           % [kg]

% ----- Cube B (panels attach here) -----
cubeB.side = 1.10;          % [m] side length
cubeB.mass = 400;           % [kg]

% ----- Connecting cylinder (thin hollow tube) -----
tube.length    = 0.08;      % [m] length along X-axis
tube.R_outer   = 0.5;      % [m] outer radius
tube.R_inner   = 0.475;      % [m] inner radius
tube.mass      = 10;        % [kg]

% ----- Solar panels (two identical, extending ±Y from Cube B) -----
panel.span     = 4.00;      % [m] extent in ±Y away from cube face
panel.length   = 1.00;      % [m] extent along X (chord on panel)
panel.mass     = 20;        % [kg] mass of each panel
panel.angle    = 0;         % [deg] dihedral: tilt about X-axis (0 = flat, >0 = tips up)
panel.angle_y  = 90;         % [deg] cant: tilt about Y-axis (0 = flat, >0 = tips toward +X)
panel.x_center = 2.28;        % [m] absolute X-coordinate of panel hinge center
                            %     leave [] to default to center of Cube B

% ----- Surface optical properties -----
% (specular + diffuse + absorptivity = 1)
cube_rho_s   = 0.10;   cube_rho_d   = 0.10;   cube_alpha = 0.80;
panel_rho_s  = 0.05;   panel_rho_d  = 0.05;   panel_alpha = 0.90;
tube_rho_s   = 0.10;   tube_rho_d   = 0.10;   tube_alpha = 0.80;

% ----- Orbit parameters -----
orbit.altitude = 800;        % [km] circular orbit altitude
orbit.incl     = 98.8;       % [deg] SSO

% ----- Magnetic residual dipole (A·m²) in body frame -----
mag.dipole_body = [0.1; 0.05; 0.02];   % [A·m²]

% ----- Simulation settings -----
sim.t_span    = [0, 100*3600];  % [s] 1000 hours
sim.plot_last = 0.5;               % [hr] plot window (last N hours)
sim.q0        = [1;0;0;0];      % initial quaternion [q0; q1; q2; q3] scalar-first
sim.omega0    = [0;0;0];        % [rad/s] initial body angular velocity

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
x_cubeA_end   = a;
x_tube_start  = a;
x_tube_end    = a + Lc;
x_cubeB_start = a + Lc;
x_cubeB_end   = a + Lc + b;

% Centroids in body frame
r_cubeA = [a/2;           0; 0];
r_tube  = [a + Lc/2;      0; 0];
r_cubeB = [a + Lc + b/2;  0; 0];

% Panel position and orientation
% Default X-center to middle of Cube B if not specified
if isempty(panel.x_center)
    panel_x_center = x_cubeB_start + b/2;
else
    panel_x_center = panel.x_center;
end

alpha_p = deg2rad(panel.angle);    % dihedral (about X-axis)
beta_p  = deg2rad(panel.angle_y);  % cant    (about Y-axis)
Sp = panel.span;

% Rotation matrices
Rx_a  = [1, 0, 0; 0, cos(alpha_p), -sin(alpha_p); 0, sin(alpha_p), cos(alpha_p)];
Rx_na = [1, 0, 0; 0, cos(alpha_p),  sin(alpha_p); 0,-sin(alpha_p), cos(alpha_p)]; % Rx(-alpha)
Ry_b  = [cos(beta_p), 0, sin(beta_p); 0, 1, 0; -sin(beta_p), 0, cos(beta_p)];

% Combined rotations (cant applied after dihedral)
%   +Y panel: dihedral tips +Z via Rx(+alpha), then cant via Ry(beta)
%   -Y panel: dihedral tips +Z via Rx(-alpha), then cant via Ry(beta)
R_panelP = Ry_b * Rx_a;    % +Y panel rotation
R_panelN = Ry_b * Rx_na;   % -Y panel rotation

% Span directions from hinge (unrotated: ±Y)
span_dirP = R_panelP * [0; 1; 0];   % +Y panel span
span_dirN = R_panelN * [0;-1; 0];   % -Y panel span

% Panel centroids: hinge at cube ±Y edge, center at hinge + (Sp/2)*span_dir
r_panelP = [panel_x_center; b/2; 0] + (Sp/2) * span_dirP;
r_panelN = [panel_x_center;-b/2; 0] + (Sp/2) * span_dirN;

% Orbit
orbit.R      = const.Re + orbit.altitude * 1e3;
orbit.n      = sqrt(const.mu / orbit.R^3);
orbit.V      = sqrt(const.mu / orbit.R);
orbit.period = 2*pi / orbit.n;

[rho_atm, ~] = atm_density_exp(orbit.altitude);

%% ==============  MASS PROPERTIES  =======================================

% --- Cube inertia about own centroid: I = m*s²/6 for each axis ---
I_cubeA = (cubeA.mass * a^2 / 6) * eye(3);
I_cubeB = (cubeB.mass * b^2 / 6) * eye(3);

% --- Hollow cylinder inertia about own centroid ---
I_tube_xx = tube.mass * (Ro^2 + Ri^2) / 2;                          % axial
I_tube_yy = tube.mass * (3*(Ro^2 + Ri^2) + Lc^2) / 12;             % transverse
I_tube = diag([I_tube_xx, I_tube_yy, I_tube_yy]);

% --- Panel inertia about own centroid (thin plate) ---
%   panel.length along local-X, panel.span along local-Y, thin in local-Z
Lp = panel.length;
I_panel_local = panel.mass / 12 * diag([Sp^2, Lp^2, Lp^2 + Sp^2]);

% Rotate to body frame using combined dihedral+cant rotation (different per panel)
I_panelP_body = R_panelP * I_panel_local * R_panelP';   % +Y panel
I_panelN_body = R_panelN * I_panel_local * R_panelN';   % -Y panel

% --- Composite mass and centroid ---
M_total = cubeA.mass + tube.mass + cubeB.mass + 2*panel.mass;
r_cm = (cubeA.mass * r_cubeA + tube.mass * r_tube + cubeB.mass * r_cubeB ...
        + panel.mass * r_panelP + panel.mass * r_panelN) / M_total;

% --- Parallel axis theorem to composite CM ---
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

%% ==============  PROPAGATE ATTITUDE DYNAMICS  ===========================
fprintf('\nPropagating attitude for %.1f hours...\n', diff(sim.t_span)/3600);

x0 = [sim.q0; sim.omega0];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',10);
[t, X] = ode15s(@(t,x) attitude_eom(t, x, I_total, orbit, const, ...
               surfaces, mag, rho_atm), sim.t_span, x0, opts);

%[t, X] = ode45(@(t,x) attitude_eom(t, x, I_total, orbit, const, ...
               %surfaces, mag, rho_atm), sim.t_span, x0, opts);

% Renormalize quaternions
qnorm = sqrt(sum(X(:,1:4).^2, 2));
X(:,1:4) = X(:,1:4) ./ qnorm;
omega = X(:,5:7);

%% ==============  COMPUTE TORQUE TIME HISTORIES  =========================
N = length(t);
T_gg   = zeros(N,3);
T_srp  = zeros(N,3);
T_aero = zeros(N,3);
T_mag  = zeros(N,3);
T_tot  = zeros(N,3);

for k = 1:N
    q = X(k,1:4)';
    [Tgg, Tsrp, Taero, Tmag] = compute_torques(t(k), q, orbit, const, ...
                                     surfaces, mag, rho_atm);
    T_gg(k,:)   = Tgg';
    T_srp(k,:)  = Tsrp';
    T_aero(k,:) = Taero';
    T_mag(k,:)  = Tmag';
    T_tot(k,:)  = (Tgg + Tsrp + Taero + Tmag)';
end

%% ==============  CUMULATIVE ANGULAR IMPULSE (REACTION WHEEL SIZING)  ====
%  H(t) = integral_0^t T(tau) dtau   [N·m·s]
%  This is the momentum that reaction wheels must absorb if attitude is
%  continuously corrected to zero error.

H_gg   = cumtrapz(t, T_gg);
H_srp  = cumtrapz(t, T_srp);
H_aero = cumtrapz(t, T_aero);
H_mag  = cumtrapz(t, T_mag);
H_tot  = cumtrapz(t, T_tot);

% Peak-to-peak momentum swing per axis (worst-case storage requirement)
H_range = max(H_tot) - min(H_tot);   % [N·m·s] per axis

% Net secular accumulation over simulation
H_net = H_tot(end,:);                 % [N·m·s] per axis
H_rate = H_net / diff(sim.t_span);   % [N·m] average secular torque

fprintf('\n====== REACTION WHEEL SIZING SUMMARY ======\n');
fprintf('Simulation duration: %.1f hours (%.1f orbits)\n', ...
    diff(sim.t_span)/3600, diff(sim.t_span)/orbit.period);
fprintf('\nPeak-to-peak momentum swing (max storage required):\n');
fprintf('  H_x = %.4e  N·m·s\n', H_range(1));
fprintf('  H_y = %.4e  N·m·s\n', H_range(2));
fprintf('  H_z = %.4e  N·m·s\n', H_range(3));
fprintf('  |H|  = %.4e  N·m·s  (RSS)\n', norm(H_range));
fprintf('\nNet secular momentum accumulated over simulation:\n');
fprintf('  H_x = %.4e  N·m·s\n', H_net(1));
fprintf('  H_y = %.4e  N·m·s\n', H_net(2));
fprintf('  H_z = %.4e  N·m·s\n', H_net(3));
fprintf('\nAverage secular torque (momentum build rate):\n');
fprintf('  T_x = %.4e  N·m\n', H_rate(1));
fprintf('  T_y = %.4e  N·m\n', H_rate(2));
fprintf('  T_z = %.4e  N·m\n', H_rate(3));
fprintf('\nPer-orbit cyclic momentum (peak-to-peak / num_orbits):\n');
n_orbits = diff(sim.t_span) / orbit.period;
fprintf('  H_x = %.4e  N·m·s/orbit\n', H_range(1) / n_orbits);
fprintf('  H_y = %.4e  N·m·s/orbit\n', H_range(2) / n_orbits);
fprintf('  H_z = %.4e  N·m·s/orbit\n', H_range(3) / n_orbits);

%% ==============  PLOTS  =================================================
t_hr = t / 3600;
t_end = t_hr(end);
mask = t_hr >= (t_end - sim.plot_last);

% --- Torque components ---
figure('Name','Perturbation Torques','Position',[100 100 1200 800]);

subplot(2,2,1);
plot(t_hr(mask), T_gg(mask,:)*1e6, 'LineWidth', 1.2);
xlabel('Time [hr]'); ylabel('Torque [\muN\cdotm]');
title('Gravity-Gradient Torque');
legend('τ_x','τ_y','τ_z'); grid on;

subplot(2,2,2);
plot(t_hr(mask), T_srp(mask,:)*1e6, 'LineWidth', 1.2);
xlabel('Time [hr]'); ylabel('Torque [\muN\cdotm]');
title('Solar Radiation Pressure Torque');
legend('τ_x','τ_y','τ_z'); grid on;

subplot(2,2,3);
plot(t_hr(mask), T_aero(mask,:)*1e6, 'LineWidth', 1.2);
xlabel('Time [hr]'); ylabel('Torque [\muN\cdotm]');
title('Aerodynamic Drag Torque');
legend('τ_x','τ_y','τ_z'); grid on;

subplot(2,2,4);
plot(t_hr(mask), T_mag(mask,:)*1e6, 'LineWidth', 1.2);
xlabel('Time [hr]'); ylabel('Torque [\muN\cdotm]');
title('Residual Magnetic Dipole Torque');
legend('τ_x','τ_y','τ_z'); grid on;

sgtitle(sprintf('Disturbance Torques  (last %.0f hr of %.0f hr sim)', ...
    sim.plot_last, diff(sim.t_span)/3600));

% --- Attitude & rates ---
figure('Name','Attitude & Angular Rates','Position',[100 100 900 800]);
tiledlayout(2,1);

nexttile;
plot(t_hr(mask), omega(mask,:) * 180/pi, 'LineWidth', 1.2);
xlabel('Time [hr]'); ylabel('Angular Rate [deg/s]');
title('Body Angular Velocity Components');
legend('\omega_x','\omega_y','\omega_z'); grid on;

euler = zeros(N,3);
for k = 1:N
    q = X(k,1:4)';
    C_bi = quat2dcm_sfirst(q);
    theta_LVLH = orbit.n * t(k);
    C_li = Rz(-theta_LVLH - pi/2) * Rx(deg2rad(orbit.incl));
    C_bl = C_bi * C_li';
    euler(k,:) = dcm2euler321(C_bl) * 180/pi;
end

nexttile;
euler_plot = euler(mask,:);
t_plot = t_hr(mask);
for col = 1:3
    jumps = find(abs(diff(euler_plot(:,col))) > 300);
    euler_plot(jumps, col) = NaN;
end
plot(t_plot, euler_plot, 'LineWidth', 1.2);
xlabel('Time [hr]'); ylabel('Angle [deg]');
title('Attitude Angles (3-2-1 Euler from LVLH)');
legend('Roll (\phi)','Pitch (\theta)','Yaw (\psi)'); grid on;
ylim([-180 180]); yticks(-180:90:180);

% --- Total torque magnitude ---
figure('Name','Total Torque Magnitude','Position',[100 100 900 400]);
plot(t_hr(mask), vecnorm(T_tot(mask,:),2,2)*1e6, 'LineWidth', 1.2, 'Color',[0.8 0.2 0.2]);
xlabel('Time [hr]'); ylabel('|T_{total}| [\muN\cdotm]');
title(sprintf('Total Disturbance Torque Magnitude  (last %.0f hr)', sim.plot_last));
grid on;

% --- Cumulative angular impulse (full simulation) ---
figure('Name','Cumulative Angular Impulse','Position',[100 100 1200 700]);
tiledlayout(2,2);

nexttile;
plot(t_hr, H_tot, 'LineWidth', 1.3);
xlabel('Time [hr]'); ylabel('H [N\cdotm\cdots]');
title('Total Cumulative Angular Impulse (All Sources)');
legend('H_x','H_y','H_z'); grid on;

nexttile;
plot(t_hr, H_gg, 'LineWidth', 1.1);
hold on;
plot(t_hr, H_srp, '--', 'LineWidth', 1.1);
xlabel('Time [hr]'); ylabel('H [N\cdotm\cdots]');
title('Cumulative Impulse: GG (solid) & SRP (dashed)');
legend('GG_x','GG_y','GG_z','SRP_x','SRP_y','SRP_z'); grid on;

nexttile;
plot(t_hr, H_aero, 'LineWidth', 1.1);
hold on;
plot(t_hr, H_mag, '--', 'LineWidth', 1.1);
xlabel('Time [hr]'); ylabel('H [N\cdotm\cdots]');
title('Cumulative Impulse: Aero (solid) & Mag (dashed)');
legend('Aero_x','Aero_y','Aero_z','Mag_x','Mag_y','Mag_z'); grid on;

nexttile;
plot(t_hr, vecnorm(H_tot,2,2), 'LineWidth', 1.5, 'Color', [0.8 0.15 0.15]);
xlabel('Time [hr]'); ylabel('|H| [N\cdotm\cdots]');
title('Total Momentum Magnitude (RW Storage Envelope)');
grid on;

sgtitle('Reaction Wheel Momentum Budget — Cumulative Angular Impulse');

% --- Per-source contribution bar chart ---
figure('Name','Momentum Budget Breakdown','Position',[100 100 700 450]);
H_final_abs = [abs(H_gg(end,:)); abs(H_srp(end,:)); abs(H_aero(end,:)); abs(H_mag(end,:))];
bar(H_final_abs);
set(gca, 'XTickLabel', {'Gravity-Gradient','SRP','Aerodynamic','Magnetic'});
ylabel('|H| at end of sim [N\cdotm\cdots]');
title('Net Angular Impulse by Source and Axis');
legend('X','Y','Z'); grid on;

fprintf('Done.\n');

%% ========================================================================
%                        LOCAL FUNCTIONS
%% ========================================================================

function xdot = attitude_eom(t, x, I, orbit, const, surfaces, mag, rho_atm)
    q = x(1:4);
    w = x(5:7);
    q = q / norm(q);

    [Tgg, Tsrp, Taero, Tmag] = compute_torques(t, q, orbit, const, ...
                                     surfaces, mag, rho_atm);
    T_total = Tgg + Tsrp + Taero + Tmag;

    wdot = I \ (T_total - cross(w, I*w));

    Omega = [ 0,    w(3), -w(2),  w(1);
             -w(3), 0,     w(1),  w(2);
              w(2), -w(1), 0,     w(3);
             -w(1), -w(2), -w(3), 0   ];
    qdot = 0.5 * Omega * q;

    xdot = [qdot; wdot];
end

function [Tgg, Tsrp, Taero, Tmag] = compute_torques(t, q, orbit, const, ...
                                          surfaces, mag, rho_atm)
    C_bi = quat2dcm_sfirst(q);

    theta = orbit.n * t;
    r_hat_i = [cos(theta); sin(theta); 0];
    r_hat_b = C_bi * r_hat_i;

    % 1. Gravity-Gradient
    I_total = surfaces.I_total;
    Tgg = 3 * orbit.n^2 * cross(r_hat_b, I_total * r_hat_b);

    % Sun direction (simplified: fixed in inertial frame)
    sun_hat_i = [1; 0; 0];
    sun_hat_b = C_bi * sun_hat_i;

    % Velocity direction
    v_hat_i = [-sin(theta); cos(theta); 0];
    v_hat_b = C_bi * v_hat_i;

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

    % ---- Solar panels (two panels, each with front & back) ----
    %  Combined rotation: dihedral (Rx) then cant (Ry)
    %  +Y panel: Ry(beta) * Rx(+alpha)     -Y panel: Ry(beta) * Rx(-alpha)
    Lp = panel.length;
    A_panel = Lp * Sp;

    ca = cos(alpha_p);  sa = sin(alpha_p);
    cb = cos(beta_p);   sb = sin(beta_p);
    Ry_b  = [cb, 0, sb; 0, 1, 0; -sb, 0, cb];
    Rx_a  = [1, 0, 0; 0, ca, -sa; 0, sa, ca];
    Rx_na = [1, 0, 0; 0, ca,  sa; 0,-sa, ca];
    R_pY = Ry_b * Rx_a;
    R_nY = Ry_b * Rx_na;

    % +Y panel: span dir, front normal, centroid
    span_pY   = R_pY * [0;1;0];
    normal_pY = R_pY * [0;0;1];     % front-face normal
    r_pp = [panel_xc; b/2; 0] + (Sp/2) * span_pY - r_cm;

    s.normal = normal_pY;   s.area = A_panel;  s.r_cp = r_pp;
    s.rho_s = p_rs;  s.rho_d = p_rd;
    list = [list, s];
    s.normal = -normal_pY;  s.area = A_panel;  s.r_cp = r_pp;
    list = [list, s];

    % -Y panel: span dir, front normal, centroid
    span_nY   = R_nY * [0;-1;0];
    normal_nY = R_nY * [0;0;1];     % front-face normal
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
    x_A0, x_B0, panel_xc)

    figure('Name','Satellite Geometry','Position',[80 80 1200 800]);
    hold on;

    col_cubeA  = [0.55 0.65 0.75];
    col_cubeB  = [0.65 0.55 0.70];
    col_tube   = [0.50 0.50 0.50];
    col_panel  = [0.15 0.20 0.55];
    col_back   = [0.75 0.75 0.70];
    col_cm     = [1.0  0.15 0.15];
    col_axes   = [0.9 0.1 0.1; 0.1 0.7 0.1; 0.1 0.2 0.9];

    % --- Draw Cube A ---
    draw_cube(x_A0, -a/2, -a/2, a, a, a, col_cubeA);

    % --- Draw Cube B ---
    draw_cube(x_B0, -b/2, -b/2, b, b, b, col_cubeB);

    % --- Draw connecting cylinder ---
    n_circ = 40;
    theta_c = linspace(0, 2*pi, n_circ);
    x_ends = [a, a + Lc];
    [Theta, Xc] = meshgrid(theta_c, x_ends);
    Yc = Ro * cos(Theta);
    Zc = Ro * sin(Theta);
    surf(Xc, Yc, Zc, 'FaceColor', col_tube, 'EdgeColor', 'none', ...
        'FaceAlpha', 0.7, 'FaceLighting', 'gouraud');

    % --- Draw solar panels (with dihedral + cant) ---
    Lp = panel.length;
    ca = cos(alpha_p);  sa = sin(alpha_p);
    cb = cos(beta_p);   sb = sin(beta_p);
    Ry_b  = [cb, 0, sb; 0, 1, 0; -sb, 0, cb];
    Rx_a  = [1, 0, 0; 0, ca, -sa; 0, sa, ca];
    Rx_na = [1, 0, 0; 0, ca,  sa; 0,-sa, ca];

    for sign = [1, -1]
        % Hinge at cube edge
        hinge = [panel_xc; sign*b/2; 0];

        % Combined rotation for this panel
        if sign > 0
            R_p = Ry_b * Rx_a;     % +Y: Ry(beta) * Rx(+alpha)
        else
            R_p = Ry_b * Rx_na;    % -Y: Ry(beta) * Rx(-alpha)
        end

        % Local panel axes (before rotation):
        %   chord along X: [-Lp/2, +Lp/2]
        %   span along sign*Y: [0, Sp]
        % 4 corners in local frame (relative to hinge)
        corners_local = [
            -Lp/2, 0, 0;        % near-hinge, -X edge
             Lp/2, 0, 0;        % near-hinge, +X edge
             Lp/2, sign*Sp, 0;  % far-tip, +X edge
            -Lp/2, sign*Sp, 0;  % far-tip, -X edge
        ]';

        % Rotate corners and translate to hinge position
        corners_body = R_p * corners_local + hinge;

        fill3(corners_body(1,:), corners_body(2,:), corners_body(3,:), col_panel, ...
            'FaceAlpha', 0.9, 'EdgeColor', [0.1 0.1 0.4], 'LineWidth', 1.5);

        % Back face (slight offset along panel normal)
        n_face = R_p * [0;0;1];
        offset = 0.005 * n_face;
        fill3(corners_body(1,:) - offset(1), ...
              corners_body(2,:) - offset(2), ...
              corners_body(3,:) - offset(3), col_back, ...
            'FaceAlpha', 0.7, 'EdgeColor', 'none');
    end

    % --- CM marker ---
    plot3(r_cm(1), r_cm(2), r_cm(3), 'o', 'Color', col_cm, ...
        'MarkerSize', 10, 'MarkerFaceColor', col_cm, 'LineWidth', 2);
    text(r_cm(1), r_cm(2), r_cm(3) - 0.35, '  CM', ...
        'Color', col_cm, 'FontSize', 11, 'FontWeight', 'bold');

    % --- Body axes ---
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

    % --- Principal axes ---
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

    % --- Dimension annotations ---
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

    % --- Lighting and formatting ---
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
    % Draw a filled cube from corner (x0,y0,z0) with dimensions (dx,dy,dz)
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

function C = quat2dcm_sfirst(q)
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    C = [q0^2+q1^2-q2^2-q3^2,  2*(q1*q2+q0*q3),      2*(q1*q3-q0*q2);
         2*(q1*q2-q0*q3),      q0^2-q1^2+q2^2-q3^2,  2*(q2*q3+q0*q1);
         2*(q1*q3+q0*q2),      2*(q2*q3-q0*q1),      q0^2-q1^2-q2^2+q3^2];
end

function eul = dcm2euler321(C)
    pitch = -asin(C(1,3));
    roll  = atan2(C(2,3), C(3,3));
    yaw   = atan2(C(1,2), C(1,1));
    eul = [roll; pitch; yaw];
end

function C = Rx(a)
    C = [1, 0, 0; 0, cos(a), sin(a); 0, -sin(a), cos(a)];
end

function C = Rz(a)
    C = [cos(a), sin(a), 0; -sin(a), cos(a), 0; 0, 0, 1];
end