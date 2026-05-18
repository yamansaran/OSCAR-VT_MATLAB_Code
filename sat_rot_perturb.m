%% ========================================================================
%  SATELLITE ROTATIONAL TENDENCY ANALYSIS
%  Cylinder + Angled Solar Panel - Perturbation Torques
%  ========================================================================
%  Computes and propagates attitude disturbance torques on a satellite
%  modeled as a uniform-density cylinder with an angled flat-plate solar
%  panel attached to one end.
%
%  Body frame:  X = along cylinder axis (from base toward panel end)
%               Y = lateral 
%               Z = completes right-hand
%
%  Perturbations modeled:
%    1. Gravity-gradient torque
%    2. Solar radiation pressure (SRP) torque
%    3. Aerodynamic drag torque (free-molecular flow)
%    4. Residual magnetic dipole torque
%
%  Author : OSCAR@VT
%           Yaman Saran
%  ========================================================================
clear; clc; close all;

%% =====================  USER-DEFINED PARAMETERS  ========================
% ----- Cylinder body -----
cyl.length   = 4.19;        % [m]  length along X-axis
cyl.diameter = 1.88;        % [m]  outer diameter
cyl.mass     = 1420;        % [kg] total cylinder mass (includes internals)
cyl.cm_frac  = 0.75;        % CM location along cylinder axis
                            %   0 = base, 0.5 = geometric center, 1 = +X end

% ----- Solar panel (flat plate) -----
panel.width  = 6.14;        % [m]  panel span (along body-Y direction)
panel.length = 2.73;        % [m]  panel extent (unfolded length)
panel.mass   = 42;         % [kg] panel mass
panel.angle  = -30;         % [deg] tilt angle of panel w.r.t. cylinder
                           %       axis (0 = panel lies along +X)

% ----- Rod / boom connecting cylinder to panel -----
rod.length    = 1.5;      % [m]  rod length (0 = panel directly on cylinder)
rod.angle_z   = 45;       % [deg] rod tilt from +X axis toward +Z
                          %       90 = rod goes straight up (+Z)
                          %       0  = rod extends along +X (no Z offset)
rod.mass      = 5.0;      % [kg]  rod mass (uniform rod assumed)

% ----- Attachment locations -----
panel.attach_frac     = 1.0;   % along cylinder length (0=base, 1=+X end)
rod.cyl_attach_angle  = 0;     % [deg] around cylinder (0=+Z, 90=+Y, 180=-Z)
panel.rod_attach_frac = 0.5;   % along panel length (0=base, 0.5=center, 1=tip)

% ----- Surface optical properties -----
% (specular reflectivity, diffuse reflectivity, absorptivity must sum to 1)
cyl.rho_s   = 0.1;   % cylinder specular reflectivity
cyl.rho_d   = 0.1;   % cylinder diffuse reflectivity
cyl.alpha_a = 0.8;   % cylinder absorptivity

panel.rho_s   = 0.05;  % panel specular reflectivity
panel.rho_d   = 0.05;  % panel diffuse reflectivity
panel.alpha_a = 0.90;  % panel absorptivity

% ----- Orbit parameters -----
orbit.altitude = 850;      % [km] circular orbit altitude
orbit.incl     = 98.8;     % [deg] SSO orbit

% ----- Magnetic residual dipole (A*m^2) in body frame -----
mag.dipole_body = [0.1; 0.05; 0.02];  % [A*m^2]

% ----- Simulation settings -----
sim.t_span    = [0, 1000*3600]; % [s] simulate for 1000 hours
sim.plot_last = 2;            % [hr] only plot this many hours from the end
sim.q0        = [1;0;0;0];    % initial quaternion [q4; q1; q2; q3] (scalar-first)
sim.omega0    = [0;0;0];      % [rad/s] initial body angular velocity

% ----- Environment constants -----
const.mu      = 3.986004418e14;  % [m^3/s^2] Earth GM
const.Re      = 6371e3;          % [m] Earth radius
const.P_sun   = 4.56e-6;         % [N/m^2] solar radiation pressure at 1 AU
const.c       = 2.998e8;         % [m/s] speed of light
const.B0      = 7.94e15;         % [T*m^3] Earth magnetic dipole moment (IGRF approx)

%% ====================  DERIVED QUANTITIES  ==============================
R_cyl = cyl.diameter / 2;
L_cyl = cyl.length;

% Orbit
orbit.R = const.Re + orbit.altitude * 1e3;   % [m] orbital radius
orbit.n = sqrt(const.mu / orbit.R^3);        % [rad/s] mean motion
orbit.V = sqrt(const.mu / orbit.R);          % [m/s] orbital speed
orbit.period = 2*pi / orbit.n;               % [s]

% Atmospheric density (exponential model for LEO)
[rho_atm, H_scale] = atm_density_exp(orbit.altitude);

%% ==============  MASS PROPERTIES (Cylinder + Panel)  ====================
% --- Cylinder inertia about its own CM ---
% Ixx (axial) is unaffected by axial CM offset
Ixx_cyl = 0.5 * cyl.mass * R_cyl^2;

% Iyy, Izz for a uniform cylinder are about the geometric center;
% shift to the declared CM using inverse parallel-axis theorem
Iyy_geo = cyl.mass * (3*R_cyl^2 + L_cyl^2) / 12;
Izz_geo = Iyy_geo;
d_ax = (cyl.cm_frac - 0.5) * L_cyl;   % axial offset: declared CM minus geo center
Iyy_cyl = Iyy_geo - cyl.mass * d_ax^2;
Izz_cyl = Izz_geo - cyl.mass * d_ax^2;

% Cylinder centroid in body frame (offset along X-axis)
r_cyl = [cyl.cm_frac * L_cyl; 0; 0];

% --- Panel inertia about its own centroid (thin plate) ---
Lp = panel.length;
Wp = panel.width;
I_panel_local = panel.mass/12 * diag([Wp^2, Lp^2, Lp^2 + Wp^2]);

% Rod geometry: base on cylinder, extends at rod.angle_z from +X toward +Z
attach_x = panel.attach_frac * L_cyl;
psi = deg2rad(rod.angle_z);    % rod tilt angle from +X
rod_dir = [cos(psi); 0; sin(psi)];  % rod unit direction in body frame

% Rod base on cylinder surface (azimuthal position set by cyl_attach_angle)
gamma = deg2rad(rod.cyl_attach_angle);  % 0=+Z, 90=+Y, 180=-Z, 270=-Y
rod_base = [attach_x; R_cyl*sin(gamma); R_cyl*cos(gamma)];
rod_tip  = rod_base + rod.length * rod_dir;

% Panel direction (unit vector along panel length in body frame)
phi = deg2rad(panel.angle);
panel_dir = [cos(phi); 0; sin(phi)];

% Panel centroid: rod_tip is at fraction panel.rod_attach_frac along panel
f = panel.rod_attach_frac;
panel_base_pt = rod_tip - f * Lp * panel_dir;       % panel start edge
panel_centroid = rod_tip + (0.5 - f) * Lp * panel_dir; % geometric center
r_panel = panel_centroid;

% Rotation matrix: panel frame -> body frame (rotation about Y by -phi)
R_panel2body = [ cos(phi), 0, sin(phi);
                 0,        1, 0;
                -sin(phi), 0, cos(phi)];

I_panel_body = R_panel2body * I_panel_local * R_panel2body';

% --- Rod inertia about its own centroid (thin uniform rod) ---
r_rod = (rod_base + rod_tip) / 2;  % rod centroid
Lr = rod.length;
u_rod = rod_dir;
if abs(u_rod(2)) < 0.9
    v_rod = cross(u_rod, [0;1;0]); v_rod = v_rod / norm(v_rod);
else
    v_rod = cross(u_rod, [1;0;0]); v_rod = v_rod / norm(v_rod);
end
w_rod = cross(u_rod, v_rod);
R_rod2body = [u_rod, v_rod, w_rod];
I_rod_local = rod.mass * Lr^2 / 12 * diag([0, 1, 1]);
I_rod_body  = R_rod2body * I_rod_local * R_rod2body';

% --- Composite centroid (cylinder + rod + panel) ---
M_total = cyl.mass + rod.mass + panel.mass;
r_cm = (cyl.mass * r_cyl + rod.mass * r_rod + panel.mass * r_panel) / M_total;

% --- Parallel axis theorem for composite inertia about CM ---
d_cyl   = r_cyl   - r_cm;
d_rod   = r_rod   - r_cm;
d_panel = r_panel - r_cm;

I_cyl_cm   = diag([Ixx_cyl, Iyy_cyl, Izz_cyl]) + ...
             cyl.mass * (dot(d_cyl,d_cyl)*eye(3) - d_cyl*d_cyl');
I_rod_cm   = I_rod_body + ...
             rod.mass * (dot(d_rod,d_rod)*eye(3) - d_rod*d_rod');
I_panel_cm = I_panel_body + ...
             panel.mass * (dot(d_panel,d_panel)*eye(3) - d_panel*d_panel');

I_total = I_cyl_cm + I_rod_cm + I_panel_cm;

fprintf('====== COMPOSITE MASS PROPERTIES ======\n');
fprintf('Total mass       : %.2f kg\n', M_total);
fprintf('CM location (body): [%.4f, %.4f, %.4f] m\n', r_cm);
fprintf('Inertia tensor (kg*m^2):\n');
disp(I_total);

%% ==============  GEOMETRY FOR TORQUE CALCULATIONS  ======================
n_cyl_panels = 12;  % azimuthal facets for cylinder side
surfaces = build_surface_model(cyl, panel, R_cyl, L_cyl, r_cm, ...
                                phi, n_cyl_panels, panel_base_pt);
surfaces.I_total = I_total;

%% ==============  3D GEOMETRY VISUALIZATION  =============================
plot_satellite_geometry(R_cyl, L_cyl, panel, phi, r_cm, I_total, rod, rod_base, rod_tip, panel_base_pt);

%% ==============  PROPAGATE ATTITUDE DYNAMICS  ===========================
fprintf('\nPropagating attitude for %.1f hours...\n', diff(sim.t_span)/3600);

% State: [q4; q1; q2; q3; wx; wy; wz]  (7 x 1)
x0 = [sim.q0; sim.omega0];

opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',10);
[t, X] = ode45(@(t,x) attitude_eom(t, x, I_total, orbit, const, ...
               surfaces, cyl, panel, mag, rho_atm), sim.t_span, x0, opts);

% Renormalize quaternions
qnorm = sqrt(sum(X(:,1:4).^2, 2));
X(:,1:4) = X(:,1:4) ./ qnorm;

omega = X(:,5:7);   % body rates [rad/s]

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
                                     surfaces, cyl, panel, mag, rho_atm);
    T_gg(k,:)   = Tgg';
    T_srp(k,:)  = Tsrp';
    T_aero(k,:) = Taero';
    T_mag(k,:)  = Tmag';
    T_tot(k,:)  = (Tgg + Tsrp + Taero + Tmag)';
end

%% ==============  PLOTS  =================================================
t_hr = t / 3600;
t_end = t_hr(end);
mask = t_hr >= (t_end - sim.plot_last);   % only plot the last N hours

figure('Name','Perturbation Torques','Position',[100 100 1200 800]);

subplot(2,2,1);
plot(t_hr(mask), T_gg(mask,:)*1e6, 'LineWidth', 1.2);
xlabel('Time [hr]'); ylabel('Torque [\muN\cdotm]');
title('Gravity-Gradient Torque');
legend('T_x','T_y','T_z'); grid on;

subplot(2,2,2);
plot(t_hr(mask), T_srp(mask,:)*1e6, 'LineWidth', 1.2);
xlabel('Time [hr]'); ylabel('Torque [\muN\cdotm]');
title('Solar Radiation Pressure Torque');
legend('T_x','T_y','T_z'); grid on;

subplot(2,2,3);
plot(t_hr(mask), T_aero(mask,:)*1e6, 'LineWidth', 1.2);
xlabel('Time [hr]'); ylabel('Torque [\muN\cdotm]');
title('Aerodynamic Drag Torque');
legend('T_x','T_y','T_z'); grid on;

subplot(2,2,4);
plot(t_hr(mask), T_mag(mask,:)*1e6, 'LineWidth', 1.2);
xlabel('Time [hr]'); ylabel('Torque [\muN\cdotm]');
title('Residual Magnetic Dipole Torque');
legend('T_x','T_y','T_z'); grid on;

sgtitle(sprintf('Disturbance Torques  (last %.0f hr of %.0f hr sim)', ...
    sim.plot_last, diff(sim.t_span)/3600));

figure('Name','Attitude & Angular Rates','Position',[100 100 900 800]);

tiledlayout(2,1);

% --- Angular velocity ---
nexttile;
plot(t_hr(mask), omega(mask,:) * 180/pi, 'LineWidth', 1.2);
xlabel('Time [hr]'); ylabel('Angular Rate [deg/s]');
title('Body Angular Velocity Components');
legend('\omega_x','\omega_y','\omega_z');
grid on;

% --- Euler angles (3-2-1 from LVLH) ---
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

% Break wrap-around lines: replace samples at discontinuities with NaN
for col = 1:3
    jumps = find(abs(diff(euler_plot(:,col))) > 300);
    euler_plot(jumps, col) = NaN;
end

plot(t_plot, euler_plot, 'LineWidth', 1.2);
xlabel('Time [hr]'); ylabel('Angle [deg]');
title('Attitude Angles (3-2-1 Euler from LVLH)');
legend('Roll (\phi)','Pitch (\theta)','Yaw (\psi)');
grid on;
ylim([-180 180]);
yticks(-180:90:180);

fprintf('Done.\n');

% --- Total torque magnitude ---
figure('Name','Total Torque Magnitude','Position',[100 100 900 400]);
plot(t_hr(mask), vecnorm(T_tot(mask,:),2,2)*1e6, 'LineWidth', 1.2, 'Color',[0.8 0.2 0.2]);
xlabel('Time [hr]'); ylabel('|T_{total}| [\muN\cdotm]');
title(sprintf('Total Disturbance Torque Magnitude  (last %.0f hr)', sim.plot_last));
grid on;


%% ========================================================================
%                        LOCAL FUNCTIONS
%% ========================================================================

function xdot = attitude_eom(t, x, I, orbit, const, surfaces, cyl, panel, mag, rho_atm)
    q = x(1:4);
    w = x(5:7);
    q = q / norm(q);

    [Tgg, Tsrp, Taero, Tmag] = compute_torques(t, q, orbit, const, ...
                                     surfaces, cyl, panel, mag, rho_atm);
    T_total = Tgg + Tsrp + Taero + Tmag;

    % Euler's equation: I*wdot = T - w x (I*w)
    wdot = I \ (T_total - cross(w, I*w));

    % Quaternion kinematics (scalar-first: q = [q4; q1; q2; q3])
    Omega = [ 0,    w(3), -w(2),  w(1);
             -w(3), 0,     w(1),  w(2);
              w(2), -w(1), 0,     w(3);
             -w(1), -w(2), -w(3), 0   ];
    qdot = 0.5 * Omega * q;

    xdot = [qdot; wdot];
end

function [Tgg, Tsrp, Taero, Tmag] = compute_torques(t, q, orbit, const, ...
                                          surfaces, cyl, panel, mag, rho_atm)
    C_bi = quat2dcm_sfirst(q);

    theta = orbit.n * t;
    r_hat_i = [cos(theta); sin(theta); 0];
    r_hat_b = C_bi * r_hat_i;

    % --- 1. Gravity-Gradient Torque ---
    I_total = surfaces.I_total;
    Tgg = 3 * orbit.n^2 * cross(r_hat_b, I_total * r_hat_b);

    % --- Sun direction (simplified: fixed in inertial frame) ---
    sun_hat_i = [1; 0; 0];
    sun_hat_b = C_bi * sun_hat_i;

    % Velocity direction in inertial frame
    v_hat_i = [-sin(theta); cos(theta); 0];
    v_hat_b = C_bi * v_hat_i;

    % --- 2 & 3. SRP and Aero torques ---
    Tsrp = [0;0;0];
    Taero = [0;0;0];

    for s = 1:length(surfaces.list)
        surf = surfaces.list(s);
        n_hat = surf.normal;
        A     = surf.area;
        r_cp  = surf.r_cp;
        rho_s = surf.rho_s;
        rho_d = surf.rho_d;

        % -- SRP --
        cos_sun = dot(sun_hat_b, n_hat);
        if cos_sun > 0
            f_srp = -const.P_sun * A * ( ...
                     (1 - rho_s)*sun_hat_b + ...
                     2*(rho_s*cos_sun + rho_d/3)*n_hat ) * cos_sun;
            Tsrp = Tsrp + cross(r_cp, f_srp);
        end

        % -- Aerodynamic drag --
        cos_v = dot(v_hat_b, n_hat);
        if cos_v > 0
            Cd = 2.2;
            f_aero = -0.5 * rho_atm * orbit.V^2 * Cd * A * cos_v * v_hat_b;
            Taero = Taero + cross(r_cp, f_aero);
        end
    end

    % --- 4. Magnetic torque ---
    B_i = earth_mag_field(r_hat_i, orbit.R, const);
    B_b = C_bi * B_i;
    Tmag = cross(mag.dipole_body, B_b);
end

function surfaces = build_surface_model(cyl, panel, R, L, r_cm, phi, n_az, panel_base_pt)
    idx = 0;
    list = [];

    % --- Cylinder side panels ---
    dtheta = 2*pi / n_az;
    for i = 1:n_az
        th = (i - 0.5) * dtheta;
        idx = idx + 1;
        s.normal = [0; cos(th); sin(th)];
        s.area   = R * dtheta * L;
        s.r_cp = [L/2; R*cos(th); R*sin(th)] - r_cm;
        s.rho_s = cyl.rho_s;
        s.rho_d = cyl.rho_d;
        list = [list, s]; %#ok<AGROW>
    end

    % --- Cylinder end caps ---
    A_cap = pi * R^2;

    idx = idx + 1;
    s.normal = [1; 0; 0];
    s.area   = A_cap;
    s.r_cp   = [L; 0; 0] - r_cm;
    s.rho_s  = cyl.rho_s;
    s.rho_d  = cyl.rho_d;
    list = [list, s];

    idx = idx + 1;
    s.normal = [-1; 0; 0];
    s.area   = A_cap;
    s.r_cp   = [0; 0; 0] - r_cm;
    s.rho_s  = cyl.rho_s;
    s.rho_d  = cyl.rho_d;
    list = [list, s];

    % --- Solar panel (two sides) ---
    A_panel = panel.length * panel.width;
    panel_dir = [cos(phi); 0; sin(phi)];
    panel_center = panel_base_pt + (panel.length/2) * panel_dir;
    r_panel_cm = panel_center - r_cm;

    n_front = [-sin(phi); 0; cos(phi)];
    n_back  = -n_front;

    s.normal = n_front;
    s.area   = A_panel;
    s.r_cp   = r_panel_cm;
    s.rho_s  = panel.rho_s;
    s.rho_d  = panel.rho_d;
    list = [list, s];

    s.normal = n_back;
    s.area   = A_panel;
    s.r_cp   = r_panel_cm;
    s.rho_s  = panel.rho_s;
    s.rho_d  = panel.rho_d;
    list = [list, s];

    surfaces.list = list;
    surfaces.I_total = [];
end

function [rho, H] = atm_density_exp(alt_km)
    alt_ref = [200, 300, 400, 500, 600, 700, 800, 900, 1000];
    rho_ref = [2.53e-10, 6.24e-12, 3.73e-13, 4.89e-14, ...
               9.15e-15, 2.39e-15, 8.19e-16, 3.58e-16, 1.84e-16];
    H_ref   = [37.5, 53.6, 58.5, 63.8, 71.8, 88.7, 124.6, 181.1, 268.0];

    if alt_km < alt_ref(1)
        rho = rho_ref(1);
        H = H_ref(1);
    elseif alt_km >= alt_ref(end)
        rho = rho_ref(end);
        H = H_ref(end);
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

function plot_satellite_geometry(R, L, panel, phi, r_cm, I_total, rod, rod_base, rod_tip, panel_base_pt)
    figure('Name','Satellite Geometry','Position',[80 80 1100 750]);

    col_cyl     = [0.65 0.70 0.75];
    col_cap     = [0.55 0.60 0.65];
    col_panel   = [0.15 0.20 0.55];
    col_back    = [0.75 0.75 0.70];
    col_cm      = [1.0  0.15 0.15];
    col_axes    = [0.9 0.1 0.1;
                   0.1 0.7 0.1;
                   0.1 0.2 0.9];

    n_circ = 60;
    n_len  = 2;
    theta_c = linspace(0, 2*pi, n_circ);
    x_c     = linspace(0, L, n_len);
    [Theta, Xc] = meshgrid(theta_c, x_c);
    Yc = R * cos(Theta);
    Zc = R * sin(Theta);

    surf(Xc, Yc, Zc, 'FaceColor', col_cyl, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.85, 'FaceLighting', 'gouraud');
    hold on;

    theta_cap = linspace(0, 2*pi, n_circ);
    r_cap     = linspace(0, R, 10);
    [Tc, Rc] = meshgrid(theta_cap, r_cap);
    Yc_cap = Rc .* cos(Tc);
    Zc_cap = Rc .* sin(Tc);

    Xc_base = zeros(size(Yc_cap));
    surf(Xc_base, Yc_cap, Zc_cap, 'FaceColor', col_cap, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.9, 'FaceLighting', 'gouraud');

    Xc_top = L * ones(size(Yc_cap));
    surf(Xc_top, Yc_cap, Zc_cap, 'FaceColor', col_cap, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.9, 'FaceLighting', 'gouraud');

    Lp = panel.length;
    Wp = panel.width;
    pb = panel_base_pt;
    rt = rod_tip;

    n_pu = 20;  n_pv = 10;
    u_p = linspace(0, Lp, n_pu);
    v_p = linspace(-Wp/2, Wp/2, n_pv);
    [Up, Vp] = meshgrid(u_p, v_p);

    Xp = pb(1) + Up * cos(phi);
    Yp = pb(2) + Vp;
    Zp = pb(3) + Up * sin(phi);

    surf(Xp, Yp, Zp, 'FaceColor', col_panel, ...
        'EdgeColor', [0.3 0.3 0.5], 'EdgeAlpha', 0.15, ...
        'FaceAlpha', 0.92, 'FaceLighting', 'gouraud');

    panel_normal = [-sin(phi); 0; cos(phi)];
    offset = 0.005;
    surf(Xp - offset*panel_normal(1), Yp, Zp - offset*panel_normal(3), ...
        'FaceColor', col_back, 'EdgeColor', 'none', ...
        'FaceAlpha', 0.7, 'FaceLighting', 'gouraud');

    bx = pb(1) + [0, Lp*cos(phi), Lp*cos(phi), 0, 0];
    by = pb(2) + [-Wp/2, -Wp/2, Wp/2, Wp/2, -Wp/2];
    bz = pb(3) + [0, Lp*sin(phi), Lp*sin(phi), 0, 0];
    plot3(bx, by, bz, 'Color', [0.1 0.1 0.4], 'LineWidth', 1.8);

    plot3(rt(1), rt(2), rt(3), 'd', 'Color', [0.8 0.5 0.0], ...
        'MarkerSize', 8, 'MarkerFaceColor', [0.9 0.6 0.1], 'LineWidth', 1.5);

    rb = rod_base;
    plot3([rb(1), rt(1)], [rb(2), rt(2)], [rb(3), rt(3)], ...
        'Color', [0.3 0.3 0.3], 'LineWidth', 3.0);
    plot3(rb(1), rb(2), rb(3), 'o', 'Color', [0.4 0.4 0.4], ...
        'MarkerSize', 5, 'MarkerFaceColor', [0.4 0.4 0.4]);

    rod_mid = (rb + rt) / 2;
    if rod.length > 0
        text(rod_mid(1) + 0.08, rod_mid(2) + 0.08, rod_mid(3) + 0.08, ...
            sprintf('rod: %.2f m', rod.length), ...
            'FontSize', 8, 'Color', [0.4 0.4 0.4]);
    end

    plot3(r_cm(1), r_cm(2), r_cm(3), 'o', 'Color', col_cm, ...
        'MarkerSize', 10, 'MarkerFaceColor', col_cm, 'LineWidth', 2);
    text(r_cm(1), r_cm(2), r_cm(3) - R*0.35, '  CM', ...
        'Color', col_cm, 'FontSize', 11, 'FontWeight', 'bold');

    ax_len = max(L, Lp) * 0.45;
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

    [V_princ, D_princ] = eig(I_total);
    I_princ = diag(D_princ);
    [I_princ, sortIdx] = sort(I_princ);
    V_princ = V_princ(:, sortIdx);
    pa_len = ax_len * 0.7;
    pa_colors = [1.0 0.5 0.0; 0.6 0.0 0.8; 0.0 0.7 0.7];
    for i = 1:3
        v = V_princ(:,i);
        quiver3(r_cm(1), r_cm(2), r_cm(3), ...
                v(1)*pa_len, v(2)*pa_len, v(3)*pa_len, 0, ...
                'Color', pa_colors(i,:), 'LineWidth', 1.5, ...
                'LineStyle', '--', 'MaxHeadSize', 0.3);
    end

    y_dim = -R - 0.15*L;
    plot3([0, L], [y_dim, y_dim], [0, 0], 'k-', 'LineWidth', 1);
    plot3([0, 0], [y_dim-0.05, y_dim+0.05], [0, 0], 'k-', 'LineWidth', 1);
    plot3([L, L], [y_dim-0.05, y_dim+0.05], [0, 0], 'k-', 'LineWidth', 1);
    text(L/2, y_dim - 0.08*L, 0, sprintf('%.2f m', L), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.3 0.3 0.3]);

    x_dim = -0.15*L;
    plot3([x_dim, x_dim], [0, 0], [-R, R], 'k-', 'LineWidth', 1);
    plot3([x_dim-0.03, x_dim+0.03], [0, 0], [-R, -R], 'k-', 'LineWidth', 1);
    plot3([x_dim-0.03, x_dim+0.03], [0, 0], [R, R], 'k-', 'LineWidth', 1);
    text(x_dim - 0.08*L, 0, 0, sprintf('\\varnothing %.2f m', 2*R), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.3 0.3 0.3], ...
        'Rotation', 90);

    p_base = pb;
    p_tip  = pb + [Lp*cos(phi); 0; Lp*sin(phi)];
    p_mid  = (p_base + p_tip) / 2;
    offset_dir = [sin(phi); 0; cos(phi)] * 0.15;
    plot3([p_base(1), p_tip(1)] + offset_dir(1), ...
          [0, 0], ...
          [p_base(3), p_tip(3)] + offset_dir(3), 'k-', 'LineWidth', 1);
    text(p_mid(1) + offset_dir(1)*2, 0, p_mid(3) + offset_dir(3)*2, ...
        sprintf('%.2f m', Lp), 'FontSize', 9, 'Color', [0.3 0.3 0.3]);

    light('Position', [2*L, -3*L, 2*L], 'Style', 'infinite');
    light('Position', [-L, 2*L, L], 'Style', 'infinite', 'Color', [0.3 0.3 0.3]);
    material([0.6 0.7 0.3 10]);

    axis equal;
    grid on;
    xlabel('X_{body} [m]'); ylabel('Y_{body} [m]'); zlabel('Z_{body} [m]');
    title(sprintf(['Satellite Geometry  |  Cyl: %.2f m \\times \\varnothing%.2f m' ...
           '  |  Panel: %.2f \\times %.2f m @ %.0f\\circ  |  Rod: %.2f m @ %.0f\\circ'], ...
           L, 2*R, Lp, Wp, rad2deg(phi), rod.length, rod.angle_z), 'FontSize', 11);

    view(135, 25);

    info_str = { ...
        sprintf('I_1 = %.2f  kg{\\cdot}m^2  (min)', I_princ(1)), ...
        sprintf('I_2 = %.2f  kg{\\cdot}m^2  (mid)', I_princ(2)), ...
        sprintf('I_3 = %.2f  kg{\\cdot}m^2  (max)', I_princ(3))};

    h_leg = gobjects(3,1);
    for i = 1:3
        h_leg(i) = plot3(NaN, NaN, NaN, '--', 'Color', pa_colors(i,:), 'LineWidth', 1.5);
    end
    legend(h_leg, info_str, 'Location', 'northeast', 'FontSize', 9);

    set(gcf, 'Color', 'w');
    rotate3d on;
end