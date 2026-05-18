%=========================================================================
% DYNAMIC SAFETY MARGIN (DSM) & VELOCITY-RANGE CORRIDOR ANALYSIS
%=========================================================================
% Computes approach risk based on chaser thruster authority. For each LVLH axis,
%  solves an LP to find the max pure-translation braking force achievable 
% (zero net torque), then builds an anisotropic velocity-range corridor and 
% evaluates DSM/RI along a simulated CW approach trajectory.
%
% Loads thruster configuration from rcs_config.mat (output of
% rcs_trade_study_independent.m) instead of hardcoding layouts.
%
% Outputs:
%   1) Per-axis max braking acceleration (from LP over thruster layout)
%   2) Velocity-range corridor plot (v_max vs range, per axis & ellipsoidal)
%   3) Simulated approach with DSM & RI time histories
%   4) Go/no-go phase diagram
%
%=========================================================================
clear; clc; close all;

%% ====================== USER CONFIG ====================================

% --- RCS configuration file (from rcs_trade_study_independent.m) ---
rcs_mat_file = 'rcs_config.mat';

% APPROACH / MISSION PARAMETERS
r_kos       = 10.0;        % [m] keep-out sphere radius
t_react     = 0.008;         % [s] GNC latency + valve response + settling
dv_remaining = 50.0;       % [m/s] remaining delta-V budget
dv_abort_reserve = 10.0;   % [m/s] reserve for abort maneuver (untouchable)
n_orbit     = 1.1068e-3;   % [rad/s] target mean motion (GEO ~ 7.292e-5,
                            %   LEO 800km ~ 1.1068e-3).

% APPROACH TRAJECTORY INITIAL CONDITIONS (LVLH, R-bar approach example)
%   State: [x_R; x_I; x_C; v_R; v_I; v_C] in target LVLH
%   R = radial, I = in-track (along-track), C = cross-track
r0 = [25; 0; 0];          % [m] start 200 m radial
v0 = [-0.05; 0.0; 0.0];   % [m/s] closing at 5 cm/s radially
t_sim = 4000;               % [s] simulation duration
dt    = 1.0;                % [s] timestep

% BODY-TO-LVLH DCM (identity = body aligned with LVLH)
%   Rows of C_BL map LVLH axes to body axes.
%   body X -> LVLH R,  body Y -> LVLH I,  body Z -> LVLH C
%   Modify if your body frame convention differs.
C_BL = eye(3);  % C_BL * v_LVLH = v_body

% RISK INDEX THRESHOLDS
RI_green  = 0.5;   % below this = comfortable
RI_yellow = 0.8;   % above this = caution
RI_red    = 1.0;   % above this = abort

%% ====================== LOAD RCS CONFIG FROM .MAT ======================
fprintf('Loading RCS configuration from: %s\n', rcs_mat_file);
if ~isfile(rcs_mat_file)
    error('RCS config file not found: %s\nRun rcs_trade_study_independent.m first.', rcs_mat_file);
end

cfg = load(rcs_mat_file);

% --- Extract Asterix config ---
ast.name     = cfg.asterix.name;
ast.mass     = cfg.asterix.mass;
ast.Lx       = cfg.asterix.Lx;
ast.Ly       = cfg.asterix.Ly;
ast.Lz       = cfg.asterix.Lz;
ast.com      = cfg.asterix.com;
ast.I        = cfg.asterix.I;
ast.thrust   = cfg.asterix.thrust;
ast.Isp      = cfg.asterix.Isp;
ast.thr_mass = cfg.asterix.thr_mass;
ast.thr_pos  = cfg.asterix.thr_pos;
ast.thr_dir  = cfg.asterix.thr_dir;
ast.ctrl_mode = cfg.asterix.ctrl_mode;

% --- Extract Obelix config ---
obx.name     = cfg.obelix.name;
obx.mass     = cfg.obelix.mass;
obx.Lx       = cfg.obelix.Lx;
obx.Ly       = cfg.obelix.Ly;
obx.Lz       = cfg.obelix.Lz;
obx.com      = cfg.obelix.com;
obx.I        = cfg.obelix.I;
obx.thrust   = cfg.obelix.thrust;
obx.Isp      = cfg.obelix.Isp;
obx.thr_mass = cfg.obelix.thr_mass;
obx.thr_pos  = cfg.obelix.thr_pos;
obx.thr_dir  = cfg.obelix.thr_dir;
obx.ctrl_mode = cfg.obelix.ctrl_mode;
obx.geom_center_global = cfg.obelix.geom_center_global;

N_THR_AST = cfg.asterix.n_thrusters;
N_THR_OBX = cfg.obelix.n_thrusters;

fprintf('  Asterix: %s mode, %d thrusters, %.1f kg\n', ast.ctrl_mode, N_THR_AST, ast.mass);
fprintf('  Obelix:  %s mode, %d thrusters, %.1f kg\n', obx.ctrl_mode, N_THR_OBX, obx.mass);
if isfield(cfg, 'timestamp')
    fprintf('  Config generated: %s\n', cfg.timestamp);
end
if isfield(cfg, 'settings')
    s = cfg.settings;
    fprintf('  Trade study weights: w_solo=%.1f, w_combined=%.1f, w_mating=%.1f\n', ...
        s.w_solo, s.w_combined, s.w_mating);
end

%% ====================== COMBINE INTO MATED STACK ========================
% Asterix body origin at global [0,0,0].  Obelix body origin offset by
% geom_center_global (from trade study).  They mate at the +X face of
% Asterix / -X face of Obelix.  Thruster directions are unchanged (bodies
% share the same orientation).
%=========================================================================
obx_offset = obx.geom_center_global;  % e.g. [+2.2, 0, 0]

% --- Combined mass ---
sc.mass = ast.mass + obx.mass;

% --- Combined COM in global frame ---
ast_com_global = ast.com;
obx_com_global = obx.com + obx_offset;
sc.com = (ast.mass * ast_com_global + obx.mass * obx_com_global) / sc.mass;

fprintf('\nAsterix COM (global): [%.3f, %.3f, %.3f] m\n', ast_com_global);
fprintf('Obelix  COM (global): [%.3f, %.3f, %.3f] m\n', obx_com_global);
fprintf('Combined COM (global): [%.3f, %.3f, %.3f] m\n', sc.com);

% --- Combined inertia about combined COM (parallel axis theorem) ---
d_ast = ast_com_global - sc.com;
d_obx = obx_com_global - sc.com;
I_shift = @(I_loc, m, d) I_loc + m * (dot(d,d)*eye(3) - d'*d);
sc.I = I_shift(ast.I, ast.mass, d_ast) + I_shift(obx.I, obx.mass, d_obx);

fprintf('Combined inertia about COM [kg*m^2]:\n');
fprintf('  Ixx=%.2f  Iyy=%.2f  Izz=%.2f\n', sc.I(1,1), sc.I(2,2), sc.I(3,3));
fprintf('  Ixy=%.2f  Ixz=%.2f  Iyz=%.2f\n', sc.I(1,2), sc.I(1,3), sc.I(2,3));

% --- Merge thruster layouts into global frame ---
%   Asterix thrusters: positions already in Asterix body frame (= global).
%   Obelix thrusters:  positions in Obelix body frame, shift by obx_offset.
%   Directions unchanged (same body orientation).
ast_thr_pos_global = ast.thr_pos;
obx_thr_pos_global = obx.thr_pos + obx_offset;

sc.thr_pos = [ast_thr_pos_global; obx_thr_pos_global];
sc.thr_dir = [ast.thr_dir;       obx.thr_dir];
sc.thrust  = ast.thrust;   % same thruster model on both

sc.name = sprintf('OSCAR (%s+%s)', ast.name, obx.name);
sc.Isp  = ast.Isp;

fprintf('\nMerged stack: %d %s + %d %s = %d thrusters\n', ...
        N_THR_AST, ast.name, N_THR_OBX, obx.name, N_THR_AST + N_THR_OBX);
fprintf('=== Dynamic Safety Margin Analysis for %s ===\n\n', sc.name);

%% ====================== BUILD FORCE/TORQUE MAP =========================
% B matrix: 6 x N_thr,  B * u = [F_body; tau_body]
% u_i in [0, 1] = duty cycle of thruster i
N_thr = size(sc.thr_pos, 1);
B = zeros(6, N_thr);
for i = 1:N_thr
    r_i = sc.thr_pos(i,:)' - sc.com';      % moment arm from COM
    f_i = sc.thrust * sc.thr_dir(i,:)';     % force vector
    B(1:3, i) = f_i;                        % force contribution
    B(4:6, i) = cross(r_i, f_i);            % torque contribution
end

fprintf('Thruster layout: %d thrusters, %.2f N each\n', N_thr, sc.thrust);
fprintf('Vehicle mass: %.1f kg\n\n', sc.mass);

%% ====================== LP: MAX BRAKING FORCE PER AXIS =================
% For each LVLH direction (±R, ±I, ±C), find max pure-translation force.
% Decision variable: u (N_thr x 1), each in [0,1]
% Maximize: d' * B_force * u   (i.e. minimize -d' * B_force * u)
% Subject to: B_torque * u = 0  (zero net torque for pure translation)
%             0 <= u <= 1

B_force  = B(1:3, :);   % 3 x N_thr (body frame forces)
B_torque = B(4:6, :);   % 3 x N_thr (body frame torques)

lb = zeros(N_thr, 1);
ub = ones(N_thr, 1);
opts = optimoptions('linprog', 'Display', 'off', 'Algorithm', 'dual-simplex');

% LVLH direction labels and unit vectors
dir_labels = {'+R','-R','+I','-I','+C','-C'};
dir_vecs_lvlh = [eye(3), -eye(3)];  % columns: +R,+I,+C,-R,-I,-C
dir_vecs_lvlh = [dir_vecs_lvlh(:,1), dir_vecs_lvlh(:,4), ...
                 dir_vecs_lvlh(:,2), dir_vecs_lvlh(:,5), ...
                 dir_vecs_lvlh(:,3), dir_vecs_lvlh(:,6)];

max_force = zeros(6,1);   % [+R,-R,+I,-I,+C,-C]
max_accel = zeros(6,1);

fprintf('--- Per-axis max pure-translation force (LP, zero torque) ---\n');
for k = 1:6
    d_lvlh = dir_vecs_lvlh(:, k);
    d_body = C_BL * d_lvlh;            % map LVLH direction to body frame
    
    % Minimize -d_body' * B_force * u
    f_obj = -(d_body' * B_force)';     % N_thr x 1
    
    % Equality constraint: B_torque * u = 0
    Aeq = B_torque;
    beq = zeros(3, 1);
    
    [u_opt, fval, exitflag] = linprog(f_obj, [], [], Aeq, beq, lb, ub, opts);
    
    if exitflag == 1
        F_vec = B_force * u_opt;       % resultant force in body frame
        max_force(k) = d_body' * F_vec; % component along desired direction
        max_accel(k) = max_force(k) / sc.mass;
    else
        % If zero-torque is infeasible, relax: allow small torque
        tau_tol = 0.5;  % [N·m] allowable residual torque
        A_ineq  = [B_torque; -B_torque];
        b_ineq  = tau_tol * ones(6, 1);
        [u_opt, fval, exitflag] = linprog(f_obj, A_ineq, b_ineq, ...
                                           [], [], lb, ub, opts);
        if exitflag == 1
            F_vec = B_force * u_opt;
            max_force(k) = d_body' * F_vec;
            max_accel(k) = max_force(k) / sc.mass;
            fprintf('  %s: %.4f N  (%.5f m/s^2) [torque relaxed]\n', ...
                    dir_labels{k}, max_force(k), max_accel(k));
        else
            fprintf('  %s: INFEASIBLE - no thrusters for this direction\n', ...
                    dir_labels{k});
            max_accel(k) = 0;
        end
        continue
    end
    fprintf('  %s: %.4f N  (%.6f m/s^2)\n', dir_labels{k}, ...
            max_force(k), max_accel(k));
end

% Extract braking acceleration per axis (deceleration = force opposite to
% velocity direction, so braking along +R uses the -R force and vice versa)
a_brake_R = max(max_accel(2), 1e-12);  % braking along -R (opposing +R vel)
a_brake_I = max(max_accel(4), 1e-12);  % braking along -I
a_brake_C = max(max_accel(6), 1e-12);  % braking along -C

% For positive velocity components, braking uses the opposite direction
a_brake = struct('Rpos', max_accel(2), 'Rneg', max_accel(1), ...
                 'Ipos', max_accel(4), 'Ineg', max_accel(3), ...
                 'Cpos', max_accel(6), 'Cneg', max_accel(5));

fprintf('\nBraking accelerations [m/s^2]:\n');
fprintf('  R-axis: +vel -> %.6f,  -vel -> %.6f\n', a_brake.Rpos, a_brake.Rneg);
fprintf('  I-axis: +vel -> %.6f,  -vel -> %.6f\n', a_brake.Ipos, a_brake.Ineg);
fprintf('  C-axis: +vel -> %.6f,  -vel -> %.6f\n', a_brake.Cpos, a_brake.Cneg);

%% ====================== CW PROPAGATION =================================
n = n_orbit;

    function Phi = CW_STM(n, dt)
        s = sin(n*dt); c = cos(n*dt);
        Phi = [
            4-3*c,      0, 0,   s/n,        2*(1-c)/n,     0;
            6*(s-n*dt), 1, 0,  -2*(1-c)/n,  (4*s-3*n*dt)/n, 0;
            0,          0, c,   0,           0,              s/n;
            3*n*s,      0, 0,   c,           2*s,            0;
            -6*n*(1-c), 0, 0,  -2*s,         4*c-3,          0;
            0,          0, -n*s, 0,           0,              c
        ];
    end

% Propagate free-drift approach trajectory
t_vec = 0:dt:t_sim;
N_steps = length(t_vec);
state = zeros(6, N_steps);
state(:,1) = [r0; v0];

for k = 2:N_steps
    Phi = CW_STM(n, dt);
    state(:,k) = Phi * state(:,k-1);
end

pos_R = state(1,:);  pos_I = state(2,:);  pos_C = state(3,:);
vel_R = state(4,:);  vel_I = state(5,:);  vel_C = state(6,:);
range_vec = sqrt(pos_R.^2 + pos_I.^2 + pos_C.^2);

%% ====================== DSM & RI COMPUTATION ============================
DSM_vec  = zeros(1, N_steps);
RI_vec   = zeros(1, N_steps);
vmax_vec = zeros(1, N_steps);
vclos_vec = zeros(1, N_steps);

for k = 1:N_steps
    r_vec = state(1:3, k);
    v_vec = state(4:6, k);
    range = norm(r_vec);
    r_hat = r_vec / max(range, 1e-10);
    
    % Closing velocity (positive = approaching target)
    v_closing = -dot(v_vec, r_hat);
    vclos_vec(k) = v_closing;
    
    if v_closing <= 0 || range <= r_kos
        if range <= r_kos
            DSM_vec(k) = 0;
            RI_vec(k)  = Inf;
        else
            DSM_vec(k) = Inf;
            RI_vec(k)  = 0;
        end
        vmax_vec(k) = NaN;
        continue
    end
    
    d_avail = range - r_kos;
    
    vR = vel_R(k);  vI = vel_I(k);  vC = vel_C(k);
    
    % Per-axis braking accel (direction-dependent)
    aR = (vR > 0) * a_brake.Rpos + (vR <= 0) * a_brake.Rneg;
    aI = (vI > 0) * a_brake.Ipos + (vI <= 0) * a_brake.Ineg;
    aC = (vC > 0) * a_brake.Cpos + (vC <= 0) * a_brake.Cneg;
    aR = max(aR, 1e-12);
    aI = max(aI, 1e-12);
    aC = max(aC, 1e-12);
    
    % Per-axis stopping distance
    d_stop_R = abs(vR)*t_react + vR^2/(2*aR);
    d_stop_I = abs(vI)*t_react + vI^2/(2*aI);
    d_stop_C = abs(vC)*t_react + vC^2/(2*aC);
    
    % Ellipsoidal RI
    RI_vec(k) = sqrt( (d_stop_R/d_avail)^2 + ...
                      (d_stop_I/d_avail)^2 + ...
                      (d_stop_C/d_avail)^2 );
    
    % Scalar DSM
    v_abs = [abs(vR); abs(vI); abs(vC)];
    a_eff = [aR; aI; aC];
    
    if v_closing > 1e-10
        w = v_abs / sum(v_abs + 1e-15);
        a_eff_scalar = dot(w, a_eff);
    else
        a_eff_scalar = min(a_eff);
    end
    
    d_stop_total = v_closing * t_react + v_closing^2 / (2*a_eff_scalar);
    DSM_vec(k) = d_avail / d_stop_total;
    
    vmax_vec(k) = -t_react*a_eff_scalar + ...
                  sqrt((t_react*a_eff_scalar)^2 + 2*a_eff_scalar*d_avail);
    
    dv_stop = v_closing;
    if dv_stop > (dv_remaining - dv_abort_reserve)
        DSM_vec(k) = min(DSM_vec(k), 0.1);
    end
end

%% ====================== VELOCITY-RANGE CORRIDOR ========================
range_plot = linspace(r_kos + 0.5, max(range_vec)*1.1, 500);

a_axes = [a_brake.Rpos, a_brake.Ipos, a_brake.Cpos];
a_worst = min(a_axes);
a_best  = max(a_axes);

vmax_corridor = zeros(3, length(range_plot));
axis_labels_short = {'R','I','C'};
for j = 1:length(range_plot)
    d = range_plot(j) - r_kos;
    for ax = 1:3
        a = a_axes(ax);
        vmax_corridor(ax, j) = -t_react*a + sqrt((t_react*a)^2 + 2*a*d);
    end
end
vmax_worst = -t_react*a_worst + sqrt((t_react*a_worst).^2 + 2*a_worst*(range_plot - r_kos));
vmax_best  = -t_react*a_best  + sqrt((t_react*a_best).^2  + 2*a_best*(range_plot - r_kos));

%% ====================== FIGURES ========================================

% --- Figure 1: Velocity-Range Corridor ---
figure('Name','Velocity-Range Corridor','Position',[100 100 900 600]);
hold on; grid on;
colors = lines(5);

fill([range_plot fliplr(range_plot)], ...
     [vmax_worst*RI_green fliplr(zeros(size(range_plot)))], ...
     [0.7 1.0 0.7], 'EdgeColor','none', 'FaceAlpha', 0.3);
fill([range_plot fliplr(range_plot)], ...
     [vmax_worst*RI_yellow fliplr(vmax_worst*RI_green)], ...
     [1.0 1.0 0.7], 'EdgeColor','none', 'FaceAlpha', 0.3);
fill([range_plot fliplr(range_plot)], ...
     [vmax_worst fliplr(vmax_worst*RI_yellow)], ...
     [1.0 0.85 0.7], 'EdgeColor','none', 'FaceAlpha', 0.3);

for ax = 1:3
    plot(range_plot, vmax_corridor(ax,:), '--', 'Color', colors(ax,:), ...
         'LineWidth', 1.2, 'DisplayName', sprintf('v_{max} %s-axis', axis_labels_short{ax}));
end
plot(range_plot, vmax_worst, 'r-', 'LineWidth', 2.5, 'DisplayName', 'v_{max} (worst axis) — ABORT LIMIT');

valid = vclos_vec > 0 & range_vec > r_kos;
scatter(range_vec(valid), vclos_vec(valid), 8, t_vec(valid), 'filled', ...
        'DisplayName', 'Approach trajectory');
cb = colorbar; ylabel(cb, 'Time [s]');

xline(r_kos, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Keep-Out Sphere');

xlabel('Range to Target [m]');
ylabel('Closing Speed [m/s]');
title(sprintf('Velocity-Range Corridor — %s (t_{react} = %.1f s)', sc.name, t_react));
legend('Location','northwest');
set(gca, 'XDir', 'reverse');
xlim([0 max(range_vec)*1.05]);

% --- Figure 2: DSM & RI Time Histories ---
figure('Name','DSM & RI vs Time','Position',[100 100 900 700]);

subplot(2,1,1);
plot(t_vec, range_vec, 'b-', 'LineWidth', 1.5); hold on; grid on;
yline(r_kos, 'r--', 'LineWidth', 1.5);
ylabel('Range [m]');
title(sprintf('Approach Timeline — %s', sc.name));
legend('Range','KOS');

subplot(2,1,2);
DSM_plot = min(DSM_vec, 10);
plot(t_vec, DSM_plot, 'b-', 'LineWidth', 1.5); hold on; grid on;
yline(1, 'r--', 'LineWidth', 1.5);
yline(2, 'Color', [0.9 0.7 0], 'LineStyle', '--');
ylabel('DSM');
legend('DSM','DSM = 1 (limit)','DSM = 2 (caution)');
ylim([0 min(max(DSM_plot)*1.1, 10)]);

% --- Figure 3: Relative Trajectory in LVLH ---
figure('Name','LVLH Trajectory','Position',[100 100 700 600]);
plot3(pos_I, pos_C, pos_R, 'b-', 'LineWidth', 1.2); hold on; grid on;
plot3(pos_I(1), pos_C(1), pos_R(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor','g');
plot3(0, 0, 0, 'r^', 'MarkerSize', 12, 'MarkerFaceColor','r');

[sx,sy,sz] = sphere(20);
surf(r_kos*sx, r_kos*sy, r_kos*sz, 'FaceAlpha', 0.15, ...
     'FaceColor','r', 'EdgeColor','none');

xlabel('In-Track [m]'); ylabel('Cross-Track [m]'); zlabel('Radial [m]');
title('Relative Trajectory (LVLH)');
legend('Trajectory','Start','Target','KOS','Location','best');
axis equal;
view(45, 30);

% --- Figure 4: Per-axis thrust authority summary ---
figure('Name','Thrust Authority','Position',[100 100 700 400]);
bar_data = [max_accel(1) max_accel(2); ...
            max_accel(3) max_accel(4); ...
            max_accel(5) max_accel(6)] * 1000;
b = bar(bar_data, 'grouped'); grid on;
b(1).FaceColor = [0.3 0.6 0.9];
b(2).FaceColor = [0.9 0.4 0.3];
set(gca, 'XTickLabel', {'R (radial)','I (in-track)','C (cross-track)'});
ylabel('Max Braking Accel [mm/s^2]');
title(sprintf('Anisotropic Thrust Authority — %s', sc.name));
legend('+direction','-direction','Location','best');

%% ====================== RISK BY RANGE (PARAMETRIC SWEEP) ================
range_sweep = linspace(r_kos + 0.1, 15, 600);
speed_sweep = linspace(0.001, 0.20, 400);

a_axes_brake = [a_brake.Rpos, a_brake.Rneg; ...
                a_brake.Ipos, a_brake.Ineg; ...
                a_brake.Cpos, a_brake.Cneg];
a_worst_brake = min(a_axes_brake(:));
a_per_axis_worst = min(a_axes_brake, [], 2);

vmax_R = zeros(size(range_sweep));
vmax_I = zeros(size(range_sweep));
vmax_C = zeros(size(range_sweep));
vmax_W = zeros(size(range_sweep));

for j = 1:length(range_sweep)
    d = range_sweep(j) - r_kos;
    for ax = 1:3
        a = a_per_axis_worst(ax);
        val = -t_react*a + sqrt((t_react*a)^2 + 2*a*d);
        switch ax
            case 1, vmax_R(j) = val;
            case 2, vmax_I(j) = val;
            case 3, vmax_C(j) = val;
        end
    end
    a = a_worst_brake;
    vmax_W(j) = -t_react*a + sqrt((t_react*a)^2 + 2*a*d);
end

[RR, VV] = meshgrid(range_sweep, speed_sweep);
DD = RR - r_kos;
D_stop = VV * t_react + VV.^2 / (2 * a_worst_brake);
RI_map = D_stop ./ DD;
RI_map(DD <= 0) = NaN;

%% --- Figure 5: Risk-by-Range Heatmap ---
figure('Name','Risk by Range','Position',[100 100 1000 700]);

%subplot(2,1,1);
%imagesc(range_sweep-10, (speed_sweep+0.05) * 100, RI_map);
imagesc(range_sweep(51:end)-10.5, (speed_sweep+0.05) * 100, RI_map(:, 51:end));
set(gca, 'YDir', 'normal', 'XDir', 'reverse');
hold on;

plot(range_sweep-10.5, (vmax_W+0.05) * 100, 'r-', 'LineWidth', 2.5);
plot(range_sweep-10.5, (vmax_R+0.05) * 100, '--', 'Color', [0.2 0.6 1], 'LineWidth', 1.5);
plot(range_sweep-10.5, (vmax_I+0.05) * 100, '--', 'Color', [1 0.5 0.1], 'LineWidth', 1.5);
plot(range_sweep-10.5, (vmax_C+0.05) * 100, '--', 'Color', [0.3 0.8 0.3], 'LineWidth', 1.5);

xline(r_kos, 'w--', 'LineWidth', 1.5);

cmap = [linspace(0.2,1,128)' linspace(0.8,0.2,128)' linspace(0.2,0.2,128)';
        linspace(1,0.6,128)' linspace(0.2,0,128)'   linspace(0.2,0,128)'];
colormap(gca, cmap);
caxis([0 2]);
cb = colorbar; ylabel(cb, 'Risk Index (RI)');

xlabel('Range to Target [m]');
ylabel('Closing Speed [cm/s]');
title(sprintf('Risk Index Map — %s  (t_{react} = %.1f s, worst-axis braking)', ...
      sc.name, t_react));
legend('v_{max} worst','v_{max} R','v_{max} I','v_{max} C', ...
       'Location','northwest');
set(gca, 'FontSize', 30, 'LineWidth', 2.25);

%% --- Figure 6: Max allowable speed vs range ---
subplot(2,1,2);
hold on; grid on;

fill([range_sweep fliplr(range_sweep)], ...
     [vmax_W*100*RI_green zeros(size(range_sweep))], ...
     [0.7 1 0.7], 'EdgeColor','none','FaceAlpha',0.35);
fill([range_sweep fliplr(range_sweep)], ...
     [vmax_W*100*RI_yellow vmax_W*100*RI_green], ...
     [1 1 0.6], 'EdgeColor','none','FaceAlpha',0.35);
fill([range_sweep fliplr(range_sweep)], ...
     [vmax_W*100 vmax_W*100*RI_yellow], ...
     [1 0.8 0.6], 'EdgeColor','none','FaceAlpha',0.35);

plot(range_sweep, vmax_R * 100, '-', 'Color', [0.2 0.6 1], 'LineWidth', 2, ...
     'DisplayName', 'R-axis limit');
plot(range_sweep, vmax_I * 100, '-', 'Color', [1 0.5 0.1], 'LineWidth', 2, ...
     'DisplayName', 'I-axis limit');
plot(range_sweep, vmax_C * 100, '-', 'Color', [0.3 0.8 0.3], 'LineWidth', 2, ...
     'DisplayName', 'C-axis limit');
plot(range_sweep, vmax_W * 100, 'r-', 'LineWidth', 2.5, ...
     'DisplayName', 'Worst-axis (ABORT)');

set(gca, 'XDir', 'reverse');
xlabel('Range to Target [m]');
ylabel('Max Allowable Closing Speed [cm/s]');
title('Velocity-Range Safe Corridor by Axis');
legend('Location','northwest');
xlim([r_kos max(range_sweep)]);

%% --- Console: tabulated risk-by-range ---
fprintf('\n============================================================\n');
fprintf('  MAX ALLOWABLE CLOSING SPEED BY RANGE\n');
fprintf('  (worst-direction braking per axis, t_react = %.1f s)\n', t_react);
fprintf('============================================================\n');
fprintf('  Range [m]  |  v_max R   v_max I   v_max C   v_max WORST\n');
fprintf('             |  [cm/s]    [cm/s]    [cm/s]    [cm/s]     \n');
fprintf('-------------------------------------------------------------\n');

table_ranges = [15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 250, 300];
table_ranges = table_ranges(table_ranges > r_kos);

for rr = table_ranges
    d = rr - r_kos;
    vmR = -t_react*a_per_axis_worst(1) + sqrt((t_react*a_per_axis_worst(1))^2 + 2*a_per_axis_worst(1)*d);
    vmI = -t_react*a_per_axis_worst(2) + sqrt((t_react*a_per_axis_worst(2))^2 + 2*a_per_axis_worst(2)*d);
    vmC = -t_react*a_per_axis_worst(3) + sqrt((t_react*a_per_axis_worst(3))^2 + 2*a_per_axis_worst(3)*d);
    vmW = -t_react*a_worst_brake + sqrt((t_react*a_worst_brake)^2 + 2*a_worst_brake*d);
    fprintf('  %7.1f    |  %7.3f   %7.3f   %7.3f   %7.3f\n', ...
            rr, vmR*100, vmI*100, vmC*100, vmW*100);
end
fprintf('-------------------------------------------------------------\n');

fprintf('\n  MINIMUM SAFE RANGE for given closing speed (worst-axis):\n');
fprintf('  Closing speed  |  Min range (RI=1)  |  Min range (RI=0.5)\n');
fprintf('  --------------------------------------------------------\n');
test_speeds = [0.01, 0.02, 0.03, 0.05, 0.10, 0.15, 0.20];
for vs = test_speeds
    d_stop = vs * t_react + vs^2 / (2 * a_worst_brake);
    r_min_abort   = d_stop + r_kos;
    r_min_comfort = d_stop / RI_green + r_kos;
    fprintf('    %5.2f cm/s    |    %7.2f m        |    %7.2f m\n', ...
            vs*100, r_min_abort, r_min_comfort);
end
fprintf('  --------------------------------------------------------\n');

%% ====================== SUMMARY TABLE ==================================
fprintf('\n========================================\n');
fprintf('  DYNAMIC SAFETY MARGIN SUMMARY\n');
fprintf('========================================\n');
fprintf('Vehicle:          %s\n', sc.name);
fprintf('Mass:             %.1f kg\n', sc.mass);
fprintf('Thrusters:        %d x %.2f N\n', N_thr, sc.thrust);
fprintf('  %s: %d thr (%s)\n', ast.name, N_THR_AST, ast.ctrl_mode);
fprintf('  %s: %d thr (%s)\n', obx.name, N_THR_OBX, obx.ctrl_mode);
fprintf('Reaction time:    %.1f s\n', t_react);
fprintf('Keep-out sphere:  %.1f m\n', r_kos);
fprintf('Orbit mean motion: %.4e rad/s\n', n);
fprintf('RCS config file:  %s\n', rcs_mat_file);
fprintf('----------------------------------------\n');
fprintf('Max braking accel [mm/s^2]:\n');
fprintf('  R+: %8.4f   R-: %8.4f\n', max_accel(1)*1e3, max_accel(2)*1e3);
fprintf('  I+: %8.4f   I-: %8.4f\n', max_accel(3)*1e3, max_accel(4)*1e3);
fprintf('  C+: %8.4f   C-: %8.4f\n', max_accel(5)*1e3, max_accel(6)*1e3);
fprintf('----------------------------------------\n');

[min_range, idx_min] = min(range_vec);
fprintf('Closest approach:  %.2f m at t = %.0f s\n', min_range, t_vec(idx_min));
fprintf('  DSM at closest:  %.3f\n', DSM_vec(idx_min));
fprintf('  RI  at closest:  %.3f\n', RI_vec(idx_min));
fprintf('  Closing speed:   %.4f m/s\n', vclos_vec(idx_min));

idx_yellow = find(RI_vec >= RI_yellow & RI_vec < Inf, 1, 'first');
idx_red    = find(RI_vec >= RI_red & RI_vec < Inf, 1, 'first');
if ~isempty(idx_yellow)
    fprintf('  RI hits CAUTION (%.1f) at t = %.0f s, range = %.1f m\n', ...
            RI_yellow, t_vec(idx_yellow), range_vec(idx_yellow));
end
if ~isempty(idx_red)
    fprintf('  RI hits ABORT   (%.1f) at t = %.0f s, range = %.1f m\n', ...
            RI_red, t_vec(idx_red), range_vec(idx_red));
end
if isempty(idx_yellow) && isempty(idx_red)
    fprintf('  RI stays GREEN throughout approach.\n');
end

fprintf('========================================\n');
fprintf('Done.\n');