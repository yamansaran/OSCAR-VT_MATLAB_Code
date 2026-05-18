%% RCS Configuration Detailed Control Authority Analysis
%  Shows exactly which thrusters fire for each maneuver direction:
%    - ±X, ±Y, ±Z translation
%    - ±Roll, ±Pitch, ±Yaw rotation
%
%  Uses constrained least-squares (lsqlin) from the Optimization Toolbox:
%    - Translation: achieve desired force, MINIMIZE parasitic torque
%    - Rotation:    achieve desired torque, MINIMIZE parasitic force
%    - Bounds:      each thruster either OFF or in [min_throttle, 1.0]
%
%  Two-pass solver handles the min-throttle dead zone:
%    Pass 1: solve relaxed (0 to 1), identify thrusters below min_throttle
%    Pass 2: lock those to zero, re-solve surviving thrusters in [min, 1]
%
%  SEPARATE ASTERIX & OBELIX CONFIGS - loaded from .mat file produced by
%  rcs_trade_study_independent.m
%
%  USAGE:
%    1. Run rcs_trade_study_independent.m  (produces .mat file)
%    2. Set MAT_FILE below to the output filename
%    3. Run this script
%
%  Yaman Saran - OSCAR@VT
%  -----------------------------------------------------------------------
clear; clc; close all;

%% ====================== SOLVER SETTINGS ================================
MIN_THROTTLE = 0.18;
MAX_THROTTLE = 1.30;
PARASITIC_TOL  = 0.01;
MAX_CORR_ITER  = 250;
target_force  = 3.0;   % [N]  desired net force for translation maneuvers
target_torque = 3.0;   % [Nm] desired net torque for rotation maneuvers

%% ====================== LOAD CONFIGURATION FROM .MAT ==================
MAT_FILE = 'rcs_config.mat';

fprintf('Loading configuration from: %s\n', MAT_FILE);
cfg = load(MAT_FILE);
fprintf('  Generated: %s\n\n', cfg.timestamp);

% --- Unpack Asterix ---
ast_thr_pos = cfg.asterix.thr_pos;
ast_thr_dir = cfg.asterix.thr_dir;

ast.name     = cfg.asterix.name;
ast.mass     = cfg.asterix.mass;
ast.Lx       = cfg.asterix.Lx;
ast.Ly       = cfg.asterix.Ly;
ast.Lz       = cfg.asterix.Lz;
ast.com      = cfg.asterix.com + [-0.3, 0, 0];
ast.I        = cfg.asterix.I;
ast.thrust   = cfg.asterix.thrust;
ast.Isp      = cfg.asterix.Isp;
ast.thr_mass = cfg.asterix.thr_mass;

% --- Unpack Obelix ---
obx_thr_pos = cfg.obelix.thr_pos;
obx_thr_dir = cfg.obelix.thr_dir;

obx.name     = cfg.obelix.name;
obx.mass     = cfg.obelix.mass + 50;
obx.Lx       = cfg.obelix.Lx;
obx.Ly       = cfg.obelix.Ly;
obx.Lz       = cfg.obelix.Lz;
obx.com      = cfg.obelix.com + [0.5, 0, 0];
obx.I        = cfg.obelix.I;
obx.thrust   = cfg.obelix.thrust;
obx.Isp      = cfg.obelix.Isp;
obx.thr_mass = cfg.obelix.thr_mass;
obx.ctrl_mode = cfg.obelix.ctrl_mode;
obx.geom_center_global = cfg.obelix.geom_center_global;

% --- Report what was loaded ---
fprintf('  %s: %d thrusters, ctrl_mode = %s\n', ast.name, size(ast_thr_pos,1), cfg.asterix.ctrl_mode);
fprintf('    sigma_min = %.4f   rank = %d   pure_fx = %s\n', ...
    cfg.asterix.sigma_min, cfg.asterix.rank, mat2str(cfg.asterix.pure_fx));
fprintf('  %s: %d thrusters, ctrl_mode = %s\n', obx.name, size(obx_thr_pos,1), obx.ctrl_mode);
if isfield(cfg.obelix, 'rot_sigma_min')
    fprintf('    rot_sigma_min = %.4f   rot_rank = %d\n', ...
        cfg.obelix.rot_sigma_min, cfg.obelix.rot_rank);
end
fprintf('  Combined: pure_ctrl = %s   sigma_min = %.4f   score = %.4f\n', ...
    mat2str(cfg.combined.pure_ctrl), cfg.combined.sigma_min, cfg.combined.score);
if isfield(cfg, 'alternatives')
    fprintf('  %d alternative configs available in .mat\n', length(cfg.alternatives));
end
fprintf('\n');

%% ====================== COMBINED PARAMETERS ============================
mating_plane_x = ast.Lx / 2;   % [m]

%% =======================================================================
% DERIVED QUANTITIES
%% =======================================================================

N_ast = size(ast_thr_pos, 1);
N_obx = size(obx_thr_pos, 1);
N_comb = N_ast + N_obx;

assert(size(ast_thr_dir,1) == N_ast, 'Asterix thr_pos/thr_dir row count mismatch');
assert(size(obx_thr_dir,1) == N_obx, 'Obelix thr_pos/thr_dir row count mismatch');

% Obelix inertia about its own COM
if ~isfield(obx, 'I') || isempty(obx.I)
    obx.Ixx = (1/12)*obx.mass*(obx.Ly^2 + obx.Lz^2);
    obx.Iyy = (1/12)*obx.mass*(obx.Lx^2 + obx.Lz^2);
    obx.Izz = (1/12)*obx.mass*(obx.Lx^2 + obx.Ly^2);
    obx.Ixx_com = obx.Ixx + obx.mass*(obx.com(2)^2 + obx.com(3)^2);
    obx.Iyy_com = obx.Iyy + obx.mass*(obx.com(1)^2 + obx.com(3)^2);
    obx.Izz_com = obx.Izz + obx.mass*(obx.com(1)^2 + obx.com(2)^2);
    obx.I = diag([obx.Ixx_com, obx.Iyy_com, obx.Izz_com]);
end

% Combined mass, COM, inertia
comb.mass = ast.mass + obx.mass + (N_ast + N_obx) * ast.thr_mass;
obx_com_global = obx.geom_center_global + obx.com;
comb.com = (ast.mass * ast.com + obx.mass * obx_com_global) / (ast.mass + obx.mass);

d_ast = comb.com - ast.com;
d_obx = comb.com - obx_com_global;
ast_Ixx = ast.I(1,1); ast_Iyy = ast.I(2,2); ast_Izz = ast.I(3,3);
obx_Ixx = obx.I(1,1); obx_Iyy = obx.I(2,2); obx_Izz = obx.I(3,3);
I_comb_xx = ast_Ixx + ast.mass*(d_ast(2)^2 + d_ast(3)^2) ...
          + obx_Ixx + obx.mass*(d_obx(2)^2 + d_obx(3)^2);
I_comb_yy = ast_Iyy + ast.mass*(d_ast(1)^2 + d_ast(3)^2) ...
          + obx_Iyy + obx.mass*(d_obx(1)^2 + d_obx(3)^2);
I_comb_zz = ast_Izz + ast.mass*(d_ast(1)^2 + d_ast(3)^2) ...
          + obx_Izz + obx.mass*(d_obx(1)^2 + d_obx(3)^2);
comb.I = diag([I_comb_xx, I_comb_yy, I_comb_zz]);

% Obelix thruster positions in global frame (for combined analysis)
obx_thr_pos_global = obx_thr_pos + obx.geom_center_global;

%% ====================== BUILD B MATRICES ===============================
% Asterix solo (about Asterix COM)
r_ast = ast_thr_pos - ast.com;
F_ast = ast.thrust * ast_thr_dir;
M_ast = cross(r_ast, F_ast, 2);
B_ast_solo = [F_ast'; M_ast'];                          % 6 × N_ast

% Obelix solo (about Obelix COM in body frame)
r_obx = obx_thr_pos - obx.com;
F_obx = obx.thrust * obx_thr_dir;
M_obx = cross(r_obx, F_obx, 2);
B_obx_solo = [F_obx'; M_obx'];                          % 6 × N_obx

% Combined (all thrusters about combined COM in global frame)
r_ca = ast_thr_pos      - comb.com;   % Asterix body frame ≡ global frame
r_co = obx_thr_pos_global - comb.com;
Fca = ast.thrust * ast_thr_dir;      Mca = cross(r_ca, Fca, 2);
Fco = obx.thrust * obx_thr_dir;      Mco = cross(r_co, Fco, 2);
B_comb = [Fca' Fco'; Mca' Mco'];                        % 6 × N_comb

%% ====================== BASIC METRICS =================================
rk_as = rank(B_ast_solo, 1e-8);  sv_as = svd(B_ast_solo);
rk_os = rank(B_obx_solo, 1e-8);  sv_os = svd(B_obx_solo);
rk_c  = rank(B_comb, 1e-8);      sv_c  = svd(B_comb);

fprintf('=============================================================\n');
fprintf('       RCS DETAILED CONTROL AUTHORITY ANALYSIS\n');
fprintf('       (Loaded from %s)\n', MAT_FILE);
fprintf('=============================================================\n\n');
fprintf('Solver: lsqlin (constrained least-squares)\n');
fprintf('Min throttle: %.0f%% (thrusters either OFF or >= %.0f%%)\n', ...
    MIN_THROTTLE*100, MIN_THROTTLE*100);
fprintf('Max throttle: %.0f%%\n\n', MAX_THROTTLE*100);

fprintf('--- %s (ctrl: %s) ---\n', ast.name, cfg.asterix.ctrl_mode);
fprintf('  Thrusters: %d     Thrust: %.2f N     Isp: %.0f s\n', N_ast, ast.thrust, ast.Isp);
fprintf('  Mass: %.1f kg    Dims: %.2f x %.2f x %.2f m\n', ast.mass, ast.Lx, ast.Ly, ast.Lz);
fprintf('  COM (body): (%+.3f, %+.3f, %+.3f) m\n', ast.com);
fprintf('  Solo: rank = %d/6   sigma_min = %.4f   cond = %.1f\n', ...
    rk_as, sv_as(end), sv_as(1)/max(sv_as(end),1e-15));

fprintf('\n--- %s (ctrl: %s) ---\n', obx.name, obx.ctrl_mode);
fprintf('  Thrusters: %d     Thrust: %.2f N     Isp: %.0f s\n', N_obx, obx.thrust, obx.Isp);
fprintf('  Mass: %.1f kg    Dims: %.2f x %.2f x %.2f m\n', obx.mass, obx.Lx, obx.Ly, obx.Lz);
fprintf('  COM (body): (%+.3f, %+.3f, %+.3f) m\n', obx.com);
fprintf('  Geom center (global): (%+.3f, %+.3f, %+.3f) m\n', obx.geom_center_global);
fprintf('  Solo: rank = %d/6   sigma_min = %.4f   cond = %.1f\n', ...
    rk_os, sv_os(end), sv_os(1)/max(sv_os(end),1e-15));
if isfield(cfg.obelix, 'rot_rank')
    B_obx_rot = B_obx_solo(4:6,:);
    sv_rot = svd(B_obx_rot);
    fprintf('  Attitude sub-B: rot_rank = %d/3   rot_sigma_min = %.4f\n', ...
        rank(B_obx_rot, 1e-8), sv_rot(end));
end

fprintf('\n--- Combined ---\n');
fprintf('  Total thrusters: %d (%d + %d)\n', N_comb, N_ast, N_obx);
fprintf('  COM (global): (%+.3f, %+.3f, %+.3f) m\n', comb.com);
fprintf('  Mass: %.1f kg\n', comb.mass);
fprintf('  rank = %d/6   sigma_min = %.4f   cond = %.1f\n', ...
    rk_c, sv_c(end), sv_c(1)/max(sv_c(end),1e-15));
fprintf('  Pure control (from trade study): %s\n\n', mat2str(cfg.combined.pure_ctrl));

%% ====================== THRUSTER TABLES ================================
half_Lx_ast = ast.Lx/2;
half_Lx_obx = obx.Lx/2;

fprintf('=============================================================\n');
fprintf('  THRUSTER LAYOUT - %s (body frame)\n', ast.name);
fprintf('=============================================================\n');
fprintf('  %3s  %24s   %14s   %-12s  Wrench [F; M]\n', ...
    'Thr', 'Position [m]', 'Direction', 'Face');
fprintf('  %s\n', repmat('-', 1, 100));
for i = 1:N_ast
    if     abs(ast_thr_pos(i,1) - half_Lx_ast) < 1e-6; face = '+X mating';
    elseif abs(ast_thr_pos(i,1) + half_Lx_ast) < 1e-6; face = '-X outer';
    else;  face = '';
    end
    w = B_ast_solo(:,i);
    wrench_str = sprintf('F(%+.2f,%+.2f,%+.2f) M(%+.3f,%+.3f,%+.3f)', w(1:3), w(4:6));
    fprintf('  T%d  (%+7.3f, %+7.3f, %+7.3f)  (%+d, %+d, %+d)   %-12s  %s\n', ...
        i, ast_thr_pos(i,:), ast_thr_dir(i,:), face, wrench_str);
end
fprintf('\n');

fprintf('=============================================================\n');
fprintf('  THRUSTER LAYOUT - %s (body frame)\n', obx.name);
fprintf('=============================================================\n');
fprintf('  %3s  %24s   %14s   %-12s  Wrench [F; M]\n', ...
    'Thr', 'Position [m]', 'Direction', 'Face');
fprintf('  %s\n', repmat('-', 1, 100));
for i = 1:N_obx
    if     abs(obx_thr_pos(i,1) + half_Lx_obx) < 1e-6; face = '-X mating';
    elseif abs(obx_thr_pos(i,1) - half_Lx_obx) < 1e-6; face = '+X outer';
    else;  face = '';
    end
    w = B_obx_solo(:,i);
    wrench_str = sprintf('F(%+.2f,%+.2f,%+.2f) M(%+.3f,%+.3f,%+.3f)', w(1:3), w(4:6));
    fprintf('  T%d  (%+7.3f, %+7.3f, %+7.3f)  (%+d, %+d, %+d)   %-12s  %s\n', ...
        i, obx_thr_pos(i,:), obx_thr_dir(i,:), face, wrench_str);
end
fprintf('\n');

%% ====================== B-MATRIX COLUMN TABLES =========================
dof_labels  = {'Fx','Fy','Fz','Mx','My','Mz'};
dof_names   = {'+X Trans','-X Trans','+Y Trans','-Y Trans', ...
               '+Z Trans','-Z Trans','+Roll','-Roll', ...
               '+Pitch','-Pitch','+Yaw','-Yaw'};

for sc_idx = 1:2
    if sc_idx == 1; B_show = B_ast_solo; n_show = N_ast; lbl = ast.name;
    else;           B_show = B_obx_solo; n_show = N_obx; lbl = obx.name;
    end
    fprintf('=============================================================\n');
    fprintf('  B-MATRIX COLUMNS (%s Solo)\n', lbl);
    fprintf('=============================================================\n');
    fprintf('  %3s', 'Thr');
    for d = 1:6; fprintf('%10s', dof_labels{d}); end
    fprintf('\n');
    fprintf('  %s\n', repmat('-', 1, 3 + 6*10));
    for i = 1:n_show
        fprintf('  T%d ', i);
        for d = 1:6
            val = B_show(d,i);
            if abs(val) < 1e-10; fprintf('%10s', '.');
            else;                 fprintf('%+10.4f', val);
            end
        end
        fprintf('\n');
    end
    fprintf('\n');
end

%% ====================== MANEUVER DEFINITIONS ==========================
maneuvers = struct();
idx = 0;
signs = [+1, -1];
for dof = 1:6
    for s = 1:2
        idx = idx + 1;
        maneuvers(idx).dof   = dof;
        maneuvers(idx).sign  = signs(s);
        maneuvers(idx).name  = dof_names{idx};
        maneuvers(idx).label = dof_labels{dof};
    end
end

%% ====================== SOLVE THRUSTER ALLOCATIONS =====================
lsq_opts = optimoptions('lsqlin', 'Display', 'off', ...
    'Algorithm', 'interior-point', ...
    'OptimalityTolerance', 1e-12, ...
    'ConstraintTolerance', 1e-10);

% --- Asterix solo ---
fprintf('=============================================================\n');
fprintf('  THRUSTER FIRING COMMANDS - %s SOLO\n', ast.name);
fprintf('=============================================================\n');
fprintf('  Objective: achieve desired DOF, minimize parasitic wrench.\n');
fprintf('  Solver: lsqlin with two-pass min-throttle = %.0f%%\n\n', MIN_THROTTLE*100);

ast_results = solve_all_maneuvers(B_ast_solo, maneuvers, N_ast, MIN_THROTTLE, MAX_THROTTLE, lsq_opts,target_force,target_torque);
print_maneuver_results(ast_results, maneuvers, N_ast, [ast.name ' Solo'], dof_labels);

% --- Obelix solo ---
fprintf('\n=============================================================\n');
fprintf('  THRUSTER FIRING COMMANDS - %s SOLO\n', obx.name);
fprintf('=============================================================\n\n');

obx_results = solve_all_maneuvers(B_obx_solo, maneuvers, N_obx, MIN_THROTTLE, MAX_THROTTLE, lsq_opts,target_force,target_torque);
print_maneuver_results(obx_results, maneuvers, N_obx, [obx.name ' Solo'], dof_labels);

% --- Combined ---
fprintf('\n=============================================================\n');
fprintf('  THRUSTER FIRING COMMANDS - COMBINED\n');
fprintf('=============================================================\n\n');

comb_results = solve_all_maneuvers(B_comb, maneuvers, N_comb, MIN_THROTTLE, MAX_THROTTLE, lsq_opts,target_force,target_torque);
print_maneuver_results(comb_results, maneuvers, N_comb, 'Combined', dof_labels);

%% ====================== QUICK-REFERENCE FIRING TABLES ==================

% ---- Asterix Solo ----
fprintf('\n=============================================================\n');
fprintf('  QUICK-REFERENCE FIRING TABLE (%s Solo)\n', ast.name);
fprintf('=============================================================\n');
fprintf('  X = full   %%% = partial (shown)   . = off\n');
fprintf('  Min throttle = %.0f%%\n\n', MIN_THROTTLE*100);
fprintf('  %-14s', 'Maneuver');
for i = 1:N_ast; fprintf('  T%d ', i); end
fprintf('  Parasitic\n');
fprintf('  %s\n', repmat('-', 1, 14 + N_ast*5 + 11));
for m = 1:12
    r = ast_results(m);
    fprintf('  %-14s', maneuvers(m).name);
    for i = 1:N_ast; fprintf('%s', fmt_thr(r.cmd(i))); end
    fprintf('  %8.4f', r.parasitic_norm);
    if ~r.min_thr_ok; fprintf('  *'); end
    fprintf('\n');
end
print_min_throttle_violations(ast_results, maneuvers);

% ---- Obelix Solo ----
fprintf('\n=============================================================\n');
fprintf('  QUICK-REFERENCE FIRING TABLE (%s Solo)\n', obx.name);
fprintf('=============================================================\n');
fprintf('  X = full   %%% = partial (shown)   . = off\n');
fprintf('  Min throttle = %.0f%%\n\n', MIN_THROTTLE*100);
fprintf('  %-14s', 'Maneuver');
for i = 1:N_obx; fprintf('  T%d ', i); end
fprintf('  Parasitic\n');
fprintf('  %s\n', repmat('-', 1, 14 + N_obx*5 + 11));
for m = 1:12
    r = obx_results(m);
    fprintf('  %-14s', maneuvers(m).name);
    for i = 1:N_obx; fprintf('%s', fmt_thr(r.cmd(i))); end
    fprintf('  %8.4f', r.parasitic_norm);
    if ~r.min_thr_ok; fprintf('  *'); end
    fprintf('\n');
end
print_min_throttle_violations(obx_results, maneuvers);

% ---- Combined ----
fprintf('\n=============================================================\n');
fprintf('  QUICK-REFERENCE FIRING TABLE (Combined)\n');
fprintf('=============================================================\n');
fprintf('  %s: A1–A%d (%d thr)    %s: O1–O%d (%d thr)\n\n', ...
    ast.name, N_ast, N_ast, obx.name, N_obx, N_obx);
fprintf('  %-14s |', 'Maneuver');
for i = 1:N_ast; fprintf(' A%d ', i); end
fprintf('|');
for i = 1:N_obx; fprintf(' O%d ', i); end
fprintf('| Parasitic\n');
fprintf('  %s\n', repmat('-', 1, 14 + 2 + N_ast*4 + 1 + N_obx*4 + 11));
for m = 1:12
    r = comb_results(m);
    fprintf('  %-14s |', maneuvers(m).name);
    for i = 1:N_ast; fprintf('%s', fmt_thr(r.cmd(i))); end
    fprintf('|');
    for i = 1:N_obx; fprintf('%s', fmt_thr(r.cmd(N_ast + i))); end
    fprintf('| %8.4f', r.parasitic_norm);
    if ~r.min_thr_ok; fprintf(' *'); end
    fprintf('\n');
end
print_min_throttle_violations(comb_results, maneuvers);

%% ====================== CORRECTED (PARASITIC-CANCELLED) SOLUTIONS ======
fprintf('\n=============================================================\n');
fprintf('  CORRECTED THRUSTER COMMANDS (Parasitic Cancellation)\n');
fprintf('=============================================================\n');
fprintf('  Iteratively cancels parasitic wrench via counter-firing.\n');
fprintf('  Tolerance = %.4f    Max iterations = %d\n', PARASITIC_TOL, MAX_CORR_ITER);
fprintf('  Min throttle = %.0f%%    Max throttle = %.0f%%\n\n', MIN_THROTTLE*100, MAX_THROTTLE*100);

ast_corrected  = solve_corrected(B_ast_solo, maneuvers, N_ast, ...
    MIN_THROTTLE, MAX_THROTTLE, PARASITIC_TOL, MAX_CORR_ITER, lsq_opts,target_force,target_torque);
obx_corrected  = solve_corrected(B_obx_solo, maneuvers, N_obx, ...
    MIN_THROTTLE, MAX_THROTTLE, PARASITIC_TOL, MAX_CORR_ITER, lsq_opts,target_force,target_torque);
comb_corrected = solve_corrected(B_comb, maneuvers, N_comb, ...
    MIN_THROTTLE, MAX_THROTTLE, PARASITIC_TOL, MAX_CORR_ITER, lsq_opts,target_force,target_torque);

% ---- Print detailed correction chains ----
for sc_idx = 1:3
    if sc_idx == 1
        lbl = [ast.name ' Solo']; corr = ast_corrected; nthr = N_ast;
    elseif sc_idx == 2
        lbl = [obx.name ' Solo']; corr = obx_corrected; nthr = N_obx;
    else
        lbl = 'Combined'; corr = comb_corrected; nthr = N_comb;
    end
    fprintf('  ---- %s Correction Details ----\n\n', lbl);
    for m = 1:12
        r = corr(m);
        fprintf('  --- %s ---\n', maneuvers(m).name);
        fprintf('  Iterations: %d    Converged: %s\n', r.n_iters, bool_str(r.converged));
        fprintf('  Parasitic: %.4f → %.4f  (%.1f%% reduction)\n', ...
            r.parasitic_initial, r.parasitic_norm, ...
            100*(1 - r.parasitic_norm/max(r.parasitic_initial, 1e-15)));
        fprintf('  Firing: ');
        for i = 1:nthr
            if r.cmd(i) > 0.01
                fprintf('T%d(%3.0f%%) ', i, r.cmd(i)*100);
            end
        end
        fprintf('\n');
        fprintf('  Final wrench: ');
        for d = 1:6
            if abs(r.w_actual(d)) > 1e-4
                fprintf('%s=%+.4f  ', dof_labels{d}, r.w_actual(d));
            end
        end
        fprintf('\n');
        if r.parasitic_norm > PARASITIC_TOL
            fprintf('  Residual parasitic: ');
            for d = 1:6
                if abs(r.parasitic(d)) > 1e-4
                    fprintf('%s=%+.4f  ', dof_labels{d}, r.parasitic(d));
                end
            end
            fprintf('  (||%.4f||)\n', r.parasitic_norm);
        else
            fprintf('  Parasitic fully cancelled\n');
        end
        if ~r.min_thr_ok
            fprintf('  [MIN THROTTLE NOT MET - lowest on = %.0f%%]\n', r.min_thr_actual*100);
        end
        fprintf('\n');
    end
end

%% ======= CORRECTED QUICK-REFERENCE FIRING TABLES ======================

% ---- Asterix Solo Corrected ----
fprintf('=============================================================\n');
fprintf('  CORRECTED QUICK-REFERENCE FIRING TABLE (%s Solo)\n', ast.name);
fprintf('=============================================================\n');
fprintf('  Total throttle after parasitic cancellation.\n');
fprintf('  X = full   %%% = partial   . = off\n\n');
fprintf('\n');
fprintf('  %-14s', 'Maneuver');
for i = 1:N_ast; fprintf('  T%d ', i); end
fprintf('  Para.Before  Para.After  Iters  Conv\n');
fprintf('  %s\n', repmat('-', 1, 14 + N_ast*5 + 48));
for m = 1:12
    r_cor = ast_corrected(m);
    fprintf('  %-14s', maneuvers(m).name);
for i = 1:N_ast; fprintf('%s', fmt_thr(r_cor.cmd(i))); end
    fprintf('  %10.4f  %10.4f  %5d  %4s', ...
        r_cor.parasitic_initial, r_cor.parasitic_norm, ...
        r_cor.n_iters, bool_str(r_cor.converged));
if ~r_cor.min_thr_ok; fprintf(' *'); end
    fprintf('\n');
end
fprintf('\n');
% ---- Obelix Solo Corrected ----
fprintf('=============================================================\n');
fprintf('  CORRECTED QUICK-REFERENCE FIRING TABLE (%s Solo)\n', obx.name);
fprintf('=============================================================\n\n');
fprintf('  %-14s', 'Maneuver');
for i = 1:N_obx; fprintf('  T%d ', i); end
fprintf('  Para.Before  Para.After  Iters  Conv\n');
fprintf('  %s\n', repmat('-', 1, 14 + N_obx*5 + 48));
for m = 1:12
    r_cor = obx_corrected(m);
    fprintf('  %-14s', maneuvers(m).name);
    for i = 1:N_obx; fprintf('%s', fmt_thr(r_cor.cmd(i))); end
    fprintf('  %10.4f  %10.4f  %5d  %4s', ...
        r_cor.parasitic_initial, r_cor.parasitic_norm, ...
        r_cor.n_iters, bool_str(r_cor.converged));
    if ~r_cor.min_thr_ok; fprintf(' *'); end
    fprintf('\n');
end
fprintf('\n');

% ---- Combined Corrected ----
fprintf('=============================================================\n');
fprintf('  CORRECTED QUICK-REFERENCE FIRING TABLE (Combined)\n');
fprintf('=============================================================\n');
fprintf('  %s: A1–A%d    %s: O1–O%d\n\n', ast.name, N_ast, obx.name, N_obx);
fprintf('  %-14s |', 'Maneuver');
for i = 1:N_ast; fprintf(' A%d ', i); end
fprintf('|');
for i = 1:N_obx; fprintf(' O%d ', i); end
fprintf('| Para.Bef  Para.Aft  It Conv\n');
hdr_len = 14 + 2 + N_ast*4 + 1 + N_obx*4 + 32;
fprintf('  %s\n', repmat('-', 1, hdr_len));
for m = 1:12
    r_cor = comb_corrected(m);
    fprintf('  %-14s |', maneuvers(m).name);
    for i = 1:N_ast; fprintf('%s', fmt_thr(r_cor.cmd(i))); end
    fprintf('|');
    for i = 1:N_obx; fprintf('%s', fmt_thr(r_cor.cmd(N_ast + i))); end
    fprintf('| %7.4f   %7.4f  %2d %4s', ...
        r_cor.parasitic_initial, r_cor.parasitic_norm, ...
        r_cor.n_iters, bool_str(r_cor.converged));
    if ~r_cor.min_thr_ok; fprintf(' *'); end
    fprintf('\n');
end
fprintf('\n');

%% ======= BEFORE/AFTER COMPARISON ======================================
fprintf('=============================================================\n');
fprintf('  BEFORE/AFTER PARASITIC COMPARISON\n');
fprintf('=============================================================\n');
fprintf('  How much parasitic wrench was eliminated by correction.\n\n');
fprintf('  %-14s  %20s  %20s  %20s\n', '', [ast.name ' Solo'], [obx.name ' Solo'], 'Combined');
fprintf('  %-14s  %9s → %9s  %9s → %9s  %9s → %9s\n', 'Maneuver', ...
    'Before', 'After', 'Before', 'After', 'Before', 'After');
fprintf('  %s\n', repmat('-', 1, 78));
for m = 1:12
    ra = ast_corrected(m);
    ro = obx_corrected(m);
    rc = comb_corrected(m);
    fprintf('  %-14s  %9.4f → %9.4f  %9.4f → %9.4f  %9.4f → %9.4f\n', ...
        maneuvers(m).name, ...
        ra.parasitic_initial, ra.parasitic_norm, ...
        ro.parasitic_initial, ro.parasitic_norm, ...
        rc.parasitic_initial, rc.parasitic_norm);
end
fprintf('\n');

%% ====================== SIGNED AUTHORITY SUMMARY =======================
fprintf('=============================================================\n');
fprintf('  SIGNED AUTHORITY SUMMARY\n');
fprintf('=============================================================\n');
fprintf('  Force/torque achieved and resulting acceleration.\n\n');
fprintf('  %-14s  %10s %12s  %10s %12s  %10s %12s  %10s %10s %10s\n', 'Maneuver', ...
    'Ast F/T', 'Ast Accel', 'Obx F/T', 'Obx Accel', 'Comb F/T', 'Comb Accel', ...
    'Ast Par.', 'Obx Par.', 'Comb Par.');
fprintf('  %s\n', repmat('-', 1, 130));
for m = 1:12
    ra = ast_results(m);
    ro = obx_results(m);
    rc = comb_results(m);
    dof = maneuvers(m).dof;

    wa = ra.w_actual(dof); wo = ro.w_actual(dof); wc = rc.w_actual(dof);

    if dof <= 3
        aa = wa / ast.mass;   ao = wo / obx.mass;   ac = wc / comb.mass;
        unit_ft = 'N';
    else
        ax = dof - 3;
        aa = wa / ast.I(ax,ax);  ao = wo / obx.I(ax,ax);  ac = wc / comb.I(ax,ax);
        unit_ft = 'Nm';
    end
    fprintf('  %-14s  %+8.4f %s %+10.6f  %+8.4f %s %+10.6f  %+8.4f %s %+10.6f   %8.4f  %8.4f  %8.4f\n', ...
        maneuvers(m).name, ...
        wa, unit_ft, aa, wo, unit_ft, ao, wc, unit_ft, ac, ...
        ra.parasitic_norm, ro.parasitic_norm, rc.parasitic_norm);
end
fprintf('\n');

%% ====================== AUTHORITY ASYMMETRY ============================
fprintf('=============================================================\n');
fprintf('  AUTHORITY ASYMMETRY (+ vs - per axis)\n');
fprintf('=============================================================\n');
fprintf('  Ratio of |+| to |-| authority.  1.0 = perfectly symmetric.\n\n');
fprintf('  %-10s  %10s %10s %10s\n', 'Axis', ast.name, obx.name, 'Combined');
fprintf('  %s\n', repmat('-', 1, 42));
axis_labels_asym = {'X Trans', 'Y Trans', 'Z Trans', 'Roll', 'Pitch', 'Yaw'};
for d = 1:6
    pos_a = abs(ast_results(2*d-1).w_actual(maneuvers(2*d-1).dof));
    neg_a = abs(ast_results(2*d).w_actual(maneuvers(2*d).dof));
    pos_o = abs(obx_results(2*d-1).w_actual(maneuvers(2*d-1).dof));
    neg_o = abs(obx_results(2*d).w_actual(maneuvers(2*d).dof));
    pos_c = abs(comb_results(2*d-1).w_actual(maneuvers(2*d-1).dof));
    neg_c = abs(comb_results(2*d).w_actual(maneuvers(2*d).dof));
    ratio_a = safe_ratio(pos_a, neg_a);
    ratio_o = safe_ratio(pos_o, neg_o);
    ratio_c = safe_ratio(pos_c, neg_c);
    fprintf('  %-10s  %10.3f %10.3f %10.3f\n', axis_labels_asym{d}, ratio_a, ratio_o, ratio_c);
end
fprintf('\n');

%% ====================== COUPLING MATRICES ==============================
short_names = {'+X','-X','+Y','-Y','+Z','-Z','+Rl','-Rl','+Pi','-Pi','+Yw','-Yw'};

for sc_idx = 1:3
    if sc_idx == 1;     res = ast_results;  lbl = [ast.name ' Solo'];
    elseif sc_idx == 2; res = obx_results;  lbl = [obx.name ' Solo'];
    else;               res = comb_results; lbl = 'Combined';
    end
    fprintf('=============================================================\n');
    fprintf('  COUPLING MATRIX (%s)\n', lbl);
    fprintf('=============================================================\n');
    fprintf('  Rows = wrench DOF produced.  Cols = commanded maneuver.\n');
    fprintf('  Normalized to commanded DOF magnitude.\n\n');
    fprintf('  Cmd DOF →  ');
    for m = 1:12; fprintf('%6s', short_names{m}); end
    fprintf('\n');
    fprintf('  %s\n', repmat('-', 1, 12 + 12*6));
    for d = 1:6
        fprintf('  %-9s  ', dof_labels{d});
        for m = 1:12
            r = res(m);
            cmd_dof_val = abs(r.w_actual(maneuvers(m).dof));
            if cmd_dof_val > 1e-10
                coupling = r.w_actual(d) / cmd_dof_val;
            else
                coupling = 0;
            end
            if abs(coupling) < 0.005
                fprintf('     .');
            else
                fprintf(' %+5.2f', coupling);
            end
        end
        fprintf('\n');
    end
    fprintf('\n');
end

%% ====================== PURE TRANSLATION PAIRS =========================
fprintf('=============================================================\n');
fprintf('  PURE TRANSLATION PAIRS (zero net torque)\n');
fprintf('=============================================================\n');
ax_names = {'X','Y','Z'};
for mode = 1:3
    if mode == 1
        fprintf('\n  %s SOLO (%d thrusters):\n', ast.name, N_ast);
        pos_t = ast_thr_pos; dir_t = ast_thr_dir; com_off = ast.com;
        nthr = N_ast; thrust_val = ast.thrust;
    elseif mode == 2
        fprintf('\n  %s SOLO (%d thrusters):\n', obx.name, N_obx);
        pos_t = obx_thr_pos; dir_t = obx_thr_dir; com_off = obx.com;
        nthr = N_obx; thrust_val = obx.thrust;
    else
        fprintf('\n  COMBINED (%d thrusters):\n', N_comb);
        pos_t = [ast_thr_pos; obx_thr_pos_global];
        dir_t = [ast_thr_dir; obx_thr_dir];
        com_off = comb.com;
        nthr = N_comb; thrust_val = ast.thrust;  % assumes same thrust
    end
    pos_rel = pos_t - com_off;
    for ax = 1:3
        for sgn = [+1.1, -1.1]
            e_ax = zeros(1,3); e_ax(ax) = sgn;
            if sgn > 0; lbl = ['+' ax_names{ax}];
            else;        lbl = ['-' ax_names{ax}];
            end
            fprintf('    %s: ', lbl);
            found = false;
            for i = 1:nthr
                F_i = thrust_val * dir_t(i,:);
                M_i = cross(pos_rel(i,:), F_i);
                if norm(M_i) < 1e-6 && dot(F_i/norm(F_i), e_ax) > 0.99
                    fprintf('T%d alone (%.3f N)  ', i, norm(F_i));
                    found = true;
                end
            end
            for i = 1:nthr
                for j = (i+1):nthr
                    F_net = thrust_val*(dir_t(i,:) + dir_t(j,:));
                    M_net = cross(pos_rel(i,:), thrust_val*dir_t(i,:)) + ...
                            cross(pos_rel(j,:), thrust_val*dir_t(j,:));
                    if norm(F_net) > 1e-6 && norm(M_net) < 1e-6
                        F_hat = F_net / norm(F_net);
                        if dot(F_hat, e_ax) > 0.99
                            fprintf('T%d+T%d (%.3f N)  ', i, j, norm(F_net));
                            found = true;
                        end
                    end
                end
            end
            if ~found; fprintf('(no pure pair - see lsqlin solution above)'); end
            fprintf('\n');
        end
    end
end

%% ====================== PURE ROTATION PAIRS ============================
fprintf('\n=============================================================\n');
fprintf('  PURE ROTATION PAIRS (zero net force)\n');
fprintf('=============================================================\n');
rot_names = {'Roll','Pitch','Yaw'};
for mode = 1:3
    if mode == 1
        fprintf('\n  %s SOLO (%d thrusters):\n', ast.name, N_ast);
        B_use = B_ast_solo; nthr = N_ast;
    elseif mode == 2
        fprintf('\n  %s SOLO (%d thrusters):\n', obx.name, N_obx);
        B_use = B_obx_solo; nthr = N_obx;
    else
        fprintf('\n  COMBINED (%d thrusters):\n', N_comb);
        B_use = B_comb; nthr = N_comb;
    end
    for ax = 1:3
        for sgn = [+1.1, -1.1]
            if sgn > 0; lbl = ['+' rot_names{ax}];
            else;        lbl = ['-' rot_names{ax}];
            end
            fprintf('    %s: ', lbl);
            found = false;
            for i = 1:nthr
                for j = (i+1):nthr
                    F_net = B_use(1:3, i) + B_use(1:3, j);
                    M_net = B_use(4:6, i) + B_use(4:6, j);
                    if norm(F_net) < 1e-6 && norm(M_net) > 1e-6
                        M_hat = M_net / norm(M_net);
                        e_ax = zeros(3,1); e_ax(ax) = sgn;
                        if dot(M_hat, e_ax) > 0.99
                            fprintf('T%d+T%d (%.4f Nm)  ', i, j, norm(M_net));
                            found = true;
                        end
                    end
                end
            end
            if ~found; fprintf('(no pure pair - see lsqlin solution above)'); end
            fprintf('\n');
        end
    end
end

%% ====================== FIGURES ========================================
fig_w = 750; fig_h = 500;

% ---- Asterix Solo 3D ----
fig1 = figure('Name', [ast.name ' Solo Layout'], 'Position', [20 350 fig_w fig_h]);
ax1 = axes(fig1); hold(ax1, 'on');
draw_box(ax1, ast.Lx, ast.Ly, ast.Lz, [0 0 0], [0.65 0.75 0.92], [0.3 0.3 0.5]);
draw_thrusters_labeled(ax1, ast_thr_pos, ast_thr_dir, max([ast.Lx ast.Ly ast.Lz])/2*0.32, ...
    [0.85 0.12 0.08]);
plot3(ax1, ast.com(1), ast.com(2), ast.com(3), 'k+', 'MarkerSize', 14, 'LineWidth', 2.5);
text(ax1, ast.com(1)+0.06, ast.com(2)+0.06, ast.com(3)+0.06, 'COM', ...
    'FontSize', 10, 'FontWeight', 'bold');
yr_m = [-ast.Ly/2 ast.Ly/2 ast.Ly/2 -ast.Ly/2];
zr_m = [-ast.Lz/2 -ast.Lz/2 ast.Lz/2 ast.Lz/2];
patch(ax1, half_Lx_ast*ones(1,4), yr_m, zr_m, ...
    'FaceColor', [1 0.85 0.3], 'FaceAlpha', 0.2, ...
    'EdgeColor', [0.8 0.6 0.0], 'LineWidth', 1.8, 'LineStyle', '--');
title(ax1, sprintf('%s Solo - %d thr  rank=%d  \\sigma_{min}=%.3f', ...
    ast.name, N_ast, rk_as, sv_as(end)), 'FontSize', 12);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
axis equal; grid on; view(135, 25);

% ---- Obelix Solo 3D ----
fig1b = figure('Name', [obx.name ' Solo Layout'], 'Position', [20+fig_w+10 350 fig_w fig_h]);
ax1b = axes(fig1b); hold(ax1b, 'on');
draw_box(ax1b, obx.Lx, obx.Ly, obx.Lz, [0 0 0], [0.6 0.9 0.7], [0.3 0.5 0.3]);
draw_thrusters_labeled(ax1b, obx_thr_pos, obx_thr_dir, max([obx.Lx obx.Ly obx.Lz])/2*0.32, ...
    [0.1 0.6 0.2]);
plot3(ax1b, obx.com(1), obx.com(2), obx.com(3), 'k+', 'MarkerSize', 14, 'LineWidth', 2.5);
text(ax1b, obx.com(1)+0.06, obx.com(2)+0.06, obx.com(3)+0.06, 'COM', ...
    'FontSize', 10, 'FontWeight', 'bold');
yr_m_o = [-obx.Ly/2 obx.Ly/2 obx.Ly/2 -obx.Ly/2];
zr_m_o = [-obx.Lz/2 -obx.Lz/2 obx.Lz/2 obx.Lz/2];
patch(ax1b, -half_Lx_obx*ones(1,4), yr_m_o, zr_m_o, ...
    'FaceColor', [1 0.85 0.3], 'FaceAlpha', 0.2, ...
    'EdgeColor', [0.8 0.6 0.0], 'LineWidth', 1.8, 'LineStyle', '--');
title(ax1b, sprintf('%s Solo (%s) - %d thr  rank=%d  \\sigma_{min}=%.3f', ...
    obx.name, obx.ctrl_mode, N_obx, rk_os, sv_os(end)), 'FontSize', 12);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
axis equal; grid on; view(135, 25);

% ---- Combined 3D ----
fig2 = figure('Name', 'Combined Layout', 'Position', [20 50 fig_w fig_h]);
ax2 = axes(fig2); hold(ax2, 'on');
draw_box(ax2, ast.Lx, ast.Ly, ast.Lz, [0 0 0], [0.65 0.75 0.92], [0.3 0.3 0.5]);
draw_box(ax2, obx.Lx, obx.Ly, obx.Lz, obx.geom_center_global, [0.6 0.9 0.7], [0.3 0.5 0.3]);
draw_thrusters_labeled(ax2, ast_thr_pos, ast_thr_dir, max([ast.Lx ast.Ly ast.Lz])/2*0.3, ...
    [0.85 0.12 0.08]);
draw_thrusters_labeled(ax2, obx_thr_pos_global, obx_thr_dir, max([obx.Lx obx.Ly obx.Lz])/2*0.3, ...
    [0.1 0.6 0.2]);
yr = [-max(ast.Ly,obx.Ly)/2 max(ast.Ly,obx.Ly)/2 max(ast.Ly,obx.Ly)/2 -max(ast.Ly,obx.Ly)/2]*1.15;
zr = [-max(ast.Lz,obx.Lz)/2 -max(ast.Lz,obx.Lz)/2 max(ast.Lz,obx.Lz)/2 max(ast.Lz,obx.Lz)/2]*1.15;
patch(ax2, mating_plane_x*ones(1,4), yr, zr, 'FaceColor', [1 0.85 0.3], ...
    'FaceAlpha', 0.15, 'EdgeColor', [0.8 0.6 0.0], 'LineWidth', 1.5, 'LineStyle', '--');
plot3(ax2, comb.com(1), comb.com(2), comb.com(3), 'kp', ...
    'MarkerSize', 16, 'MarkerFaceColor', 'k');
text(ax2, comb.com(1)+0.05, 0.05, 0.05, 'COM_{comb}', 'FontWeight','bold');
pure_str = ''; if cfg.combined.pure_ctrl; pure_str = '  PURE CTRL OK'; end
title(ax2, sprintf('Combined - %d+%d thr  rank=%d  \\sigma_{min}=%.3f%s', ...
    N_ast, N_obx, rk_c, sv_c(end), pure_str), 'FontSize', 12);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
axis equal; grid on; view(135, 25);

% ---- Firing heatmaps ----
for sc_idx = 1:3
    if sc_idx == 1
        res = ast_corrected; nthr = N_ast; lbl = [ast.name ' Solo (Corrected)'];
        xtick_fn = @(i) sprintf('T%d', i);
    elseif sc_idx == 2
        res = obx_corrected; nthr = N_obx; lbl = [obx.name ' Solo (Corrected)'];
        xtick_fn = @(i) sprintf('T%d', i);
    else
        res = comb_corrected; nthr = N_comb; lbl = 'Combined (Corrected)';
    end
    fig_hm = figure('Name', [lbl ' Heatmap'], ...
        'Position', [20 + (sc_idx-1)*350, 20, max(400, nthr*35 + 120), 440]);
    cmd_mat = zeros(12, nthr);
    for m = 1:12; cmd_mat(m,:) = res(m).cmd'; end
    imagesc(cmd_mat, [0 1]); colormap(hot);
    cb = colorbar; ylabel(cb, 'Command (0=off, 1=full)');
    set(gca, 'YTick', 1:12, 'YTickLabel', {maneuvers.name});
    xtick_labels = cell(1, nthr);
    for i = 1:nthr; xtick_labels{i} = xtick_fn(i); end
    set(gca, 'XTick', 1:nthr, 'XTickLabel', xtick_labels);
    if nthr > 12; set(gca, 'XTickLabelRotation', 45); end
    xlabel('Thruster'); ylabel('Maneuver');
    title(sprintf('%s Thruster Commands (min thr = %.0f%%)', lbl, MIN_THROTTLE*100), 'FontSize', 11);
    % Separator line for combined
    if sc_idx == 3
        hold on;
        plot([N_ast+0.5, N_ast+0.5], [0.5, 12.5], 'c-', 'LineWidth', 2);
    end
    % Value labels
    for m = 1:12
        for i = 1:nthr
            val = cmd_mat(m,i);
            if val > 0.01
                if val > 0.6; clr = 'k'; else; clr = 'w'; end
                text(i, m, sprintf('%.0f%%', val*100), 'HorizontalAlignment', 'center', ...
                    'FontSize', max(6, 8 - floor(nthr/10)), 'Color', clr, 'FontWeight', 'bold');
            end
        end
    end
end

% ---- +/- Authority bars ----
fig5 = figure('Name', 'Authority Comparison', 'Position', [20 20 fig_w*2 350]);
for sc_idx = 1:3
    if sc_idx == 1;     res = ast_results;  lbl = ast.name;
    elseif sc_idx == 2; res = obx_results;  lbl = obx.name;
    else;               res = comb_results; lbl = 'Combined';
    end
    auth_pos = zeros(1,6); auth_neg = zeros(1,6);
    for d = 1:6
        auth_pos(d) = abs(res(2*d-1).w_actual(d));
        auth_neg(d) = abs(res(2*d).w_actual(d));
    end
    subplot(1,3,sc_idx);
    bar_data = [auth_pos; auth_neg]';
    b = bar(bar_data, 'grouped');
    b(1).FaceColor = [0.25 0.45 0.75]; b(2).FaceColor = [0.8 0.25 0.25];
    set(gca, 'XTickLabel', dof_labels); ylabel('[N] or [N*m]');
    title(sprintf('%s: + vs -', lbl), 'FontSize', 11);
    legend('+', '-', 'Location', 'best'); grid on;
end

% ---- Parasitic wrench comparison bar ----
fig6 = figure('Name', 'Parasitic Wrench', 'Position', [20 20 700 300]);
para_a = zeros(1,12); para_o = zeros(1,12); para_c = zeros(1,12);
for m = 1:12
    para_a(m) = ast_results(m).parasitic_norm;
    para_o(m) = obx_results(m).parasitic_norm;
    para_c(m) = comb_results(m).parasitic_norm;
end
bar_data_p = [para_a; para_o; para_c]';
b = bar(bar_data_p, 'grouped');
b(1).FaceColor = [0.5 0.6 0.85]; b(2).FaceColor = [0.5 0.8 0.6]; b(3).FaceColor = [0.7 0.6 0.8];
set(gca, 'XTick', 1:12, 'XTickLabel', short_names, 'XTickLabelRotation', 45);
ylabel('||Parasitic wrench||');
title(sprintf('Parasitic Wrench per Maneuver (min thr = %.0f%%)', MIN_THROTTLE*100), 'FontSize', 11);
legend(ast.name, obx.name, 'Combined', 'Location', 'best'); grid on;

% ---- Before/After parasitic comparison bar chart ----
fig7 = figure('Name', 'Parasitic Before/After', 'Position', [20 20 1000 350]);
for sc_idx = 1:3
    if sc_idx == 1;     corr = ast_corrected;  lbl = ast.name;
    elseif sc_idx == 2; corr = obx_corrected;  lbl = obx.name;
    else;               corr = comb_corrected;  lbl = 'Combined';
    end
    para_bef = zeros(1,12); para_aft = zeros(1,12);
    for m = 1:12
        para_bef(m) = corr(m).parasitic_initial;
        para_aft(m) = corr(m).parasitic_norm;
    end
    subplot(1,3,sc_idx);
    bar_ba = [para_bef; para_aft]';
    b = bar(bar_ba, 'grouped');
    b(1).FaceColor = [0.8 0.4 0.4]; b(2).FaceColor = [0.3 0.7 0.4];
    set(gca, 'XTick', 1:12, 'XTickLabel', short_names, 'XTickLabelRotation', 45);
    ylabel('||Parasitic||'); title(sprintf('%s: Before vs After', lbl), 'FontSize', 11);
    legend('Before', 'After', 'Location', 'best'); grid on;
end

fprintf('\n=== Analysis Complete ===\n');


%% ====================== SOLVER FUNCTION ================================
function results = solve_all_maneuvers(B, maneuvers, n_thr, min_throttle, max_throttle, opts, target_force,target_torque)
    MAX_ITER = 5;
    n_man = length(maneuvers);
    results = struct();

    for m = 1:n_man
        dof = maneuvers(m).dof;
        sgn = maneuvers(m).sign;
        cmd_row   = dof;
        para_rows = setdiff(1:6, cmd_row);

        C     = B(para_rows, :);
        d_vec = zeros(size(C,1), 1);
        Aeq = B(cmd_row, :);
        if dof <= 3
            beq = sgn * target_force;
        else
            beq = sgn * target_torque;
        end
        lb = zeros(n_thr, 1);
        ub = max_throttle * ones(n_thr, 1);

        if sgn > 0; feasible = any(B(cmd_row,:) > 1e-10);
        else;       feasible = any(B(cmd_row,:) < -1e-10);
        end

        if ~feasible
            results(m).cmd            = zeros(n_thr, 1);
            results(m).w_actual       = zeros(6, 1);
            results(m).parasitic      = zeros(6, 1);
            results(m).parasitic_norm = 0;
            results(m).n_active       = 0;
            results(m).solver_flag    = -99;
            results(m).min_thr_ok     = true;
            results(m).min_thr_actual = 0;
            continue;
        end

        % L2 regularization: prefer sparse (minimum-throttle) solutions
        reg_w = 1e-3;
        C     = [C; reg_w * eye(n_thr)];
        d_vec = [d_vec; zeros(n_thr, 1)];

        [u_sol, ~, ~, flag1] = lsqlin(C, d_vec, [], [], Aeq, beq, lb, ub, [], opts);

        if flag1 <= 0
            w_des = zeros(6,1); w_des(dof) = sgn;
            [u_sol, ~] = lsqnonneg(B, w_des);
            if max(u_sol) > 1e-12; u_sol = u_sol / max(u_sol); end
        end

        min_thr_enforced = true;
        if min_throttle > 0 && flag1 > 0
            u_best = u_sol;
            for iter = 1:MAX_ITER
                in_dead_zone = (u_sol > 1e-6) & (u_sol < min_throttle);
                if ~any(in_dead_zone); break; end
                active = (u_sol >= min_throttle);
                lb_c = zeros(n_thr, 1); ub_c = max_throttle * ones(n_thr, 1);
                lb_c(active) = min_throttle; lb_c(~active) = 0; ub_c(~active) = 0;
                [u_new, ~, ~, flag_c] = lsqlin(C, d_vec, [], [], Aeq, beq, lb_c, ub_c, [], opts);
                if flag_c > 0; u_sol = u_new;
                else; u_sol = u_best; min_thr_enforced = false; break;
                end
            end
            still_bad = (u_sol > 1e-6) & (u_sol < min_throttle);
            if any(still_bad); min_thr_enforced = false; end
        end

        w_actual = B * u_sol;
        w_parasitic = w_actual; w_parasitic(cmd_row) = 0;
        nonzero_cmds = u_sol(u_sol > 1e-6);
        if isempty(nonzero_cmds); actual_min = 0; else; actual_min = min(nonzero_cmds); end

        results(m).cmd            = u_sol;
        results(m).w_actual       = w_actual;
        results(m).parasitic      = w_parasitic;
        results(m).parasitic_norm = norm(w_parasitic);
        results(m).n_active       = sum(u_sol > 1e-6);
        results(m).solver_flag    = flag1;
        results(m).min_thr_ok     = min_thr_enforced;
        results(m).min_thr_actual = actual_min;
    end
end

%% ====================== DISPLAY FUNCTIONS ==============================
function print_maneuver_results(results, maneuvers, n_thr, mode_label, dof_labels)
    for m = 1:length(maneuvers)
        r = results(m);
        fprintf('  --- %s: %s ---\n', mode_label, maneuvers(m).name);
        fprintf('  Active thrusters: %d/%d', r.n_active, n_thr);
        if r.solver_flag <= 0; fprintf('  [solver flag = %d]', r.solver_flag); end
        if ~r.min_thr_ok; fprintf('  [MIN THROTTLE NOT MET - lowest on = %.0f%%]', r.min_thr_actual*100); end
        fprintf('\n');
        fprintf('  Firing: ');
        any_firing = false;
        for i = 1:n_thr
            if r.cmd(i) > 0.01; fprintf('T%d(%3.0f%%) ', i, r.cmd(i)*100); any_firing = true; end
        end
        if ~any_firing; fprintf('(none)'); end
        fprintf('\n');
        fprintf('  Wrench: ');
        for d = 1:6
            if abs(r.w_actual(d)) > 1e-6; fprintf('%s=%+.4f  ', dof_labels{d}, r.w_actual(d)); end
        end
        fprintf('\n');
        if r.parasitic_norm > 0.001
            fprintf('  Parasitic: ');
            for d = 1:6
                if abs(r.parasitic(d)) > 0.001; fprintf('%s=%+.4f  ', dof_labels{d}, r.parasitic(d)); end
            end
            fprintf('  (||%.4f||)\n', r.parasitic_norm);
        else
            fprintf('  Clean - no parasitic coupling\n');
        end
        fprintf('\n');
    end
end

function print_min_throttle_violations(results, maneuvers)
    any_viol = any(~[results.min_thr_ok]);
    if any_viol
        fprintf('\n  * = min throttle could NOT be enforced.\n');
        fprintf('      Violating maneuvers:\n');
        for m = 1:12
            if ~results(m).min_thr_ok
                fprintf('        %-14s  lowest on = %4.0f%%\n', ...
                    maneuvers(m).name, results(m).min_thr_actual*100);
            end
        end
    end
    fprintf('\n');
end

function draw_box(ax, Lx, Ly, Lz, center, fc, ec)
    a = Lx/2; b = Ly/2; c = Lz/2;
    v = center + [-a -b -c; a -b -c; a b -c; -a b -c;
                  -a -b  c; a -b  c; a b  c; -a b  c];
    f = [1 2 3 4; 5 6 7 8; 1 2 6 5; 3 4 8 7; 1 4 8 5; 2 3 7 6];
    patch(ax, 'Vertices', v, 'Faces', f, 'FaceColor', fc, 'FaceAlpha', 0.15, ...
          'EdgeColor', ec, 'LineWidth', 1.2);
end

function draw_thrusters_labeled(ax, pos, dir, scale, clr)
    for i = 1:size(pos,1)
        p = pos(i,:); d = dir(i,:);
        quiver3(ax, p(1),p(2),p(3), d(1)*scale, d(2)*scale, d(3)*scale, 0, ...
                'Color', clr, 'LineWidth', 2.2, 'MaxHeadSize', 0.45);
        plot3(ax, p(1), p(2), p(3), 'o', 'MarkerSize', 7, ...
              'MarkerFaceColor', clr, 'Color', clr);
        text(ax, p(1)+d(1)*scale*1.2, p(2)+d(2)*scale*1.2, p(3)+d(3)*scale*1.2, ...
            sprintf('T%d', i), 'FontSize', 12, 'Color', 'k', 'FontWeight', 'bold');
    end
end

function s = bool_str(b)
    if b; s = 'yes'; else; s = 'NO'; end
end

function s = fmt_thr(val)
    if val < 0.01;                    s = '   . ';
    elseif val > 0.99 && val < 1.01;  s = '   X ';
    else;                              s = sprintf(' %3.0f%%', val*100);
    end
end

function r = safe_ratio(a, b)
    if b > 1e-10; r = a/b; else; r = Inf; end
end

function lbl = comb_thr_label(i, N_ast)
    if i <= N_ast; lbl = sprintf('A.T%d', i);
    else;          lbl = sprintf('O.T%d', i - N_ast);
    end
end

function results = solve_corrected(B, maneuvers, n_thr, min_throttle, max_throttle, tol, max_iter, opts,target_force,target_torque)
    corr_opts = opts;
    n_man = length(maneuvers);
    results = struct();

    for m = 1:n_man
        dof = maneuvers(m).dof;
        sgn = maneuvers(m).sign;
        cmd_row   = dof;
        para_rows = setdiff(1:6, cmd_row);

        C_prim = B(para_rows, :);
        d_prim = zeros(size(C_prim,1), 1);
        Aeq = B(cmd_row, :);
        if dof <= 3
            beq = sgn * target_force;
        else
            beq = sgn * target_torque;
        end
        lb = zeros(n_thr, 1); ub = max_throttle * ones(n_thr, 1);

        if sgn > 0; feasible = any(B(cmd_row,:) > 1e-10);
        else;       feasible = any(B(cmd_row,:) < -1e-10);
        end

        if ~feasible
            results(m).cmd              = zeros(n_thr, 1);
            results(m).w_actual         = zeros(6, 1);
            results(m).parasitic        = zeros(6, 1);
            results(m).parasitic_norm   = 0;
            results(m).parasitic_initial = 0;
            results(m).n_active         = 0;
            results(m).n_iters          = 0;
            results(m).converged        = true;
            results(m).min_thr_ok       = true;
            results(m).min_thr_actual   = 0;
            continue;
        end

        reg_w  = 1e-3;
        C_prim = [C_prim; reg_w * eye(n_thr)];
        d_prim = [d_prim; zeros(n_thr, 1)];

        [u_primary, ~, ~, flag1] = lsqlin(C_prim, d_prim, [], [], Aeq, beq, lb, ub, [], opts);
        if flag1 <= 0
            w_des = zeros(6,1); w_des(dof) = sgn;
            [u_primary, ~] = lsqnonneg(B, w_des);
            if max(u_primary) > 1e-12; u_primary = u_primary / max(u_primary); end
        end

        u_total = u_primary;
        w_actual = B * u_total;
        w_para = w_actual; w_para(cmd_row) = 0;
        parasitic_initial = norm(w_para);

        n_iters = 0;
        converged = (norm(w_para) < tol);

        min_thr_enforced = true;
        for iter = 1:max_iter
            if converged; break; end

            in_dead = (u_total > 1e-6) & (u_total < min_throttle);
            if ~any(in_dead) && n_iters > 0; break; end

            lb_iter = zeros(n_thr, 1);
            ub_iter = max_throttle * ones(n_thr, 1);
            if min_throttle > 0
                active = (u_total >= min_throttle);
                lb_iter(active) = min_throttle;
                ub_iter(in_dead) = 0;
                lb_iter(in_dead) = 0;
            end

            [u_new, ~, ~, flag_iter] = lsqlin(C_prim, d_prim, [], [], ...
                                                Aeq, beq, lb_iter, ub_iter, [], corr_opts);
            if flag_iter <= 0; break; end

            w_new = B * u_new;
            w_para_new = w_new; w_para_new(cmd_row) = 0;

            if n_iters > 0 && norm(w_para_new) >= norm(w_para) - 1e-8
                break;
            end

            u_total = u_new; w_actual = w_new; w_para = w_para_new;
            n_iters = n_iters + 1;
            converged = (norm(w_para) < tol);
        end

        if min_throttle > 0
            still_bad = (u_total > 1e-6) & (u_total < min_throttle);
            if any(still_bad); min_thr_enforced = false; end
        end

        w_actual = B * u_total;
        w_parasitic = w_actual; w_parasitic(cmd_row) = 0;
        nonzero_cmds = u_total(u_total > 1e-6);
        if isempty(nonzero_cmds); actual_min = 0; else; actual_min = min(nonzero_cmds); end

        results(m).cmd              = u_total;
        results(m).w_actual         = w_actual;
        results(m).parasitic        = w_parasitic;
        results(m).parasitic_norm   = norm(w_parasitic);
        results(m).parasitic_initial = parasitic_initial;
        results(m).n_active         = sum(u_total > 1e-6);
        results(m).n_iters          = n_iters;
        results(m).converged        = converged;
        results(m).min_thr_ok       = min_thr_enforced;
        results(m).min_thr_actual   = actual_min;
    end
end

function plot_maneuver_3d(maneuver_idx, mode, result, maneuver, ...
    thr_pos_a, thr_dir_a, N_ast_thr, ast_sc, ...
    thr_pos_o_global, thr_dir_o, mating_plane_x, max_throttle)
%PLOT_MANEUVER_3D  3D visualization of thruster firing for a single maneuver.
    dof = maneuver.dof;
    sgn = maneuver.sign;

    if contains(mode, 'ast')
        pos_all = thr_pos_a; dir_all = thr_dir_a;
        com = ast_sc.com; n_total = N_ast_thr;
        mode_str = [ast_sc.name ' Solo'];
        Lx = ast_sc.Lx; Ly = ast_sc.Ly; Lz = ast_sc.Lz;
    else
        pos_all = [thr_pos_a; thr_pos_o_global];
        dir_all = [thr_dir_a; thr_dir_o];
        com = [mating_plane_x, 0, 0];
        n_total = size(pos_all, 1);
        mode_str = 'Combined';
        Lx = ast_sc.Lx; Ly = ast_sc.Ly; Lz = ast_sc.Lz;
    end

    cmd = result.cmd;
    arrow_scale = max([Lx, Ly, Lz])/2 * 0.35;

    fig = figure('Name', sprintf('%s %s', mode_str, maneuver.name), ...
        'Position', [50 + mod(maneuver_idx-1, 3)*380, ...
                     400 - floor((maneuver_idx-1)/3)*350, 370, 330]);
    ax = axes(fig); hold(ax, 'on');

    draw_box_fn(ax, Lx, Ly, Lz, [0 0 0], [0.7 0.8 0.95], [0.4 0.4 0.6]);

    for i = 1:n_total
        p = pos_all(i,:); d = dir_all(i,:);
        u_i = cmd(i);
        if u_i < 0.01
            plot3(ax, p(1), p(2), p(3), 'o', 'MarkerSize', 4, ...
                'MarkerFaceColor', [0.7 0.7 0.7], 'Color', [0.5 0.5 0.5]);
            text(ax, p(1), p(2), p(3)+0.06, sprintf('T%d', i), ...
                'FontSize', 6, 'Color', [0.6 0.6 0.6], 'HorizontalAlignment', 'center');
        else
            frac = min(u_i / max_throttle, 1.0);
            clr = throttle_color(frac);
            len = arrow_scale * (0.3 + 0.7 * u_i);
            quiver3(ax, p(1), p(2), p(3), d(1)*len, d(2)*len, d(3)*len, 0, ...
                'Color', clr, 'LineWidth', 2.5, 'MaxHeadSize', 0.4);
            plot3(ax, p(1), p(2), p(3), 'o', 'MarkerSize', 7, ...
                'MarkerFaceColor', clr, 'Color', clr*0.7);
            lbl_pos = p + d * len * 1.15;
            text(ax, lbl_pos(1), lbl_pos(2), lbl_pos(3), ...
                sprintf('T%d\n%.0f%%', i, u_i*100), ...
                'FontSize', 7, 'Color', clr*0.6, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center');
        end
    end

    plot3(ax, com(1), com(2), com(3), 'k+', 'MarkerSize', 10, 'LineWidth', 2);

    des_scale = arrow_scale * 1.5;
    if dof <= 3
        maneuver_arrow = zeros(1,3); maneuver_arrow(dof) = sgn;
        quiver3(ax, com(1), com(2), com(3), ...
            maneuver_arrow(1)*des_scale, maneuver_arrow(2)*des_scale, maneuver_arrow(3)*des_scale, 0, ...
            'Color', [0.0 0.7 0.0], 'LineWidth', 3.5, 'MaxHeadSize', 0.3);
    else
        draw_rotation_arrow(ax, com, dof-3, sgn, des_scale*0.6, [0.0 0.7 0.0]);
    end

    w_para = result.parasitic;
    para_force = w_para(1:3); para_torque = w_para(4:6);
    para_scale = arrow_scale * 0.8;
    if norm(para_force) > 0.01
        pf_dir = para_force / max(norm(para_force), 1e-10);
        pf_len = para_scale * norm(para_force);
        quiver3(ax, com(1), com(2), com(3), pf_dir(1)*pf_len, pf_dir(2)*pf_len, pf_dir(3)*pf_len, 0, ...
            'Color', [0.8 0.2 0.2], 'LineWidth', 1.5, 'LineStyle', '--', 'MaxHeadSize', 0.3);
    end
    if norm(para_torque) > 0.01
        pt_dir = para_torque / max(norm(para_torque), 1e-10);
        pt_len = para_scale * norm(para_torque) * 0.5;
        quiver3(ax, com(1), com(2), com(3), pt_dir(1)*pt_len, pt_dir(2)*pt_len, pt_dir(3)*pt_len, 0, ...
            'Color', [0.8 0.4 0.1], 'LineWidth', 1.5, 'LineStyle', ':', 'MaxHeadSize', 0.3);
    end

    title(ax, sprintf('%s: %s\nParasitic: %.4f → %.4f  (%d iters)', ...
        mode_str, maneuver.name, result.parasitic_initial, result.parasitic_norm, result.n_iters), ...
        'FontSize', 10);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis(ax, 'equal'); grid(ax, 'on'); view(ax, 135, 25);
end

function clr = throttle_color(frac)
    if frac < 0.5
        t = frac * 2;
        clr = [t, t, 1-t];
    else
        t = (frac - 0.5) * 2;
        clr = [1, 1-t, 0];
    end
end

function draw_rotation_arrow(ax, center, rot_axis, rot_sign, radius, clr)
    n_pts = 20;
    theta = linspace(0, rot_sign * pi * 0.75, n_pts);
    pts = zeros(n_pts, 3);
    for k = 1:n_pts
        switch rot_axis
            case 1; pts(k,:) = center + radius * [0, cos(theta(k)), sin(theta(k))];
            case 2; pts(k,:) = center + radius * [cos(theta(k)), 0, sin(theta(k))];
            case 3; pts(k,:) = center + radius * [cos(theta(k)), sin(theta(k)), 0];
        end
    end
    plot3(ax, pts(:,1), pts(:,2), pts(:,3), '-', 'Color', clr, 'LineWidth', 3);
    dp = pts(end,:) - pts(end-1,:);
    dp = dp / norm(dp) * radius * 0.25;
    quiver3(ax, pts(end,1), pts(end,2), pts(end,3), dp(1), dp(2), dp(3), 0, ...
        'Color', clr, 'LineWidth', 3, 'MaxHeadSize', 1.5);
end

function draw_box_fn(ax, Lx, Ly, Lz, center, fc, ec)
    a = Lx/2; b = Ly/2; c = Lz/2;
    v = center + [-a -b -c; a -b -c; a b -c; -a b -c;
                  -a -b  c; a -b  c; a b  c; -a b  c];
    f = [1 2 3 4; 5 6 7 8; 1 2 6 5; 3 4 8 7; 1 4 8 5; 2 3 7 6];
    patch(ax, 'Vertices', v, 'Faces', f, 'FaceColor', fc, 'FaceAlpha', 0.12, ...
          'EdgeColor', ec, 'LineWidth', 1.0);
end