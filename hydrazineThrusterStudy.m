%% RCS Vertex Thruster Placement Trade Study for OSCAR (Asterix & Obelix)
%  INDEPENDENT THRUSTER COUNTS - Asterix and Obelix can have different
%  numbers of thrusters.  Solo configs are enumerated independently, then
%  the top K from each are paired for combined evaluation.
%
%  CONTROL REQUIREMENTS (configurable):
%    Asterix solo  : ctrl_mode + configurable pure control axes
%    Obelix solo   : ctrl_mode + configurable pure control axes
%    Combined      : full 6-DOF + configurable pure control axes
%
%  Yaman Saran - OSCAR@VT
%  -----------------------------------------------------------------------
clear; clc; close all;

%% ====================== ASTERIX INPUTS =================================
ast.name = 'Asterix';
ast.Lx   = 2.2;      ast.Ly = 1.1;      ast.Lz = 1.1;
ast.mass  = 382;
ast.Ixx   = [];  ast.Iyy = [];  ast.Izz = [];
ast.com_offset = [0.8, 0, 0];

%% ====================== OBELIX INPUTS ==================================
obx.name = 'Obelix';
obx.Lx   = 2.2;      obx.Ly = 1.1;      obx.Lz = 1.1;
obx.mass  = 581;
obx.Ixx   = [];  obx.Iyy = [];  obx.Izz = [];
obx.com_offset = [-0.8, 0, 0];

%% ====================== THRUSTER INPUTS ================================
thr.thrust = 1.1;     thr.Isp = 235;
thr.mass   = 0.330;   thr.length = 0.15;
N_THR_AST  = 11;       % Number of thrusters on Asterix
N_THR_OBX  = 10;       % Number of thrusters on Obelix

%% ====================== CONTROL MODE & PURE CONTROL REQUIREMENTS =======
%  ctrl_mode: 'full6dof' -> requires rank 6 of full B (6xN)
%             'attitude' -> requires rank 3 of rotational sub-B (rows 4:6)
%
%  pure_req: Nx2 matrix specifying pure-control (lsqnonneg u>=0) requirements.
%    Column 1: DOF index (1=Fx, 2=Fy, 3=Fz, 4=Mx, 5=My, 6=Mz)
%    Column 2: sign requirement (+1, -1, or 0 = both +/-)
%    Use [] for no pure-control requirements.
%
%  For 'attitude' mode, only DOFs 4-6 (Mx,My,Mz) are valid in pure_req.
%  Bidirectionality in all 3 rotation axes is always enforced in attitude
%  mode as a fundamental controllability requirement.
%
%  EXAMPLES:
%    [1, 0]                   -> pure +/-Fx only
%    [1,+1; 2,0; 4,-1]       -> pure +Fx, pure +/-Fy, pure -Mx
%    [(1:6)', zeros(6,1)]     -> all 6 axes, both signs (12 wrenches)
%    []                       -> no pure-control requirements
% -----------------------------------------------------------------------

ast_ctrl_mode = 'full6dof';
ast_pure_req  = [1, 0];                    % pure +/-Fx

obx_ctrl_mode = 'full6dof';
obx_pure_req  = [1,0];                        % pure +/-Fx

comb_pure_req = [(1:6)', zeros(6,1)];      % all 6 axes, both signs

%% ====================== OPTIMIZATION WEIGHTS ===========================
w_solo     = 1.0;     % weight on each spacecraft's solo sigma_min
w_combined = 1.0;     % weight on combined sigma_min
w_mating   = 0.8;     % penalty per mating-face thruster

K_TOP      = 300;     % top solo configs from each SC to pair for combined eval
N_PURE_TOP = 200;     % how many top-scored combined pairs to LP-check for pure ctrl

%% ====================== DERIVED GEOMETRY: ASTERIX ======================
ast.a = ast.Lx/2;  ast.b = ast.Ly/2;  ast.c = ast.Lz/2;
ast.vertices = [ ast.a  ast.b  ast.c;
                 ast.a -ast.b  ast.c;
                -ast.a -ast.b  ast.c;
                -ast.a  ast.b  ast.c;
                 ast.a  ast.b -ast.c;
                 ast.a -ast.b -ast.c;
                -ast.a -ast.b -ast.c;
                -ast.a  ast.b -ast.c];
ast.mating_verts = [1, 2, 5, 6];   % +x face
ast.outer_verts  = [3, 4, 7, 8];
ast.com = ast.com_offset;

if isempty(ast.Ixx); ast.Ixx = (1/12)*ast.mass*(ast.Ly^2 + ast.Lz^2); end
if isempty(ast.Iyy); ast.Iyy = (1/12)*ast.mass*(ast.Lx^2 + ast.Lz^2); end
if isempty(ast.Izz); ast.Izz = (1/12)*ast.mass*(ast.Lx^2 + ast.Lz^2); end
ast.Ixx_com = ast.Ixx + ast.mass*(ast.com_offset(2)^2 + ast.com_offset(3)^2);
ast.Iyy_com = ast.Iyy + ast.mass*(ast.com_offset(1)^2 + ast.com_offset(3)^2);
ast.Izz_com = ast.Izz + ast.mass*(ast.com_offset(1)^2 + ast.com_offset(2)^2);
ast.I = diag([ast.Ixx_com, ast.Iyy_com, ast.Izz_com]);

%% ====================== DERIVED GEOMETRY: OBELIX =======================
obx.a = obx.Lx/2;  obx.b = obx.Ly/2;  obx.c = obx.Lz/2;
obx.geom_center_global = [ast.a + obx.a, 0, 0];   % Obelix sits beside Asterix
obx.com_global = obx.geom_center_global + obx.com_offset;
obx.vertices_local = [ obx.a  obx.b  obx.c;
                        obx.a -obx.b  obx.c;
                       -obx.a -obx.b  obx.c;
                       -obx.a  obx.b  obx.c;
                        obx.a  obx.b -obx.c;
                        obx.a -obx.b -obx.c;
                       -obx.a -obx.b -obx.c;
                       -obx.a  obx.b -obx.c];
obx.mating_verts = [3, 4, 7, 8];   % -x face (in Obelix body frame)
obx.outer_verts  = [1, 2, 5, 6];
obx.vertices_global = obx.vertices_local + obx.geom_center_global;

if isempty(obx.Ixx); obx.Ixx = (1/12)*obx.mass*(obx.Ly^2 + obx.Lz^2); end
if isempty(obx.Iyy); obx.Iyy = (1/12)*obx.mass*(obx.Lx^2 + obx.Lz^2); end
if isempty(obx.Izz); obx.Izz = (1/12)*obx.mass*(obx.Lx^2 + obx.Lz^2); end
obx.Ixx_com = obx.Ixx + obx.mass*(obx.com_offset(2)^2 + obx.com_offset(3)^2);
obx.Iyy_com = obx.Iyy + obx.mass*(obx.com_offset(1)^2 + obx.com_offset(3)^2);
obx.Izz_com = obx.Izz + obx.mass*(obx.com_offset(1)^2 + obx.com_offset(2)^2);
obx.I = diag([obx.Ixx_com, obx.Iyy_com, obx.Izz_com]);

%% ====================== COMBINED VEHICLE ===============================
combined_mass = ast.mass + obx.mass + (N_THR_AST + N_THR_OBX)*thr.mass;
combined_com  = (ast.mass * ast.com + obx.mass * obx.com_global) / (ast.mass + obx.mass);

d_ast = combined_com - ast.com;
d_obx = combined_com - obx.com_global;
I_comb_xx = ast.Ixx_com + ast.mass*(d_ast(2)^2 + d_ast(3)^2) ...
          + obx.Ixx_com + obx.mass*(d_obx(2)^2 + d_obx(3)^2);
I_comb_yy = ast.Iyy_com + ast.mass*(d_ast(1)^2 + d_ast(3)^2) ...
          + obx.Iyy_com + obx.mass*(d_obx(1)^2 + d_obx(3)^2);
I_comb_zz = ast.Izz_com + ast.mass*(d_ast(1)^2 + d_ast(3)^2) ...
          + obx.Izz_com + obx.mass*(d_obx(1)^2 + d_obx(3)^2);
I_combined = diag([I_comb_xx, I_comb_yy, I_comb_zz]);

%% ====================== EXPAND & VALIDATE PURE REQUIREMENTS ===========
dof_labels_all = {'Fx','Fy','Fz','Mx','My','Mz'};

ast_pure_exp  = expand_pure_req(ast_pure_req);
obx_pure_exp  = expand_pure_req(obx_pure_req);
comb_pure_exp = expand_pure_req(comb_pure_req);

% Validate: attitude mode can only require DOFs 4-6
if strcmp(ast_ctrl_mode, 'attitude') && ~isempty(ast_pure_exp)
    assert(all(ast_pure_exp(:,1) >= 4 & ast_pure_exp(:,1) <= 6), ...
        'Asterix attitude mode: pure_req DOFs must be 4-6 (Mx,My,Mz)');
end
if strcmp(obx_ctrl_mode, 'attitude') && ~isempty(obx_pure_exp)
    assert(all(obx_pure_exp(:,1) >= 4 & obx_pure_exp(:,1) <= 6), ...
        'Obelix attitude mode: pure_req DOFs must be 4-6 (Mx,My,Mz)');
end

fprintf('=== OSCAR Dual-Spacecraft RCS Trade Study (Configurable Pure Control) ===\n');
fprintf('\n--- Control Requirements ---\n');
fprintf('  Asterix:  mode=%s', ast_ctrl_mode);
if isempty(ast_pure_exp); fprintf(',  pure_req: NONE\n');
else; fprintf(',  pure_req: %s\n', format_pure_req(ast_pure_exp, dof_labels_all)); end
fprintf('  Obelix:   mode=%s', obx_ctrl_mode);
if isempty(obx_pure_exp); fprintf(',  pure_req: NONE\n');
else; fprintf(',  pure_req: %s\n', format_pure_req(obx_pure_exp, dof_labels_all)); end
fprintf('  Combined: full 6-DOF');
if isempty(comb_pure_exp); fprintf(',  pure_req: NONE\n');
else; fprintf(',  pure_req: %s\n', format_pure_req(comb_pure_exp, dof_labels_all)); end

fprintf('\n--- %s ---\n', ast.name);
fprintf('  Dimensions: %.2f x %.2f x %.2f m,  mass = %.1f kg,  N_thr = %d\n', ast.Lx, ast.Ly, ast.Lz, ast.mass, N_THR_AST);
fprintf('  COM in body frame: (%.3f, %.3f, %.3f) m\n', ast.com);
fprintf('  Inertia (about COM): Ixx=%.2f  Iyy=%.2f  Izz=%.2f  [kg*m^2]\n', ast.Ixx_com, ast.Iyy_com, ast.Izz_com);
fprintf('\n--- %s ---\n', obx.name);
fprintf('  Dimensions: %.2f x %.2f x %.2f m,  mass = %.1f kg,  N_thr = %d\n', obx.Lx, obx.Ly, obx.Lz, obx.mass, N_THR_OBX);
fprintf('  COM in global frame: (%.3f, %.3f, %.3f) m\n', obx.com_global);
fprintf('  Inertia (about COM): Ixx=%.2f  Iyy=%.2f  Izz=%.2f  [kg*m^2]\n', obx.Ixx_com, obx.Iyy_com, obx.Izz_com);
fprintf('\n--- Combined ---\n');
fprintf('  COM: (%.3f, %.3f, %.3f) m\n', combined_com);
fprintf('  Inertia: Ixx=%.2f  Iyy=%.2f  Izz=%.2f  [kg*m^2]\n', I_comb_xx, I_comb_yy, I_comb_zz);
fprintf('\n--- Thrusters ---\n');
fprintf('  F=%.2f N,  Isp=%.0f s,  m=%.3f kg\n\n', thr.thrust, thr.Isp, thr.mass);

%% ========== PRECOMPUTE MASTER B MATRICES ===============================
nv = 8;
N_PLACEMENTS = nv * 3;   % 24

master_pos  = zeros(N_PLACEMENTS, 3);
master_dir  = zeros(N_PLACEMENTS, 3);
master_vtx  = zeros(N_PLACEMENTS, 1);
master_didx = zeros(N_PLACEMENTS, 1);

for v = 1:nv
    vtx = ast.vertices(v,:);
    sx = sign(vtx(1));  sy = sign(vtx(2));  sz = sign(vtx(3));
    dirs_v = [-sx 0 0; 0 -sy 0; 0 0 -sz];
    for d = 1:3
        p = 3*(v-1) + d;
        master_pos(p,:)  = vtx;
        master_dir(p,:)  = dirs_v(d,:);
        master_vtx(p)    = v;
        master_didx(p)   = d;
    end
end

% --- Asterix solo B (about Asterix COM) ---
r_ast_solo = master_pos - ast.com;
F_ast_solo = thr.thrust * master_dir;
M_ast_solo = cross(r_ast_solo, F_ast_solo, 2);
B_master_ast_solo = [F_ast_solo'; M_ast_solo'];

% --- Obelix solo B (about Obelix COM in body frame) ---
r_obx_solo = master_pos - obx.com_offset;
F_obx_solo = thr.thrust * master_dir;
M_obx_solo = cross(r_obx_solo, F_obx_solo, 2);
B_master_obx_solo = [F_obx_solo'; M_obx_solo'];

% --- Combined B: Asterix placements about combined COM ---
r_comb_ast = master_pos - combined_com;
F_comb_ast = thr.thrust * master_dir;
M_comb_ast = cross(r_comb_ast, F_comb_ast, 2);
B_master_comb_ast = [F_comb_ast'; M_comb_ast'];

% --- Combined B: Obelix placements in global frame about combined COM ---
master_pos_obx_global = master_pos + obx.geom_center_global;
r_comb_obx = master_pos_obx_global - combined_com;
F_comb_obx = thr.thrust * master_dir;
M_comb_obx = cross(r_comb_obx, F_comb_obx, 2);
B_master_comb_obx = [F_comb_obx'; M_comb_obx'];

B_master_comb = [B_master_comb_ast, B_master_comb_obx];  % 6x48

% Mating penalty lookup per placement
mate_penalty_ast = zeros(N_PLACEMENTS, 1);
mate_penalty_obx = zeros(N_PLACEMENTS, 1);
for p = 1:N_PLACEMENTS
    v = master_vtx(p);
    if ismember(v, ast.mating_verts)
        if abs(master_dir(p,1)) > 1e-10; mate_penalty_ast(p) = 1.0;
        else;                              mate_penalty_ast(p) = 0.5; end
    end
    if ismember(v, obx.mating_verts)
        if abs(master_dir(p,1)) > 1e-10; mate_penalty_obx(p) = 1.0;
        else;                              mate_penalty_obx(p) = 0.5; end
    end
end

fprintf('Precomputed master B matrices: 24 cols each\n');

%% ====================== ENUMERATION HELPER =============================
subsets_for_k = cell(4,1);
for k = 0:3
    masks = [];
    for code = 0:7
        bits = bitget(code, 1:3);
        if sum(bits) == k
            masks = [masks; bits]; %#ok<AGROW>
        end
    end
    subsets_for_k{k+1} = masks;
end

%% ====================== ENUMERATE & EVALUATE ASTERIX SOLO ==============
fprintf('\n===== Enumerating %s solo configs (%d thrusters, %s) =====\n', ...
    ast.name, N_THR_AST, upper(ast_ctrl_mode));
[ast_configs, ast_solo_sigma, ast_solo_rank, ast_mate_pen, ast_n_mate, ast_n_outer, ast_mult, ast_pure_flag] = ...
    enumerate_and_evaluate_solo(N_THR_AST, nv, subsets_for_k, B_master_ast_solo, ...
                                 mate_penalty_ast, ast.mating_verts, ast.outer_verts, ...
                                 ast_ctrl_mode, ast_pure_exp);

if strcmp(ast_ctrl_mode, 'full6dof')
    req_rank = 6;
else
    req_rank = 3;
end
n_ast_ctrl = sum(ast_solo_rank == req_rank & ast_solo_sigma > 0);
n_ast_pure = sum(ast_pure_flag);
fprintf('  Total configs: %d,  Controllable (rank %d): %d (%.1f%%)\n', ...
    size(ast_configs,1), req_rank, n_ast_ctrl, 100*n_ast_ctrl/size(ast_configs,1));
if ~isempty(ast_pure_exp)
    fprintf('  With pure-ctrl requirements met: %d (%.1f%% of controllable)\n', ...
        n_ast_pure, 100*n_ast_pure/max(n_ast_ctrl,1));
end

%% ====================== ENUMERATE & EVALUATE OBELIX SOLO ===============
fprintf('\n===== Enumerating %s solo configs (%d thrusters, %s) =====\n', ...
    obx.name, N_THR_OBX, upper(obx_ctrl_mode));
[obx_configs, obx_solo_sigma, obx_solo_rank, obx_mate_pen, obx_n_mate, obx_n_outer, obx_mult, obx_pure_flag] = ...
    enumerate_and_evaluate_solo(N_THR_OBX, nv, subsets_for_k, B_master_obx_solo, ...
                                 mate_penalty_obx, obx.mating_verts, obx.outer_verts, ...
                                 obx_ctrl_mode, obx_pure_exp);

if strcmp(obx_ctrl_mode, 'full6dof')
    req_rank_obx = 6;
else
    req_rank_obx = 3;
end
n_obx_ctrl = sum(obx_solo_rank == req_rank_obx & obx_solo_sigma > 0);
n_obx_pure = sum(obx_pure_flag);
fprintf('  Total configs: %d,  Controllable (rank %d): %d (%.1f%%)\n', ...
    size(obx_configs,1), req_rank_obx, n_obx_ctrl, 100*n_obx_ctrl/size(obx_configs,1));
if ~isempty(obx_pure_exp)
    fprintf('  With pure-ctrl requirements met: %d (%.1f%% of controllable)\n', ...
        n_obx_pure, 100*n_obx_pure/max(n_obx_ctrl,1));
end

%% ====================== SELECT TOP-K FROM EACH FOR PAIRING ============
% Asterix: must pass pure-ctrl requirements (if any)
ast_solo_score = ast_solo_sigma - w_mating * ast_mate_pen;
if ~isempty(ast_pure_exp)
    ast_solo_score(~ast_pure_flag) = -Inf;   % REQUIRE pure-ctrl
    n_ast_eligible = n_ast_pure;
else
    ast_solo_score(ast_solo_rank ~= req_rank | ast_solo_sigma == 0) = -Inf;
    n_ast_eligible = n_ast_ctrl;
end
[~, ast_sort] = sort(ast_solo_score, 'descend');
K_ast = min(K_TOP, n_ast_eligible);
ast_top_idx = ast_sort(1:K_ast);

% Obelix: must pass pure-ctrl requirements (if any), else controllable
obx_solo_score = obx_solo_sigma - w_mating * obx_mate_pen;
if ~isempty(obx_pure_exp)
    obx_solo_score(~obx_pure_flag) = -Inf;
    n_obx_eligible = n_obx_pure;
else
    obx_solo_score(obx_solo_rank ~= req_rank_obx | obx_solo_sigma == 0) = -Inf;
    n_obx_eligible = n_obx_ctrl;
end
[~, obx_sort] = sort(obx_solo_score, 'descend');
K_obx = min(K_TOP, n_obx_eligible);
obx_top_idx = obx_sort(1:K_obx);

fprintf('\nPairing top %d %s x top %d %s = %d combined evaluations\n', ...
    K_ast, ast.name, K_obx, obx.name, K_ast * K_obx);

%% ====================== JOINT COMBINED SEARCH ==========================
% Derive combined bidirectional DOFs from comb_pure_req
comb_bidir_dofs = [];
if ~isempty(comb_pure_exp)
    comb_bidir_dofs = unique(comb_pure_exp(:,1))';
end
has_comb_pure = ~isempty(comb_pure_exp);

fprintf('Evaluating combined pairs (rank 6');
if has_comb_pure; fprintf(' + pure ctrl check'); end
fprintf(')...\n');
tic;

n_pairs = K_ast * K_obx;
pair_score    = -Inf(n_pairs, 1);
pair_sig_comb = zeros(n_pairs, 1);
pair_rank_c   = zeros(n_pairs, 1);
pair_ast_idx  = zeros(n_pairs, 1);
pair_obx_idx  = zeros(n_pairs, 1);
pair_bidir_ok = false(n_pairs, 1);

n_rank6 = 0;
n_bidir_pass = 0;
pi_count = 0;

for ia = 1:K_ast
    ai = ast_top_idx(ia);
    cols_ast = ast_configs(ai,:);
    sig_ast  = ast_solo_sigma(ai);
    pen_ast  = ast_mate_pen(ai);

    for io = 1:K_obx
        pi_count = pi_count + 1;
        oi = obx_top_idx(io);
        cols_obx = obx_configs(oi,:);
        sig_obx  = obx_solo_sigma(oi);
        pen_obx  = obx_mate_pen(oi);

        pair_ast_idx(pi_count) = ai;
        pair_obx_idx(pi_count) = oi;

        % Combined B: Asterix columns from first 24, Obelix from last 24
        col_idx_comb = [cols_ast, cols_obx + N_PLACEMENTS];
        B_comb = B_master_comb(:, col_idx_comb);
        rk_c = rank(B_comb, 1e-8);
        pair_rank_c(pi_count) = rk_c;

        if rk_c < 6; continue; end
        n_rank6 = n_rank6 + 1;

        % Bidirectional pre-filter on required DOFs
        if has_comb_pure
            bidir_ok = true;
            for bd = comb_bidir_dofs
                if ~(any(B_comb(bd,:) > 1e-10) && any(B_comb(bd,:) < -1e-10))
                    bidir_ok = false; break;
                end
            end
            pair_bidir_ok(pi_count) = bidir_ok;
            if ~bidir_ok; continue; end
        else
            pair_bidir_ok(pi_count) = true;  % no bidir requirement
        end
        n_bidir_pass = n_bidir_pass + 1;

        sv_c = svd(B_comb);
        pair_sig_comb(pi_count) = sv_c(end);

        pair_score(pi_count) = w_solo * sig_ast + w_solo * sig_obx ...
                             + w_combined * sv_c(end) ...
                             - w_mating * (pen_ast + pen_obx);
    end
end

elapsed = toc;
fprintf('Evaluated: %d pairs  (%.2f sec)\n', pi_count, elapsed);
fprintf('Combined rank-6: %d', n_rank6);
if has_comb_pure; fprintf('  |  Bidir pre-filter pass: %d', n_bidir_pass); end
fprintf('\n');

%% ====================== PURE CONTROL LP CHECK ON TOP CANDIDATES ========
pair_pure_ctrl = false(pi_count, 1);
n_checked = 0;
n_pure_ok = 0;

if has_comb_pure
    fprintf('\nChecking pure control feasibility (lsqnonneg) on top %d candidates...\n', ...
        min(N_PURE_TOP, n_bidir_pass));

    [~, pair_sort_prelim] = sort(pair_score(1:pi_count), 'descend');

    for kk = 1:min(N_PURE_TOP * 3, pi_count)
        pi = pair_sort_prelim(kk);
        if pair_score(pi) == -Inf; break; end
        n_checked = n_checked + 1;

        ai = pair_ast_idx(pi);  oi = pair_obx_idx(pi);
        col_idx_comb = [ast_configs(ai,:), obx_configs(oi,:) + N_PLACEMENTS];
        B_comb = B_master_comb(:, col_idx_comb);

        if check_pure_control_general(B_comb, comb_pure_exp)
            pair_pure_ctrl(pi) = true;
            n_pure_ok = n_pure_ok + 1;
        end

        if n_pure_ok >= N_PURE_TOP; break; end
    end
    fprintf('  Checked %d pairs, %d pass combined pure-ctrl requirements\n', n_checked, n_pure_ok);
else
    % No combined pure requirements: all rank-6 pairs are valid
    pair_pure_ctrl(1:pi_count) = pair_bidir_ok(1:pi_count) & (pair_score(1:pi_count) > -Inf);
    n_pure_ok = sum(pair_pure_ctrl(1:pi_count));
    fprintf('\nNo combined pure-ctrl requirements. %d rank-6 pairs available.\n', n_pure_ok);
end

% Sort: pure-ctrl pairs first (by score), then others
pair_sort_key = pair_score(1:pi_count);
if has_comb_pure
    pair_sort_key(pair_bidir_ok(1:pi_count) & ~pair_pure_ctrl(1:pi_count)) = ...
        pair_sort_key(pair_bidir_ok(1:pi_count) & ~pair_pure_ctrl(1:pi_count)) - 1000;
end
[~, pair_sort] = sort(pair_sort_key, 'descend');

%% ====================== DISTRIBUTION BREAKDOWN =========================
fprintf('\n========== DISTRIBUTION BREAKDOWN ==========\n');
fprintf('(Asterix mating / outer) x (Obelix mating / outer)\n');
valid_pairs = pair_pure_ctrl(1:pi_count);
if ~any(valid_pairs)
    fprintf('  (No pairs with verified pure control - showing bidir-pass pairs instead)\n');
    valid_pairs = pair_bidir_ok(1:pi_count);
end
for nma = 0:N_THR_AST
    noa = N_THR_AST - nma;
    for nmo = 0:N_THR_OBX
        noo = N_THR_OBX - nmo;
        mask_d = valid_pairs;
        for pi = find(valid_pairs)'
            ai = pair_ast_idx(pi); oi = pair_obx_idx(pi);
            if ast_n_mate(ai) ~= nma || obx_n_mate(oi) ~= nmo
                mask_d(pi) = false;
            end
        end
        if ~any(mask_d); continue; end
        idx_d = find(mask_d);
        [best_s, ~] = max(pair_score(idx_d));
        fprintf('  A(%d/%d) O(%d/%d):  %4d pairs  |  best score=%.4f\n', ...
            nma, noa, nmo, noo, sum(mask_d), best_s);
    end
end

%% ====================== TOP COMBINED CONFIGS ===========================
n_show = min(15, n_pure_ok);
if n_show == 0; n_show = min(15, n_bidir_pass); end
if n_show == 0; n_show = min(15, n_rank6); end

fprintf('\n========== TOP %d COMBINED CONFIGURATIONS ==========\n', n_show);
fprintf('%4s  %10s  %10s  %8s %8s %8s %8s  Pure  Ast_mult  Obx_mult\n', ...
    'Rank', 'A(mate/out)', 'O(mate/out)', 'sig_A', 'sig_O', 'sig_C', 'Score');
fprintf('%s\n', repmat('-', 1, 115));
for k = 1:n_show
    pi = pair_sort(k);
    ai = pair_ast_idx(pi);  oi = pair_obx_idx(pi);
    if has_comb_pure
        if pair_pure_ctrl(pi); pure_str = 'YES';
        elseif pair_bidir_ok(pi); pure_str = 'no';
        else; pure_str = '?'; end
    else
        pure_str = 'N/A';
    end
    fprintf('%4d    %d/%d          %d/%d       %8.4f %8.4f %8.4f %8.4f  %3s   %s  %s\n', ...
        k, ast_n_mate(ai), ast_n_outer(ai), obx_n_mate(oi), obx_n_outer(oi), ...
        ast_solo_sigma(ai), obx_solo_sigma(oi), pair_sig_comb(pi), pair_score(pi), ...
        pure_str, mat2str(ast_mult(ai,:)), mat2str(obx_mult(oi,:)));
end

%% ====================== BEST CONFIG DETAILED ANALYSIS ==================
best_pi   = pair_sort(1);
best_ai   = pair_ast_idx(best_pi);
best_oi   = pair_obx_idx(best_pi);

% Asterix best config
best_cols_ast = ast_configs(best_ai,:);
best_pos_ast  = master_pos(best_cols_ast,:);
best_dir_ast  = master_dir(best_cols_ast,:);

% Obelix best config
best_cols_obx = obx_configs(best_oi,:);
best_pos_obx_body   = master_pos(best_cols_obx,:);
best_dir_obx_body   = master_dir(best_cols_obx,:);
best_pos_obx_global = best_pos_obx_body + obx.geom_center_global;

% Reconstruct B matrices
B_ast_solo = B_master_ast_solo(:, best_cols_ast);
B_obx_solo = B_master_obx_solo(:, best_cols_obx);
B_obx_rot  = B_obx_solo(4:6,:);
col_idx_comb = [best_cols_ast, best_cols_obx + N_PLACEMENTS];
B_comb = B_master_comb(:, col_idx_comb);

if strcmp(ast_ctrl_mode, 'full6dof')
    met_ast = evaluate_config(B_ast_solo, ast.I, N_THR_AST, thr.mass);
else
    met_ast = evaluate_config_rot(B_ast_solo, ast.I, N_THR_AST, thr.mass);
end
if strcmp(obx_ctrl_mode, 'full6dof')
    met_obx = evaluate_config(B_obx_solo, obx.I, N_THR_OBX, thr.mass);
else
    met_obx = evaluate_config_rot(B_obx_solo, obx.I, N_THR_OBX, thr.mass);
end
met_comb = evaluate_config(B_comb, I_combined, N_THR_AST + N_THR_OBX, thr.mass);

fprintf('\n--- Best Configuration Detail ---\n');
fprintf('%s (%s): %d mating / %d outer,  mult = %s\n', ...
    ast.name, ast_ctrl_mode, ast_n_mate(best_ai), ast_n_outer(best_ai), mat2str(ast_mult(best_ai,:)));
fprintf('%s (%s): %d mating / %d outer,  mult = %s\n', ...
    obx.name, obx_ctrl_mode, obx_n_mate(best_oi), obx_n_outer(best_oi), mat2str(obx_mult(best_oi,:)));

if strcmp(ast_ctrl_mode, 'full6dof')
    fprintf('%s solo:  rank=%d  sigma_min=%.4f  sigma_max=%.4f  cond=%.2f  pure=%s\n', ...
        ast.name, met_ast.rank, met_ast.sigma_min, met_ast.sigma_max, met_ast.cond, ...
        mat2str(ast_pure_flag(best_ai)));
else
    fprintf('%s solo (ROT):  rot_rank=%d  rot_sigma_min=%.4f  pure=%s\n', ...
        ast.name, met_ast.rot_rank, met_ast.rot_sigma_min, mat2str(ast_pure_flag(best_ai)));
end
if strcmp(obx_ctrl_mode, 'full6dof')
    fprintf('%s solo:  rank=%d  sigma_min=%.4f  sigma_max=%.4f  cond=%.2f  pure=%s\n', ...
        obx.name, met_obx.rank, met_obx.sigma_min, met_obx.sigma_max, met_obx.cond, ...
        mat2str(obx_pure_flag(best_oi)));
else
    fprintf('%s solo (ROT):  rot_rank=%d  rot_sigma_min=%.4f  pure=%s\n', ...
        obx.name, met_obx.rot_rank, met_obx.rot_sigma_min, mat2str(obx_pure_flag(best_oi)));
end
fprintf('Combined:      rank=%d  sigma_min=%.4f  sigma_max=%.4f  cond=%.2f  pure_ctrl=%s\n', ...
    met_comb.rank, met_comb.sigma_min, met_comb.sigma_max, met_comb.cond, ...
    mat2str(pair_pure_ctrl(best_pi)));

%% --- Pure Control Feasibility (Best Config, all 12 wrenches for reference) ---
fprintf('\n--- Pure Control Feasibility (Best Config, all 12 wrenches) ---\n');
for d = 1:6
    for s = [1, -1]
        e_d = zeros(6,1); e_d(d) = s;
        [u_sol, resnorm] = lsqnonneg(B_comb, e_d);
        sign_str = '+'; if s < 0; sign_str = '-'; end

        % Check if this wrench was required
        is_required = false;
        if has_comb_pure
            for rr = 1:size(comb_pure_exp, 1)
                if comb_pure_exp(rr,1) == d && comb_pure_exp(rr,2) == s
                    is_required = true; break;
                end
            end
        end
        req_tag = ''; if is_required; req_tag = ' [REQUIRED]'; end

        if resnorm < 1e-6
            n_active = sum(u_sol > 1e-8);
            fprintf('  %s%s: FEASIBLE  (%d thrusters active, residual=%.2e)%s\n', ...
                sign_str, dof_labels_all{d}, n_active, resnorm, req_tag);
        else
            fprintf('  %s%s: INFEASIBLE  (residual=%.4f)%s\n', ...
                sign_str, dof_labels_all{d}, resnorm, req_tag);
        end
    end
end

%% --- Asterix solo pure control detail ---
fprintf('\n--- %s Solo Pure Control Detail ---\n', ast.name);
if strcmp(ast_ctrl_mode, 'full6dof')
    B_ast_check = B_ast_solo;
    n_rows_ast = 6;
else
    B_ast_check = B_obx_rot;  % not used for ast, just for structure
    B_ast_check = B_ast_solo(4:6,:);
    n_rows_ast = 3;
end
for d = 1:6
    for s = [1, -1]
        % For attitude mode, only show rotation DOFs
        if strcmp(ast_ctrl_mode, 'attitude') && d < 4; continue; end

        if strcmp(ast_ctrl_mode, 'attitude')
            e_d = zeros(3,1); e_d(d-3) = s;
            [u_sol, resnorm] = lsqnonneg(B_ast_solo(4:6,:), e_d);
        else
            e_d = zeros(6,1); e_d(d) = s;
            [u_sol, resnorm] = lsqnonneg(B_ast_solo, e_d);
        end
        sign_str = '+'; if s < 0; sign_str = '-'; end

        % Check if required
        is_required = false;
        if ~isempty(ast_pure_exp)
            for rr = 1:size(ast_pure_exp, 1)
                if ast_pure_exp(rr,1) == d && ast_pure_exp(rr,2) == s
                    is_required = true; break;
                end
            end
        end
        req_tag = ''; if is_required; req_tag = ' [REQUIRED]'; end

        if resnorm < 1e-6
            active = find(u_sol > 1e-8);
            fprintf('  %s%s: FEASIBLE  thrusters: ', sign_str, dof_labels_all{d});
            for ii = 1:length(active)
                fprintf('T%d(%.3f) ', active(ii), u_sol(active(ii)));
            end
            fprintf('%s\n', req_tag);
        else
            fprintf('  %s%s: INFEASIBLE  (residual=%.4f)%s\n', ...
                sign_str, dof_labels_all{d}, resnorm, req_tag);
        end
    end
end

%% --- Obelix solo pure control detail ---
fprintf('\n--- %s Solo Pure Control Detail ---\n', obx.name);
for d = 1:6
    for s = [1, -1]
        if strcmp(obx_ctrl_mode, 'attitude') && d < 4; continue; end

        if strcmp(obx_ctrl_mode, 'attitude')
            e_d = zeros(3,1); e_d(d-3) = s;
            [u_sol, resnorm] = lsqnonneg(B_obx_solo(4:6,:), e_d);
        else
            e_d = zeros(6,1); e_d(d) = s;
            [u_sol, resnorm] = lsqnonneg(B_obx_solo, e_d);
        end
        sign_str = '+'; if s < 0; sign_str = '-'; end

        is_required = false;
        if ~isempty(obx_pure_exp)
            for rr = 1:size(obx_pure_exp, 1)
                if obx_pure_exp(rr,1) == d && obx_pure_exp(rr,2) == s
                    is_required = true; break;
                end
            end
        end
        req_tag = ''; if is_required; req_tag = ' [REQUIRED]'; end

        if resnorm < 1e-6
            active = find(u_sol > 1e-8);
            fprintf('  %s%s: FEASIBLE  thrusters: ', sign_str, dof_labels_all{d});
            for ii = 1:length(active)
                fprintf('T%d(%.3f) ', active(ii), u_sol(active(ii)));
            end
            fprintf('%s\n', req_tag);
        else
            fprintf('  %s%s: INFEASIBLE  (residual=%.4f)%s\n', ...
                sign_str, dof_labels_all{d}, resnorm, req_tag);
        end
    end
end

fprintf('\n%s thruster layout (body frame):\n', ast.name);
for i = 1:N_THR_AST
    v = master_vtx(best_cols_ast(i));
    mate_str = '';
    if ismember(v, ast.mating_verts); mate_str = ' [MATING]'; end
    fprintf('  T%d @ V%d (%+.2f,%+.2f,%+.2f) -> dir (%+d,%+d,%+d)%s\n', ...
        i, v, best_pos_ast(i,:), best_dir_ast(i,:), mate_str);
end

fprintf('\n%s thruster layout (body frame):\n', obx.name);
for i = 1:N_THR_OBX
    v = master_vtx(best_cols_obx(i));
    mate_str = '';
    if ismember(v, obx.mating_verts); mate_str = ' [MATING]'; end
    fprintf('  T%d @ V%d (%+.2f,%+.2f,%+.2f) -> dir (%+d,%+d,%+d)%s\n', ...
        i, v, best_pos_obx_body(i,:), best_dir_obx_body(i,:), mate_str);
end

fprintf('\n%s solo B matrix:\n', ast.name);
disp(B_ast_solo);
fprintf('%s solo SVs: ', ast.name);  fprintf('%.4f  ', svd(B_ast_solo)); fprintf('\n');

fprintf('\n%s solo B matrix (full 6-row, attitude rows 4:6 active):\n', obx.name);
disp(B_obx_solo);
fprintf('%s rotational SVs: ', obx.name);  fprintf('%.4f  ', svd(B_obx_rot)); fprintf('\n');

fprintf('\nCombined B matrix:\n');
disp(B_comb);
fprintf('Combined SVs: ');  fprintf('%.4f  ', svd(B_comb)); fprintf('\n');

%% ====================== ANGULAR ACCELERATION ==========================
fprintf('\n--- Angular Acceleration Authority [rad/s^2] ---\n');
fprintf('%-15s %10s %10s %10s\n', 'Mode', 'alpha_x', 'alpha_y', 'alpha_z');
if strcmp(ast_ctrl_mode, 'full6dof')
    fprintf('%-15s %10.4f %10.4f %10.4f\n', [ast.name ' solo'], met_ast.alpha_max);
else
    fprintf('%-15s %10.4f %10.4f %10.4f\n', [ast.name ' att'], met_ast.alpha_max);
end
if strcmp(obx_ctrl_mode, 'full6dof')
    fprintf('%-15s %10.4f %10.4f %10.4f\n', [obx.name ' solo'], met_obx.alpha_max);
else
    fprintf('%-15s %10.4f %10.4f %10.4f\n', [obx.name ' att'], met_obx.alpha_max);
end
fprintf('%-15s %10.4f %10.4f %10.4f\n', 'Combined', met_comb.alpha_max);

%% ====================== PER-AXIS AUTHORITY =============================
fprintf('\n--- Per-Axis Control Authority ---\n');
dof_labels_full = {'Fx','Fy','Fz','Mx','My','Mz'};
dof_labels_rot  = {'Mx','My','Mz'};
fprintf('%-15s', 'Mode');
for d = 1:6; fprintf('%10s', dof_labels_full{d}); end
fprintf('%10s\n', 'Weakest');
fprintf('%s\n', repmat('-', 1, 85));

% Asterix solo
if strcmp(ast_ctrl_mode, 'full6dof')
    Bp = pinv(B_ast_solo);
    auth = zeros(1,6);
    for d = 1:6
        e_d = zeros(6,1); e_d(d) = 1;
        auth(d) = 1.0 / norm(Bp * e_d);
    end
    [~, worst_idx] = min(auth);
    fprintf('%-15s', [ast.name ' solo']);
    for d = 1:6; fprintf('%10.4f', auth(d)); end
    fprintf('%10s\n', dof_labels_full{worst_idx});
else
    fprintf('%-15s', [ast.name ' att']);
    Bp_rot = pinv(B_ast_solo(4:6,:));
    auth_rot = zeros(1,3);
    for d = 1:3
        e_d = zeros(3,1); e_d(d) = 1;
        auth_rot(d) = 1.0 / norm(Bp_rot * e_d);
    end
    [~, worst_rot] = min(auth_rot);
    fprintf('%10s%10s%10s', '---', '---', '---');
    for d = 1:3; fprintf('%10.4f', auth_rot(d)); end
    fprintf('%10s\n', dof_labels_rot{worst_rot});
end

% Obelix solo
if strcmp(obx_ctrl_mode, 'full6dof')
    Bp = pinv(B_obx_solo);
    auth = zeros(1,6);
    for d = 1:6
        e_d = zeros(6,1); e_d(d) = 1;
        auth(d) = 1.0 / norm(Bp * e_d);
    end
    [~, worst_idx] = min(auth);
    fprintf('%-15s', [obx.name ' solo']);
    for d = 1:6; fprintf('%10.4f', auth(d)); end
    fprintf('%10s\n', dof_labels_full{worst_idx});
else
    fprintf('%-15s', [obx.name ' att']);
    Bp_rot = pinv(B_obx_rot);
    auth_rot = zeros(1,3);
    for d = 1:3
        e_d = zeros(3,1); e_d(d) = 1;
        auth_rot(d) = 1.0 / norm(Bp_rot * e_d);
    end
    [~, worst_rot] = min(auth_rot);
    fprintf('%10s%10s%10s', '---', '---', '---');
    for d = 1:3; fprintf('%10.4f', auth_rot(d)); end
    fprintf('%10s\n', dof_labels_rot{worst_rot});
end

% Combined (always full 6-DOF)
Bp = pinv(B_comb);
auth = zeros(1,6);
for d = 1:6
    e_d = zeros(6,1); e_d(d) = 1;
    auth(d) = 1.0 / norm(Bp * e_d);
end
[~, worst_idx] = min(auth);
fprintf('%-15s', 'Combined');
for d = 1:6; fprintf('%10.4f', auth(d)); end
fprintf('%10s\n', dof_labels_full{worst_idx});

%% ====================== SIGNED AUTHORITY ===============================
fprintf('\n--- Signed Control Authority ---\n');
fprintf('%-15s', 'Mode');
for d = 1:6; fprintf('  %5s+  %5s-', dof_labels_full{d}, dof_labels_full{d}); end
fprintf('\n');
fprintf('%s\n', repmat('-', 1, 15 + 6*14));
for mode = 1:3
    if mode == 1; B_use = B_ast_solo; lbl = [ast.name ' solo'];
    elseif mode == 2; B_use = B_obx_solo; lbl = [obx.name ' (full B)'];
    else;            B_use = B_comb; lbl = 'Combined';
    end
    fprintf('%-15s', lbl);
    for d = 1:6
        pos_sum = sum(B_use(d, B_use(d,:) > 0));
        neg_sum = sum(B_use(d, B_use(d,:) < 0));
        fprintf('  %+7.3f %+7.3f', pos_sum, neg_sum);
    end
    fprintf('\n');
end
fprintf('\n');

%% ====================== TRANSLATION PAIR CHECK =========================
fprintf('\n--- Translation Pair Check ---\n');
ax_names = {'X','Y','Z'};
for mode = 1:3
    if mode == 1
        fprintf('\n  %s SOLO (%d thrusters, %s):\n', ast.name, N_THR_AST, ast_ctrl_mode);
        pos_t = best_pos_ast; dir_t = best_dir_ast; com_off = ast.com;
    elseif mode == 2
        fprintf('\n  %s SOLO (%d thrusters, %s; showing for reference):\n', obx.name, N_THR_OBX, obx_ctrl_mode);
        pos_t = best_pos_obx_body; dir_t = best_dir_obx_body; com_off = obx.com_offset;
    else
        fprintf('\n  COMBINED (%d thrusters):\n', N_THR_AST + N_THR_OBX);
        pos_t = [best_pos_ast; best_pos_obx_global];
        dir_t = [best_dir_ast; best_dir_obx_body];
        com_off = combined_com;
    end
    pos_rel = pos_t - com_off;
    nthr = size(pos_t, 1);
    for ax = 1:3
        e_ax = zeros(1,3); e_ax(ax) = 1;
        fprintf('    %s: ', ax_names{ax});
        found = false;
        for i = 1:nthr
            for j = (i+1):nthr
                F_net = thr.thrust*(dir_t(i,:) + dir_t(j,:));
                M_net = cross(pos_rel(i,:), thr.thrust*dir_t(i,:)) + ...
                        cross(pos_rel(j,:), thr.thrust*dir_t(j,:));
                if norm(F_net) > 1e-6
                    F_hat = F_net / norm(F_net);
                    if abs(dot(F_hat, e_ax)) > 0.99 && norm(M_net) < 1e-6
                        fprintf('T%d+T%d ', i, j);
                        found = true;
                    end
                end
            end
        end
        if ~found; fprintf('(no pure pair)'); end
        fprintf('\n');
    end
end

%% ========== COPY-PASTE OUTPUT: ASTERIX =================================
fprintf('\n');
fprintf('%%=========================================================================\n');
fprintf('%% ASTERIX CONFIG (%s)\n', ast_ctrl_mode);
fprintf('%%=========================================================================\n');

fprintf('ast_thr_pos = [ ...\n');
for i = 1:N_THR_AST
    fprintf('  %+8.4f %+8.4f %+8.4f', best_pos_ast(i,:));
    if i < N_THR_AST; fprintf('; ...  %% T%d @ V%d\n', i, master_vtx(best_cols_ast(i)));
    else;              fprintf('  ...  %% T%d @ V%d\n', i, master_vtx(best_cols_ast(i)));
    end
end
fprintf('];\n\n');

fprintf('ast_thr_dir = [ ...\n');
for i = 1:N_THR_AST
    fprintf('  %+d %+d %+d', best_dir_ast(i,:));
    if i < N_THR_AST; fprintf('; ...\n'); else; fprintf('  ...\n'); end
end
fprintf('];\n\n');

fprintf('ast.name     = ''%s'';\n', ast.name);
fprintf('ast.mass     = %.1f;      %% [kg]\n', ast.mass);
fprintf('ast.Lx       = %.4f;    %% [m]\n', ast.Lx);
fprintf('ast.Ly       = %.4f;    %% [m]\n', ast.Ly);
fprintf('ast.Lz       = %.4f;    %% [m]\n', ast.Lz);
fprintf('ast.com      = [%+.4f %+.4f %+.4f];  %% [m] COM in body frame\n', ast.com);
fprintf('ast.I        = diag([%.4f, %.4f, %.4f]);  %% [kg*m^2] about COM\n', ...
    ast.Ixx_com, ast.Iyy_com, ast.Izz_com);
fprintf('ast.thrust   = %.4f;    %% [N] per thruster\n', thr.thrust);
fprintf('ast.Isp      = %.1f;     %% [s]\n', thr.Isp);
fprintf('ast.thr_mass = %.4f;    %% [kg] per thruster\n', thr.mass);
fprintf('ast.ctrl_mode = ''%s'';\n', ast_ctrl_mode);
if ~isempty(ast_pure_req)
    fprintf('ast.pure_req = %s;  %% [dof, sign]\n', mat2str(ast_pure_req));
end
fprintf('%%=========================================================================\n\n');

%% ========== COPY-PASTE OUTPUT: OBELIX =================================
fprintf('%%=========================================================================\n');
fprintf('%% OBELIX CONFIG (%s)\n', obx_ctrl_mode);
fprintf('%%=========================================================================\n');

fprintf('obx_thr_pos = [ ...\n');
for i = 1:N_THR_OBX
    fprintf('  %+8.4f %+8.4f %+8.4f', best_pos_obx_body(i,:));
    if i < N_THR_OBX; fprintf('; ...  %% T%d @ V%d\n', i, master_vtx(best_cols_obx(i)));
    else;              fprintf('  ...  %% T%d @ V%d\n', i, master_vtx(best_cols_obx(i)));
    end
end
fprintf('];\n\n');

fprintf('obx_thr_dir = [ ...\n');
for i = 1:N_THR_OBX
    fprintf('  %+d %+d %+d', best_dir_obx_body(i,:));
    if i < N_THR_OBX; fprintf('; ...\n'); else; fprintf('  ...\n'); end
end
fprintf('];\n\n');

fprintf('obx.name     = ''%s'';\n', obx.name);
fprintf('obx.mass     = %.1f;      %% [kg]\n', obx.mass);
fprintf('obx.Lx       = %.4f;    %% [m]\n', obx.Lx);
fprintf('obx.Ly       = %.4f;    %% [m]\n', obx.Ly);
fprintf('obx.Lz       = %.4f;    %% [m]\n', obx.Lz);
fprintf('obx.com      = [%+.4f %+.4f %+.4f];  %% [m] COM in Obelix body frame\n', obx.com_offset);
fprintf('obx.I        = diag([%.4f, %.4f, %.4f]);  %% [kg*m^2] about COM\n', ...
    obx.Ixx_com, obx.Iyy_com, obx.Izz_com);
fprintf('obx.thrust   = %.4f;    %% [N] per thruster\n', thr.thrust);
fprintf('obx.Isp      = %.1f;     %% [s]\n', thr.Isp);
fprintf('obx.thr_mass = %.4f;    %% [kg] per thruster\n', thr.mass);
fprintf('obx.ctrl_mode = ''%s'';\n', obx_ctrl_mode);
if ~isempty(obx_pure_req)
    fprintf('obx.pure_req = %s;  %% [dof, sign]\n', mat2str(obx_pure_req));
end
fprintf('obx.geom_center_global = [%+.4f %+.4f %+.4f];\n', obx.geom_center_global);
fprintf('%%=========================================================================\n\n');

%% ========== TOP 5 ALTERNATIVES =========================================
n_alt = min(5, n_show);
fprintf('%% --- Top %d config pairs (alternatives) ---\n', n_alt);
for k = 1:n_alt
    pi = pair_sort(k);
    ai = pair_ast_idx(pi); oi = pair_obx_idx(pi);
    cols_a = ast_configs(ai,:); cols_o = obx_configs(oi,:);
    if has_comb_pure
        if pair_pure_ctrl(pi); pure_str = 'YES';
        elseif pair_bidir_ok(pi); pure_str = 'no';
        else; pure_str = '?'; end
    else
        pure_str = 'N/A';
    end
    fprintf('%% #%d  score=%.4f  sig_A=%.4f  sig_O=%.4f  sig_C=%.4f  A(%d/%d) O(%d/%d) pure=%s\n', ...
        k, pair_score(pi), ast_solo_sigma(ai), obx_solo_sigma(oi), pair_sig_comb(pi), ...
        ast_n_mate(ai), ast_n_outer(ai), obx_n_mate(oi), obx_n_outer(oi), pure_str);
    fprintf('%%   ast_thr_pos = [');
    for i = 1:N_THR_AST
        fprintf('%+.4f %+.4f %+.4f', master_pos(cols_a(i),:));
        if i < N_THR_AST; fprintf('; '); end
    end
    fprintf('];\n');
    fprintf('%%   ast_thr_dir = [');
    for i = 1:N_THR_AST
        fprintf('%+d %+d %+d', master_dir(cols_a(i),:));
        if i < N_THR_AST; fprintf('; '); end
    end
    fprintf('];\n');
    fprintf('%%   obx_thr_pos = [');
    for i = 1:N_THR_OBX
        fprintf('%+.4f %+.4f %+.4f', master_pos(cols_o(i),:));
        if i < N_THR_OBX; fprintf('; '); end
    end
    fprintf('];\n');
    fprintf('%%   obx_thr_dir = [');
    for i = 1:N_THR_OBX
        fprintf('%+d %+d %+d', master_dir(cols_o(i),:));
        if i < N_THR_OBX; fprintf('; '); end
    end
    fprintf('];\n');
end

%% ====================== .MAT OUTPUT ====================================
mat_out.timestamp = datestr(now, 'yyyy-mm-dd_HH:MM:SS');

% --- Asterix config ---
mat_out.asterix.name       = ast.name;
mat_out.asterix.ctrl_mode  = ast_ctrl_mode;
mat_out.asterix.mass       = ast.mass;
mat_out.asterix.Lx         = ast.Lx;
mat_out.asterix.Ly         = ast.Ly;
mat_out.asterix.Lz         = ast.Lz;
mat_out.asterix.com        = ast.com;
mat_out.asterix.I          = ast.I;
mat_out.asterix.thrust     = thr.thrust;
mat_out.asterix.Isp        = thr.Isp;
mat_out.asterix.thr_mass   = thr.mass;
mat_out.asterix.n_thrusters = N_THR_AST;
mat_out.asterix.thr_pos    = best_pos_ast;
mat_out.asterix.thr_dir    = best_dir_ast;
mat_out.asterix.thr_vertex = master_vtx(best_cols_ast);
mat_out.asterix.mating_verts = ast.mating_verts;
mat_out.asterix.outer_verts  = ast.outer_verts;
mat_out.asterix.n_mating   = ast_n_mate(best_ai);
mat_out.asterix.n_outer    = ast_n_outer(best_ai);
mat_out.asterix.multiplicity = ast_mult(best_ai,:);
mat_out.asterix.B_solo     = B_ast_solo;
mat_out.asterix.svd_solo   = svd(B_ast_solo);
mat_out.asterix.pure_req   = ast_pure_req;
mat_out.asterix.pure_flag  = ast_pure_flag(best_ai);
mat_out.asterix.pure_fx    = ast_pure_flag(best_ai);  % backward compat alias
if strcmp(ast_ctrl_mode, 'full6dof')
    mat_out.asterix.rank       = met_ast.rank;
    mat_out.asterix.sigma_min  = met_ast.sigma_min;
    mat_out.asterix.sigma_max  = met_ast.sigma_max;
    mat_out.asterix.cond       = met_ast.cond;
    mat_out.asterix.max_force  = met_ast.max_force;
    mat_out.asterix.max_torque = met_ast.max_torque;
    mat_out.asterix.alpha_max  = met_ast.alpha_max;
else
    mat_out.asterix.rot_rank      = met_ast.rot_rank;
    mat_out.asterix.rot_sigma_min = met_ast.rot_sigma_min;
    mat_out.asterix.rot_sigma_max = met_ast.rot_sigma_max;
    mat_out.asterix.rot_cond      = met_ast.rot_cond;
    mat_out.asterix.max_torque    = met_ast.max_torque;
    mat_out.asterix.alpha_max     = met_ast.alpha_max;
    % Provide full-B fields too for backward compat (computed from full 6-row B)
    sv_full = svd(B_ast_solo);
    mat_out.asterix.rank       = rank(B_ast_solo, 1e-8);
    mat_out.asterix.sigma_min  = sv_full(end);
    mat_out.asterix.sigma_max  = sv_full(1);
    mat_out.asterix.cond       = sv_full(1) / max(sv_full(end), 1e-15);
    for ax = 1:3
        mat_out.asterix.max_force(ax) = max(abs(B_ast_solo(ax,:)));
    end
end

% --- Obelix config ---
mat_out.obelix.name       = obx.name;
mat_out.obelix.ctrl_mode  = obx_ctrl_mode;
mat_out.obelix.mass       = obx.mass;
mat_out.obelix.Lx         = obx.Lx;
mat_out.obelix.Ly         = obx.Ly;
mat_out.obelix.Lz         = obx.Lz;
mat_out.obelix.com        = obx.com_offset;
mat_out.obelix.I          = obx.I;
mat_out.obelix.thrust     = thr.thrust;
mat_out.obelix.Isp        = thr.Isp;
mat_out.obelix.thr_mass   = thr.mass;
mat_out.obelix.n_thrusters = N_THR_OBX;
mat_out.obelix.thr_pos    = best_pos_obx_body;
mat_out.obelix.thr_dir    = best_dir_obx_body;
mat_out.obelix.thr_vertex = master_vtx(best_cols_obx);
mat_out.obelix.mating_verts = obx.mating_verts;
mat_out.obelix.outer_verts  = obx.outer_verts;
mat_out.obelix.n_mating   = obx_n_mate(best_oi);
mat_out.obelix.n_outer    = obx_n_outer(best_oi);
mat_out.obelix.multiplicity = obx_mult(best_oi,:);
mat_out.obelix.B_solo     = B_obx_solo;
mat_out.obelix.B_rot      = B_obx_rot;
mat_out.obelix.svd_rot    = svd(B_obx_rot);
mat_out.obelix.pure_req   = obx_pure_req;
mat_out.obelix.pure_flag  = obx_pure_flag(best_oi);
mat_out.obelix.pure_fx    = obx_pure_flag(best_oi);  % backward compat alias
if strcmp(obx_ctrl_mode, 'full6dof')
    mat_out.obelix.rank       = met_obx.rank;
    mat_out.obelix.sigma_min  = met_obx.sigma_min;
    mat_out.obelix.sigma_max  = met_obx.sigma_max;
    mat_out.obelix.cond       = met_obx.cond;
    mat_out.obelix.max_force  = met_obx.max_force;
    mat_out.obelix.max_torque = met_obx.max_torque;
    mat_out.obelix.alpha_max  = met_obx.alpha_max;
    % Also populate rot fields for backward compat
    sv_rot = svd(B_obx_rot);
    mat_out.obelix.rot_rank      = rank(B_obx_rot, 1e-8);
    mat_out.obelix.rot_sigma_min = sv_rot(end);
    mat_out.obelix.rot_sigma_max = sv_rot(1);
    mat_out.obelix.rot_cond      = sv_rot(1) / max(sv_rot(end), 1e-15);
else
    mat_out.obelix.rot_rank      = met_obx.rot_rank;
    mat_out.obelix.rot_sigma_min = met_obx.rot_sigma_min;
    mat_out.obelix.rot_sigma_max = met_obx.rot_sigma_max;
    mat_out.obelix.rot_cond      = met_obx.rot_cond;
    mat_out.obelix.max_torque    = met_obx.max_torque;
    mat_out.obelix.alpha_max     = met_obx.alpha_max;
    % Also populate full-B fields for backward compat
    sv_full = svd(B_obx_solo);
    mat_out.obelix.rank       = rank(B_obx_solo, 1e-8);
    mat_out.obelix.sigma_min  = sv_full(end);
    mat_out.obelix.sigma_max  = sv_full(1);
    mat_out.obelix.cond       = sv_full(1) / max(sv_full(end), 1e-15);
end
mat_out.obelix.geom_center_global = obx.geom_center_global;

% --- Combined config ---
mat_out.combined.mass       = combined_mass;
mat_out.combined.com        = combined_com;
mat_out.combined.I          = I_combined;
mat_out.combined.n_thrusters = N_THR_AST + N_THR_OBX;
mat_out.combined.B          = B_comb;
mat_out.combined.svd        = svd(B_comb);
mat_out.combined.rank       = met_comb.rank;
mat_out.combined.sigma_min  = met_comb.sigma_min;
mat_out.combined.sigma_max  = met_comb.sigma_max;
mat_out.combined.cond       = met_comb.cond;
mat_out.combined.pure_req   = comb_pure_req;
mat_out.combined.pure_ctrl  = pair_pure_ctrl(best_pi);
mat_out.combined.max_force  = met_comb.max_force;
mat_out.combined.max_torque = met_comb.max_torque;
mat_out.combined.alpha_max  = met_comb.alpha_max;
mat_out.combined.score      = pair_score(best_pi);

% --- Top alternatives ---
n_alt_save = min(5, n_pure_ok);
if n_alt_save == 0; n_alt_save = min(5, n_bidir_pass); end
if n_alt_save == 0; n_alt_save = min(5, n_rank6); end
for k = 1:n_alt_save
    pi = pair_sort(k);
    ai = pair_ast_idx(pi);  oi = pair_obx_idx(pi);
    cols_a = ast_configs(ai,:);  cols_o = obx_configs(oi,:);
    alt.rank           = k;
    alt.score          = pair_score(pi);
    alt.pure_ctrl      = pair_pure_ctrl(pi);
    alt.sig_ast        = ast_solo_sigma(ai);
    alt.sig_obx        = obx_solo_sigma(oi);
    alt.sig_comb       = pair_sig_comb(pi);
    alt.ast_n_mate     = ast_n_mate(ai);
    alt.ast_n_outer    = ast_n_outer(ai);
    alt.obx_n_mate     = obx_n_mate(oi);
    alt.obx_n_outer    = obx_n_outer(oi);
    alt.ast_thr_pos    = master_pos(cols_a,:);
    alt.ast_thr_dir    = master_dir(cols_a,:);
    alt.obx_thr_pos    = master_pos(cols_o,:);
    alt.obx_thr_dir    = master_dir(cols_o,:);
    alt.ast_mult       = ast_mult(ai,:);
    alt.obx_mult       = obx_mult(oi,:);
    mat_out.alternatives(k) = alt;
end

% --- Weights / settings used ---
mat_out.settings.w_solo     = w_solo;
mat_out.settings.w_combined = w_combined;
mat_out.settings.w_mating   = w_mating;
mat_out.settings.K_TOP      = K_TOP;
mat_out.settings.N_PURE_TOP = N_PURE_TOP;
mat_out.settings.ast_ctrl_mode = ast_ctrl_mode;
mat_out.settings.obx_ctrl_mode = obx_ctrl_mode;
mat_out.settings.ast_pure_req  = ast_pure_req;
mat_out.settings.obx_pure_req  = obx_pure_req;
mat_out.settings.comb_pure_req = comb_pure_req;

mat_filename = 'rcs_config.mat';
save(mat_filename, '-struct', 'mat_out');
fprintf('\nSaved configuration to: %s\n', mat_filename);

%% ====================== FIGURES ========================================
fig_w = 700;  fig_h = 480;  bar_h = 210;

% --- Asterix Solo 3D ---
fig1 = figure('Name', [ast.name ' Solo (' upper(ast_ctrl_mode) ')'], 'Position', [20 300 fig_w fig_h]);
ax1 = axes(fig1); hold(ax1, 'on');
draw_spacecraft(ax1, ast.vertices, best_pos_ast, best_dir_ast, ast.a, ast.b, ast.c, ...
    [0.6 0.7 0.9], [0.3 0.3 0.5], [0.9 0.15 0.1], 'V', ast.mating_verts);
mfv = ast.vertices(ast.mating_verts,:);
[~, ord] = sort(atan2(mfv(:,3), mfv(:,2)));
patch(ax1, 'Vertices', mfv(ord,:), 'Faces', 1:4, ...
      'FaceColor', [1 0.85 0.3], 'FaceAlpha', 0.3, ...
      'EdgeColor', [0.8 0.6 0.0], 'LineWidth', 2, 'LineStyle', '--');
plot3(ax1, ast.com(1),ast.com(2),ast.com(3), 'k+', 'MarkerSize', 12, 'LineWidth', 2);
if strcmp(ast_ctrl_mode, 'full6dof')
    title(ax1, sprintf('%s Solo (6-DOF) - %d thr  rank=%d  \\sigma_{min}=%.3f', ...
        ast.name, N_THR_AST, met_ast.rank, met_ast.sigma_min));
else
    title(ax1, sprintf('%s Solo (ATT) - %d thr  rot\\_rank=%d  rot\\_\\sigma_{min}=%.3f', ...
        ast.name, N_THR_AST, met_ast.rot_rank, met_ast.rot_sigma_min));
end
xlabel(ax1,'X'); ylabel(ax1,'Y'); zlabel(ax1,'Z');
axis(ax1,'equal'); grid(ax1,'on'); view(ax1,135,25);

% --- Obelix Solo 3D ---
fig1b = figure('Name', [obx.name ' Solo (' upper(obx_ctrl_mode) ')'], 'Position', [20+fig_w+10 300 fig_w fig_h]);
ax1b = axes(fig1b); hold(ax1b, 'on');
draw_spacecraft(ax1b, obx.vertices_local, best_pos_obx_body, best_dir_obx_body, ...
    obx.a, obx.b, obx.c, [0.6 0.9 0.7], [0.3 0.5 0.3], [0.1 0.6 0.2], 'V', obx.mating_verts);
mfv_o = obx.vertices_local(obx.mating_verts,:);
[~, ord_o] = sort(atan2(mfv_o(:,3), mfv_o(:,2)));
patch(ax1b, 'Vertices', mfv_o(ord_o,:), 'Faces', 1:4, ...
      'FaceColor', [1 0.85 0.3], 'FaceAlpha', 0.3, ...
      'EdgeColor', [0.8 0.6 0.0], 'LineWidth', 2, 'LineStyle', '--');
plot3(ax1b, obx.com_offset(1),obx.com_offset(2),obx.com_offset(3), 'k+', 'MarkerSize', 12, 'LineWidth', 2);
if strcmp(obx_ctrl_mode, 'full6dof')
    title(ax1b, sprintf('%s Solo (6-DOF) - %d thr  rank=%d  \\sigma_{min}=%.3f', ...
        obx.name, N_THR_OBX, met_obx.rank, met_obx.sigma_min));
else
    title(ax1b, sprintf('%s Solo (ATT) - %d thr  rot\\_rank=%d  rot\\_\\sigma_{min}=%.3f', ...
        obx.name, N_THR_OBX, met_obx.rot_rank, met_obx.rot_sigma_min));
end
xlabel(ax1b,'X'); ylabel(ax1b,'Y'); zlabel(ax1b,'Z');
axis(ax1b,'equal'); grid(ax1b,'on'); view(ax1b,135,25);

% --- Combined 3D ---
fig2 = figure('Name','Combined', 'Position', [20 50 fig_w fig_h]);
ax2 = axes(fig2); hold(ax2, 'on');
draw_spacecraft(ax2, ast.vertices, best_pos_ast, best_dir_ast, ast.a, ast.b, ast.c, ...
    [0.6 0.7 0.9], [0.3 0.3 0.5], [0.9 0.15 0.1], 'A.', ast.mating_verts);
draw_spacecraft(ax2, obx.vertices_global, best_pos_obx_global, best_dir_obx_body, ...
    obx.a, obx.b, obx.c, [0.6 0.9 0.7], [0.3 0.5 0.3], [0.1 0.6 0.2], 'O.', obx.mating_verts);
yr = [-max(ast.b,obx.b) max(ast.b,obx.b) max(ast.b,obx.b) -max(ast.b,obx.b)]*1.15;
zr = [-max(ast.c,obx.c) -max(ast.c,obx.c) max(ast.c,obx.c) max(ast.c,obx.c)]*1.15;
patch(ax2, ast.a*ones(1,4), yr, zr, 'FaceColor', [1 0.85 0.3], ...
    'FaceAlpha', 0.2, 'EdgeColor', [0.8 0.6 0.0], 'LineWidth', 1.5, 'LineStyle', '--');
plot3(ax2, combined_com(1),combined_com(2),combined_com(3), 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'k');
pure_title = '';
if has_comb_pure && pair_pure_ctrl(best_pi); pure_title = '  PURE CTRL OK'; end
title(ax2, sprintf('%s(%d) + %s(%d) - rank=%d  \\sigma_{min}=%.3f%s', ...
    ast.name, N_THR_AST, obx.name, N_THR_OBX, met_comb.rank, met_comb.sigma_min, pure_title));
xlabel(ax2,'X'); ylabel(ax2,'Y'); zlabel(ax2,'Z');
axis(ax2,'equal'); grid(ax2,'on'); view(ax2,135,25);

% --- Bar charts ---
fig3 = figure('Name', 'Authority', 'Position', [20+fig_w+10 50 fig_w*2+10 bar_h]);
if strcmp(ast_ctrl_mode, 'full6dof')
    subplot(2,3,1); bar(met_ast.max_torque,'FaceColor',[0.5 0.6 0.85]);
    set(gca,'XTickLabel',{'X','Y','Z'}); ylabel('Torque [N·m]'); title([ast.name ' Torque']); grid on;
    subplot(2,3,4); bar(met_ast.max_force,'FaceColor',[0.5 0.6 0.85]);
    set(gca,'XTickLabel',{'X','Y','Z'}); ylabel('Force [N]'); title([ast.name ' Force']); grid on;
else
    subplot(2,3,1); bar(met_ast.max_torque,'FaceColor',[0.5 0.6 0.85]);
    set(gca,'XTickLabel',{'X','Y','Z'}); ylabel('Torque [N·m]'); title([ast.name ' Torque (att)']); grid on;
    subplot(2,3,4); bar(zeros(1,3),'FaceColor',[0.5 0.6 0.85]);
    set(gca,'XTickLabel',{'X','Y','Z'}); ylabel('Force [N]'); title([ast.name ' Force (N/A)']); grid on;
end
if strcmp(obx_ctrl_mode, 'full6dof')
    subplot(2,3,2); bar(met_obx.max_torque,'FaceColor',[0.5 0.8 0.6]);
    set(gca,'XTickLabel',{'X','Y','Z'}); ylabel('Torque [N·m]'); title([obx.name ' Torque']); grid on;
    subplot(2,3,5); bar(met_obx.max_force,'FaceColor',[0.5 0.8 0.6]);
    set(gca,'XTickLabel',{'X','Y','Z'}); ylabel('Force [N]'); title([obx.name ' Force']); grid on;
else
    subplot(2,3,2); bar(met_obx.max_torque,'FaceColor',[0.5 0.8 0.6]);
    set(gca,'XTickLabel',{'X','Y','Z'}); ylabel('Torque [N·m]'); title([obx.name ' Torque (att)']); grid on;
    subplot(2,3,5); bar(zeros(1,3),'FaceColor',[0.5 0.8 0.6]);
    set(gca,'XTickLabel',{'X','Y','Z'}); ylabel('Force [N]'); title([obx.name ' Force (N/A)']); grid on;
end
subplot(2,3,3); bar(met_comb.max_torque,'FaceColor',[0.7 0.6 0.8]);
set(gca,'XTickLabel',{'X','Y','Z'}); ylabel('Torque [N·m]'); title('Comb Torque'); grid on;
subplot(2,3,6); bar(met_comb.max_force,'FaceColor',[0.7 0.6 0.8]);
set(gca,'XTickLabel',{'X','Y','Z'}); ylabel('Force [N]'); title('Comb Force'); grid on;

% --- Pareto scatter ---
fig4 = figure('Name','Trade Space', 'Position', [20 750 fig_w fig_h]);
valid_pi = find(pair_bidir_ok(1:pi_count) & pair_score(1:pi_count) > -Inf);
total_mate = zeros(length(valid_pi), 1);
sig_a_plot = zeros(length(valid_pi), 1);
sig_c_plot = zeros(length(valid_pi), 1);
pure_plot  = false(length(valid_pi), 1);
for idx = 1:length(valid_pi)
    pi = valid_pi(idx);
    total_mate(idx) = ast_n_mate(pair_ast_idx(pi)) + obx_n_mate(pair_obx_idx(pi));
    sig_a_plot(idx) = ast_solo_sigma(pair_ast_idx(pi)) + obx_solo_sigma(pair_obx_idx(pi));
    sig_c_plot(idx) = pair_sig_comb(pi);
    pure_plot(idx)  = pair_pure_ctrl(pi);
end
scatter(sig_a_plot(~pure_plot), sig_c_plot(~pure_plot), 12, total_mate(~pure_plot), ...
    'filled', 'MarkerFaceAlpha', 0.3, 'Marker', 'o');
hold on;
if any(pure_plot)
    scatter(sig_a_plot(pure_plot), sig_c_plot(pure_plot), 25, total_mate(pure_plot), ...
        'filled', 'MarkerFaceAlpha', 0.7, 'Marker', 'd');
end
for k = 1:min(3, n_show)
    pi = pair_sort(k);
    sa = ast_solo_sigma(pair_ast_idx(pi)) + obx_solo_sigma(pair_obx_idx(pi));
    sc = pair_sig_comb(pi);
    plot(sa, sc, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    text(sa+0.005, sc+0.005, sprintf('#%d', k), ...
        'FontSize', 8, 'Color', 'r', 'FontWeight', 'bold');
end
cb = colorbar; ylabel(cb, 'Total # on mating faces');
colormap(flipud(autumn));
xlabel('\sigma_{min,A} + \sigma_{min,O}  (solo sum)');
ylabel('\sigma_{min} Combined');
if has_comb_pure
    title('Trade Space (diamonds = pure ctrl OK)'); grid on;
    legend({'bidir only', 'pure ctrl OK', 'top picks'}, 'Location', 'best');
else
    title('Trade Space'); grid on;
    legend({'configs', 'top picks'}, 'Location', 'best');
end

fprintf('\n=== Done. ===\n');

%% ====================== HELPER FUNCTIONS ===============================

function req_expanded = expand_pure_req(pure_req)
% Expand pure_req matrix: sign=0 becomes two rows (+1 and -1).
%   pure_req : Nx2 matrix [dof, sign] or []
%   Returns  : Mx2 matrix [dof, sign] with only +1/-1 in column 2
    if isempty(pure_req)
        req_expanded = zeros(0, 2);
        return;
    end
    req_expanded = zeros(0, 2);
    for i = 1:size(pure_req, 1)
        if pure_req(i, 2) == 0
            req_expanded = [req_expanded; pure_req(i,1), +1; pure_req(i,1), -1]; %#ok<AGROW>
        else
            req_expanded = [req_expanded; pure_req(i,:)]; %#ok<AGROW>
        end
    end
end

function str = format_pure_req(pure_req_exp, dof_labels)
% Format expanded pure_req into a readable string like "+Fx -Fx +My"
    parts = cell(1, size(pure_req_exp, 1));
    for i = 1:size(pure_req_exp, 1)
        if pure_req_exp(i,2) > 0; s = '+'; else; s = '-'; end
        parts{i} = [s, dof_labels{pure_req_exp(i,1)}];
    end
    str = strjoin(parts, ' ');
end

function feasible = check_pure_control_general(B, pure_req_expanded)
% Check lsqnonneg feasibility for each (dof, sign) pair in pure_req_expanded.
%   B                  : n_dof x n_thr matrix
%   pure_req_expanded  : Mx2 matrix [dof, sign] (sign is +1 or -1 only)
    tol = 1e-6;
    feasible = true;
    n_dof = size(B, 1);
    for i = 1:size(pure_req_expanded, 1)
        d = pure_req_expanded(i, 1);
        s = pure_req_expanded(i, 2);
        e_d = zeros(n_dof, 1);
        e_d(d) = s;
        [~, resnorm] = lsqnonneg(B, e_d);
        if resnorm > tol
            feasible = false;
            return;
        end
    end
end

function [configs, solo_sigma, solo_rank, mate_pen, n_mate, n_outer, mult_out, pure_flag] = ...
    enumerate_and_evaluate_solo(N_THR, nv, subsets_for_k, B_master_solo, ...
                                 mate_penalty_per_placement, mating_verts, outer_verts, ...
                                 ctrl_mode, pure_req_expanded)
% Enumerate all N_THR-thruster configs on 8 vertices x 3 edge-tangent directions.
%
% ctrl_mode:
%   'full6dof'  - rank 6 of full B
%   'attitude'  - rank 3 of B(4:6,:), bidirectional in all 3 rot axes always enforced
%
% pure_req_expanded:
%   Mx2 matrix [dof, sign] (expanded, sign is +1/-1 only).
%   For attitude mode, DOFs must be 4-6 (mapped internally to rows 1-3 of B_rot).
%   Empty [] means no pure-control requirements.
%
% Returns:
%   configs    - M x N_THR matrix of placement indices (1..24)
%   solo_sigma - M x 1 minimum singular values (of relevant sub-matrix)
%   solo_rank  - M x 1 rank (6 for full, 3 for attitude)
%   mate_pen   - M x 1 total mating penalty
%   n_mate     - M x 1 number of thrusters on mating-face vertices
%   n_outer    - M x 1 number on outer vertices
%   mult_out   - M x nv multiplicity vectors
%   pure_flag  - M x 1 logical: true if pure control requirements passed

    is_attitude = strcmp(ctrl_mode, 'attitude');
    has_pure = ~isempty(pure_req_expanded);

    % For attitude mode, remap DOF indices 4-6 -> 1-3 for B_rot
    if is_attitude && has_pure
        pure_req_rot = pure_req_expanded;
        pure_req_rot(:,1) = pure_req_rot(:,1) - 3;
    elseif is_attitude
        pure_req_rot = zeros(0, 2);
    else
        pure_req_rot = []; % unused
    end

    % Derive bidirectional DOFs for pre-filtering (from pure_req)
    if has_pure
        if is_attitude
            bidir_dofs_from_pure = unique(pure_req_rot(:,1))';
        else
            bidir_dofs_from_pure = unique(pure_req_expanded(:,1))';
        end
    else
        bidir_dofs_from_pure = [];
    end

    % Enumerate multiplicity vectors
    mult_vecs = zeros(0, nv);
    for code = 0:(4^nv - 1)
        kvec = zeros(1, nv);
        temp = code;
        for v = 1:nv
            kvec(v) = mod(temp, 4);
            temp = floor(temp / 4);
        end
        if sum(kvec) == N_THR
            mult_vecs = [mult_vecs; kvec]; %#ok<AGROW>
        end
    end
    n_mult = size(mult_vecs, 1);

    % Count total configs
    total_configs = 0;
    for m = 1:n_mult
        nc = 1;
        for v = 1:nv
            nc = nc * nchoosek(3, mult_vecs(m,v));
        end
        total_configs = total_configs + nc;
    end
    fprintf('  %d multiplicity vectors,  %d total configs\n', n_mult, total_configs);

    % Preallocate
    configs    = zeros(total_configs, N_THR);
    solo_sigma = zeros(total_configs, 1);
    solo_rank  = zeros(total_configs, 1);
    mate_pen   = zeros(total_configs, 1);
    n_mate     = zeros(total_configs, 1);
    n_outer    = zeros(total_configs, 1);
    mult_out   = zeros(total_configs, nv);
    pure_flag  = false(total_configs, 1);

    cfg_count = 0;
    n_rank_pass = 0;

    for m = 1:n_mult
        kvec = mult_vecs(m,:);
        n_mate_fixed  = sum(kvec(mating_verts));
        n_outer_fixed = sum(kvec(outer_verts));

        subset_options = cell(nv, 1);
        radices = zeros(1, nv);
        for v = 1:nv
            subset_options{v} = subsets_for_k{kvec(v)+1};
            radices(v) = size(subset_options{v}, 1);
        end
        n_combos_m = prod(radices);

        for combo = 0:(n_combos_m - 1)
            cfg_count = cfg_count + 1;

            temp = combo;
            col_idx = zeros(1, N_THR);
            ci_pos = 0;
            for v = 1:nv
                si = mod(temp, radices(v)) + 1;
                temp = floor(temp / radices(v));
                mask = subset_options{v}(si,:);
                for d = 1:3
                    if mask(d)
                        ci_pos = ci_pos + 1;
                        col_idx(ci_pos) = 3*(v-1) + d;
                    end
                end
            end

            configs(cfg_count,:)  = col_idx;
            mult_out(cfg_count,:) = kvec;
            n_mate(cfg_count)     = n_mate_fixed;
            n_outer(cfg_count)    = n_outer_fixed;
            mate_pen(cfg_count)   = sum(mate_penalty_per_placement(col_idx));

            B_full = B_master_solo(:, col_idx);

            if is_attitude
                % ---- ATTITUDE MODE ----
                B_rot = B_full(4:6, :);
                rk = rank(B_rot, 1e-8);
                solo_rank(cfg_count) = rk;
                if rk < 3; continue; end

                % Always require bidirectional in all 3 rotation axes
                bidir_ok = true;
                for row = 1:3
                    if ~(any(B_rot(row,:) > 1e-10) && any(B_rot(row,:) < -1e-10))
                        bidir_ok = false; break;
                    end
                end
                if ~bidir_ok; solo_rank(cfg_count) = 0; continue; end

                sv = svd(B_rot);
                solo_sigma(cfg_count) = sv(end);
                n_rank_pass = n_rank_pass + 1;

                % Pure control check (on B_rot with remapped DOFs)
                if ~has_pure
                    pure_flag(cfg_count) = true;  % no requirements beyond bidirectionality
                else
                    pure_flag(cfg_count) = check_pure_control_general(B_rot, pure_req_rot);
                end
            else
                % ---- FULL 6-DOF MODE ----
                rk = rank(B_full, 1e-8);
                solo_rank(cfg_count) = rk;
                if rk < 6; continue; end

                % Bidirectional pre-filter on pure_req DOFs
                bidir_ok = true;
                for bd = bidir_dofs_from_pure
                    if ~(any(B_full(bd,:) > 1e-10) && any(B_full(bd,:) < -1e-10))
                        bidir_ok = false; break;
                    end
                end
                if ~bidir_ok; continue; end

                sv = svd(B_full);
                solo_sigma(cfg_count) = sv(end);
                n_rank_pass = n_rank_pass + 1;

                % Pure control check
                if ~has_pure
                    pure_flag(cfg_count) = true;  % no requirements
                else
                    pure_flag(cfg_count) = check_pure_control_general(B_full, pure_req_expanded);
                end
            end
        end
    end

    fprintf('  Rank-pass: %d,  Pure-ctrl pass: %d\n', n_rank_pass, sum(pure_flag));

    % Trim
    configs    = configs(1:cfg_count,:);
    solo_sigma = solo_sigma(1:cfg_count);
    solo_rank  = solo_rank(1:cfg_count);
    mate_pen   = mate_pen(1:cfg_count);
    n_mate     = n_mate(1:cfg_count);
    n_outer    = n_outer(1:cfg_count);
    mult_out   = mult_out(1:cfg_count,:);
    pure_flag  = pure_flag(1:cfg_count);
end

function metrics = evaluate_config(B, I_mat, n_thr, thr_mass)
    metrics.n_thrusters    = n_thr;
    metrics.total_thr_mass = n_thr * thr_mass;
    metrics.rank           = rank(B, 1e-8);
    metrics.controllable   = (metrics.rank == 6);
    sv = svd(B);
    metrics.sigma_min = sv(end);
    metrics.sigma_max = sv(1);
    metrics.cond      = sv(1) / max(sv(end), 1e-15);
    for ax = 1:3
        metrics.max_force(ax)  = max(abs(B(ax,:)));
        metrics.max_torque(ax) = max(abs(B(3+ax,:)));
    end
    metrics.alpha_max = I_mat \ metrics.max_torque';
end

function metrics = evaluate_config_rot(B_full, I_mat, n_thr, thr_mass)
% Evaluate attitude-only metrics using rotational sub-matrix (rows 4:6).
    B_rot = B_full(4:6, :);
    metrics.n_thrusters    = n_thr;
    metrics.total_thr_mass = n_thr * thr_mass;
    metrics.rot_rank       = rank(B_rot, 1e-8);
    sv_rot = svd(B_rot);
    metrics.rot_sigma_min  = sv_rot(end);
    metrics.rot_sigma_max  = sv_rot(1);
    metrics.rot_cond       = sv_rot(1) / max(sv_rot(end), 1e-15);
    for ax = 1:3
        metrics.max_torque(ax) = max(abs(B_rot(ax,:)));
    end
    metrics.alpha_max = I_mat \ metrics.max_torque';
end

function draw_spacecraft(ax, verts, thr_pos, thr_dir, a, b, c, ...
                         clr_face, clr_edge, clr_arrow, label_prefix, mating_v)
    faces_def = [1 2 3 4; 5 6 7 8; 1 2 6 5; 3 4 8 7; 1 4 8 5; 2 3 7 6];
    patch(ax, 'Vertices', verts, 'Faces', faces_def, ...
          'FaceColor', clr_face, 'FaceAlpha', 0.2, ...
          'EdgeColor', clr_edge, 'LineWidth', 1.3);
    scale = max([a b c]) * 0.30;
    nthr = size(thr_pos, 1);
    for i = 1:nthr
        p = thr_pos(i,:);  d = thr_dir(i,:);
        quiver3(ax, p(1),p(2),p(3), d(1)*scale, d(2)*scale, d(3)*scale, 0, ...
                'Color', clr_arrow, 'LineWidth', 2, 'MaxHeadSize', 0.5);
        plot3(ax, p(1), p(2), p(3), 'o', 'MarkerSize', 6, ...
              'MarkerFaceColor', clr_arrow, 'Color', clr_arrow);
    end
    nv_l = size(verts, 1);
    for v = 1:nv_l
        if ismember(v, mating_v); lbl_clr = [0.7 0.5 0.0];
        else;                     lbl_clr = [0.2 0.2 0.6]; end
        text(ax, verts(v,1)*1.08, verts(v,2)*1.08, verts(v,3)*1.08, ...
            sprintf('%s%d', label_prefix, v), 'FontSize', 8, 'Color', lbl_clr);
    end
end