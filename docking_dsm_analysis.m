%=========================================================================
% COOPERATIVE DOCKING DSM / RISK INDEX ANALYSIS
%=========================================================================
%
% PURPOSE:
%   Simulates and evaluates the safety of a cooperative docking approach
%   between two spacecraft: Asterix (the chaser) and Obelix (the target).
%   The code computes whether the chaser can always stop or abort before
%   violating spatial/angular constraints, producing a time-varying
%   "Risk Index" (RI) that flight controllers monitor during prox-ops.
%
% HIGH-LEVEL WORKFLOW:
%   1) Load thruster configurations for both vehicles from rcs_config.mat
%   2) Build force/torque "B-matrices" mapping thruster duty cycles to
%      body-frame wrenches (forces + torques)
%   3) Solve Linear Programs (LPs) to find max pure-translation braking
%      force and max pure-rotation braking torque along each axis, for
%      each vehicle independently, then sum for cooperative authority
%   4) Propagate a guided approach trajectory using the Clohessy-Wiltshire
%      (CW) relative-motion equations with a PD guidance controller
%   5) At each timestep, compute a multi-channel Risk Index:
%        RI = max(RI_axial, RI_lateral, RI_angular)
%      where each channel compares stopping distance to available margin
%   6) Project the state at docking contact to check capture envelope
%      feasibility (lateral offset, speed, attitude within tolerances)
%   7) Generate velocity-range corridor plots, risk heatmaps, and a
%      console summary with safety margins
%
% PHYSICS MODELS:
%   1) Clohessy-Wiltshire (CW) Relative Motion:
%      The CW equations linearizes the two-body problem about a circular
%      reference orbit. In the LVLH frame, the relative motion of a
%      chaser w.r.t. a target on a circular orbit of mean motion n is:
%        x_ddot - 2ny_dot - 3n²x = a_x     (radial, R)
%        y_ddot + 2nx_dot       = a_y     (in-track, I)
%        zz_ddot + n²z       = a_z    (cross-track, C)
%      where (x,y,z) = (R,I,C) and a = control acceleration.
%      The coupling between R and I arises from Coriolis and gravity-
%      gradient effects in the rotating frame. Cross-track motion (C)
%      is decoupled and is a simple harmonic oscillator at frequency n.
%      The state transition matrix Φ(Δt) gives closed-form propagation.
%
%   2) Thruster Force/Torque Model:
%      Each thruster i at position rᵢ (relative to CoM) with thrust
%      direction d̂ᵢ and magnitude F produces:
%        Force:  Fᵢ = F · d̂ᵢ              (body frame)
%        Torque: τᵢ = (rᵢ - rCoM) × Fᵢ    (body frame)
%      The B-matrix stacks these for all N thrusters:
%        [F_total; τ_total] = B · u,  u ∈ [0,1]^N
%      This is a linear map, enabling LP-based optimization.
%
%   3) LP for Maximum Braking Authority:
%      To find the maximum deceleration along direction d̂ with zero
%      parasitic torque (pure translation):
%        max  d̂^T · B_force · u
%        s.t. B_torque · u = 0    (no net torque)
%             0 ≤ u ≤ 1           (thruster duty cycle limits)
%      If infeasible (can't achieve zero torque), the constraint is
%      relaxed to |B_torque · u| ≤ τ_tol. Similarly for pure-rotation
%      torque LP with zero-force constraint.
%
%   4) Cooperative Braking:
%      During cooperative docking, both vehicles can fire thrusters to
%      arrest relative motion. The relative deceleration is:
%        a_rel = F_ast/m_ast + F_obx/m_obx
%      because each vehicle decelerates its own mass independently.
%      The chaser brakes in the approach direction; the target brakes
%      in the opposite direction (resisting the approach). This roughly
%      doubles available braking authority compared to chaser-only.
%
%   5) Stopping Distance Model:
%      Given current speed v, reaction time t_react, and braking
%      acceleration a_brake, the stopping distance is:
%        d_stop = v · t_react + v² / (2 · a_brake)
%      The first term is the coast distance during GNC latency (sensor
%      processing + valve actuation + thruster settling). The second
%      is the kinematic braking distance assuming constant deceleration.
%      This is a conservative (worst-case) model: it assumes the
%      vehicle coasts at full speed during t_react before braking begins.
%
%   6) Risk Index Channels:
%      - AXIAL:   RI = d_stop_axial / range_to_target
%                 Can the chaser stop before hitting the target?
%      - LATERAL: RI = d_stop_lateral / (corridor_wall - lateral_offset)
%                 Can the chaser stop before hitting the corridor wall?
%                 The corridor narrows linearly from far to near.
%      - ANGULAR: RI = angular_stopping_distance / (att_tolerance - att_error)
%                 Can the relative attitude be arrested before exceeding
%                 the capture envelope's attitude tolerance?
%      Total RI = max(all channels). This is conservative: the worst
%      channel governs overall safety.
%
%   7) Approach Corridor (Funnel):
%      A rectangular corridor that narrows linearly from a wide opening
%      at range_far to a tight aperture at contact. This models the
%      physical constraint that the chaser must enter the docking port
%      (small cross-section) from a wider acquisition region. The
%      half-width at range r is:
%        hw(r) = hw_close + (r/r_far) · (hw_far - hw_close)
%
%   8) Speed Profile:
%      A piecewise-linear schedule of commanded closing speed vs. range.
%      As the chaser gets closer, it slows down - this is standard
%      practice in proximity operations to maintain DSM > 1 with
%      decreasing available margin. The PD controller tracks this
%      profile, adjusting thrust to match the commanded speed.
%
%   9) Capture Envelope:
%      At the moment of physical contact (range ≈ 0), the relative
%      state must be within tight tolerances for the docking mechanism
%      to engage successfully:
%        - Lateral offset     < lat_tol      (alignment with dock port)
%        - Closing speed      < axial_vel_max (gentle contact)
%        - Lateral speed      < lat_vel_max   (no sideways scraping)
%        - Attitude error     < att_tol       (angular alignment)
%        - Attitude rate      < att_rate_max  (no tumbling at contact)
%      The code projects current state forward to contact time and
%      checks feasibility continuously.
%
%   10) Attitude Dynamics (Simplified):
%       The code uses a linearized Euler integration for relative
%       attitude (small-angle approximation):
%         ω(t+dt) = ω(t) + α_cmd · dt
%         θ(t+dt) = θ(t) + ω(t)·dt + 0.5·α_cmd·dt²
%       This is valid for the small misalignment angles (< few degrees)
%       expected during a well-controlled docking approach. Full
%       quaternion propagation (as used in the related attitude
%       dynamics simulation) is unnecessary here.
%
% INPUTS:
%   rcs_config.mat - thruster geometry, spacecraft properties for both
%   vehicles, generated by rcs_trade_study_independent.m
%
% OUTPUTS:
%   8 figures (velocity-range corridors, risk heatmaps, trajectory,
%   guidance telemetry, capture feasibility, control effort)
%   Console summary with safety margins, key events, and risk statistics
%
% RELATED SCRIPTS:
%   rcs_trade_study_independent.m  - generates rcs_config.mat
%   rcs_analyze_config_independent.m - detailed B-matrix/SVD analysis
%
%=========================================================================
clear; clc; close all;

%% ====================== USER CONFIG ====================================
% All tunable parameters are collected here for easy modification.

% --- RCS configuration file ---
rcs_mat_file = 'rcs_config.mat';

% DOCKING CORRIDOR GEOMETRY
%   The approach corridor is modeled as a rectangular funnel that narrows
%   linearly from a wide far end to a tight aperture at the docking port.
%   The chaser must remain within this corridor at all times.
%   Physically, this represents the field-of-regard of the docking sensor
%   (e.g., lidar, camera) and the geometric constraint of the docking
%   mechanism's capture latch range.
%
%   half_width_far:   lateral half-width at the far end (range = range_far)
%   half_width_close: lateral half-width at contact (range = 0), roughly
%                     matching the docking mechanism's capture cross-section
%   range_far:        distance at which the corridor begins (approach gate)
corridor.half_width_far   = 2.0;    % [m]
corridor.half_width_close = 0.10;   % [m]
corridor.range_far        = 50.0;   % [m]

% CAPTURE ENVELOPE (at contact, must be within these tolerances)
%   These are the maximum allowable relative states at the moment the
%   docking mechanism engages. Exceeding any of these causes a failed
%   capture (miss, bounce-off, or mechanism damage).
%   Values are typical for small-satellite soft-capture mechanisms.
capture.lat_tol       = 0.08;    % [m] max lateral offset at contact
capture.axial_vel_max = 0.03;    % [m/s] max closing speed at contact
capture.lat_vel_max   = 0.01;    % [m/s] max lateral rate at contact
capture.att_tol       = 2.0;     % [deg] max attitude misalignment
capture.att_rate_max  = 0.5;     % [deg/s] max relative attitude rate

% APPROACH / MISSION PARAMETERS
%   t_react: total GNC latency from sensing an off-nominal condition to
%            the thrusters producing useful force. Includes:
%              - sensor measurement + processing delay
%              - GNC algorithm execution time
%              - valve command latency + valve opening time
%              - thruster rise time to steady-state
%            This is the critical parameter in the stopping distance model:
%            during t_react, the vehicle coasts at full speed.
%
%   dv_remaining: remaining propellant budget (as delta-V). The chaser
%                 cannot expend more than (dv_remaining - dv_abort_reserve).
%
%   dv_abort_reserve: delta-V held in reserve for an abort maneuver.
%                     This is "untouchable" during normal approach - if the
%                     remaining budget drops below this, the approach is
%                     declared budget-limited and RI is artificially raised.
%
%   n_orbit: mean motion of the target orbit, n = sqrt(μ/a³).
%            For a circular orbit at altitude h above Earth:
%              a = R_earth + h, n = sqrt(μ_earth / a³)
%            At ~800 km LEO: n ≈ 1.107e-3 rad/s, period ≈ 5680 s ≈ 95 min.
%            This sets the CW dynamics timescale.
t_react     = 2.0;          % [s]
dv_remaining = 50.0;        % [m/s]
dv_abort_reserve = 10.0;    % [m/s]
n_orbit     = 1.1068e-3;    % [rad/s] (~800 km circular LEO)

% APPROACH TRAJECTORY INITIAL CONDITIONS (LVLH)
%   The chaser starts at r0 in the LVLH frame of the target (at origin).
%   V-bar approach: the chaser is "behind" the target along the in-track
%   axis and closes along +I toward the target.
%
%   r0 = [R; I; C] = [0; 50; 0] m - 50 m behind along V-bar, on-axis
%   v0 = [vR; vI; vC] = [0; -0.03; 0] m/s - closing at 3 cm/s along I
%
%   Note: LVLH I is positive in the velocity direction. The target is at
%   the origin, so the chaser at I=50 is 50 m "behind" (in the orbit
%   sense) and must close in the -I direction. The closing speed is
%   v_close = -vI = +0.03 m/s (positive = approaching).
r0 = [0; 50; 0];            % [m] LVLH position
v0 = [0; -0.03; 0];         % [m/s] LVLH velocity
t_sim = 3000;                % [s] simulation duration
dt    = 0.5;                 % [s] integration timestep

% INITIAL ATTITUDE STATE (relative: chaser w.r.t. target)
%   Small initial misalignments that the attitude controller must null
%   before docking. These represent the error remaining after the coarse
%   alignment phase that precedes the final approach.
att0      = [0.5; -0.3; 0.2];   % [deg] initial roll/pitch/yaw error
att_rate0 = [0.02; -0.01; 0.01]; % [deg/s] initial relative attitude rates

% APPROACH SPEED PROFILE
%   A piecewise-linear schedule mapping range-to-target → commanded
%   closing speed. This is the core guidance law: the chaser slows down
%   as it gets closer, maintaining a safe stopping margin at all times.
%   The profile is designed so that at each range, the stopping distance
%   (at the commanded speed) is a safe fraction of the remaining range.
%
%   At 50 m: 3.0 cm/s (acquisition speed)
%   At 30 m: 2.0 cm/s (approach speed)
%   At 15 m: 1.0 cm/s (fine approach)
%   At  5 m: 0.5 cm/s (final approach)
%   At  1 m: 0.3 cm/s (creep)
%   At  0 m: 0.2 cm/s (contact speed - within capture.axial_vel_max)
speed_profile.range = [50,  30,   15,    5,    1,    0];  % [m]
speed_profile.v_cmd = [0.03, 0.02, 0.01, 0.005, 0.003, 0.002];  % [m/s]

% GUIDANCE GAINS (PD controller)
%   Translational PD: a_cmd = -Kp·(x - x_ref) - Kd·(v - v_ref)
%     Kp [1/s²] sets the position-error-to-acceleration gain.
%     Kd [1/s]  sets the velocity-error-to-acceleration (damping) gain.
%     Larger Kp → stiffer tracking, risk of overshoot.
%     Larger Kd → more damping, slower convergence.
%     Rule of thumb: Kd ≈ 2·sqrt(Kp) for critical damping (ζ = 1).
%     Here gains are conservative (underdamped would cause oscillation
%     with the thruster quantization and saturation).
%
%   Attitude PD: α_cmd = -Kp_att·θ_err - Kd_att·ω_err
%     Same structure but in angular quantities. Units are such that
%     deg input → deg/s² output (the gains absorb unit conversion).
guid.Kp_axial  = 0.002;   % [1/s²] axial position gain
guid.Kd_axial  = 0.05;    % [1/s]   axial velocity gain
guid.Kp_lat    = 0.003;   % [1/s²] lateral position gain
guid.Kd_lat    = 0.06;    % [1/s]   lateral velocity gain
guid.Kp_att    = 0.5;     % [1/s²] attitude proportional gain
guid.Kd_att    = 1.5;     % [1/s]   attitude derivative gain

% BODY-TO-LVLH DIRECTION COSINE MATRIX (DCM)
%   C_BL transforms a vector from LVLH frame to body frame:
%     v_body = C_BL · v_lvlh
%
%   For a V-bar approach, the docking axis (body +X, where the mating
%   face is) must point toward the target. Since the chaser starts at
%   I = +50 m and the target is at the origin, the chaser's +X axis
%   must point in the LVLH -I direction.
%
%   Convention chosen here:
%     body X → LVLH -I  (docking axis toward target)
%     body Y → LVLH +R  (radial)
%     body Z → LVLH +C  (cross-track)
%
%   This makes C_BL the matrix whose rows are the body axes expressed
%   in LVLH coordinates:
%     Row 1 (body X in LVLH): [0, -1, 0]
%     Row 2 (body Y in LVLH): [1,  0, 0]
%     Row 3 (body Z in LVLH): [0,  0, 1]
C_BL = [0  -1  0;
        1   0  0;
        0   0  1];

% Docking axis unit vector in LVLH (from chaser toward target)
%   The chaser is at I=+50 and the target is at origin, so the direction
%   from chaser to target is -I = [0; -1; 0].
dock_axis_lvlh = [0; -1; 0];

% RISK INDEX THRESHOLDS
%   These define the green/yellow/red safety zones for flight-controller
%   displays. Standard practice in proximity operations.
%   GREEN:  RI < 0.5 → comfortable margin, nominal operations
%   YELLOW: 0.5 ≤ RI < 0.8 → caution, monitor closely
%   RED:    RI ≥ 0.8 → consider abort, approach limits of safe braking
%   RI ≥ 1.0: stopping distance exceeds available margin (abort trigger)
RI_green  = 0.5;
RI_yellow = 0.8;
RI_red    = 1.0;

%% ====================== LOAD RCS CONFIG ================================
% Load the thruster configuration .mat file generated by the RCS trade
% study. This contains for each vehicle:
%   - mass, inertia tensor I [3×3], dimensions (Lx, Ly, Lz)
%   - center of mass (com) position
%   - thruster positions (thr_pos, N×3), thrust directions (thr_dir, N×3)
%   - thrust magnitude, Isp, control mode ('6DOF' or 'att_only')
%   - number of thrusters
fprintf('Loading RCS configuration from: %s\n', rcs_mat_file);
if ~isfile(rcs_mat_file)
    error('RCS config file not found: %s\nRun rcs_trade_study_independent.m first.', rcs_mat_file);
end
cfg = load(rcs_mat_file);

% --- Asterix (chaser) ---
% The chaser is the active vehicle that maneuvers toward the target.
ast.name      = cfg.asterix.name;
ast.mass      = cfg.asterix.mass;         % [kg] spacecraft wet mass
ast.Lx        = cfg.asterix.Lx;           % [m] body X dimension
ast.Ly        = cfg.asterix.Ly;           % [m] body Y dimension
ast.Lz        = cfg.asterix.Lz;           % [m] body Z dimension
ast.com       = cfg.asterix.com;          % [m] center of mass in body frame
ast.I         = cfg.asterix.I;            % [kg·m²] 3×3 inertia tensor about CoM
ast.thrust    = cfg.asterix.thrust;       % [N] single thruster force magnitude
ast.Isp       = cfg.asterix.Isp;          % [s] specific impulse
ast.thr_pos   = cfg.asterix.thr_pos;      % [m] N×3 thruster positions (body frame)
ast.thr_dir   = cfg.asterix.thr_dir;      % [-] N×3 thruster unit direction vectors
ast.ctrl_mode = cfg.asterix.ctrl_mode;    % '6DOF' or 'att_only'
ast.n_thr     = cfg.asterix.n_thrusters;  % number of thrusters

% --- Obelix (target, cooperative) ---
% The target is station-keeping at the LVLH origin. In a cooperative
% docking scenario, the target also fires thrusters to help brake the
% relative motion (unlike a passive target where only the chaser brakes).
obx.name      = cfg.obelix.name;
obx.mass      = cfg.obelix.mass;
obx.Lx        = cfg.obelix.Lx;
obx.Ly        = cfg.obelix.Ly;
obx.Lz        = cfg.obelix.Lz;
obx.com       = cfg.obelix.com;
obx.I         = cfg.obelix.I;
obx.thrust    = cfg.obelix.thrust;
obx.Isp       = cfg.obelix.Isp;
obx.thr_pos   = cfg.obelix.thr_pos;
obx.thr_dir   = cfg.obelix.thr_dir;
obx.ctrl_mode = cfg.obelix.ctrl_mode;
obx.n_thr     = cfg.obelix.n_thrusters;
obx.geom_center_global = cfg.obelix.geom_center_global;

fprintf('  Asterix (chaser): %s mode, %d thr, %.1f kg, [%.2f x %.2f x %.2f] m\n', ...
        ast.ctrl_mode, ast.n_thr, ast.mass, ast.Lx, ast.Ly, ast.Lz);
fprintf('  Obelix (target):  %s mode, %d thr, %.1f kg, [%.2f x %.2f x %.2f] m\n', ...
        obx.ctrl_mode, obx.n_thr, obx.mass, obx.Lx, obx.Ly, obx.Lz);
if isfield(cfg, 'timestamp')
    fprintf('  Config generated: %s\n', cfg.timestamp);
end

fprintf('\n=== Cooperative Docking DSM Analysis: %s -> %s ===\n\n', ...
        ast.name, obx.name);

%% ====================== BUILD B-MATRICES (EACH VEHICLE) ================
% The B-matrix (also called the "thruster influence matrix" or "control
% effectiveness matrix") maps thruster duty cycles to body-frame wrenches.
%
% For N thrusters, each with force magnitude F and direction d̂ᵢ at
% position rᵢ relative to the center of mass:
%
%   B = [ F·d̂₁  F·d̂₂  ...  F·d̂ₙ ]   ← 3×N force partition
%       [ r₁×(F·d̂₁)  ...  rₙ×(F·d̂ₙ) ] ← 3×N torque partition
%
% So B is 6×N, and the wrench is:
%   [F_total; τ_total] = B · u,   u ∈ [0,1]^N
%
% Each column of B represents the wrench produced by one thruster at
% full duty cycle. The top 3 rows give force, bottom 3 give torque.
% Off-diagonal inertia terms in the torque don't appear here because
% the torque is computed geometrically (r × F), not dynamically -
% the inertia tensor is used later when converting torque to angular
% acceleration.

    function B = build_B_matrix(thr_pos, thr_dir, com, thrust_mag)
        % Inputs:
        %   thr_pos:    N×3 thruster positions in body frame [m]
        %   thr_dir:    N×3 unit thrust direction vectors (direction of
        %               force on spacecraft, not exhaust direction)
        %   com:        1×3 center of mass in body frame [m]
        %   thrust_mag: scalar thrust magnitude per thruster [N]
        %
        % Output:
        %   B: 6×N control effectiveness matrix
        N = size(thr_pos, 1);
        B = zeros(6, N);
        for i = 1:N
            r_i = thr_pos(i,:)' - com';       % moment arm: thruster pos relative to CoM
            f_i = thrust_mag * thr_dir(i,:)';  % force vector from this thruster
            B(1:3, i) = f_i;                   % force contribution
            B(4:6, i) = cross(r_i, f_i);       % torque contribution (r × F)
        end
    end

B_ast = build_B_matrix(ast.thr_pos, ast.thr_dir, ast.com, ast.thrust);
B_obx = build_B_matrix(obx.thr_pos, obx.thr_dir, obx.com, obx.thrust);

fprintf('Asterix B-matrix: %d thrusters, %.2f N each\n', ast.n_thr, ast.thrust);
fprintf('Obelix  B-matrix: %d thrusters, %.2f N each\n', obx.n_thr, obx.thrust);

%% ====================== LP SOLVER (FORCE-ONLY AND TORQUE-ONLY) =========
% These functions solve linear programs to find the maximum force (or
% torque) achievable along a specified direction, subject to:
%   1) Thruster duty cycles u ∈ [0, 1] (can't fire below 0 or above 100%)
%   2) Zero parasitic wrench in the complementary partition:
%      - For force LP: zero net torque (pure translation, no tumbling)
%      - For torque LP: zero net force (pure rotation, no drift)
%
% The LP formulation (for force):
%   maximize   d̂ᵀ · B_force · u       (force projection along desired direction)
%   subject to B_torque · u = 0        (zero torque)
%              0 ≤ u ≤ 1               (duty cycle bounds)
%
% In linprog, we minimize the negative: min -(d̂ᵀ · B_force)ᵀ · u
%
% If the equality constraint is infeasible (the thruster geometry cannot
% produce the desired force with exactly zero torque), the constraint is
% relaxed to an inequality: |B_torque · u| ≤ τ_tol. This trades a small
% residual torque for usable translational authority - the attitude
% controller can compensate for small parasitic torques.

    function [f_max, a_max, exitflag] = solve_force_lp(B, mass, d_body)
        % Find maximum pure-translation force along direction d_body.
        %
        % Inputs:
        %   B:      6×N B-matrix for this vehicle
        %   mass:   vehicle mass [kg]
        %   d_body: 3×1 desired force direction (body frame, unit vector)
        %
        % Outputs:
        %   f_max:    max force along d_body [N]
        %   a_max:    corresponding acceleration [m/s²]
        %   exitflag: linprog exit status (1 = success)
        
        N = size(B, 2);          % number of thrusters
        B_f = B(1:3, :);         % force partition of B-matrix
        B_t = B(4:6, :);         % torque partition of B-matrix
        
        % Objective: maximize d̂ᵀ·B_f·u → minimize -(d̂ᵀ·B_f)ᵀ·u
        f_obj = -(d_body' * B_f)';
        lb = zeros(N, 1);        % lower bound: can't fire negative
        ub = ones(N, 1);         % upper bound: 100% duty cycle max
        opts = optimoptions('linprog', 'Display', 'off', 'Algorithm', 'dual-simplex');
        
        % First attempt: strict zero-torque equality constraint
        % B_torque · u = 0 (no parasitic torque at all)
        [u_opt, ~, exitflag] = linprog(f_obj, [], [], B_t, zeros(3,1), lb, ub, opts);
        
        if exitflag ~= 1
            % Infeasible with zero torque - relax to bounded torque.
            % Allow up to τ_tol of residual torque in each axis.
            % This is formulated as: -τ_tol ≤ B_t·u ≤ τ_tol
            % → [B_t; -B_t]·u ≤ [τ_tol; τ_tol]
            tau_tol = 0.5;  % [N·m] allowable residual torque
            A_ineq = [B_t; -B_t];
            b_ineq = tau_tol * ones(6, 1);
            [u_opt, ~, exitflag] = linprog(f_obj, A_ineq, b_ineq, [], [], lb, ub, opts);
        end
        
        if exitflag == 1
            F_vec = B_f * u_opt;           % resulting force vector [N]
            f_max = d_body' * F_vec;       % projection along desired direction
            a_max = f_max / mass;          % acceleration = F/m
        else
            f_max = 0;  % no feasible solution found
            a_max = 0;
        end
    end

    function [tau_max, alpha_max, exitflag] = solve_torque_lp(B, I_tensor, d_body)
        % Find maximum pure-rotation torque about axis d_body.
        %
        % Inputs:
        %   B:        6×N B-matrix
        %   I_tensor: 3×3 inertia tensor [kg·m²]
        %   d_body:   3×1 desired torque axis (body frame, unit vector)
        %
        % Outputs:
        %   tau_max:   max torque about d_body [N·m]
        %   alpha_max: corresponding angular acceleration [rad/s²]
        %   exitflag:  linprog exit status
        
        N = size(B, 2);
        B_f = B(1:3, :);
        B_t = B(4:6, :);
        
        % Objective: maximize d̂ᵀ·B_t·u
        f_obj = -(d_body' * B_t)';
        lb = zeros(N, 1);
        ub = ones(N, 1);
        opts = optimoptions('linprog', 'Display', 'off', 'Algorithm', 'dual-simplex');
        
        % Constraint: zero net force (B_f · u = 0) for pure rotation
        [u_opt, ~, exitflag] = linprog(f_obj, [], [], B_f, zeros(3,1), lb, ub, opts);
        
        if exitflag ~= 1
            % Relax: allow small residual force
            f_tol = 0.05;  % [N] allowable residual force
            A_ineq = [B_f; -B_f];
            b_ineq = f_tol * ones(6, 1);
            [u_opt, ~, exitflag] = linprog(f_obj, A_ineq, b_ineq, [], [], lb, ub, opts);
        end
        
        if exitflag == 1
            tau_vec = B_t * u_opt;              % resulting torque vector [N·m]
            tau_max = d_body' * tau_vec;         % projection along desired axis
            % Angular acceleration: τ = I·α → α = τ/I_axis
            % For rotation about axis d̂, the effective moment of inertia is:
            %   I_axis = d̂ᵀ · I · d̂
            % This is the scalar moment of inertia about that axis.
            I_axis = d_body' * I_tensor * d_body;
            alpha_max = tau_max / I_axis;        % [rad/s²]
        else
            tau_max = 0;
            alpha_max = 0;
        end
    end

%% ====================== TRANSLATIONAL BRAKING (DUAL-VEHICLE LP) ========
% For each of the 6 LVLH directions (±R, ±I, ±C), find the maximum
% braking acceleration available when both vehicles fire cooperatively.
%
% Procedure:
%   1) Convert the LVLH direction to body frame using C_BL
%   2) Solve the force LP for Asterix in that body-frame direction
%   3) Solve the force LP for Obelix in the OPPOSITE body-frame
%      direction (Obelix brakes against the approach, so it fires
%      thrusters in the opposite sense)
%   4) Sum the accelerations: a_rel = a_ast + a_obx
%      (each vehicle decelerates its own mass, so relative deceleration
%      is the sum of individual decelerations)
%
% The result is a table of max braking accelerations per LVLH axis,
% which is then used in the RI computation.

dir_labels = {'+R','-R','+I','-I','+C','-C'};
dir_vecs_lvlh = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1]';

fprintf('\n--- Per-axis max pure-translation braking force (dual-vehicle LP) ---\n');
fprintf('  Dir   |  Asterix [N]  Obelix [N]  Combined [N]  Accel [mm/s^2]\n');
fprintf('  ------|------------------------------------------------------------\n');

max_accel_ast = zeros(6,1);
max_accel_obx = zeros(6,1);
max_accel_combined = zeros(6,1);

% Masses for acceleration computation
m_ast = ast.mass;
m_obx = obx.mass;

for k = 1:6
    d_lvlh = dir_vecs_lvlh(:, k);
    
    % Asterix: convert desired LVLH braking direction to body frame
    d_body_ast = C_BL * d_lvlh;
    [f_ast, ~, ef1] = solve_force_lp(B_ast, m_ast, d_body_ast);
    a_ast = f_ast / m_ast;
    
    % Obelix: brakes in the OPPOSITE LVLH direction.
    % If the relative velocity is in the +I direction, Asterix brakes
    % by thrusting in +I (opposing its own motion), while Obelix brakes
    % by thrusting in -I (opposing the chaser's approach from its
    % perspective). Both contribute to reducing relative velocity.
    d_body_obx = C_BL * (-d_lvlh);
    [f_obx, ~, ef2] = solve_force_lp(B_obx, m_obx, d_body_obx);
    a_obx = f_obx / m_obx;
    
    % Combined relative deceleration (both positive contributions)
    max_accel_ast(k) = max(a_ast, 0);
    max_accel_obx(k) = max(a_obx, 0);
    max_accel_combined(k) = max_accel_ast(k) + max_accel_obx(k);
    
    fprintf('  %s   |  %8.4f     %8.4f     %8.4f      %8.4f\n', ...
            dir_labels{k}, f_ast, f_obx, ...
            f_ast + f_obx, max_accel_combined(k)*1e3);
end

% Map to a struct for easy lookup during simulation.
% Convention: "Rpos" means braking authority for positive R velocity
% (i.e., if moving in +R, we need to brake in -R, so we use the -R
% LP result, which is index 2 in our direction list).
a_brake = struct( ...
    'Rpos', max_accel_combined(2), 'Rneg', max_accel_combined(1), ...
    'Ipos', max_accel_combined(4), 'Ineg', max_accel_combined(3), ...
    'Cpos', max_accel_combined(6), 'Cneg', max_accel_combined(5));

% Chaser-only braking (for abort analysis: if target becomes
% non-cooperative, only Asterix can brake)
a_brake_chaser = struct( ...
    'Rpos', max_accel_ast(2), 'Rneg', max_accel_ast(1), ...
    'Ipos', max_accel_ast(4), 'Ineg', max_accel_ast(3), ...
    'Cpos', max_accel_ast(6), 'Cneg', max_accel_ast(5));

fprintf('\nCooperative braking accelerations [mm/s^2]:\n');
fprintf('  R-axis: +vel -> %.4f,  -vel -> %.4f\n', a_brake.Rpos*1e3, a_brake.Rneg*1e3);
fprintf('  I-axis: +vel -> %.4f,  -vel -> %.4f\n', a_brake.Ipos*1e3, a_brake.Ineg*1e3);
fprintf('  C-axis: +vel -> %.4f,  -vel -> %.4f\n', a_brake.Cpos*1e3, a_brake.Cneg*1e3);

fprintf('\nChaser-only braking accelerations [mm/s^2] (abort/non-cooperative):\n');
fprintf('  R-axis: +vel -> %.4f,  -vel -> %.4f\n', a_brake_chaser.Rpos*1e3, a_brake_chaser.Rneg*1e3);
fprintf('  I-axis: +vel -> %.4f,  -vel -> %.4f\n', a_brake_chaser.Ipos*1e3, a_brake_chaser.Ineg*1e3);
fprintf('  C-axis: +vel -> %.4f,  -vel -> %.4f\n', a_brake_chaser.Cpos*1e3, a_brake_chaser.Cneg*1e3);

%% ====================== ANGULAR BRAKING (TORQUE-ONLY LP) ===============
% Same approach as translational, but for rotational axes.
% For each body axis (X=roll, Y=pitch, Z=yaw), find the maximum pure-
% rotation angular acceleration (with zero net force) for each vehicle.
%
% Both vehicles contribute to relative attitude correction: if the
% chaser has a +2° roll error, both vehicles can fire thrusters to
% produce roll torque (chaser torques itself one way, target torques
% itself the other way, both reducing the relative error).
%
% The relative angular acceleration is:
%   α_rel = τ_ast/I_ast + τ_obx/I_obx
% (analogous to translational cooperative braking)

body_axes = eye(3);  % columns: body X, Y, Z
axis_names = {'Roll (X)', 'Pitch (Y)', 'Yaw (Z)'};

alpha_max_ast = zeros(3,2);       % [+direction, -direction] per axis
alpha_max_obx = zeros(3,2);
alpha_max_combined = zeros(3,2);

fprintf('\n--- Per-axis max pure-rotation angular accel (dual-vehicle LP) ---\n');
fprintf('  Axis      | Asterix [deg/s^2]  Obelix [deg/s^2]  Combined [deg/s^2]\n');
fprintf('  ----------|-----------------------------------------------------------\n');

for ax = 1:3
    for sgn = 1:2
        d = body_axes(:, ax) * (3 - 2*sgn);  % sgn=1 → +1, sgn=2 → -1
        
        [~, alpha_a, ~] = solve_torque_lp(B_ast, ast.I, d);
        [~, alpha_o, ~] = solve_torque_lp(B_obx, obx.I, d);
        
        alpha_max_ast(ax, sgn) = max(alpha_a, 0);
        alpha_max_obx(ax, sgn) = max(alpha_o, 0);
        alpha_max_combined(ax, sgn) = alpha_max_ast(ax, sgn) + alpha_max_obx(ax, sgn);
    end
    
    fprintf('  %s+ |  %10.4f          %10.4f          %10.4f\n', ...
            axis_names{ax}, rad2deg(alpha_max_ast(ax,1)), ...
            rad2deg(alpha_max_obx(ax,1)), rad2deg(alpha_max_combined(ax,1)));
    fprintf('  %s- |  %10.4f          %10.4f          %10.4f\n', ...
            axis_names{ax}, rad2deg(alpha_max_ast(ax,2)), ...
            rad2deg(alpha_max_obx(ax,2)), rad2deg(alpha_max_combined(ax,2)));
end

% Worst-direction angular braking per axis:
% Take the minimum of +/- authority for each axis, since the RI
% computation must assume the worst case (we don't know which direction
% the error will be when we need to brake).
alpha_brake = zeros(3,1);
for ax = 1:3
    alpha_brake(ax) = min(alpha_max_combined(ax,:));
    alpha_brake(ax) = max(alpha_brake(ax), 1e-12);  % avoid divide-by-zero
end

fprintf('\nWorst-case angular braking [deg/s^2]: Roll=%.4f  Pitch=%.4f  Yaw=%.4f\n', ...
        rad2deg(alpha_brake(1)), rad2deg(alpha_brake(2)), rad2deg(alpha_brake(3)));

%% ====================== GUIDED CW PROPAGATION ===========================
% This section simulates the closed-loop approach trajectory.
%
% The CW (Clohessy-Wiltshire) State Transition Matrix propagates the
% relative state [r; v] forward in time under the linearized dynamics
% of two bodies in nearby circular orbits. The STM accounts for:
%   - Gravity gradient (tidal) forces (3n²x term in radial equation)
%   - Coriolis coupling between radial and in-track motion
%   - Simple harmonic oscillation in cross-track
%
% On top of the natural CW dynamics, the PD guidance controller
% computes a commanded acceleration, which is applied as an impulsive
% velocity change (dv = a·dt) each timestep. This is a good
% approximation when dt is small compared to the thruster on-time
% and the orbital period.

n = n_orbit;

% Clohessy-Wiltshire State Transition Matrix Φ(Δt)
% Maps state x(t) to x(t+Δt) under free-drift (no thrust):
%   x(t+Δt) = Φ(Δt) · x(t)
% where x = [R; I; C; vR; vI; vC]
%
% The CW STM is derived by solving the linearized Hill equations:
%   ẍ - 2nẏ - 3n²x = 0   (radial)
%   ÿ + 2nẋ       = 0   (in-track)
%   z̈ + n²z       = 0   (cross-track)
% The solution involves sin(nt), cos(nt) terms and their derivatives.
    function Phi = CW_STM(n, dt)
        s = sin(n*dt); c = cos(n*dt);
        Phi = [
            % Position rows (top 3): how position at t maps to position at t+dt
            4-3*c,      0, 0,   s/n,        2*(1-c)/n,     0;       % R(t+dt)
            6*(s-n*dt), 1, 0,  -2*(1-c)/n,  (4*s-3*n*dt)/n, 0;      % I(t+dt): includes secular drift (n*dt term)
            0,          0, c,   0,           0,              s/n;    % C(t+dt): decoupled SHM
            % Velocity rows (bottom 3):
            3*n*s,      0, 0,   c,           2*s,            0;      % vR(t+dt)
            -6*n*(1-c), 0, 0,  -2*s,         4*c-3,          0;     % vI(t+dt): note the -3 term (secular drift)
            0,          0, -n*s, 0,           0,              c      % vC(t+dt): decoupled SHM
        ];
    end

% Speed profile interpolator
%   Given the current axial range, return the commanded closing speed
%   from the piecewise-linear speed schedule. Clamps to the schedule
%   endpoints if range is outside the defined breakpoints.
    function v_cmd = speed_cmd(range, sp)
        v_cmd = interp1(sp.range, sp.v_cmd, ...
                        min(max(range, sp.range(end)), sp.range(1)), ...
                        'linear', sp.v_cmd(1));
    end

% Acceleration saturation
%   Clips commanded acceleration to the available thruster authority.
%   Positive and negative limits may differ (asymmetric thruster layout).
    function a_sat = saturate_accel(a_cmd, a_max_pos, a_max_neg)
        if a_cmd > 0
            a_sat = min(a_cmd, a_max_pos);
        else
            a_sat = max(a_cmd, -a_max_neg);
        end
    end

% --- Initialize state arrays ---
t_vec = 0:dt:t_sim;
N_steps = length(t_vec);
state = zeros(6, N_steps);       % [R; I; C; vR; vI; vC] in LVLH
state(:,1) = [r0; v0];

% Attitude states (small-angle, linearized)
att_state  = zeros(3, N_steps);   % [roll; pitch; yaw] misalignment [deg]
rate_state = zeros(3, N_steps);   % [p; q; r] relative attitude rates [deg/s]
att_state(:,1)  = att0;
rate_state(:,1) = att_rate0;

% Telemetry logging
accel_cmd   = zeros(3, N_steps);   % commanded translational accel [m/s²]
alpha_cmd   = zeros(3, N_steps);   % commanded angular accel [deg/s²]
dv_used     = zeros(1, N_steps);   % cumulative delta-V expended [m/s]
v_cmd_log   = zeros(1, N_steps);   % commanded axial closing speed [m/s]

% Cooperative translational authority struct (from LP results)
a_max_trans = struct( ...
    'Rpos', a_brake.Rpos, 'Rneg', a_brake.Rneg, ...
    'Ipos', a_brake.Ipos, 'Ineg', a_brake.Ineg, ...
    'Cpos', a_brake.Cpos, 'Cneg', a_brake.Cneg);

% Angular authority in deg/s² (from torque LP)
alpha_max = rad2deg(alpha_brake);  % [deg/s²] per body axis

docked = false;
dock_time = NaN;

% === MAIN SIMULATION LOOP ===
% At each timestep:
%   1) Decompose state into axial (along docking axis) and lateral components
%   2) Check if docking contact has occurred
%   3) Compute translational guidance commands (PD on speed profile + lateral centering)
%   4) Compute attitude guidance commands (PD to null misalignment)
%   5) Propagate translation via CW STM + impulsive delta-V
%   6) Propagate attitude via Euler integration
%   7) Accumulate delta-V expenditure

for k = 2:N_steps
    % --- Decompose current state ---
    r_vec = state(1:3, k-1);    % position in LVLH [m]
    v_vec = state(4:6, k-1);    % velocity in LVLH [m/s]
    d_hat = dock_axis_lvlh / norm(dock_axis_lvlh);  % docking axis unit vec
    
    % Axial range: projection of position onto docking axis.
    % The chaser starts at I=+50 and the target is at origin.
    % dock_axis_lvlh = [0;-1;0], so:
    %   range_now = -dot([0;50;0], [0;-1;0]) = -(-50) = +50 (correct)
    range_now = -dot(r_vec, d_hat);
    
    % Closing speed: rate of decrease of range.
    % v_close > 0 means the chaser is approaching the target.
    v_close   = -dot(v_vec, d_hat);
    
    % Lateral offset: position component perpendicular to docking axis.
    % This should stay near zero for a centered approach.
    r_lat = r_vec - dot(r_vec, d_hat) * d_hat;
    v_lat_vec = v_vec - dot(v_vec, d_hat) * d_hat;
    
    % --- Check for docking contact ---
    % Contact occurs when range drops below 0.05 m and the chaser is
    % approaching (not receding) at less than 1.5× the capture speed limit.
    if range_now <= 0.05 && v_close >= 0 && v_close <= capture.axial_vel_max * 1.5
        if ~docked
            docked = true;
            dock_time = t_vec(k-1);
            fprintf('  DOCKING at t = %.0f s, range = %.3f m, v = %.4f m/s\n', ...
                    dock_time, range_now, v_close);
        end
        % After docking, hold state constant (no further dynamics)
        state(:,k) = state(:,k-1);
        att_state(:,k) = att_state(:,k-1);
        rate_state(:,k) = rate_state(:,k-1);
        dv_used(k) = dv_used(k-1);
        continue
    end
    
    % === TRANSLATIONAL GUIDANCE (PD controller) ===
    
    % AXIAL CONTROL: Track the speed profile.
    % Look up the commanded closing speed for the current range.
    v_target = speed_cmd(range_now, speed_profile);
    v_cmd_log(k-1) = v_target;
    
    % Axial guidance: velocity-tracking PD
    %   a_cmd = Kd · (v_target - v_close)
    % No position term - axial guidance follows the speed schedule,
    % not a position reference. This is a "speed-scheduled approach":
    % the position evolves naturally from integrating the speed profile.
    a_cmd_axial = guid.Kd_axial * (v_target - v_close);
    
    % Map axial command to LVLH: closing is along -d_hat direction,
    % so a positive axial command (speed up closing) maps to -d_hat.
    a_lvlh_axial = -a_cmd_axial * d_hat;
    
    % LATERAL CONTROL: PD to drive lateral offset and rate to zero.
    % Reference is the corridor centerline (zero lateral offset, zero
    % lateral rate). This keeps the chaser on the docking axis.
    %   a_cmd = -Kp·r_lat - Kd·v_lat
    a_lvlh_lat = -guid.Kp_lat * r_lat - guid.Kd_lat * v_lat_vec;
    
    % Total translational acceleration command (LVLH)
    a_cmd_lvlh = a_lvlh_axial + a_lvlh_lat;
    
    % Saturate to available thruster authority (from LP results).
    % Each LVLH axis has independent positive/negative limits because
    % the thruster geometry may be asymmetric.
    a_cmd_lvlh(1) = saturate_accel(a_cmd_lvlh(1), a_max_trans.Rpos, a_max_trans.Rneg);
    a_cmd_lvlh(2) = saturate_accel(a_cmd_lvlh(2), a_max_trans.Ipos, a_max_trans.Ineg);
    a_cmd_lvlh(3) = saturate_accel(a_cmd_lvlh(3), a_max_trans.Cpos, a_max_trans.Cneg);
    
    accel_cmd(:, k-1) = a_cmd_lvlh;
    
    % === ATTITUDE GUIDANCE (PD controller) ===
    % Drive relative attitude error and rates to zero.
    % Target attitude = zero misalignment (perfect alignment for docking).
    %   α_cmd = -Kp_att · θ_err - Kd_att · ω_err
    att_err  = att_state(:, k-1);    % [deg] current attitude error
    att_rate = rate_state(:, k-1);   % [deg/s] current attitude rates
    
    alpha_cmd_raw = -guid.Kp_att * att_err - guid.Kd_att * att_rate;
    
    % Saturate per-axis angular acceleration to available torque authority
    for ax = 1:3
        alpha_cmd_raw(ax) = max(min(alpha_cmd_raw(ax), alpha_max(ax)), -alpha_max(ax));
    end
    alpha_cmd(:, k-1) = alpha_cmd_raw;
    
    % === PROPAGATE TRANSLATION ===
    % Two-step process:
    %   1) CW free-drift: x_drift = Φ(dt) · x(k-1)
    %      This accounts for natural orbital mechanics (gravity gradient,
    %      Coriolis coupling) during the timestep.
    %   2) Apply control impulse: Δv = a_cmd · dt
    %      Added to the velocity after free-drift propagation.
    %      This is equivalent to a constant-acceleration-over-dt
    %      approximation (impulse at the end of the step).
    Phi = CW_STM(n, dt);
    state_drift = Phi * state(:, k-1);    % free-drift state
    dv_step = a_cmd_lvlh * dt;            % impulsive Δv from controller
    state(:, k) = state_drift;
    state(4:6, k) = state(4:6, k) + dv_step;  % add control Δv to velocity
    
    % === PROPAGATE ATTITUDE ===
    % Simple kinematic integration (small-angle, linearized):
    %   ω(k) = ω(k-1) + α·dt               (angular rate update)
    %   θ(k) = θ(k-1) + ω(k-1)·dt + ½α·dt²  (angle update, 2nd-order)
    % This is a second-order Taylor expansion of the rotation.
    % Valid for small angles (< ~5°) and short timesteps.
    rate_state(:, k) = rate_state(:, k-1) + alpha_cmd_raw * dt;
    att_state(:, k)  = att_state(:, k-1) + rate_state(:, k-1) * dt + ...
                       0.5 * alpha_cmd_raw * dt^2;
    
    % === ACCUMULATE DELTA-V ===
    % Total Δv expended so far. This tracks propellant consumption.
    % norm(dv_step) gives the magnitude regardless of direction.
    dv_used(k) = dv_used(k-1) + norm(dv_step);
end

% Fill last entry of telemetry logs
v_cmd_log(end) = v_cmd_log(end-1);
accel_cmd(:, end) = accel_cmd(:, end-1);

% Extract position and velocity component time histories for plotting
pos_R = state(1,:);  pos_I = state(2,:);  pos_C = state(3,:);
vel_R = state(4,:);  vel_I = state(5,:);  vel_C = state(6,:);

fprintf('\nGuided approach simulation complete.\n');
fprintf('  Duration: %.0f s,  Delta-V used: %.4f m/s\n', t_sim, dv_used(end));
if docked
    fprintf('  DOCKED at t = %.0f s\n', dock_time);
else
    fprintf('  Did not reach docking port in %.0f s\n', t_sim);
end

%% ====================== APPROACH CORRIDOR GEOMETRY =====================
% Compute corridor half-width as a function of axial range.
% The corridor narrows linearly: wide far away, tight at the docking port.

    function hw = corridor_halfwidth(range, corridor)
        % Linear interpolation between far and close half-widths.
        % frac = 0 at contact (range=0), frac = 1 at range_far.
        frac = min(max(range / corridor.range_far, 0), 1);
        hw = corridor.half_width_close + ...
             frac * (corridor.half_width_far - corridor.half_width_close);
    end

% Decompose the full trajectory into axial and lateral components
% relative to the docking axis.
d_hat = dock_axis_lvlh / norm(dock_axis_lvlh);

range_axial = zeros(1, N_steps);
offset_lat  = zeros(1, N_steps);
vel_axial   = zeros(1, N_steps);
vel_lat     = zeros(1, N_steps);

for k = 1:N_steps
    r_vec = state(1:3, k);
    v_vec = state(4:6, k);
    
    range_axial(k) = -dot(r_vec, d_hat);   % positive = chaser behind target
    r_lat = r_vec - dot(r_vec, d_hat) * d_hat;
    offset_lat(k)  = norm(r_lat);           % magnitude of lateral offset
    
    vel_axial(k) = -dot(v_vec, d_hat);      % positive = closing
    v_lat = v_vec - dot(v_vec, d_hat) * d_hat;
    vel_lat(k) = norm(v_lat);               % magnitude of lateral velocity
end

%% ====================== MULTI-CHANNEL DSM & RI =========================
% Compute the Risk Index at each timestep across three safety channels.
%
% CHANNEL 1 - AXIAL (along docking axis):
%   RI_axial = d_stop_axial / range_to_target
%   "Can the chaser stop before hitting the target?"
%   d_stop = v_close · t_react + v_close² / (2·a_brake)
%
% CHANNEL 2 - LATERAL (perpendicular to docking axis):
%   RI_lateral = d_stop_lateral / lateral_margin
%   "Can the chaser stop before hitting the corridor wall?"
%   lateral_margin = corridor_halfwidth - current_lateral_offset
%   Lateral stopping distance uses the RSS (root sum square) of R and C
%   components, each with their own braking authority.
%
% CHANNEL 3 - ATTITUDE (rotational):
%   RI_angular = angular_stopping_distance / angular_margin
%   "Can the relative attitude be arrested before exceeding tolerance?"
%   angular_margin = attitude_tolerance - current_attitude_error
%   Angular stopping distance: same kinematic model but with angular
%   quantities (ω instead of v, α instead of a, θ instead of d).
%
% TOTAL RI = max(RI_axial, RI_lateral, RI_angular)
%   The worst channel governs. This is conservative: the system is only
%   as safe as its most-constrained degree of freedom.
%
% DSM = 1/RI (Dynamic Safety Margin, the traditional metric).
%   DSM > 1: can stop in time. DSM < 1: cannot stop - abort required.

RI_axial   = zeros(1, N_steps);
RI_lateral = zeros(1, N_steps);
RI_angular = zeros(1, N_steps);
RI_total   = zeros(1, N_steps);
DSM_total  = zeros(1, N_steps);
hw_vec     = zeros(1, N_steps);

% Map docking axis (LVLH I for V-bar approach) to the appropriate
% braking authority from the LP results.
a_axial_pos  = a_brake.Ipos;    % braking for positive closing velocity
a_axial_neg  = a_brake.Ineg;
a_lat_R_pos  = a_brake.Rpos;    % lateral braking in R axis
a_lat_R_neg  = a_brake.Rneg;
a_lat_C_pos  = a_brake.Cpos;    % lateral braking in C axis
a_lat_C_neg  = a_brake.Cneg;

for k = 1:N_steps
    rng = range_axial(k);
    
    if rng <= 0
        % Past the target - RI is infinite (undefined/unsafe)
        RI_axial(k)  = Inf;
        RI_lateral(k) = Inf;
        RI_angular(k) = Inf;
        RI_total(k)  = Inf;
        DSM_total(k) = 0;
        continue
    end
    
    % --- Channel 1: AXIAL RI ---
    v_close = vel_axial(k);
    if v_close > 0
        % Chaser is approaching: compute stopping distance
        % d_stop = v·t_react + v²/(2a)
        %   First term: coast during reaction time (worst case)
        %   Second term: kinematic braking distance at constant deceleration
        a_ax = max(a_axial_pos, 1e-12);  % avoid division by zero
        d_stop_ax = v_close * t_react + v_close^2 / (2 * a_ax);
        RI_axial(k) = d_stop_ax / rng;   % ratio of stopping dist to range
    else
        RI_axial(k) = 0;  % receding from target - no axial risk
    end
    
    % --- Channel 2: LATERAL RI ---
    % The available lateral margin is the corridor half-width minus the
    % current lateral offset. This shrinks as the chaser approaches
    % (corridor narrows) and/or as lateral offset increases.
    hw = corridor_halfwidth(rng, corridor);
    hw_vec(k) = hw;
    lat_margin = hw - offset_lat(k);
    lat_margin = max(lat_margin, 1e-6);  % clamp to avoid division by zero
    
    % Compute lateral stopping distance per LVLH axis (R and C),
    % then take the RSS to get the total lateral stopping distance.
    vR = vel_R(k);  vC = vel_C(k);
    aR = max((vR > 0) * a_lat_R_pos + (vR <= 0) * a_lat_R_neg, 1e-12);
    aC = max((vC > 0) * a_lat_C_pos + (vC <= 0) * a_lat_C_neg, 1e-12);
    
    d_stop_R = abs(vR) * t_react + vR^2 / (2 * aR);
    d_stop_C = abs(vC) * t_react + vC^2 / (2 * aC);
    d_stop_lat = sqrt(d_stop_R^2 + d_stop_C^2);  % RSS of R and C components
    
    RI_lateral(k) = d_stop_lat / lat_margin;
    
    % --- Channel 3: ATTITUDE RI ---
    % Angular stopping distance: ω·t_react + ω²/(2α)
    % Computed per body axis, then RSS'd.
    att_err   = att_state(:, k);     % [deg] current misalignment
    att_rates = rate_state(:, k);    % [deg/s] current attitude rates
    
    att_tol_rad = capture.att_tol;   % [deg] allowed misalignment
    
    ang_stop = zeros(3,1);
    for ax = 1:3
        omega = abs(att_rates(ax)) * pi/180;  % convert to rad/s for computation
        alpha = alpha_brake(ax);               % [rad/s²] braking authority
        ang_stop_rad = omega * t_react + omega^2 / (2 * alpha);
        ang_stop(ax) = ang_stop_rad * 180/pi;  % convert back to degrees
    end
    
    % Angular margin: tolerance minus current error magnitude
    att_err_mag = norm(att_err);
    att_margin = max(att_tol_rad - att_err_mag, 1e-6);
    ang_stop_total = norm(ang_stop);  % RSS of per-axis angular stopping distances
    
    RI_angular(k) = ang_stop_total / att_margin;
    
    % --- Combined RI: worst channel governs ---
    RI_total(k) = max([RI_axial(k), RI_lateral(k), RI_angular(k)]);
    DSM_total(k) = 1 / max(RI_total(k), 1e-12);  % DSM = 1/RI
    
    % --- Delta-V budget check ---
    % If the delta-V needed to stop exceeds the remaining budget
    % (minus abort reserve), artificially inflate RI to signal
    % a budget-limited situation. This forces an abort before
    % propellant is depleted.
    dv_stop_est = v_close;  % rough estimate: dv to stop ≈ current speed
    if dv_stop_est > (dv_remaining - dv_abort_reserve)
        DSM_total(k) = min(DSM_total(k), 0.1);
        RI_total(k) = max(RI_total(k), 10);
    end
end

%% ====================== CAPTURE ENVELOPE PROJECTION ====================
% At each timestep, project the current state forward to the moment of
% contact (assuming current rates persist unchanged) and check whether
% the projected state falls within the capture envelope.
%
% This is a simple linear extrapolation:
%   t_contact = range / v_close
%   projected_lateral = current_lateral + v_lateral · t_contact
%   projected_attitude = current_att + att_rate · t_contact
%
% Each capture criterion is normalized to its tolerance (0 = perfect,
% 1 = at the limit). The worst criterion determines overall feasibility.
% Feasible = 1 if all criteria ≤ 1, infeasible = 0 otherwise.

capture_feasible = zeros(1, N_steps);
capture_margin   = zeros(1, N_steps);

for k = 1:N_steps
    rng = range_axial(k);
    v_close = vel_axial(k);
    
    if v_close <= 0 || rng <= 0
        % Not approaching or already past - projection is meaningless
        capture_feasible(k) = NaN;
        capture_margin(k)   = NaN;
        continue
    end
    
    % Time to contact at current closing rate (linear estimate)
    t_contact = rng / v_close;
    
    % Project lateral offset and attitude to contact time
    proj_lat = offset_lat(k) + vel_lat(k) * t_contact;
    proj_att = norm(att_state(:, k) + rate_state(:, k) * t_contact);
    proj_att_rate = norm(rate_state(:, k));
    
    % Normalize each criterion to its tolerance
    c1 = proj_lat / capture.lat_tol;          % lateral offset / tolerance
    c2 = v_close / capture.axial_vel_max;     % closing speed / max allowed
    c3 = vel_lat(k) / capture.lat_vel_max;    % lateral speed / max allowed
    c4 = proj_att / capture.att_tol;           % attitude error / tolerance
    c5 = proj_att_rate / capture.att_rate_max; % attitude rate / max allowed
    
    % Worst criterion determines feasibility
    capture_margin(k) = max([c1, c2, c3, c4, c5]);
    capture_feasible(k) = double(capture_margin(k) <= 1.0);
end

%% ====================== VELOCITY-RANGE CORRIDORS =======================
% Compute the maximum safe closing speed as a function of range.
% This is the velocity at which the stopping distance equals the range:
%   d_stop(v_max) = range  →  v·t_react + v²/(2a) = r
%   → v²/(2a) + v·t_react - r = 0
%   → v = -a·t_react + sqrt((a·t_react)² + 2·a·r)   (quadratic formula)
%
% This defines the "velocity corridor": the region in velocity-range
% space where RI < 1. Anything above this curve means the chaser cannot
% stop before hitting the target (or corridor wall).
%
% Two corridors are computed:
%   1) Cooperative (both vehicles brake): wider corridor, more margin
%   2) Chaser-only (abort case): narrower corridor

range_plot = linspace(0.5, corridor.range_far * 1.2, 500);

% Axial corridor using the quadratic formula above
a_axial = max(a_axial_pos, 1e-12);
a_axial_chaser_only = max(a_brake_chaser.Ipos, 1e-12);
vmax_axial_coop   = -t_react*a_axial + sqrt((t_react*a_axial).^2 + 2*a_axial*range_plot);
vmax_axial_chaser = -t_react*a_axial_chaser_only + ...
                    sqrt((t_react*a_axial_chaser_only).^2 + 2*a_axial_chaser_only*range_plot);
vmax_axial_coop   = max(vmax_axial_coop, 0);
vmax_axial_chaser = max(vmax_axial_chaser, 0);

% Lateral corridor: max lateral speed before hitting the corridor wall.
% Uses the worst-case lateral braking authority (min of R and C axes)
% and the corridor half-width as the available margin.
a_lat_axes = [min(a_brake.Rpos, a_brake.Rneg), ...
              min(a_brake.Cpos, a_brake.Cneg)];
a_lat_worst = min(a_lat_axes);

vmax_lat = zeros(1, length(range_plot));
hw_plot  = zeros(1, length(range_plot));
for j = 1:length(range_plot)
    hw_plot(j) = corridor_halfwidth(range_plot(j), corridor);
    d = hw_plot(j);  % worst case: chaser is centered, full half-width available
    vmax_lat(j) = -t_react*a_lat_worst + sqrt((t_react*a_lat_worst)^2 + 2*a_lat_worst*d);
end
vmax_lat = max(vmax_lat, 0);

%% ====================== FIGURES ========================================
% The following 8 figures provide a comprehensive view of approach safety.

% --- Figure 1: Axial Velocity-Range Corridor ---
% Shows the maximum safe closing speed vs. range for cooperative and
% chaser-only braking. The approach trajectory is overlaid as a
% scatter plot colored by time. Green/yellow/red zones indicate RI levels.
figure('Name','Axial Velocity-Range Corridor','Position',[100 100 900 600]);
hold on; grid on;

% Background fill for risk zones
fill([range_plot fliplr(range_plot)], ...
     [vmax_axial_coop*100*RI_green zeros(size(range_plot))]*0, ...
     [0.7 1.0 0.7], 'EdgeColor','none', 'FaceAlpha', 0.25, ...
     'HandleVisibility','off');

fill([range_plot fliplr(range_plot)], ...
     [vmax_axial_coop*100*RI_green zeros(size(range_plot))], ...
     [0.7 1.0 0.7], 'EdgeColor','none', 'FaceAlpha', 0.3);
fill([range_plot fliplr(range_plot)], ...
     [vmax_axial_coop*100*RI_yellow vmax_axial_coop*100*RI_green], ...
     [1.0 1.0 0.6], 'EdgeColor','none', 'FaceAlpha', 0.3);
fill([range_plot fliplr(range_plot)], ...
     [vmax_axial_coop*100 vmax_axial_coop*100*RI_yellow], ...
     [1.0 0.85 0.6], 'EdgeColor','none', 'FaceAlpha', 0.3);

plot(range_plot, vmax_axial_coop*100, 'b-', 'LineWidth', 2.5, ...
     'DisplayName', 'Cooperative limit (Asterix+Obelix)');
plot(range_plot, vmax_axial_chaser*100, 'r--', 'LineWidth', 2, ...
     'DisplayName', 'Chaser-only limit (abort case)');

valid = vel_axial > 0 & range_axial > 0;
scatter(range_axial(valid), vel_axial(valid)*100, 10, t_vec(valid), 'filled', ...
        'DisplayName', 'Approach trajectory');
cb = colorbar; ylabel(cb, 'Time [s]');

plot(speed_profile.range, speed_profile.v_cmd*100, 'k--', 'LineWidth', 2, ...
     'DisplayName', 'Speed profile (commanded)');

xlabel('Axial Range to Docking Port [m]');
ylabel('Axial Closing Speed [cm/s]');
title(sprintf('Axial Velocity-Range Corridor - %s -> %s', ast.name, obx.name));
legend('Location','northwest');
set(gca, 'XDir', 'reverse');  % range decreases left-to-right (approaching)
xlim([0 max(range_axial)*1.05]);

% --- Figure 2: Lateral Corridor ---
% Top: lateral offset vs. range with corridor walls
% Bottom: lateral speed vs. range with max-speed corridor
figure('Name','Lateral Corridor','Position',[100 100 900 500]);

subplot(2,1,1);
hold on; grid on;
plot(range_axial(valid), offset_lat(valid)*100, 'b-', 'LineWidth', 1.5, ...
     'DisplayName', 'Lateral offset');
hw_traj = arrayfun(@(r) corridor_halfwidth(r, corridor), range_axial(valid));
plot(range_axial(valid), hw_traj*100, 'r--', 'LineWidth', 2, ...
     'DisplayName', 'Corridor wall');
plot(range_axial(valid), -hw_traj*100, 'r--', 'LineWidth', 2, ...
     'HandleVisibility','off');
fill([range_axial(valid) fliplr(range_axial(valid))], ...
     [hw_traj*100 fliplr(-hw_traj*100)], ...
     [0.9 0.95 1], 'EdgeColor','none', 'FaceAlpha', 0.3, 'HandleVisibility','off');
set(gca, 'XDir', 'reverse');
xlabel('Axial Range [m]');
ylabel('Lateral Offset [cm]');
title('Lateral Position in Approach Corridor');
legend('Location','northwest');

subplot(2,1,2);
hold on; grid on;
plot(range_plot, vmax_lat*100, 'b-', 'LineWidth', 2, ...
     'DisplayName', 'Max lateral speed (corridor-centered)');
scatter(range_axial(valid), vel_lat(valid)*100, 8, t_vec(valid), 'filled', ...
        'DisplayName', 'Actual lateral speed');
cb = colorbar; ylabel(cb, 'Time [s]');
set(gca, 'XDir', 'reverse');
xlabel('Axial Range [m]');
ylabel('Lateral Speed [cm/s]');
title('Lateral Velocity Corridor');
legend('Location','northwest');

% --- Figure 3: Multi-Channel RI Time Histories ---
% Shows RI per channel and combined over time, with green/yellow/red zones.
figure('Name','Multi-Channel RI','Position',[100 100 1000 800]);

subplot(4,1,1);
plot(t_vec, range_axial, 'b-', 'LineWidth', 1.5); hold on; grid on;
ylabel('Axial Range [m]');
title(sprintf('Docking Approach: %s -> %s', ast.name, obx.name));

subplot(4,1,2);
RI_ax_plot = min(RI_axial, 5);     % clamp for plotting
RI_lt_plot = min(RI_lateral, 5);
RI_an_plot = min(RI_angular, 5);
plot(t_vec, RI_ax_plot, 'b-', 'LineWidth', 1.5); hold on; grid on;
plot(t_vec, RI_lt_plot, 'Color', [0.9 0.5 0.1], 'LineWidth', 1.5);
plot(t_vec, RI_an_plot, 'Color', [0.5 0.2 0.8], 'LineWidth', 1.5);
yline(RI_red, 'r--', 'LineWidth', 1.5);
yline(RI_yellow, 'Color', [0.9 0.7 0], 'LineStyle', '--');
ylabel('Channel RI');
legend('Axial','Lateral','Attitude','Abort','Caution','Location','best');
ylim([0 3]);

subplot(4,1,3);
RI_tot_plot = min(RI_total, 5);
plot(t_vec, RI_tot_plot, 'k-', 'LineWidth', 2); hold on; grid on;
yline(RI_red, 'r--', 'LineWidth', 1.5);
yline(RI_yellow, 'Color', [0.9 0.7 0], 'LineStyle', '--');
yline(RI_green, 'Color', [0.3 0.8 0.3], 'LineStyle', '--');
ylabel('Total RI');
legend('RI (max channel)','Abort','Caution','Comfortable','Location','best');
ylim([0 3]);

subplot(4,1,4);
DSM_plot = min(DSM_total, 10);
plot(t_vec, DSM_plot, 'b-', 'LineWidth', 1.5); hold on; grid on;
yline(1, 'r--', 'LineWidth', 1.5);
yline(2, 'Color', [0.9 0.7 0], 'LineStyle', '--');
ylabel('DSM');
xlabel('Time [s]');
legend('DSM','Limit (DSM=1)','Caution (DSM=2)','Location','best');
ylim([0 min(max(DSM_plot)*1.1, 10)]);

% --- Figure 4: Capture Envelope Feasibility ---
% Top: projected capture margin over time (1 = at limit)
% Bottom: attitude error magnitude vs. tolerance
figure('Name','Capture Envelope','Position',[100 100 900 500]);

subplot(2,1,1);
valid_cap = ~isnan(capture_margin);
scatter(t_vec(valid_cap), capture_margin(valid_cap), 12, ...
        capture_feasible(valid_cap), 'filled');
hold on; grid on;
yline(1.0, 'r--', 'LineWidth', 2);
colormap(gca, [1 0.3 0.3; 0.3 0.8 0.3]);  % red = infeasible, green = feasible
ylabel('Capture Margin (1 = limit)');
title('Projected Capture Envelope Feasibility');
legend('Projected state','Envelope limit','Location','best');
ylim([0 max(capture_margin(valid_cap))*1.1]);

subplot(2,1,2);
att_err_mag = vecnorm(att_state);  % ||[roll;pitch;yaw]|| at each timestep
plot(t_vec, att_err_mag, 'b-', 'LineWidth', 1.5); hold on; grid on;
yline(capture.att_tol, 'r--', 'LineWidth', 1.5);
ylabel('Attitude Error [deg]');
xlabel('Time [s]');
title('Relative Attitude Misalignment');
legend('|att error|','Capture tolerance','Location','best');

% --- Figure 5: Thrust Authority Comparison ---
% Bar chart comparing chaser-only vs. target-only vs. combined braking
% acceleration for each LVLH axis.
figure('Name','Thrust Authority Comparison','Position',[100 100 800 500]);

bar_data = [max_accel_ast(2)*1e3 max_accel_obx(2)*1e3 max_accel_combined(2)*1e3;   % -R (braking +R vel)
            max_accel_ast(4)*1e3 max_accel_obx(4)*1e3 max_accel_combined(4)*1e3;   % -I (braking +I vel)
            max_accel_ast(6)*1e3 max_accel_obx(6)*1e3 max_accel_combined(6)*1e3];  % -C (braking +C vel)
b = bar(bar_data, 'grouped'); grid on;
b(1).FaceColor = [0.3 0.6 0.9];
b(2).FaceColor = [0.9 0.5 0.2];
b(3).FaceColor = [0.3 0.8 0.3];
set(gca, 'XTickLabel', {'R (radial)','I (along-track)','C (cross-track)'});
ylabel('Braking Accel [mm/s^2]');
title('Braking Authority: Chaser vs Target vs Combined');
legend('Asterix (chaser)','Obelix (target)','Combined','Location','best');

% --- Figure 6: LVLH Trajectory with Corridor ---
% 2D projection of the approach trajectory (I vs R plane) with the
% funnel-shaped corridor overlaid. Shows the path the chaser takes
% from start to docking.
figure('Name','LVLH Trajectory + Corridor','Position',[100 100 800 600]);
hold on; grid on;

hw_near = corridor.half_width_close;
hw_far  = corridor.half_width_far;
rng_far = corridor.range_far;
corridor_I = [0, 0, -rng_far, -rng_far, 0];   % I coordinates of corridor corners
corridor_R = [hw_near, -hw_near, -hw_far, hw_far, hw_near];  % R coordinates
fill(corridor_I, corridor_R, [0.85 0.9 1], 'EdgeColor', [0.3 0.5 0.8], ...
     'FaceAlpha', 0.3, 'LineWidth', 1.5, 'DisplayName', 'Approach corridor');

plot(pos_I, pos_R, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Trajectory');
plot(pos_I(1), pos_R(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', ...
     'DisplayName', 'Start');
plot(0, 0, 'r^', 'MarkerSize', 12, 'MarkerFaceColor', 'r', ...
     'DisplayName', 'Target dock port');

xlabel('In-Track (I) [m]');
ylabel('Radial (R) [m]');
title('V-bar Approach Trajectory (LVLH I-R Plane)');
legend('Location', 'best');
axis equal;

% --- Figure 7: Guided Approach Telemetry ---
% Comprehensive time-history: range, speed tracking, attitude, delta-V
figure('Name','Guidance Telemetry','Position',[100 100 1000 900]);

subplot(4,1,1);
plot(t_vec, range_axial, 'b-', 'LineWidth', 1.5); hold on; grid on;
ylabel('Axial Range [m]');
title(sprintf('Guided Approach: %s -> %s', ast.name, obx.name));
if docked
    xline(dock_time, 'g-', 'LineWidth', 2);
    legend('Range','Docked','Location','best');
end

subplot(4,1,2);
plot(t_vec, vel_axial*100, 'b-', 'LineWidth', 1.5); hold on; grid on;
plot(t_vec, v_cmd_log*100, 'r--', 'LineWidth', 1.5);
ylabel('Axial Speed [cm/s]');
legend('Actual','Commanded','Location','best');

subplot(4,1,3);
plot(t_vec, att_state(1,:), '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.2); hold on; grid on;
plot(t_vec, att_state(2,:), '-', 'Color', [0.2 0.6 0.9], 'LineWidth', 1.2);
plot(t_vec, att_state(3,:), '-', 'Color', [0.3 0.8 0.3], 'LineWidth', 1.2);
yline(capture.att_tol, 'r--', 'LineWidth', 1);
yline(-capture.att_tol, 'r--', 'LineWidth', 1);
ylabel('Attitude Error [deg]');
legend('Roll','Pitch','Yaw','Tolerance','Location','best');

subplot(4,1,4);
plot(t_vec, dv_used*100, 'b-', 'LineWidth', 1.5); hold on; grid on;
yline((dv_remaining - dv_abort_reserve)*100, 'r--', 'LineWidth', 1.5);
ylabel('\DeltaV Used [cm/s]');
xlabel('Time [s]');
legend('\DeltaV expended','Budget (excl. abort reserve)','Location','best');

% --- Figure 8: Control Effort ---
% Translational and angular acceleration commands over time.
% Shows how hard the thrusters are working during the approach.
figure('Name','Control Effort','Position',[100 100 1000 600]);

subplot(2,1,1);
plot(t_vec, accel_cmd(1,:)*1e3, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 1); hold on; grid on;
plot(t_vec, accel_cmd(2,:)*1e3, '-', 'Color', [0.2 0.6 0.9], 'LineWidth', 1);
plot(t_vec, accel_cmd(3,:)*1e3, '-', 'Color', [0.3 0.8 0.3], 'LineWidth', 1);
ylabel('Translational Accel [mm/s^2]');
title('Control Commands');
legend('R (radial)','I (along-track)','C (cross-track)','Location','best');

subplot(2,1,2);
plot(t_vec, alpha_cmd(1,:), '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 1); hold on; grid on;
plot(t_vec, alpha_cmd(2,:), '-', 'Color', [0.2 0.6 0.9], 'LineWidth', 1);
plot(t_vec, alpha_cmd(3,:), '-', 'Color', [0.3 0.8 0.3], 'LineWidth', 1);
ylabel('Angular Accel Cmd [deg/s^2]');
xlabel('Time [s]');
legend('Roll','Pitch','Yaw','Location','best');

%% ====================== RISK-BY-RANGE HEATMAP ==========================
% 2D heatmap of RI as a function of range and closing speed.
% This is the velocity-range plane colored by risk level.
%
% For each (range, speed) point:
%   RI = d_stop(speed) / range = (v·t_react + v²/(2a)) / r
%
% Two panels:
%   Top: cooperative braking (both vehicles)
%   Bottom: chaser-only braking (abort scenario)
%
% The white curve is the RI=1 contour (velocity corridor boundary).
% Above the curve: RI > 1 (unsafe). Below: RI < 1 (safe).

range_sweep = linspace(0.5, corridor.range_far, 500);
speed_sweep = linspace(0.001, 0.10, 400);

[RR, VV] = meshgrid(range_sweep, speed_sweep);
D_stop_coop   = VV * t_react + VV.^2 / (2 * max(a_axial, 1e-12));
D_stop_chaser = VV * t_react + VV.^2 / (2 * max(a_axial_chaser_only, 1e-12));
RI_map_coop   = D_stop_coop ./ RR;
RI_map_chaser = D_stop_chaser ./ RR;

figure('Name','Docking Risk Heatmap','Position',[100 100 1000 800]);

subplot(2,1,1);
imagesc(range_sweep, speed_sweep*100, RI_map_coop);
set(gca, 'YDir', 'normal', 'XDir', 'reverse');
hold on;
plot(range_plot, vmax_axial_coop*100, 'w-', 'LineWidth', 2.5);  % RI=1 contour
cmap = [linspace(0.2,1,128)' linspace(0.8,0.2,128)' linspace(0.2,0.2,128)';
        linspace(1,0.6,128)' linspace(0.2,0,128)'   linspace(0.2,0,128)'];
colormap(gca, cmap);  % green (safe) → red (dangerous)
caxis([0 2]);
cb = colorbar; ylabel(cb, 'Risk Index');
xlabel('Axial Range [m]'); ylabel('Closing Speed [cm/s]');
title(sprintf('Cooperative Docking Risk Map - %s + %s braking', ast.name, obx.name));

subplot(2,1,2);
imagesc(range_sweep, speed_sweep*100, RI_map_chaser);
set(gca, 'YDir', 'normal', 'XDir', 'reverse');
hold on;
plot(range_plot, vmax_axial_chaser*100, 'w-', 'LineWidth', 2.5);
colormap(gca, cmap);
caxis([0 2]);
cb = colorbar; ylabel(cb, 'Risk Index');
xlabel('Axial Range [m]'); ylabel('Closing Speed [cm/s]');
title(sprintf('Abort-case Risk Map - %s only braking', ast.name));

%% ====================== CONSOLE SUMMARY ================================
% Print a comprehensive text summary of all key results.
fprintf('\n================================================================\n');
fprintf('  COOPERATIVE DOCKING DSM SUMMARY: %s -> %s\n', ast.name, obx.name);
fprintf('================================================================\n');
fprintf('Chaser: %s (%.1f kg, %d thr, %s mode)\n', ...
        ast.name, ast.mass, ast.n_thr, ast.ctrl_mode);
fprintf('Target: %s (%.1f kg, %d thr, %s mode)\n', ...
        obx.name, obx.mass, obx.n_thr, obx.ctrl_mode);
fprintf('Reaction time:     %.1f s\n', t_react);
fprintf('Orbit mean motion: %.4e rad/s\n', n);
fprintf('Approach axis:     V-bar (LVLH I)\n');
fprintf('Corridor:          %.2f m -> %.2f m half-width over %.0f m\n', ...
        corridor.half_width_far, corridor.half_width_close, corridor.range_far);
fprintf('----------------------------------------------------------------\n');
fprintf('COOPERATIVE braking [mm/s^2]:\n');
fprintf('  Axial (I):  %.4f\n', max(a_axial_pos, a_axial_neg)*1e3);
fprintf('  Lateral R:  %.4f\n', min(a_brake.Rpos, a_brake.Rneg)*1e3);
fprintf('  Lateral C:  %.4f\n', min(a_brake.Cpos, a_brake.Cneg)*1e3);
fprintf('CHASER-ONLY braking [mm/s^2]:\n');
fprintf('  Axial (I):  %.4f\n', max(a_brake_chaser.Ipos, a_brake_chaser.Ineg)*1e3);
fprintf('  Lateral R:  %.4f\n', min(a_brake_chaser.Rpos, a_brake_chaser.Rneg)*1e3);
fprintf('  Lateral C:  %.4f\n', min(a_brake_chaser.Cpos, a_brake_chaser.Cneg)*1e3);
fprintf('Angular braking [deg/s^2]:\n');
fprintf('  Roll:  %.4f   Pitch: %.4f   Yaw: %.4f\n', ...
        rad2deg(alpha_brake(1)), rad2deg(alpha_brake(2)), rad2deg(alpha_brake(3)));
fprintf('----------------------------------------------------------------\n');
fprintf('CAPTURE ENVELOPE:\n');
fprintf('  Lateral tolerance:    %.2f m\n', capture.lat_tol);
fprintf('  Axial speed limit:    %.3f m/s\n', capture.axial_vel_max);
fprintf('  Lateral speed limit:  %.3f m/s\n', capture.lat_vel_max);
fprintf('  Attitude tolerance:   %.1f deg\n', capture.att_tol);
fprintf('  Attitude rate limit:  %.1f deg/s\n', capture.att_rate_max);
fprintf('----------------------------------------------------------------\n');
fprintf('GUIDANCE PERFORMANCE:\n');
fprintf('  Speed profile: [');
for j = 1:length(speed_profile.range)
    fprintf('%.0fm:%.0fcm/s', speed_profile.range(j), speed_profile.v_cmd(j)*100);
    if j < length(speed_profile.range), fprintf(', '); end
end
fprintf(']\n');
fprintf('  Simulation time:      %.0f s\n', t_sim);
if docked
    fprintf('  DOCKING ACHIEVED at t = %.0f s\n', dock_time);
else
    fprintf('  Docking NOT reached in simulation window\n');
end
fprintf('  Total delta-V used:   %.4f m/s  (%.1f%% of budget)\n', ...
        dv_used(end), dv_used(end) / dv_remaining * 100);
fprintf('  Delta-V remaining:    %.4f m/s  (abort reserve: %.1f m/s)\n', ...
        dv_remaining - dv_used(end), dv_abort_reserve);
fprintf('----------------------------------------------------------------\n');

% Risk index statistics over the approach
valid_ri = RI_total < Inf & range_axial > 0;
if any(valid_ri)
    RI_max  = max(RI_total(valid_ri));
    RI_mean = mean(RI_total(valid_ri));
    [~, idx_ri_max] = max(RI_total .* valid_ri);
    
    fprintf('RISK INDEX SUMMARY (guided approach):\n');
    fprintf('  Peak RI:     %.3f at t = %.0f s, range = %.1f m\n', ...
            RI_max, t_vec(idx_ri_max), range_axial(idx_ri_max));
    
    % Identify which channel drove the peak RI
    ri_at_peak = [RI_axial(idx_ri_max), RI_lateral(idx_ri_max), RI_angular(idx_ri_max)];
    [~, ch_peak] = max(ri_at_peak);
    ch_names = {'AXIAL','LATERAL','ATTITUDE'};
    fprintf('    Peak driven by:  %s channel (%.3f)\n', ch_names{ch_peak}, ri_at_peak(ch_peak));
    fprintf('  Mean RI:     %.3f\n', RI_mean);
    fprintf('  %% time GREEN  (RI < %.1f):  %.1f%%\n', RI_green, ...
            100*sum(RI_total(valid_ri) < RI_green) / sum(valid_ri));
    fprintf('  %% time YELLOW (%.1f ≤ RI < %.1f):  %.1f%%\n', RI_green, RI_yellow, ...
            100*sum(RI_total(valid_ri) >= RI_green & RI_total(valid_ri) < RI_yellow) / sum(valid_ri));
    fprintf('  %% time RED    (RI ≥ %.1f):  %.1f%%\n', RI_yellow, ...
            100*sum(RI_total(valid_ri) >= RI_yellow) / sum(valid_ri));
end

% Find first caution/abort events
idx_yellow = find(RI_total >= RI_yellow & RI_total < Inf, 1, 'first');
idx_red    = find(RI_total >= RI_red & RI_total < Inf, 1, 'first');
if ~isempty(idx_yellow)
    fprintf('  RI hits CAUTION (%.1f) at t = %.0f s, range = %.1f m\n', ...
            RI_yellow, t_vec(idx_yellow), range_axial(idx_yellow));
    ri_channels = [RI_axial(idx_yellow), RI_lateral(idx_yellow), RI_angular(idx_yellow)];
    [~, ch] = max(ri_channels);
    fprintf('    Triggered by: %s channel (RI = %.3f)\n', ch_names{ch}, ri_channels(ch));
end
if ~isempty(idx_red)
    fprintf('  RI hits ABORT (%.1f) at t = %.0f s, range = %.1f m\n', ...
            RI_red, t_vec(idx_red), range_axial(idx_red));
    ri_channels = [RI_axial(idx_red), RI_lateral(idx_red), RI_angular(idx_red)];
    [~, ch] = max(ri_channels);
    fprintf('    Triggered by: %s channel (RI = %.3f)\n', ch_names{ch}, ri_channels(ch));
end
if isempty(idx_yellow) && isempty(idx_red)
    fprintf('  *** RI stays GREEN throughout entire approach. ***\n');
end

fprintf('----------------------------------------------------------------\n');
fprintf('ATTITUDE PERFORMANCE:\n');
fprintf('  Initial misalignment:  [%.2f, %.2f, %.2f] deg\n', att0(1), att0(2), att0(3));
fprintf('  Initial rates:         [%.3f, %.3f, %.3f] deg/s\n', att_rate0(1), att_rate0(2), att_rate0(3));
att_final = att_state(:, end);
rate_final = rate_state(:, end);
fprintf('  Final misalignment:    [%.4f, %.4f, %.4f] deg  (|err| = %.4f deg)\n', ...
        att_final(1), att_final(2), att_final(3), norm(att_final));
fprintf('  Final rates:           [%.5f, %.5f, %.5f] deg/s\n', ...
        rate_final(1), rate_final(2), rate_final(3));
fprintf('  Capture att tolerance: %.1f deg  ->  margin: %.2fx\n', ...
        capture.att_tol, capture.att_tol / max(norm(att_final), 1e-12));

% Capture feasibility summary
idx_first_feasible = find(capture_feasible == 1, 1, 'first');
if ~isempty(idx_first_feasible)
    fprintf('----------------------------------------------------------------\n');
    fprintf('CAPTURE ENVELOPE:\n');
    fprintf('  First feasible at t = %.0f s, range = %.1f m\n', ...
            t_vec(idx_first_feasible), range_axial(idx_first_feasible));
    pct_feasible = 100 * sum(capture_feasible == 1) / sum(~isnan(capture_feasible));
    fprintf('  Feasible for %.1f%% of approach\n', pct_feasible);
else
    fprintf('----------------------------------------------------------------\n');
    fprintf('CAPTURE ENVELOPE:\n');
    fprintf('  NEVER feasible during simulation.\n');
    fprintf('  (Projected contact state exceeds tolerances - check speed profile)\n');
end

fprintf('================================================================\n');

%% ====================== TABULATED CORRIDOR DATA ========================
% Print a table showing the maximum safe closing speed at selected ranges,
% comparing cooperative vs. chaser-only braking, and the percentage
% benefit from cooperative braking.
fprintf('\n  AXIAL VELOCITY CORRIDOR (cooperative vs chaser-only)\n');
fprintf('  Range [m]  |  v_max coop [cm/s]  v_max chaser [cm/s]  Benefit\n');
fprintf('  ------------------------------------------------------------ \n');

table_ranges = [2, 5, 10, 15, 20, 25, 30, 40, 50];
for rr = table_ranges
    % Quadratic formula: v_max = -a·t_react + sqrt((a·t_react)² + 2·a·r)
    vc = -t_react*a_axial + sqrt((t_react*a_axial)^2 + 2*a_axial*rr);
    vch = -t_react*a_axial_chaser_only + ...
          sqrt((t_react*a_axial_chaser_only)^2 + 2*a_axial_chaser_only*rr);
    vc = max(vc, 0);
    vch = max(vch, 0);
    if vch > 0
        benefit = (vc/vch - 1)*100;  % percentage increase from cooperation
    else
        benefit = Inf;
    end
    fprintf('  %7.1f    |     %7.3f            %7.3f           +%.0f%%\n', ...
            rr, vc*100, vch*100, benefit);
end
fprintf('  ------------------------------------------------------------ \n');

fprintf('\nDone.\n');