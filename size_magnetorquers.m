function size_magnetorquers(tau_req, alt_km, inc_deg, varargin)
% SIZE_MAGNETORQUERS  Compute required magnetic dipole moments for 
%   magnetorquer sizing given torque authority requirements and orbit.
%
%   results = size_magnetorquers(tau_req, alt_km, inc_deg)
%   results = size_magnetorquers(tau_req, alt_km, inc_deg, Name, Value, ...)
%
%   INPUTS (required):
%     tau_req  - [3x1] Required torque authority in body frame [N·m]
%                [tau_x; tau_y; tau_z] (e.g., from disturbance torque budget)
%     alt_km   - Orbital altitude [km]
%     inc_deg  - Orbital inclination [deg]
%
%   NAME-VALUE PAIRS (optional):
%     'RAAN'          - Right ascension of ascending node [deg] (default: 0)
%     'NumPoints'     - Number of orbit sample points (default: 360)
%     'Margin'        - Design margin factor, e.g. 1.3 = 30% (default: 1.3)
%     'DutyCycle'     - Fraction of orbit actuators are active (default: 1.0)
%     'Percentile'    - Percentile of B_perp to use for sizing (default: 5)
%                       0 = absolute min (most conservative), 
%                       5 = 5th percentile (robust), 50 = median
%     'DipoleTilt'    - Geomagnetic dipole tilt angle [deg] (default: 11.5)
%     'DipoleRADeg'   - RA of geomagnetic dipole axis [deg] (default: -72.6)
%     'Plot'          - Generate diagnostic plots (default: true)
%
%   OUTPUT:
%     results - struct with fields:
%       .m_req          - [3x1] Required magnetic dipole moment [A·m^2]
%       .m_req_margin   - [3x1] With margin and duty cycle applied [A·m^2]
%       .B_min_perp     - [3x1] Min perpendicular B per axis [T]
%       .B_percentile   - [3x1] Percentile-based B_perp per axis [T]
%       .B_orbit_LVLH   - [3xN] B field in LVLH over orbit [T]
%       .B_magnitude    - [1xN] |B| over orbit [T]
%       .true_anom_deg  - [1xN] True anomaly sample points [deg]
%       .tau_available  - [3xN] Max available torque per axis over orbit [N·m]
%                         (given the sized dipole moments)
%       .orbit          - struct with orbit parameters used
%
%   COORDINATE FRAMES:
%     LVLH (Local-Vertical-Local-Horizontal):
%       x = along-track (velocity direction for circular orbit)
%       y = cross-track (orbit normal for right-handed frame)
%       z = radial (nadir-pointing, toward Earth center)
%     The magnetorquer rods are assumed aligned with LVLH body axes.
%
%   GEOMAGNETIC MODEL:
%     Uses a tilted dipole model (first-order spherical harmonic).
%     For higher fidelity, replace the dipole model section with IGRF
%     lookup (e.g., using igrfmagm from Aerospace Toolbox if available).
%
%   EXAMPLE:
%     % 850 km SSO, need 1e-4 N·m authority in each axis
%     tau = [1e-4; 1e-4; 1e-4];
%     results = size_magnetorquers(tau, 850, 98.8);
%
%     % More conservative: 10th percentile, 40% margin, 70% duty cycle
%     results = size_magnetorquers(tau, 850, 98.8, ...
%         'Percentile', 10, 'Margin', 1.4, 'DutyCycle', 0.7);
%
%   Yaman — OSCAR@VT ADCS Sizing Tool
%   Reference: Wertz, "Space Mission Engineering: The New SMAD", Ch. 19

    %% Parse inputs
    p = inputParser;
    addRequired(p, 'tau_req', @(x) isnumeric(x) && numel(x)==3);
    addRequired(p, 'alt_km',  @(x) isnumeric(x) && isscalar(x) && x>0);
    addRequired(p, 'inc_deg', @(x) isnumeric(x) && isscalar(x));

    addParameter(p, 'RAAN',        0,      @isnumeric);
    addParameter(p, 'NumPoints',   360,    @isnumeric);
    addParameter(p, 'Margin',      1.3,    @isnumeric);
    addParameter(p, 'DutyCycle',   1.0,    @isnumeric);
    addParameter(p, 'Percentile',  5,      @isnumeric);
    addParameter(p, 'DipoleTilt',  11.5,   @isnumeric);
    addParameter(p, 'DipoleRADeg', -72.6,  @isnumeric);
    addParameter(p, 'Plot',        true,   @islogical);

    parse(p, tau_req, alt_km, inc_deg, varargin{:});
    opts = p.Results;

    tau_req = tau_req(:);  % ensure column

    %% Constants
    RE    = 6371.2e3;                % Earth mean radius [km]
    mu_E  = 7.94e22;               % Earth magnetic dipole moment [A·m^2]
    mu0   = 4*pi*1e-7;             % permeability of free space [T·m/A]
    B0    = mu0*mu_E/(4*pi*RE^3);  % equatorial surface field ~3.12e-5 T
    % Note: B0 = mu0*mu_E/(4*pi) gives dipole strength, then /r^3

    r_orbit = RE + opts.alt_km*1e3;    % orbit radius [km]
    inc     = deg2rad(opts.inc_deg);
    RAAN    = deg2rad(opts.RAAN);
    N       = opts.NumPoints;

    % Geomagnetic dipole tilt (geographic -> geomagnetic)
    tilt_angle = deg2rad(opts.DipoleTilt);
    tilt_RA    = deg2rad(opts.DipoleRADeg);

    %% Sample true anomaly around orbit
    nu = linspace(0, 2*pi, N+1);
    nu = nu(1:end-1);  % N evenly spaced points, avoid duplicate at 2*pi

    %% Compute B field at each orbit point
    % Strategy:
    %   1. Compute satellite position in ECI (assuming circular orbit)
    %   2. Transform to geocentric coords
    %   3. Evaluate tilted dipole field in geocentric frame
    %   4. Rotate B into LVLH frame

    B_LVLH = zeros(3, N);

    for k = 1:N
        % --- Satellite position in ECI (perifocal -> ECI) ---
        % For circular orbit, argument of perigee is irrelevant; 
        % use nu as the argument of latitude (omega + nu, set omega=0)
        u = nu(k);  % argument of latitude

        % Position in ECI
        r_ECI = r_orbit * [cos(RAAN)*cos(u) - sin(RAAN)*sin(u)*cos(inc);
                           sin(RAAN)*cos(u) + cos(RAAN)*sin(u)*cos(inc);
                           sin(u)*sin(inc)];

        % --- Geocentric spherical coords ---
        r_mag = norm(r_ECI);
        lat_gc = asin(r_ECI(3)/r_mag);     % geocentric latitude
        lon_gc = atan2(r_ECI(2), r_ECI(1)); % geocentric longitude

        % --- Geomagnetic colatitude (angle from dipole axis) ---
        % Dipole axis unit vector in ECI
        m_hat = [cos(tilt_RA)*sin(tilt_angle);
                 sin(tilt_RA)*sin(tilt_angle);
                 cos(tilt_angle)];

        r_hat = r_ECI / r_mag;
        cos_theta_m = dot(r_hat, m_hat);  % cos of geomagnetic colatitude

        % --- Dipole field in ECI ---
        % B = (B0*(RE/r)^3) * [3*(m_hat·r_hat)*r_hat - m_hat]
        scale = B0 * (RE / r_mag)^3;
        B_ECI = scale * (3 * cos_theta_m * r_hat - m_hat);

        % --- Build LVLH rotation matrix ---
        % z_LVLH = -r_hat (nadir)
        % velocity direction for circular orbit:
        %   v_ECI = d/dt(r_ECI), direction = orbit tangent
        z_lvlh = -r_hat;  % nadir

        % Orbit normal (angular momentum direction)
        % h = r x v; for classical elements: h_hat depends on inc, RAAN
        h_hat = [sin(RAAN)*sin(inc);
                -cos(RAAN)*sin(inc);
                 cos(inc)];

        % Along-track (velocity direction) = h x (-z_lvlh) = h x r_hat
        x_lvlh = cross(h_hat, -z_lvlh);
        x_lvlh = x_lvlh / norm(x_lvlh);

        % Cross-track completes right-handed frame
        y_lvlh = cross(z_lvlh, x_lvlh);
        y_lvlh = y_lvlh / norm(y_lvlh);

        % DCM: ECI -> LVLH
        R_ECI2LVLH = [x_lvlh'; y_lvlh'; z_lvlh'];

        B_LVLH(:, k) = R_ECI2LVLH * B_ECI;
    end

    B_mag = vecnorm(B_LVLH, 2, 1);

    %% Compute perpendicular B component for each magnetorquer axis
    % A magnetorquer rod along axis i produces torque: tau = m_i * (e_i x B)
    % |tau| = m_i * |e_i x B| = m_i * |B_perp_to_i|
    % B_perp_to_i = sqrt(B_total^2 - B_i^2)

    B_perp = zeros(3, N);
    for ax = 1:3
        B_perp(ax, :) = sqrt(B_mag.^2 - B_LVLH(ax,:).^2);
    end

    %% Size using percentile of B_perp distribution
    B_perp_sizing = zeros(3, 1);
    B_perp_min    = zeros(3, 1);
    for ax = 1:3
        B_perp_sizing(ax) = prctile(B_perp(ax,:), opts.Percentile);
        B_perp_min(ax)    = min(B_perp(ax,:));
    end

    %% Compute required dipole moments
    % Basic: m = tau / B_perp
    m_req_basic = tau_req ./ B_perp_sizing;

    % With margin and duty cycle
    m_req_final = m_req_basic * opts.Margin / opts.DutyCycle;

    %% Compute available torque over orbit (for verification)
    tau_avail = zeros(3, N);
    for ax = 1:3
        tau_avail(ax, :) = m_req_final(ax) * B_perp(ax, :);
    end

    %% Pack output struct
    results.m_req           = m_req_basic;
    results.m_req_margin    = m_req_final;
    results.B_min_perp      = B_perp_min;
    results.B_percentile    = B_perp_sizing;
    results.B_orbit_LVLH    = B_LVLH;
    results.B_magnitude     = B_mag;
    results.B_perp_all      = B_perp;
    results.true_anom_deg   = rad2deg(nu);
    results.tau_available   = tau_avail;
    results.orbit.alt_km    = opts.alt_km;
    results.orbit.inc_deg   = opts.inc_deg;
    results.orbit.RAAN_deg  = opts.RAAN;
    results.orbit.r_km      = r_orbit;
    results.params.margin   = opts.Margin;
    results.params.duty     = opts.DutyCycle;
    results.params.pctile   = opts.Percentile;

    %% Print summary
    fprintf('\n========================================\n');
    fprintf('  MAGNETORQUER SIZING RESULTS\n');
    fprintf('========================================\n');
    fprintf('Orbit:  %.0f km altitude, %.1f° inclination (SSO)\n', ...
        opts.alt_km, opts.inc_deg);
    fprintf('        r = %.1f km, RAAN = %.1f°\n', r_orbit, opts.RAAN);
    fprintf('----------------------------------------\n');
    fprintf('Design: %.0f%% margin, %.0f%% duty cycle\n', ...
        (opts.Margin-1)*100, opts.DutyCycle*100);
    fprintf('        %dth percentile B_perp for sizing\n', opts.Percentile);
    fprintf('----------------------------------------\n');
    fprintf('  Axis   tau_req [N·m]   B_perp_%d%% [T]   m_basic [A·m²]   m_final [A·m²]\n', ...
        opts.Percentile);
    fprintf('  ----   ------------   -------------   --------------   --------------\n');
    labels = {'X (along)', 'Y (cross)', 'Z (nadir)'};
    for ax = 1:3
        fprintf('  %s   %10.2e     %10.2e      %10.2f       %10.2f\n', ...
            labels{ax}, tau_req(ax), B_perp_sizing(ax), ...
            m_req_basic(ax), m_req_final(ax));
    end
    fprintf('----------------------------------------\n');
    fprintf('|B| range over orbit: [%.2e, %.2e] T\n', min(B_mag), max(B_mag));
    fprintf('========================================\n\n');

    %% Diagnostic plots
    if opts.Plot
        nu_deg = rad2deg(nu);

        figure('Name', 'Magnetorquer Sizing', 'Position', [100 100 1200 800]);

        % --- Subplot 1: B field components in LVLH ---
        subplot(2,2,1);
        plot(nu_deg, B_LVLH(1,:)*1e6, 'r-', 'LineWidth', 1.2); hold on;
        plot(nu_deg, B_LVLH(2,:)*1e6, 'g-', 'LineWidth', 1.2);
        plot(nu_deg, B_LVLH(3,:)*1e6, 'b-', 'LineWidth', 1.2);
        plot(nu_deg, B_mag*1e6, 'k--', 'LineWidth', 1.0);
        xlabel('Argument of Latitude [deg]');
        ylabel('B [\muT]');
        title('Geomagnetic Field in LVLH Frame');
        legend('B_x (along)', 'B_y (cross)', 'B_z (nadir)', '|B|', ...
            'Location', 'best');
        grid on; xlim([0 360]);

        % --- Subplot 2: Perpendicular B per rod axis ---
        subplot(2,2,2);
        plot(nu_deg, B_perp(1,:)*1e6, 'r-', 'LineWidth', 1.2); hold on;
        plot(nu_deg, B_perp(2,:)*1e6, 'g-', 'LineWidth', 1.2);
        plot(nu_deg, B_perp(3,:)*1e6, 'b-', 'LineWidth', 1.2);
        for ax = 1:3
            colors = {'r','g','b'};
            yline(B_perp_sizing(ax)*1e6, [colors{ax}, '--'], ...
                sprintf('%dth %%ile', opts.Percentile), 'LineWidth', 1.0);
        end
        xlabel('Argument of Latitude [deg]');
        ylabel('B_\perp [\muT]');
        title('Perpendicular B Component per Rod Axis');
        legend('Rod X', 'Rod Y', 'Rod Z', 'Location', 'best');
        grid on; xlim([0 360]);

        % --- Subplot 3: Available torque vs required ---
        subplot(2,2,3);
        plot(nu_deg, tau_avail(1,:)*1e6, 'r-', 'LineWidth', 1.2); hold on;
        plot(nu_deg, tau_avail(2,:)*1e6, 'g-', 'LineWidth', 1.2);
        plot(nu_deg, tau_avail(3,:)*1e6, 'b-', 'LineWidth', 1.2);
        for ax = 1:3
            yline(tau_req(ax)*1e6, [colors{ax}, '--'], 'LineWidth', 1.0);
        end
        xlabel('Argument of Latitude [deg]');
        ylabel('\tau [\muN·m]');
        title('Available Torque vs Required (with margin)');
        legend('\tau_x avail', '\tau_y avail', '\tau_z avail', ...
            '\tau_x req', '\tau_y req', '\tau_z req', 'Location', 'best');
        grid on; xlim([0 360]);

        % --- Subplot 4: Summary bar chart ---
        subplot(2,2,4);
        b = bar(m_req_basic);
        b.FaceColor = [0.3 0.6 0.9];
        set(gca, 'XTickLabel', {'X (along)', 'Y (cross)', 'Z (nadir)'});
        ylabel('Dipole Moment [A·m^2]');
        title('Required Magnetic Dipole Moment');
        grid on;
 
        sgtitle(sprintf('Magnetorquer Sizing — %0.f km, %.1f° SSO', ...
            opts.alt_km, opts.inc_deg), 'FontWeight', 'bold');
    end
end