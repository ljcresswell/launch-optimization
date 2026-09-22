%% =========================================================
%  Rocket Trajectory Test Script
%  Validates ODE, runs full integration, plots everything
%% =========================================================
clc; clear; close all;
addpath(genpath(pwd));

%% =========================================================
%  STEP 1 — Load and derive params
%% =========================================================
fprintf('=== STEP 1: Loading Parameters ===\n\n');

params = readParamsFile("UserInputValues.txt");
disp('Parameters read from file.');

Mpi  = params.solidPropMassIG;
Mp2i = params.liquidPropMassIG;

params.solidPropMassInitial  = Mpi;
params.liquidPropMassInitial = Mp2i;
params.massInitial = (params.massFrac1 + 1)*Mpi  + ...
                     (params.massFrac2 + 1)*Mp2i  + ...
                     params.payloadMass            + ...
                     params.engineMass1            + ...
                     params.engineMass2;

fprintf('  massInitial : %.2f kg\n', params.massInitial);

validateParams(params);
params.stage = 1;

%% =========================================================
%  STEP 2 — Build X0 and sanity-check ODE at t = 0
%% =========================================================
fprintf('\n=== STEP 2: ODE Sanity Check at t = 0 ===\n\n');

r0       = params.radiusEarth;
gamma_0  = params.launchAngleIG;
V0       = params.gravityTurnVelocity;
v_r0     = V0 * sin(gamma_0);
v_theta0 = V0 * cos(gamma_0);
Mi       = params.massInitial;

X0 = [r0; 0; v_r0; v_theta0; Mi];

n_nodes     = 100;
t_check     = linspace(0, 600, n_nodes);
delta_check = zeros(1, n_nodes); % approxximate solution
params.deltaInterp = @(t) interp1(t_check, delta_check, t, 'linear', 'extrap');

dXdt0 = rocketODE_rTheta(X0, params, 0);

fprintf('  dr/dt       = %+.4f m/s      (expect ~%.4f)\n', dXdt0(1), v_r0);
fprintf('  dtheta/dt   = %+.6f rad/s\n',                   dXdt0(2));
fprintf('  dr_dot/dt   = %+.4f m/s^2   (expect negative)\n', dXdt0(3));
fprintf('  dv_theta/dt = %+.6f m/s^2\n',                   dXdt0(4));
fprintf('  dm/dt       = %+.4f kg/s    (expect negative)\n', dXdt0(5));

assert(dXdt0(1) > 0,   'dr/dt should be positive at launch');
assert(dXdt0(3) > 0,   'Rocket should be going upward at launch');
assert(dXdt0(5) < 0,   'Mass should decrease during burn');
fprintf('\n  All ODE sanity checks passed.\n');

%% =========================================================
%  STEP 3 — Full integration
%% =========================================================
fprintf('\n=== STEP 3: Full Integration ===\n\n');

t_max       = 10000;
n_nodes     = 200;
t_nodes     = linspace(0, t_max, n_nodes);
delta_nodes = [0.125*ones(1, floor(n_nodes/50)), 0.0*ones(1, ceil(49*n_nodes/50))];

[X, t, aux] = rocketIntegrator(X0, t_max, t_nodes, delta_nodes, params, true);

% Unpack state
r       = X(:,1);
theta   = X(:,2);
r_dot   = X(:,3);
v_theta = X(:,4);
m       = X(:,5);

% Derived quantities
V              = sqrt(r_dot.^2 + v_theta.^2);
altitude       = r - params.radiusEarth;
flightPathAngle = atan2(r_dot, v_theta);
specificEnergy = 0.5*V.^2 - params.gravConstant*params.massEarth ./ r;
angMomentum    = r .* v_theta;

fprintf('  Steps   : %d\n',        length(t));
fprintf('  t_end   : %.1f s\n',    t(end));
fprintf('  Alt end : %.2f km\n',   altitude(end)/1e3);
fprintf('  V end   : %.2f m/s\n',  V(end));
fprintf('  m end   : %.2f kg\n',   m(end));

%% =========================================================
%  STEP 4 — Find stage-change times
%           (detect discontinuities in mass / stage proxy)
%% =========================================================

% Mass drop at engine separation → find the largest single-step mass jump
dm      = diff(m);
massDrop = find(dm < -0.5*(params.engineMass1*0.5));  % rough threshold
stageTimes = t(massDrop);                              % engine-drop moment(s)

% Detect thrust-off → thrust-on transitions via m_dot proxy
% (mass is flat in coast, falling in burns)
mDotProxy = [0; dm ./ max(diff(t), 1e-10)];
burnOn  = find(diff(mDotProxy < -0.5) ==  1);   % coast→burn
burnOff = find(diff(mDotProxy < -0.5) == -1);   % burn→coast

% Collect all transition times for vertical lines
%stageChangeTimes = unique([stageTimes; t(burnOn); t(burnOff)]);
stageChangeTimes = t(find(diff(aux.stage) ~= 0));

%% =========================================================
%  STEP 5 — Plots
%% =========================================================
LW = 1.5;

%% Helper: draw stage-change lines on an axes
    function markStages(ax, stageChangeTimes)
        yl = ylim(ax);
        hold(ax, 'on');
        for k = 1:numel(stageChangeTimes)
            xline(ax, stageChangeTimes(k), '--', ...
                'Color', [0.6 0.6 0.6], 'LineWidth', 1, ...
                'Label', sprintf('S%d', k+1), ...
                'LabelVerticalAlignment', 'top', ...
                'LabelHorizontalAlignment', 'right', ...
                'FontSize', 7);
        end
        ylim(ax, yl);
    end

%% Pull aux fields
T_thrust = aux.T;        % thrust magnitude (N)
D_drag   = aux.D;        % drag magnitude (N)
gamma    = aux.gamma;    % flight path angle (rad)
delta    = aux.delta;    % nozzle deflection (rad)
V        = aux.V;        % total speed (m/s)

% Decompose thrust into radial / tangential
% Thrust acts along velocity vector + nozzle deflection
T_r     =  T_thrust .* sin(gamma + delta);
T_theta =  T_thrust .* cos(gamma + delta);

% Decompose drag into radial / tangential (opposes velocity)
D_r     = -D_drag .* sin(gamma);
D_theta = -D_drag .* cos(gamma);

% Gravity (purely radial, inward)
mu     = params.gravConstant * params.massEarth;
F_grav = mu * m ./ r.^2;

% Net force components
F_net_r     = T_r     + D_r     - F_grav;
F_net_theta = T_theta + D_theta;

%% Create main figure and tab group
fig = figure('Name', 'Rocket Trajectory Results', 'NumberTitle', 'off', ...
             'Position', [50 50 1300 820]);
tg  = uitabgroup(fig);

%% ── TAB 1: Core State Variables ──────────────────────────────────────────
tab1 = uitab(tg, 'Title', 'State Variables');
tl1  = tiledlayout(tab1, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl1, 'Core State Variables', 'FontSize', 13, 'FontWeight', 'bold');

ax = nexttile(tl1);
plot(t, altitude/1e3, 'b', 'LineWidth', LW);
ylabel('Altitude (km)'); xlabel('Time (s)'); title('Altitude'); grid on;
markStages(ax, stageChangeTimes);

ax = nexttile(tl1);
plot(t, V, 'r', 'LineWidth', LW);
ylabel('Speed (m/s)'); xlabel('Time (s)'); title('Total Speed'); grid on;
markStages(ax, stageChangeTimes);

ax = nexttile(tl1);
plot(t, m, 'Color', [0 0.7 0.7], 'LineWidth', LW);
ylabel('Mass (kg)'); xlabel('Time (s)'); title('Mass'); grid on;
markStages(ax, stageChangeTimes);

ax = nexttile(tl1);
plot(t, rad2deg(gamma), 'm', 'LineWidth', LW);
yline(0, 'k--', 'Horizontal', 'LabelHorizontalAlignment', 'left');
ylabel('FPA (deg)'); xlabel('Time (s)'); title('Flight Path Angle'); grid on;
markStages(ax, stageChangeTimes);

ax = nexttile(tl1);
plot(t, rad2deg(theta), 'g', 'LineWidth', LW);
ylabel('\theta (deg)'); xlabel('Time (s)'); title('Angular Position'); grid on;
markStages(ax, stageChangeTimes);

ax = nexttile(tl1);
plot(t, aux.mdot, 'Color', [0.7 0 0.7], 'LineWidth', LW);
ylabel('dm/dt (kg/s)'); xlabel('Time (s)'); title('Mass Flow Rate'); grid on;
markStages(ax, stageChangeTimes);

%% ── TAB 2: Control Plot ────────────────────────────────────────
tab6 = uitab(tg, 'Title', 'Control');
tl4  = tiledlayout(tab6, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl4, 'Control Variable', 'FontSize', 13, 'FontWeight', 'bold');
 
ax = nexttile(tl4);
plot(t, rad2deg(delta), 'b', 'LineWidth', LW);
yline(0, 'k--', 'HandleVisibility', 'off');
ylabel('\delta (deg)'); xlabel('Time (s)'); title('Nozzle Deflection Angle  \delta'); grid on;
markStages(ax, stageChangeTimes);

%% ── TAB 3: Velocity Decomposition ────────────────────────────────────────
tab2 = uitab(tg, 'Title', 'Velocity Decomposition');
tl2  = tiledlayout(tab2, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl2, 'Velocity Vector Decomposition', 'FontSize', 13, 'FontWeight', 'bold');

ax = nexttile(tl2);
plot(t, V, 'cyan', 'LineWidth', LW);
ylabel('|V| (m/s)'); xlabel('Time (s)'); title('Total Speed  |V|'); grid on;
markStages(ax, stageChangeTimes);

ax = nexttile(tl2);
plot(t, r_dot, 'Color', [1 0.4 0], 'LineWidth', LW);
yline(0, 'k--');
ylabel('v_r (m/s)'); xlabel('Time (s)'); title('Radial Velocity  v_r'); grid on;
markStages(ax, stageChangeTimes);

ax = nexttile(tl2);
plot(t, v_theta, 'Color', [0.2 0.6 0.2], 'LineWidth', LW);
ylabel('v_\theta (m/s)'); xlabel('Time (s)'); title('Tangential Velocity  v_\theta'); grid on;
markStages(ax, stageChangeTimes);

ax = nexttile(tl2);
plot(t, rad2deg(gamma), 'm', 'LineWidth', LW);
yline(0, 'k--', 'Horizontal');
ylabel('FPA (deg)'); xlabel('Time (s)'); title('Flight Path Angle  \gamma'); grid on;
markStages(ax, stageChangeTimes);

%% ── TAB 4: Force Decomposition ───────────────────────────────────────────
tab3 = uitab(tg, 'Title', 'Force Decomposition');
tl3  = tiledlayout(tab3, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl3, 'Force Decomposition — Radial & Tangential', 'FontSize', 13, 'FontWeight', 'bold');

% Thrust magnitude
ax = nexttile(tl3);
plot(t, T_thrust/1e3, 'b', 'LineWidth', LW);
ylabel('Force (kN)'); xlabel('Time (s)'); title('Thrust — Magnitude'); grid on;
markStages(ax, stageChangeTimes);

% Thrust r/theta
ax = nexttile(tl3);
h1 = plot(t, T_r/1e3,     'b',  'LineWidth', LW, 'DisplayName', 'T_r');     hold on;
h2 = plot(t, T_theta/1e3, 'b--','LineWidth', LW, 'DisplayName', 'T_\theta');
legend([h1 h2], 'Location', 'best');
ylabel('Force (kN)'); xlabel('Time (s)'); title('Thrust — Radial / Tangential'); grid on;
markStages(ax, stageChangeTimes);

% Drag magnitude
ax = nexttile(tl3);
plot(t, D_drag/1e3, 'r', 'LineWidth', LW);
ylabel('Force (kN)'); xlabel('Time (s)'); title('Drag — Magnitude'); grid on;
markStages(ax, stageChangeTimes);

% Drag r/theta
ax = nexttile(tl3);
h1 = plot(t, D_r/1e3,     'r',  'LineWidth', LW, 'DisplayName', 'D_r');     hold on;
h2 = plot(t, D_theta/1e3, 'r--','LineWidth', LW, 'DisplayName', 'D_\theta');
legend([h1 h2], 'Location', 'best');
ylabel('Force (kN)'); xlabel('Time (s)'); title('Drag — Radial / Tangential'); grid on;
markStages(ax, stageChangeTimes);

% Gravity
ax = nexttile(tl3);
plot(t, F_grav/1e3, 'Color', [0.5 0.3 0], 'LineWidth', LW);
ylabel('Force (kN)'); xlabel('Time (s)'); title('Gravity (inward radial)'); grid on;
markStages(ax, stageChangeTimes);

% Net r/theta
ax = nexttile(tl3);
h1 = plot(t, F_net_r/1e3,     'cyan',  'LineWidth', LW, 'DisplayName', 'F_{net,r}');     hold on;
h2 = plot(t, F_net_theta/1e3, 'k--','LineWidth', LW, 'DisplayName', 'F_{net,\theta}');
yline(0, 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
legend([h1 h2], 'Location', 'best');
ylabel('Force (kN)'); xlabel('Time (s)'); title('Net Force — Radial / Tangential'); grid on;
markStages(ax, stageChangeTimes);

%% ── TAB 5: Auxiliary Variables ───────────────────────────────────────────
% Plot any aux fields not already covered above
skipFields = {'T', 'D', 'mdot', 'gamma', 'delta', 'V', 'stage'};
auxFields  = fieldnames(aux);
plotFields = auxFields(~ismember(auxFields, skipFields));
nAux  = numel(plotFields);

if nAux > 0
    tab4 = uitab(tg, 'Title', 'Auxiliary');
    nCols = 2;
    nRows = ceil(nAux / nCols);
    tl4   = tiledlayout(tab4, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl4, 'Auxiliary Variables', 'FontSize', 13, 'FontWeight', 'bold');

    for k = 1:nAux
        ax  = nexttile(tl4);
        fld = plotFields{k};
        dat = aux.(fld);
        if numel(dat) == numel(t)
            plot(t, dat, 'LineWidth', LW);
            xlabel('Time (s)');
        else
            plot(dat, 'LineWidth', LW);
            xlabel('Index');
        end
        title(strrep(fld, '_', '\_')); grid on;
        markStages(ax, stageChangeTimes);
    end
end

%% ── TAB 6: Inertial Trajectory ───────────────────────────────────────────
tab5 = uitab(tg, 'Title', 'Trajectory');
ax5  = axes(tab5);

x_traj = r .* cos(theta) / 1e3;
y_traj = r .* sin(theta) / 1e3;

th_circle = linspace(0, 2*pi, 500);
x_earth   = params.radiusEarth/1e3 * cos(th_circle);
y_earth   = params.radiusEarth/1e3 * sin(th_circle);
x_atmo    = params.rEndAtmo/1e3    * cos(th_circle);
y_atmo    = params.rEndAtmo/1e3    * sin(th_circle);

hold(ax5, 'on');
fill(ax5, x_earth, y_earth, [0.3 0.6 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
plot(ax5, x_atmo, y_atmo, '--', 'Color', [0.5 0.8 1], 'LineWidth', 1);
scatter(ax5, x_traj, y_traj, 4, V, 'filled');
cb = colorbar(ax5); cb.Label.String = 'Speed (m/s)';
colormap(ax5, jet);

for k = 1:numel(stageChangeTimes)
    [~, idx] = min(abs(t - stageChangeTimes(k)));
    plot(ax5, x_traj(idx), y_traj(idx), 'kv', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
    text(ax5, x_traj(idx), y_traj(idx), sprintf('  S%d', k+1), 'FontSize', 8);
end

x_final = x_traj(end);
y_final = y_traj(end);
plot(ax5, [0 x_final], [0 y_final], 'r--', 'LineWidth', 1);

axis(ax5); grid(ax5, 'on');
xlabel(ax5, 'x (km)'); ylabel(ax5, 'y (km)');
title(ax5, 'Inertial Trajectory (coloured by speed)');
legend(ax5, {'Earth', 'Atmosphere boundary', 'Trajectory', 'End Asymptote'}, 'Location', 'best');
hold(ax5, 'off');

fprintf('\nAll plots generated.\n');