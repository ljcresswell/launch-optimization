function u0 = generateInitialGuess(params, N_nodes, t_max)

fprintf('\n------ Generating IG ------\n\n');

mu = params.gravConstant * params.massEarth;
g0 = 9.81;

%% Engine characterization
m_dot2 = m_dotFunc(params.specificHeatRatio2, params.At2, ...
                   params.totalPressure2, params.totalTemp2, params.gasConstant2);
T2     = thrust(params.Ae2, params.At2, params.specificHeatRatio2, ...
                params.totalTemp2, params.totalPressure2, params.gasConstant2, ...
                params.radiusEarth, params.a);
ve2 = T2 / m_dot2;

m_dot1 = m_dotFunc(params.specificHeatRatio1, params.At1, ...
                   params.totalPressure1, params.totalTemp1, params.gasConstant1);
T1     = thrust(params.Ae1, params.At1, params.specificHeatRatio1, ...
                params.totalTemp1, params.totalPressure1, params.gasConstant1, ...
                params.radiusEarth, params.radiusEarth + 1000);

fprintf('  Solid  engine: T=%.1f kN  mdot=%.1f kg/s\n', T1/1000, m_dot1);
fprintf('  Liquid engine: Isp=%.1f s  ve=%.1f m/s\n\n', ve2/g0, ve2);

t_nodes = linspace(0, 1, N_nodes);
delta   = zeros(N_nodes, 1);

%% Characteristic scales — everything normalized by these
% These define what "order 1" means for each variable and constraint.
% Choose physically meaningful references:
scale.solid  = 1e4;          % kg  — O(10,000 kg) solid prop
scale.liquid = 1e4;          % kg  — O(10,000 kg) liquid prop
scale.angle  = deg2rad(45);  % rad — O(45 deg) launch angle
scale.r      = params.a;     % m   — normalize altitude by target orbit radius
scale.v      = sqrt(mu/params.a);  % m/s — circular velocity at target
scale.dv     = sqrt(mu/params.a);  % m/s — same reference for dv

fprintf('  Scales:\n');
fprintf('    solid  : %.0f kg\n',  scale.solid);
fprintf('    liquid : %.0f kg\n',  scale.liquid);
fprintf('    angle  : %.1f deg\n', rad2deg(scale.angle));
fprintf('    r      : %.0f km\n',  scale.r/1000);
fprintf('    v_circ : %.0f m/s\n', scale.v);
fprintf('\n');

%% Initial guess (in scaled space)
dv_est    = 3500;
R_est     = exp(dv_est / ve2);
m_dry_est = params.payloadMass + params.engineMass2;
liq_guess = m_dry_est * (R_est - 1) / (1 - params.massFrac2*(R_est - 1));

m_fixed   = params.payloadMass + params.engineMass1 + params.engineMass2 ...
          + (1 + params.massFrac2) * liq_guess;
sol_guess = (T1 / (1.5 * g0) - m_fixed) / (1 + params.massFrac1);
sol_guess = max(sol_guess, 5000);
ang_guess = deg2rad(65);

% Scale the initial guess
u_guess_s = [sol_guess / scale.solid; ...
             liq_guess / scale.liquid; ...
             ang_guess / scale.angle];

fprintf('  Initial guess (unscaled):\n');
fprintf('    solid  = %.0f kg\n',  sol_guess);
fprintf('    liquid = %.0f kg\n',  liq_guess);
fprintf('    angle  = %.1f deg\n', rad2deg(ang_guess));
m0 = params.payloadMass + params.engineMass1 + params.engineMass2 ...
   + (1+params.massFrac1)*sol_guess + (1+params.massFrac2)*liq_guess;
fprintf('    TWR    = %.2f\n\n',   T1 / (m0 * g0));

fprintf('  Initial guess (scaled):\n');
fprintf('    solid_s  = %.3f\n', u_guess_s(1));
fprintf('    liquid_s = %.3f\n', u_guess_s(2));
fprintf('    angle_s  = %.3f\n', u_guess_s(3));
fprintf('\n');

%% Bounds (in scaled space)
m_no_solid = params.payloadMass + params.engineMass1 + params.engineMass2 ...
           + (1 + params.massFrac2) * liq_guess;
solid_max  = (T1 / (1.05 * g0) - m_no_solid) / (1 + params.massFrac1);

lb_s = [1e3      / scale.solid; ...
        1e3      / scale.liquid; ...
        deg2rad(30) / scale.angle];

ub_s = [solid_max / scale.solid; ...
        2e5       / scale.liquid; ...
        deg2rad(85) / scale.angle];

fprintf('  Bounds (scaled):\n');
fprintf('    solid_s : [%.3f, %.3f]\n', lb_s(1), ub_s(1));
fprintf('    liquid_s: [%.3f, %.3f]\n', lb_s(2), ub_s(2));
fprintf('    angle_s : [%.3f, %.3f]\n', lb_s(3), ub_s(3));
fprintf('\n');

%% fmincon
opts = optimoptions('fmincon', ...
    'Algorithm',                'sqp', ...
    'Display',                  'iter', ...
    'MaxFunctionEvaluations',   500, ...
    'MaxIterations',            100, ...
    'ConstraintTolerance',      1e-3, ...
    'OptimalityTolerance',      1e-4, ...
    'StepTolerance',            1e-6, ...
    'FiniteDifferenceStepSize', 1e-2, ...   % step in SCALED space
    'FiniteDifferenceType',     'central');

objFun = @(u_s) igObjective(u_s, scale);
conFun = @(u_s) igConstraints(u_s, scale, params, t_nodes, delta, ve2, mu, t_max);

warning('off', 'all');
[u_sol_s, ~, exitflag, output] = fmincon(objFun, u_guess_s, [], [], [], [], ...
                                          lb_s, ub_s, conFun, opts);
warning('on', 'all');

% Unscale solution
solidMass   = u_sol_s(1) * scale.solid;
liquidMass  = u_sol_s(2) * scale.liquid;
launchAngle = u_sol_s(3) * scale.angle;

fprintf('\n  Exit flag  : %d\n', exitflag);
fprintf('  Iterations : %d\n', output.iterations);
fprintf('  Func evals : %d\n', output.funcCount);

%% Re-simulate to get delta guess
params_ig = buildParams(params, solidMass, liquidMass);
X0 = buildInitialState(launchAngle, params_ig);

[X_all, T_all, ~] = rocketIntegrator(X0, t_max, t_nodes, delta, params_ig, false);

% Apoapsis — still needed for reporting
X_apo = coastApoapsis(T_all, X_all);

% Find state at atmosphere exit (start of vacuum burn)
r_all = X_all(:, 1);
idx_atmo_exit = find(r_all >= params.rEndAtmo, 1, 'first');

if isempty(idx_atmo_exit)
    warning('IG: never exited atmosphere, falling back to apoapsis gamma for delta guess.');
    X_ref = X_apo;
else
    X_ref = X_all(idx_atmo_exit, :);
end

% Taper delta from flight path angle at atmo exit down to zero
gamma_exit = atan2(X_ref(3), X_ref(4));
delta_ig = gamma_exit * linspace(1, 0, N_nodes)';

u0 = [solidMass; liquidMass; launchAngle; delta_ig];

% Report still uses X_apo
m0_final = buildParams(params, solidMass, liquidMass).massInitial;
fprintf('\n=== IG Result ===\n');
fprintf('  Solid  : %.1f kg  (%.2f x scale)\n', solidMass,   u_sol_s(1));
fprintf('  Liquid : %.1f kg  (%.2f x scale)\n', liquidMass,  u_sol_s(2));
fprintf('  Angle  : %.2f deg (%.2f x scale)\n', rad2deg(launchAngle), u_sol_s(3));
fprintf('  TWR    : %.2f\n',   T1 / (m0_final * g0));
fprintf('  r_apo  : %.1f km\n', (X_apo(1)-params.radiusEarth)/1000);
fprintf('  gamma_exit : %.2f deg\n', rad2deg(gamma_exit));
fprintf('  delta_ig   : [%.3f ... %.3f] rad\n', delta_ig(1), delta_ig(end));
fprintf('=================\n\n');

end


%% -------------------------------------------------------------------------
function J = igObjective(u_s, scale)
% Minimize total propellant mass, normalized to O(1).
% J = (solid + liquid) / (scale.solid + scale.liquid)
% At the initial guess this is roughly 1.0.

solid_kg  = u_s(1) * scale.solid;
liquid_kg = u_s(2) * scale.liquid;

J = (solid_kg + liquid_kg) / (scale.solid + scale.liquid);

end


%% -------------------------------------------------------------------------
function [c, ceq] = igConstraints(u_s, scale, params, t_nodes, delta, ve2, mu, t_max)
% All constraints normalized to O(1).
%
% Inequality c <= 0:
%   c(1): apoapsis altitude must reach target          — normalized by target r
%   c(2): liquid mass must cover circularization dv    — normalized by liq_needed
%
% Equality ceq == 0:
%   ceq(1): rdot at apoapsis ~ 0                       — normalized by v_circ

FAIL_VAL = 10.0;   % returned on crash — large O(1) violation, not 1e6

solidMass   = u_s(1) * scale.solid;
liquidMass  = u_s(2) * scale.liquid;
launchAngle = u_s(3) * scale.angle;

c   = FAIL_VAL * ones(3,1);
ceq = FAIL_VAL * ones(1,1);

params_ig = buildParams(params, solidMass, liquidMass);
X0 = buildInitialState(launchAngle, params_ig);

try
    [X_all, T_all, ~] = rocketIntegrator(X0, t_max, t_nodes, ...
                                          delta, params_ig, false);
catch
    return;
end

if isempty(X_all) || size(X_all,2) < 5 || max(X_all(:,1)) < params.rEndAtmo
    return;
end

X_apo  = coastApoapsis(T_all, X_all);
r_apo  = X_apo(1);
vt_apo = X_apo(4);
rd_apo = X_apo(3);

v_circ_apo = sqrt(mu / r_apo);
dv_needed  = abs(v_circ_apo - vt_apo);

R           = exp(dv_needed / ve2);
m_dry_fixed = params.payloadMass + params.engineMass2;
denom       = 1 - params.massFrac2 * (R - 1);

if denom <= 0 || ~isfinite(R)
    return;   % dv too large for engine — stay at FAIL_VAL
end

liq_needed = m_dry_fixed * (R - 1) / denom;

% c(1): apoapsis must reach target radius
%   = (r_target - r_apo) / r_target   →  0 when r_apo = r_target, >0 when short
c(1) = (params.a - r_apo) / params.a;

% c(2): liquid prop must be sufficient
%   = (liq_needed - liquidMass) / liq_needed  →  0 when exactly sufficient
if liq_needed > 0
    c(2) = (liq_needed - liquidMass) / liq_needed;
else
    c(2) = -1;   % no liquid needed, always satisfied
end

% ceq(1): rdot at apoapsis should be zero
%   normalized by circular velocity at that radius
ceq(1) = rd_apo / scale.v;

% c(3): gamma at atmo exit must be reasonable
r_all = X_all(:,1);
idx_exit = find(r_all >= params.rEndAtmo, 1, 'first');

    if isempty(idx_exit)
        c(3) = FAIL_VAL;
    else
        X_exit    = X_all(idx_exit,:);
        gamma_exit = atan2(X_exit(3), X_exit(4));
        % Require gamma < 35 deg at atmo exit, normalized to O(1)
        c(3) = (gamma_exit - deg2rad(35)) / deg2rad(35);
    end

end


%% -------------------------------------------------------------------------
function X_apo = coastApoapsis(T_all, X_all)
% Extract apoapsis of the coast arc (first rdot + -> - crossing).
% Clips trajectory before vacuum burn pulls rdot positive again.

rdot = X_all(:,3);

first_neg = find(rdot < 0, 1, 'first');
if ~isempty(first_neg)
    back_pos = find(rdot(first_neg:end) > 0, 1, 'first');
    if ~isempty(back_pos)
        cut   = first_neg + back_pos - 2;
        T_all = T_all(1:cut);
        X_all = X_all(1:cut,:);
        rdot  = X_all(:,3);
    end
end

idx = find(rdot(1:end-1) > 0 & rdot(2:end) <= 0, 1, 'first');

if isempty(idx)
    [~, idx] = max(X_all(:,1));
    X_apo = X_all(idx,:);
    return;
end

rdot1 = rdot(idx);
rdot2 = rdot(idx+1);
alpha = max(0, min(1, rdot1 / (rdot1 - rdot2)));
X_apo = X_all(idx,:) + alpha * (X_all(idx+1,:) - X_all(idx,:));

end

