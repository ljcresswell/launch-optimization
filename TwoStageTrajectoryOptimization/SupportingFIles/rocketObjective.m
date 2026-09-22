function [f, c, ceq] = rocketObjective(u_scaled, params, N_nodes, t_max, scale)

%% Unscale
u = u_scaled .* scale;
solidPropMass  = u(1);
liquidPropMass = u(2);
launchAngle    = u(3);
delta_nodes    = u(4:end);

%% Failure return values
FAIL_F   = 10.0;
FAIL_CEQ = 10.0 * ones(2,1);
c        = [];

%% Update params
params.solidPropMassInitial  = solidPropMass;
params.liquidPropMassInitial = liquidPropMass;
params.massInitial = params.payloadMass  + params.engineMass1 ...
                   + params.engineMass2  ...
                   + (1 + params.massFrac1) * solidPropMass ...
                   + (1 + params.massFrac2) * liquidPropMass;

m_dot2 = m_dotFunc(params.specificHeatRatio2, params.At2, ...
                   params.totalPressure2, params.totalTemp2, ...
                   params.gasConstant2);
params.liquidBurnDuration = liquidPropMass * (1 - params.massFrac2) / m_dot2;

%% t_nodes must be [0,1] — Stage 3 ODE interpolates delta over burn fraction
t_nodes = linspace(0, 1, N_nodes);

%% Initial state
X0 = buildInitialState(launchAngle, params);

%% Integrate
try
    [X_all, ~, aux] = rocketIntegrator(X0, t_max, t_nodes, delta_nodes, params, true);
catch ME
    fprintf('  Integrator error: %s\n', ME.message);
    f = FAIL_F; ceq = FAIL_CEQ; return;
end

if isempty(X_all) || ~isfield(aux, 'X_stage3_end')
    f = FAIL_F; ceq = FAIL_CEQ; return;
end

%% Orbital elements at end of Stage 3
[a_f, e_f] = stateToOrbitalElements(aux.X_stage3_end, params);
if ~isfinite(a_f) || ~isfinite(e_f)
    f = FAIL_F; ceq = FAIL_CEQ; return;
end

%% Objective
prop_scale = scale(1) + scale(2);
f = (solidPropMass + liquidPropMass) / prop_scale;

%% Equality constraints
ceq(1) = (a_f - params.a) / params.a;
ceq(2) = (e_f - params.e) /0.1;
c = [];
end