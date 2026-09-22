function [X_all, T_all, aux] = rocketIntegrator(X0, t_max, t_nodes, delta_nodes, params, returnAux)
%
% The initial state should be
%           r0          = radius of earth + 1000
%           theta0      = 0
%           vr          = 0.1;
%           v_theta   = calculated from gamma_0
%           Mi          = payload mass + engines + structures + propellants
%

if nargin < 6
    returnAux = false;
end


% Interpolate delta
params.deltaInterp = @(t) interp1(t_nodes, delta_nodes, t, 'linear', 'extrap');

T_all = [];
X_all = [];
stage_all = [];
aux = [];
t0 = 0;


%% Stage 1: Solid Burn ----------------------------------------------------
params.stage = 1;
options = odeset('RelTol',1e-8,'AbsTol',1e-8, ...
    'Events', @(t,X) rocketEvents(t,X,params));
[t1,X1,~,Xe,ie] = ode45(@(t,X) rocketODE_rTheta(X,params,t), [t0 t_max], X0, options);
T_all = [T_all; t1];
X_all = [X_all; X1];
stage_all = [stage_all; ones(length(t1),1)];

if isempty(ie)
    error('Stage 1: no event triggered before t_max.');

elseif ie(end) == 1
    % Normal: solid burnout inside atmosphere
    X_sep = Xe(end,:)';
    X_sep(5) = X_sep(5) - params.engineMass1 ...
                         - params.massFrac1*params.solidPropMassInitial;
    params.massInitial = params.massInitial - params.engineMass1 ...
                       - params.massFrac1*params.solidPropMassInitial;
    t0 = t1(end);

elseif ie(end) == 2
    % Exited atmosphere while still burning continue as Stage 1b (5)
    t0 = t1(end);
    params.stage = 5;
    options_1b = odeset('RelTol',1e-8,'AbsTol',1e-8, ...
        'Events', @(t,X) rocketEvents(t,X,params));
    [t1b,X1b,~,Xe1b,ie1b] = ode45(@(t,X) rocketODE_rTheta(X,params,t), ...
                                    [t0 t_max], Xe(end,:)', options_1b);
    T_all     = [T_all;     t1b(2:end)];
    X_all     = [X_all;     X1b(2:end,:)];
    stage_all = [stage_all; ones(length(t1b)-1,1)];
    t0 = t1b(end);

    if isempty(ie1b) || ie1b(end) ~= 1
        error('Stage 1b vacuum burn did not end on burnout. Got event: %s', num2str(ie1b));
    end

    X_sep = Xe1b(end,:)';
    X_sep(5) = X_sep(5) - params.engineMass1 ...
                         - params.massFrac1*params.solidPropMassInitial;
    params.massInitial = params.massInitial - params.engineMass1 ...
                       - params.massFrac1*params.solidPropMassInitial;

elseif ie(end) == 4
    error('Stage 1: rocket hit ground before burnout.');

else
    error('Stage 1: unexpected event %d.', ie(end));
end

%% Stage 2: Coast ---------------------------------------------------------
if X_sep(1) >= params.rEndAtmo % Skip coast phase if the solid burned through all atmo
    X_liq = X_sep;
else
    params.stage = 2;
    options = odeset('RelTol',1e-8,'AbsTol',1e-8, ...
    'Events', @(t,X) rocketEvents(t,X,params));   % rebuild with updated params
    
    [t2,X2,~,Xe,ie] = ode45(@(t,X) rocketODE_rTheta(X,params,t), [t0 t_max], X_sep, options);
    T_all = [T_all; t2(2:end)];
    X_all = [X_all; X2(2:end, :)];
    stage_all = [stage_all; 2*ones(length(t2)-1, 1)];
    t0 = t2(end);
end

if isempty(ie)
    warning('Stage 2 coast reached t_max without any event.');
    X_liq = X2(end,:)';
elseif ie(end) == 2        % exited atmosphere, vacuum burn
    X_liq = Xe(end,:)';
elseif ie(end) == 4        % rocket hit ground
    warning('Rocket crashed during coast. Terminating.');
    if returnAux, aux = computeAuxData(T_all, X_all, params, stage_all); end
    return;
elseif ie(end) == 5        % full orbit before vacuum burn
    warning('Full orbit completed during coast. Terminating.');
    if returnAux, aux = computeAuxData(T_all, X_all, params, stage_all); end
    return;
else
    warning('Unexpected event %d during Stage 2 coast.', ie(end));
    if returnAux, aux = computeAuxData(T_all, X_all, params, stage_all); end
    return;
end

%% Stage 3: Vacuum burn ---------------------------------------------------
params.t_stage3_start = t0;
params.stage = 3;
aux.X_stage3_start = X_liq;
options = odeset('RelTol',1e-8,'AbsTol',1e-8, ... % Must rebuild options for
    'Events', @(t,X) rocketEvents(t,X,params));   % new stage
[t3,X3,~,Xe,ie] = ode45(@(t,X) rocketODE_rTheta(X,params,t), [t0 t_max], X_liq, options);
T_all = [T_all; t3(2:end)];
X_all = [X_all; X3(2:end,:)];
stage_all = [stage_all; 3*ones(length(t3)-1, 1)];
t0 = t3(end);

% Stage 3 exit handling
if isempty(ie)
    X_fin = X3(end,:)';
elseif ie(end) == 3
    X_fin = Xe(end,:)';
    X_fin(5) = X_fin(5) - params.engineMass2 ...
                - params.massFrac2*params.liquidPropMassInitial;
    params.massInitial = params.massInitial - params.engineMass2 ...
                - params.massFrac2*params.liquidPropMassInitial;
elseif ie(end) == 4
    X_fin = Xe(end,:)';   
elseif ie(end) == 5
    X_fin = Xe(end,:)';
else
    X_fin = X3(end,:)';
end

% need stage 3 save before rocket fails in vacuum burn
aux.X_stage3_end = X_fin;

% all sorts of ways things could go wrong
if ie(end) == 4
    warning('Rocket re-entered during vacuum burn. Terminating.');
    if returnAux, aux = computeAuxData(T_all, X_all, params, stage_all); end
    aux.X_stage3_end = X_fin;   % computeAuxData may overwrite aux, restore it
    aux.X_stage3_start = X_liq; 
    return;
elseif ie(end) == 5
    warning('Full orbit during vacuum burn. Terminating.');
    if returnAux, aux = computeAuxData(T_all, X_all, params, stage_all); end
    aux.X_stage3_end = X_fin;
    aux.X_stage3_start = X_liq; 
    return;
end


%% Stage 4: Finish integrating till hit atmo again or full rev ------------
params.stage = 4;
options = odeset('RelTol',1e-8,'AbsTol',1e-8, ...
    'Events', @(t,X) rocketEvents(t,X,params));
[t4,X4,~,~,ie4] = ode45(@(t,X) rocketODE_rTheta(X,params,t), [t0 t_max], X_fin, options);
T_all = [T_all; t4(2:end)];
X_all = [X_all; X4(2:end,:)];
stage_all = [stage_all; 4*ones(length(t4)-1, 1)];

if ~isempty(ie4)
    %fprintf('Stage 4 ended on event %d at t = %.2f s\n', ie4(end), t4(end));
end

%% Return Aux Data --------------------------------------------------------
%% Return Aux Data
if returnAux
    aux = computeAuxData(T_all, X_all, params, stage_all);
    aux.X_stage3_end   = X_fin;
    aux.X_stage3_start = X_liq;
    
    % Stage transition times
    aux.t_stage1_end = T_all(find(stage_all == 1, 1, 'last'));
    aux.t_stage2_end = T_all(find(stage_all == 2, 1, 'last'));
    aux.t_stage3_end = T_all(find(stage_all == 3, 1, 'last'));
    aux.t_stage3_start = T_all(find(stage_all == 3, 1, 'first'));
    
    % Also store indices for inertial plot
    aux.idx_stage1_end   = find(stage_all == 1, 1, 'last');
    aux.idx_stage2_end   = find(stage_all == 2, 1, 'last');
    aux.idx_stage3_end   = find(stage_all == 3, 1, 'last');
    aux.idx_stage3_start = find(stage_all == 3, 1, 'first');
end

end