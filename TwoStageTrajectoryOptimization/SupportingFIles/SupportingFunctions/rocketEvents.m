function [value, isterminal, direction] = rocketEvents(t, X, params)

    % Terminal: 
    %           1 - stop integration
    %           0 - continue integrating
    % direction:
    %           1 - value goes from - to +
    %           0 - any crossing
    %           -1 - value goes from + to -

    r = X(1);
    m = X(5);
    Mi   = params.massInitial;
    Mp1  = params.solidPropMassInitial;
    Mp2  = params.liquidPropMassInitial;
    stage = params.stage;

    % Initialize events
    value      = ones(5,1);
    isterminal = zeros(5,1);
    direction  = zeros(5,1);

    % Event 1: end solid burn
    if stage == 1 || stage == 5
        value(1)      = m - (Mi - Mp1);
        isterminal(1) = 1;
        direction(1)  = -1;
    end

    % Event 2: exit atmosphere
    if stage <= 2 && stage ~= 5 % Stage 5 being solids burning exoatmo
        value(2) = r - params.rEndAtmo;
        isterminal(2) = 1;
        direction(2) = 1;
    end

    % Event 3: end liquid burn
    if stage == 3
        value(3)      = m - (Mi - Mp1 - Mp2);
        isterminal(3) = 1;
        direction(3)  = -1;
    end

    % Event 4: hit Earth / re-enter atmo
    value(4) = r - params.radiusEarth; % 40 km alt. (to keep drag from exploding)
    isterminal(4) = 1;
    direction(4) = -1;

    % Event 5: one full revolution
    value(5) = X(2) - 2*pi;
    isterminal(5) = 1;
    direction(5) = 1;
end