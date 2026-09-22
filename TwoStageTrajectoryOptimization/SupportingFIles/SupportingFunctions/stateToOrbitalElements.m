function [a, e] = stateToOrbitalElements(X_end, params)
    r = X_end(1);
    r_dot = X_end(3);
    v_theta = X_end(4);
    mu = params.gravConstant * params.massEarth;

    h = r * v_theta;
    v_sq = r_dot^2 + v_theta^2;
    epsilon = v_sq/2 - mu/r;

    if epsilon >= 0
        % Hyperbolic or parabolic — not a closed orbit
        a = Inf;
        e = Inf;
        return;
    end

    a = -mu / (2*epsilon);
    e = sqrt(1 - h^2/(mu*a));
end