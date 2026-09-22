function [f, c, ceq] = rocketSharedCache(u_scaled, params, N, t_max, scale)
    persistent last_x last_f last_c last_ceq

    if isempty(last_x) || ~isequal(u_scaled, last_x)
        [last_f, last_c, last_ceq] = rocketObjective(u_scaled, params, N, t_max, scale);
        last_x = u_scaled;
    end

    f   = last_f;
    c   = last_c;
    ceq = last_ceq;
end