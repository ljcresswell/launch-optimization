function [c, ceq] = rocketConCached(u_scaled, params, N, t_max, scale)
    [~, c, ceq] = rocketSharedCache(u_scaled, params, N, t_max, scale);
end