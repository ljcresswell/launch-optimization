function f = rocketObjCached(u_scaled, params, N, t_max, scale)
    [f, ~, ~] = rocketSharedCache(u_scaled, params, N, t_max, scale);
end