function plotTraj(T_all, X_all, aux, params, deltaOpt, t_nodes)

    r      = X_all(:,1);
    theta  = X_all(:,2);
    rdot   = X_all(:,3);
    vtheta = X_all(:,4);
    mass   = X_all(:,5);
    
    alt    = (r - params.radiusEarth) / 1000;
    v      = sqrt(rdot.^2 + vtheta.^2);
    gamma  = rad2deg(atan2(rdot, vtheta));
    t_min  = T_all / 60;
    
    % Stage transition times
    stage_times = [];
    stage_labels = {};
    if isfield(aux, 't_stage1_end') && ~isempty(aux.t_stage1_end)
        stage_times(end+1)  = aux.t_stage1_end / 60;
        stage_labels{end+1} = 'S1 end';
    end
    if isfield(aux, 't_stage3_start') && ~isempty(aux.t_stage3_start)
        stage_times(end+1)  = aux.t_stage3_start / 60;
        stage_labels{end+1} = 'S3 start';
    end
    if isfield(aux, 't_stage3_end') && ~isempty(aux.t_stage3_end)
        stage_times(end+1)  = aux.t_stage3_end / 60;
        stage_labels{end+1} = 'S3 end';
    end

    stage_colors = {[1.0 0.6 0.0], [0.4 0.9 0.4], [0.9 0.4 0.9]}; 
    
    mu     = params.gravConstant * params.massEarth;
    [a_f, e_f] = stateToOrbitalElements(aux.X_stage3_end, params);
    
    x_traj = r .* cos(theta) / 1000;
    y_traj = r .* sin(theta) / 1000;
    
    % Compute forces along trajectory
    n = length(T_all);
    F_r = zeros(n,1);
    F_t = zeros(n,1);
    D_r = zeros(n,1);
    D_t = zeros(n,1);
    G = zeros(n,1);
    
    for i = 1:n
        ri = r(i);
        rdi = rdot(i);
        vti = vtheta(i);
        mi = mass(i);
        Vi = sqrt(rdi^2 + vti^2);
        gam_i = atan2(rdi, vti);
        g_i = mu / ri^2;
        G(i) = g_i * mi;
    
        inAtmo = ri < params.rEndAtmo;
        if inAtmo
            Di = dragFunc(ri, params.radiusEarth, Vi, ...
                          params.dragCoefficient, params.referenceArea);
        else
            Di = 0;
        end
        D_r(i) =  sin(gam_i) * Di;   % drag radial
        D_t(i) = -cos(gam_i) * Di;   % drag tangential 
        end
    
    % Figures
    fig = figure('Name','Mission Analysis','Position',[100 50 1400 900]);
    tg  = uitabgroup(fig);
    

% TAB 1 — Overview subplots
    tab1 = uitab(tg, 'Title', 'Overview');
    
    ax = axes('Parent', tab1);  % need individual axes per subplot in a uitab
    
    positions = {[0.05 0.55 0.28 0.38], [0.37 0.55 0.28 0.38], [0.69 0.55 0.28 0.38], ...
                 [0.05 0.08 0.28 0.38], [0.37 0.08 0.28 0.38], [0.69 0.08 0.28 0.38]};
    delete(ax);
    
    % Altitude
    ax1 = axes('Parent', tab1, 'Position', positions{1});
    plot(ax1, t_min, alt, 'b', 'LineWidth', 1.5);
    yline(ax1, (params.a - params.radiusEarth)/1000,       'r--', 'Target');
    yline(ax1, (params.rEndAtmo - params.radiusEarth)/1000, 'k:', 'Atmo exit');
    xlabel(ax1,'Time (min)'); ylabel(ax1,'Altitude (km)');
    title(ax1,'Altitude vs Time'); grid(ax1,'on');
    addStageLines(ax1, stage_times, stage_labels, stage_colors);

        
    % Speed
    ax2 = axes('Parent', tab1, 'Position', positions{2});
    plot(ax2, t_min, v/1000, 'b', 'LineWidth', 1.5);
    xlabel(ax2,'Time (min)'); ylabel(ax2,'Speed (km/s)');
    title(ax2,'Speed vs Time'); grid(ax2,'on');
    addStageLines(ax2, stage_times, stage_labels, stage_colors);
    
    % Flight path angle
    ax3 = axes('Parent', tab1, 'Position', positions{3});
    plot(ax3, t_min, gamma, 'b', 'LineWidth', 1.5);
    yline(ax3, 0, 'r--', 'Horizontal');
    xlabel(ax3,'Time (min)'); ylabel(ax3,'\gamma (deg)');
    title(ax3,'Flight Path Angle'); grid(ax3,'on');
    addStageLines(ax3, stage_times, stage_labels, stage_colors);
    
    % Mass
    ax4 = axes('Parent', tab1, 'Position', positions{4});
    plot(ax4, t_min, mass/1000, 'b', 'LineWidth', 1.5);
    xlabel(ax4,'Time (min)'); ylabel(ax4,'Mass (tonnes)');
    title(ax4,'Vehicle Mass'); grid(ax4,'on');
    addStageLines(ax4, stage_times, stage_labels, stage_colors);
    
    % Delta profile
    ax5 = axes('Parent', tab1, 'Position', positions{5});
    plot(ax5, t_nodes, rad2deg(deltaOpt), 'bo-', 'LineWidth',1.5, 'MarkerFaceColor','b');
    xlabel(ax5,'Burn fraction'); ylabel(ax5,'\delta (deg)');
    title(ax5,'Thrust Angle (Stage 3)'); grid(ax5,'on');
    
    % Radial velocity
    ax6 = axes('Parent', tab1, 'Position', positions{6});
    plot(ax6, t_min, rdot, 'b', 'LineWidth', 1.5); hold(ax6,'on');
    plot(ax6, t_min, vtheta, 'r', 'LineWidth', 1.5);
    yline(ax6, 0, 'k--');
    legend(ax6, 'r\_dot','v\_\theta','Location','best');
    xlabel(ax6,'Time (min)'); ylabel(ax6,'m/s');
    title(ax6,'Velocity Components'); grid(ax6,'on');
    addStageLines(ax6, stage_times, stage_labels, stage_colors);
    
    annotation(tab1,'textbox',[0.2 0.96 0.6 0.03], ...
        'String', sprintf('a = %.1f km    e = %.5f    Total prop = (see report)', ...
        (a_f-params.radiusEarth)/1000, e_f), ...
        'HorizontalAlignment','center','EdgeColor','none','FontSize',12);

% TAB 2 — Inertial trajectory

    tab2 = uitab(tg, 'Title', 'Inertial Trajectory');
    ax_i = axes('Parent', tab2, 'Position', [0.1 0.08 0.82 0.85]);
    hold(ax_i, 'on');
    
    th = linspace(0, 2*pi, 500);
    
    % Earth
    fill(ax_i, params.radiusEarth/1000 * cos(th), ...
              params.radiusEarth/1000 * sin(th), ...
              [0.3 0.6 1.0], 'EdgeColor','k','LineWidth',1.5);
    
    % Target orbit
    plot(ax_i, params.a/1000 * cos(th), params.a/1000 * sin(th), ...
         'r--', 'LineWidth', 1.5);
    
    % Powered ascent
    plot(ax_i, x_traj, y_traj, 'b', 'LineWidth', 2);
    
    % Final orbit
    X_end  = aux.X_stage3_end;
    r_f    = X_end(1);  th_f = X_end(2);  vt_f = X_end(4);
    T_orb  = 2*pi*sqrt(r_f^3 / mu);
    omega  = vt_f / r_f;
    th_ext = th_f + omega * linspace(0, T_orb, 1000);
    plot(ax_i, r_f/1000 * cos(th_ext), r_f/1000 * sin(th_ext), ...
         'g--', 'LineWidth', 1.5);
    
    % Launch marker
    plot(ax_i, x_traj(1), y_traj(1), 'ko', 'MarkerFaceColor','k','MarkerSize',8);
    
    % Stage markers — must be inside hold, before legend
    stage_idx   = {aux.idx_stage1_end, aux.idx_stage3_start, aux.idx_stage3_end};
    marker_syms = {'o', 's', '^'};
    for k = 1:length(stage_idx)
        idx = stage_idx{k};
        if ~isempty(idx) && idx <= length(x_traj)
            c = stage_colors{k};
            plot(ax_i, x_traj(idx), y_traj(idx), marker_syms{k}, ...
                 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 9, ...
                 'DisplayName', stage_labels{k});
        end
    end
    
    axis(ax_i, 'equal');
    ax_i.DataAspectRatioMode    = 'manual';  ax_i.DataAspectRatio    = [1 1 1];
    ax_i.PlotBoxAspectRatioMode = 'manual';  ax_i.PlotBoxAspectRatio = [1 1 1];
    addlistener(ax_i, 'XLim', 'PostSet', @(~,~) enforceLimits(ax_i));
    addlistener(ax_i, 'YLim', 'PostSet', @(~,~) enforceLimits(ax_i));
    
    grid(ax_i, 'on');
    xlabel(ax_i, 'x (km)'); ylabel(ax_i, 'y (km)');
    title(ax_i, sprintf('Inertial Trajectory — a=%.1f km  e=%.5f', ...
        (a_f-params.radiusEarth)/1000, e_f), 'FontSize', 13);
    
    % Legend after everything is plotted
    legend(ax_i, 'Earth', 'Target orbit', 'Powered ascent', 'Final orbit', 'Launch', ...
           stage_labels{:}, 'Location', 'best');
    
    hold(ax_i, 'off');

% TAB 3 — Forces and velocity decomposition

    tab3 = uitab(tg, 'Title', 'Forces & Velocities');
    
    pos3 = {[0.05 0.55 0.42 0.38], [0.55 0.55 0.42 0.38], ...
            [0.05 0.08 0.42 0.38], [0.55 0.08 0.42 0.38]};
    
    % Velocity decomposition
    ax31 = axes('Parent', tab3, 'Position', pos3{1});
    plot(ax31, t_min, rdot,   'b',  'LineWidth', 1.5); hold(ax31,'on');
    plot(ax31, t_min, vtheta, 'r',  'LineWidth', 1.5);
    plot(ax31, t_min, v,      'k--','LineWidth', 1.2);
    yline(ax31, 0, 'k:');
    legend(ax31, 'v_r (radial)','v_\theta (tangential)','|v| total','Location','best');
    xlabel(ax31,'Time (min)'); ylabel(ax31,'Velocity (m/s)');
    title(ax31,'Velocity Components'); grid(ax31,'on');
    addStageLines(ax31, stage_times, stage_labels, stage_colors);
    
    % Drag decomposition
    ax32 = axes('Parent', tab3, 'Position', pos3{2});
    plot(ax32, t_min, D_r/1000, 'b', 'LineWidth', 1.5); hold(ax32,'on');
    plot(ax32, t_min, D_t/1000, 'r', 'LineWidth', 1.5);
    legend(ax32, 'D_r (radial)','D_t (tangential)','Location','best');
    xlabel(ax32,'Time (min)'); ylabel(ax32,'Drag force (kN)');
    title(ax32,'Drag Components'); grid(ax32,'on');
    addStageLines(ax32, stage_times, stage_labels, stage_colors);
    
    % Thrust decomposition 
    ax33 = axes('Parent', tab3, 'Position', pos3{3});
    g_vec    = mu ./ r.^2;
    accel_r3 = gradient(rdot,   T_all);
    accel_t3 = gradient(vtheta, T_all);
    T_r = mass .* (accel_r3 - vtheta.^2./r + g_vec) + D_r;
    T_t = mass .* (accel_t3 - rdot.*vtheta./r)       + D_t;
    plot(ax33, t_min, T_r/1000, 'b', 'LineWidth', 1.5); hold(ax33,'on');
    plot(ax33, t_min, T_t/1000, 'r', 'LineWidth', 1.5);
    yline(ax33, 0, 'k:');
    legend(ax33, 'T_r (radial)','T_\theta (tangential)','Location','best');
    xlabel(ax33,'Time (min)'); ylabel(ax33,'Thrust (kN)');
    title(ax33,'Thrust Components'); grid(ax33,'on');
    addStageLines(ax33, stage_times, stage_labels, stage_colors);
    
    % Net radial and tangential acceleration
    ax34 = axes('Parent', tab3, 'Position', pos3{4});
    accel_r = gradient(rdot,   T_all); % Technically finite diff
    accel_t = gradient(vtheta, T_all);
    plot(ax34, t_min, accel_r, 'b', 'LineWidth', 1.5); hold(ax34,'on');
    plot(ax34, t_min, accel_t, 'r', 'LineWidth', 1.5);
    yline(ax34, 0, 'k:');
    legend(ax34, 'a_r (radial)','a_\theta (tangential)','Location','best');
    xlabel(ax34,'Time (min)'); ylabel(ax34,'Acceleration (m/s²)');
    title(ax34,'Net Acceleration Components'); grid(ax34,'on');
    addStageLines(ax34, stage_times, stage_labels, stage_colors);
    
    annotation(tab3,'textbox',[0.2 0.96 0.6 0.03], ...
        'String','Forces & velocity decomposition along trajectory', ...
        'HorizontalAlignment','center','EdgeColor','none','FontSize',12);
    
% TAB 4 — Ascent detail (cropped to 95% apo height) 
    tab4 = uitab(tg, 'Title', 'Ascent Detail');
    
    % Dark background
    tab4.BackgroundColor = [0.12 0.12 0.12];
    ax_z = axes('Parent', tab4, 'Position', [0.08 0.08 0.84 0.85]);
    hold(ax_z, 'on');
    
    % Style
    ax_z.Color = [0.12 0.12 0.12];
    ax_z.XColor = [0.85 0.85 0.85];
    ax_z.YColor = [0.85 0.85 0.85];
    ax_z.GridColor = [0.35 0.35 0.35];
    ax_z.GridAlpha = 0.5;
    ax_z.MinorGridColor = [0.25 0.25 0.25];
    
    % Crop to 90% of apoapsis altitude
    r_max = max(r);
    r_target = params.radiusEarth + 0.95 * (r_max - params.radiusEarth);
    idx_crop = find(r >= r_target, 1, 'first');
    if isempty(idx_crop)
        idx_crop = length(r);
    end
    
    r_crop = r(1:idx_crop);
    th_crop = theta(1:idx_crop);
    x_crop = r_crop .* cos(th_crop) / 1000;
    y_crop = r_crop .* sin(th_crop) / 1000;
    
    % Rotate so launch is at bottom
    rot_angle = -(theta(1) + pi/2);
    R_mat = [cos(rot_angle) -sin(rot_angle);
                 sin(rot_angle)  cos(rot_angle)];
    
    xy_rot = R_mat * [x_crop'; y_crop'];
    
    % Earth arc
    th_earth     = linspace(theta(1) - deg2rad(10), th_crop(end) + deg2rad(10), 400);
    xy_earth_rot = R_mat * [params.radiusEarth/1000 * cos(th_earth);
                             params.radiusEarth/1000 * sin(th_earth)];
    
    % Target altitude arc
    th_arc     = linspace(theta(1) - deg2rad(5), th_crop(end) + deg2rad(5), 400);
    xy_arc_rot = R_mat * [params.a/1000 * cos(th_arc);
                           params.a/1000 * sin(th_arc)];
    
    % 90% crop endpoint marker
    crop_rot = R_mat * [r_crop(end)*cos(th_crop(end))/1000;
                         r_crop(end)*sin(th_crop(end))/1000];
    
    % Plot
    plot(ax_z, xy_earth_rot(1,:), xy_earth_rot(2,:), ...
         'Color', [0.30 0.65 1.00], 'LineWidth', 3);
    plot(ax_z, xy_arc_rot(1,:),   xy_arc_rot(2,:), ...
         '--', 'Color', [1.00 0.40 0.40], 'LineWidth', 1.5);
    plot(ax_z, xy_rot(1,:),       xy_rot(2,:), ...
         'Color', [0.40 1.00 0.60], 'LineWidth', 2);
    plot(ax_z, xy_rot(1,1),       xy_rot(2,1), ...
         'o', 'Color', [1.00 0.80 0.20], 'MarkerFaceColor', [1.00 0.80 0.20], 'MarkerSize', 8);
    plot(ax_z, crop_rot(1),       crop_rot(2), ...
         's', 'Color', [1.00 0.40 0.40], 'MarkerFaceColor', [1.00 0.40 0.40], 'MarkerSize', 8);
    
    axis(ax_z, 'equal');
    ax_z.DataAspectRatioMode    = 'manual';  ax_z.DataAspectRatio    = [1 1 1];
    ax_z.PlotBoxAspectRatioMode = 'manual';  ax_z.PlotBoxAspectRatio = [1 1 1];
    addlistener(ax_z, 'XLim', 'PostSet', @(~,~) enforceLimits(ax_z));
    addlistener(ax_z, 'YLim', 'PostSet', @(~,~) enforceLimits(ax_z));
    
    grid(ax_z, 'on');
    xlabel(ax_z, 'x_{rot} (km)', 'Color', [0.85 0.85 0.85]);
    ylabel(ax_z, 'y_{rot} (km)', 'Color', [0.85 0.85 0.85]);
    title(ax_z, sprintf('Ascent Detail — shown to 90%% of apoapsis (%.0f km)', ...
        (r_max - params.radiusEarth)/1000), 'FontSize', 12, 'Color', [0.95 0.95 0.95]);
    
    t_crop_min = T_all(idx_crop) / 60;
    
    lg = legend(ax_z, 'Earth surface', 'Target altitude', 'Trajectory', ...
                'Launch', sprintf('95%% apoapsis (t = %.1f min)', t_crop_min), ...
                'Location', 'best');
    lg.Color     = [0.18 0.18 0.18];
    lg.TextColor = [0.85 0.85 0.85];
    lg.EdgeColor = [0.40 0.40 0.40];
    
    hold(ax_z, 'off');

end

% ----------------------------------------------------------------
function enforceLimits(ax)
    persistent busy
    if ~isempty(busy) && busy, return; end
    busy = true;
    xl   = ax.XLim;  yl = ax.YLim;
    xc   = mean(xl); yc = mean(yl);
    half = max(diff(xl), diff(yl)) / 2;
    ax.XLim = [xc-half, xc+half];
    ax.YLim = [yc-half, yc+half];
    busy = false;
end

% -------------------------------------------------------------------------
function addStageLines(ax, stage_times, stage_labels, stage_colors)
% Adds vertical lines at stage transitions to a time-series axes
    yl = ylim(ax);
    for k = 1:length(stage_times)
        c = stage_colors{mod(k-1, length(stage_colors)) + 1};
        xline(ax, stage_times(k), '--', stage_labels{k}, ...
              'Color', c, 'LineWidth', 1.2, ...
              'LabelVerticalAlignment', 'bottom', ...
              'FontSize', 8);
    end
end