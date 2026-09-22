function drag = dragFunc(r, radiusEarth, V, Cd, Aref)
    % drag  Compute aerodynamic drag force.
    %
    %   drag = drag(r, radiusEarth, V, Cd, Aref)
    %
    %   Calculates drag using the standard drag equation based on atmospheric
    %   density at altitude, velocity, drag coefficient, and reference area.
    %
    %   Inputs:
    %       r            - Distance from Earth's center
    %       radiusEarth  - Earth radius
    %       V            - Velocity
    %       Cd           - Drag coefficient
    %       Aref         - Reference area
    %
    %   Output:
    %       drag         - Aerodynamic drag force
    %
    %   Requires: atmosphere

    [~, rho] = atmosphere(r-radiusEarth);
    drag = 0.5*rho*V^2*Cd*Aref;
end