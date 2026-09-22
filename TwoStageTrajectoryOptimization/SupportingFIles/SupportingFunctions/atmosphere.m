function [Patm, rho] = atmosphere(h)
    % ATMOSPHERE Compute atmospheric temperature, pressure, and density
    %   Inputs:
    %       h   - altitude [m]
    %   Outputs:
    %       Tatm - temperature [K]
    %       Patm - pressure [Pa]
    %       rho  - density [kg/m^3]
    
    % Constants
    Tsl = 288.15;    % sea level temperature [K]
    Psl = 101325;    % sea level pressure [Pa]
    L = 0.0065;      % temperature lapse rate [K/m]
    Rair = 287.0;    % specific gas constant for air [J/(kg-K)]
    g0 = 9.80665;    % gravitational acceleration [m/s^2]
    
    % Preallocate outputs
    Tatm = zeros(size(h));
    Patm = zeros(size(h));
    rho = zeros(size(h));
    
    % Compute atmospheric properties
    for i = 1:length(h)
        if h(i) <= 11000
            Tatm(i) = Tsl - L*h(i);
            Patm(i) = Psl * (Tatm(i)/Tsl)^(g0/(Rair*L));
        else
            Tatm(i) = 216.65;  % constant in stratosphere
            P11 = Psl * (216.65/Tsl)^(g0/(Rair*L));  % pressure at 11 km
            Patm(i) = P11 * exp(-g0*(h(i)-11000)/(Rair*Tatm(i)));
        end
        rho(i) = Patm(i)/(Rair*Tatm(i));
    end
end