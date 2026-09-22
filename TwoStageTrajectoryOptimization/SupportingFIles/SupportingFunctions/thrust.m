function thrust = thrust(Ae, At, gamma, T0, P0, R, radiusEarth, r)
    % thrust  Compute rocket nozzle thrust using isentropic flow relations.
    %
    %   thrust = thrust(Ae, At, gamma, T0, P0, R, radiusEarth, r)
    %
    %   Calculates thrust as the sum of momentum and pressure forces based on
    %   chamber conditions, nozzle geometry, and ambient pressure at altitude.
    %
    %   Inputs:
    %       Ae           - Nozzle exit area
    %       At           - Nozzle throat area
    %       gamma        - Specific heat ratio
    %       T0           - Combustion Chamber Total temperature
    %       P0           - Combustion Chamber Total pressure
    %       R            - Specific gas constant
    %       radiusEarth  - Earth radius
    %       r            - Distance from Earth's center
    %
    %   Output:
    %       thrust       - Total thrust produced
    %
    %   Requires: AreaRatio2Mach, atmosphere
    exitAreaRatio = Ae/At;
    
    m_dot = m_dotFunc(gamma, At, P0, T0, R);

    Me = AreaRatio2Mach(exitAreaRatio, gamma);
    T_ratio = (1+(gamma-1)/2*Me^2)^-1;

    Pe = P0*T_ratio^(gamma/(gamma-1));
    Te = T0*T_ratio;
    Ve = Me*sqrt(gamma*R*Te);
    [Patm, ~] = atmosphere(r-radiusEarth);

    thrust = m_dot*Ve + (Pe - Patm)*Ae;
end