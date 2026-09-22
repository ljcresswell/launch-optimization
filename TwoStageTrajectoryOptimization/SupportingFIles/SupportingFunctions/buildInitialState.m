function X0 = buildInitialState(launchAngle, params)
% buildInitialState  Construct initial state vector from launch angle.
%
%   launchAngle - flight path angle from horizontal [rad]
%   Returns X0 = [r0; theta0; vr0; vtheta0; m0]

    v0 = params.gravityTurnVelocity;   % m/s gravity turn trigger speed

    r0 = params.radiusEarth + 1000;
    theta0  = 0;
    vr0 = v0 * sin(launchAngle);
    vtheta0 = v0 * cos(launchAngle);

    m0 = params.payloadMass ...
       + params.engineMass1 ...
       + params.engineMass2 ...
       + (1 + params.massFrac1) * params.solidPropMassInitial ...
       + (1 + params.massFrac2) * params.liquidPropMassInitial;

    params.massInitial = m0;
    
    X0 = [r0; theta0; vr0; vtheta0; m0];
end