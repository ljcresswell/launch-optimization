function validateParams(params)
% validateParams  Check physical consistency of params before integration.
%
%   Throws descriptive errors for any condition that would cause the
%   integrator to fail or produce nonsense results.

    fprintf('=== Pre-flight Parameter Validation ===\n\n');
    passCount = 0;
    
    % --- TWR check ---
    T_surface = thrust(params.Ae1, params.At1, params.specificHeatRatio1, ...
                       params.totalTemp1, params.totalPressure1, params.gasConstant1, ...
                       params.radiusEarth, params.radiusEarth);
    g0        = gravAccel(params.massEarth, params.gravConstant, params.radiusEarth);
    TWR       = T_surface / (params.massInitial * g0);
    fprintf('TWR (stage 1 at launch): %.3f\n', TWR);
    if TWR < 1.0
        error(['TWR = %.3f < 1.0 — rocket cannot lift off.\n' ...
               'Need thrust > %.1f N, currently %.1f N.\n' ...
               'Increase totalPressure1 or At1.'], ...
               TWR, params.massInitial * g0, T_surface);
    elseif TWR < 1.3
        warning('TWR = %.3f is low (< 1.3) — ascent will be sluggish.', TWR);
    else
        fprintf('  PASS: TWR = %.3f\n', TWR);
        passCount = passCount + 1;
    end

    % --- Nozzle geometry: exit area must be larger than throat ---
    if params.Ae1 <= params.At1
        error('Stage 1 nozzle: Ae1 (%.4f) must be greater than At1 (%.4f).', ...
              params.Ae1, params.At1);
    else
        fprintf('  PASS: Stage 1 nozzle geometry (Ae1 > At1)\n');
        passCount = passCount + 1;
    end
    if params.Ae2 <= params.At2
        error('Stage 2 nozzle: Ae2 (%.4f) must be greater than At2 (%.4f).', ...
              params.Ae2, params.At2);
    else
        fprintf('  PASS: Stage 2 nozzle geometry (Ae2 > At2)\n');
        passCount = passCount + 1;
    end

    % --- Propellant mass must be positive ---
    if params.solidPropMassInitial <= 0
        error('solidPropMassInitial must be positive. Got: %.2f', ...
              params.solidPropMassInitial);
    else
        fprintf('  PASS: Stage 1 propellant mass positive\n');
        passCount = passCount + 1;
    end
    if params.liquidPropMassInitial <= 0
        error('liquidPropMassInitial must be positive. Got: %.2f', ...
              params.liquidPropMassInitial);
    else
        fprintf('  PASS: Stage 2 propellant mass positive\n');
        passCount = passCount + 1;
    end

    % --- Propellant mass must not exceed total mass ---
    totalProp = params.solidPropMassInitial + params.liquidPropMassInitial;
    if totalProp >= params.massInitial
        error(['Total propellant (%.2f kg) >= massInitial (%.2f kg).\n' ...
               'No mass left for structure or payload.'], ...
               totalProp, params.massInitial);
    else
        fprintf('  PASS: Propellant mass < massInitial\n');
        passCount = passCount + 1;
    end

    % --- Payload and structural masses must be positive ---
    if params.payloadMass <= 0
        error('payloadMass must be positive. Got: %.2f', params.payloadMass);
    end
    if params.engineMass1 <= 0
        error('engineMass1 must be positive. Got: %.2f', params.engineMass1);
    end
    if params.engineMass2 <= 0
        error('engineMass2 must be positive. Got: %.2f', params.engineMass2);
    end
    fprintf('  PASS: Payload and engine masses positive\n');
    passCount = passCount + 1;

    % --- Mass fractions must be non-negative ---
    if params.massFrac1 < 0
        error('massFrac1 must be >= 0. Got: %.4f', params.massFrac1);
    end
    if params.massFrac2 < 0
        error('massFrac2 must be >= 0. Got: %.4f', params.massFrac2);
    end
    fprintf('  PASS: Mass fractions non-negative\n');
    passCount = passCount + 1;

    % --- Target orbit must be above atmosphere ---
    r_periapsis = params.a * (1 - params.e);
    if r_periapsis <= params.rEndAtmo
        error(['Orbit periapsis (%.2f m) is inside atmosphere (rEndAtmo = %.2f m).\n' ...
               'Increase semi-major axis or reduce eccentricity.'], ...
               r_periapsis, params.rEndAtmo);
    else
        fprintf('  PASS: Orbit periapsis above atmosphere\n');
        passCount = passCount + 1;
    end

    % --- Eccentricity must be physical ---
    if params.e < 0 || params.e >= 1
        error('Eccentricity must be in [0, 1). Got: %.4f', params.e);
    else
        fprintf('  PASS: Eccentricity in valid range\n');
        passCount = passCount + 1;
    end

    % --- Launch angle must be between 0 and pi/2 ---
    if params.launchAngleIG <= 0 || params.launchAngleIG > pi/2
        error('launchAngleIG must be in (0, pi/2]. Got: %.4f rad', ...
              params.launchAngleIG);
    else
        fprintf('  PASS: Launch angle in valid range\n');
        passCount = passCount + 1;
    end
    if params.gravityTurnVelocity <= 0
        error('gravityTurnVelocity must be positive. Got: %.2f', params.gravityTurnVelocity);
    end
    if params.gravityTurnVelocity > 500
        warning('gravityTurnVelocity = %.1f m/s is high — gravity turn may start very late.', ...
                params.gravityTurnVelocity);
    end
    fprintf('  PASS: gravityTurnVelocity = %.1f m/s\n', params.gravityTurnVelocity);
    passCount = passCount + 1;

    % --- Stage 2 burn duration check: enough propellant to matter ---
    mdot2    = m_dotFunc(params.specificHeatRatio2, params.At2, params.totalPressure2, ...
                         params.totalTemp2, params.gasConstant2);
    t_burn2  = params.liquidPropMassInitial / mdot2;
    if t_burn2 < 1.0
        warning('Stage 2 burn duration is only %.2f s — consider increasing liquidPropMassIG.', ...
                t_burn2);
    else
        fprintf('  PASS: Stage 2 burn duration = %.1f s\n', t_burn2);
        passCount = passCount + 1;
    end

    fprintf('\nValidation complete: %d checks passed.\n\n', passCount);
end