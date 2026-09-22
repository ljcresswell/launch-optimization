function params_ig = buildParams(params, solidMass, liquidMass)

params_ig = params;
params_ig.solidPropMassInitial  = solidMass;
params_ig.liquidPropMassInitial = liquidMass;
params_ig.massInitial = params.payloadMass  + params.engineMass1 ...
                      + params.engineMass2  ...
                      + (1 + params.massFrac1) * solidMass ...
                      + (1 + params.massFrac2) * liquidMass;

% Required by rocketODE Stage 3 even when delta=0
m_dot2 = m_dotFunc(params.specificHeatRatio2, params.At2, ...
                   params.totalPressure2, params.totalTemp2, params.gasConstant2);
params_ig.liquidBurnDuration = liquidMass * (1 - params.massFrac2) / m_dot2;

end