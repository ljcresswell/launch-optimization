% Target Orbit Parameters
a                    		 % Semi-major axis of the orbit (distance from Earth's center)
e                   		 % Orbit eccentricity (0 = circular orbit)

% Planet Parameters
radiusEarth         	   % Radius of Earth
massEarth          	   % Mass of Earth
gravConstant       	   % Universal gravitational constant
atmoEndPercent         % Fraction of surface density where atmosphere effectively ends

% Engine Parameters - Stage 1
totalTemp1           	 % Combustion temperature
totalPressure1       	 % Combustion chamber pressure
specificHeatRatio1   % Ratio of specific heats (gamma)
gasConstant1       	  % Gas constant for exhaust gases
Ae1               	  	  % Nozzle exit area
At1               		   % Nozzle throat area

% Engine Parameters - Stage 2
totalTemp2           	  % Combustion temperature
totalPressure2   	  % Combustion chamber pressure
specificHeatRatio2    % Ratio of specific heats (gamma)
gasConstant2             % Gas constant for exhaust gases
Ae2                  	   % Nozzle exit area
At2                   	  % Nozzle throat area

% Rocket Parameters
dragCoefficient       % Aerodynamic drag coefficient
referenceArea         % Reference area for drag calculations
payloadMass           % Mass of the payload
engineMass1           % Mass of stage 1 engine
engineMass2           % Mass of stage 2 engine
massFrac1              % Structural mass fraction for stage 1
massFrac2              % Structural mass fraction for stage 2

% Optimizer Setup
maxThrustAngle          	% Maximum thrust vector angle (radians)
gravityTurnVelocity   	% Velocity at which gravity turn begins
numberOfControlNodes  % Number of control points for optimization
maxSimTime                   % Maximum simulation time
massLowerBoundFactor  % Lower bound factor for mass constraints
massUpperBoundFactor   % Upper bound factor for mass constraints