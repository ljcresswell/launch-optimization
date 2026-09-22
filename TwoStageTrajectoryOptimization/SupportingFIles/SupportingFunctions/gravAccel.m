function g = gravAccel(massEarth, GravConstant, r)
    % gravAccel  Compute gravitational acceleration at a distance from Earth.
    %
    %   g = gravAccel(massEarth, GravConstant, r)
    %
    %   Calculates gravitational acceleration using Newton's law of gravitation
    %   as a function of distance from the Earth's center.
    %
    %   Inputs:
    %       massEarth    - Mass of Earth
    %       GravConstant - Gravitational constant
    %       r            - Distance from Earth's center
    %
    %   Output:
    %       g            - Gravitational acceleration
    
    g = massEarth*GravConstant/r^2;
end