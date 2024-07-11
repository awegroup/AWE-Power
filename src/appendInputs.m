function [inputs] = appendInputs(inputs)

% Set up selection for wind profiles
if ischar(inputs.windProfile_h)
    inputs.windProfile_h = inputs.(inputs.windProfile_h);
end
if ischar(inputs.windProfile_vw)
    inputs.windProfile_vw = inputs.(inputs.windProfile_vw);
end

inputs.b                = sqrt(inputs.AR * inputs.S); % Wing span [m]

inputs.etaGen.v_max     = inputs.v_d_max; % Maximum generator speed [m/s]

% Optimisation Initialisation
inputs.nx               = ones(1, inputs.numDeltaLelems); % For discretized stroke length

% Initial Guess
inputs.x0 = [200, deg2rad(30), deg2rad(5), 5*inputs.b, ...   % stroke length, pattern elevation, cone angle, initial turning radius
    inputs.v_d_max * inputs.nx, ...                          % reel-in speed
    inputs.Cl_maxAirfoil * inputs.Cl_eff_F * inputs.nx, ...  % reel-in C_L
    0.8 * inputs.nx, 90 * inputs.nx, ...                     % reel-out speed, kinematic ratio
    inputs.Cl_maxAirfoil * inputs.Cl_eff_F * inputs.nx];     % reel-out C_L

% Bounds for the initial guess
inputs.lb = [50, deg2rad(1), deg2rad(1), 5*inputs.b, ...
    1 * inputs.nx, 0.1 * inputs.nx, ...
    0.8 * inputs.nx, 1 * inputs.nx, ...
    0.1 * inputs.nx];

inputs.ub = [500, deg2rad(90), deg2rad(60), 10*inputs.b, ...
    inputs.v_d_max * inputs.nx, ...
    inputs.Cl_maxAirfoil * inputs.Cl_eff_F * inputs.nx, ...
    inputs.v_d_max * inputs.nx, 200 * inputs.nx, ...
    inputs.Cl_maxAirfoil * inputs.Cl_eff_F * inputs.nx];

end