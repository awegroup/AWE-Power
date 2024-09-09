function [inputs] = appendInputs(inputs)
  % APPENDINPUTS Appends additional fields and parameters to the 'inputs' structure
  %
  %   This function sets up wind profile parameters, calculates various
  %   aerodynamic and optimization-related parameters, and establishes
  %   initial guesses and bounds for an optimization problem related to
  %   flight dynamics or a similar application.
  %
  %   INPUTS:
  %       inputs - a structure containing relevant fields such as wind profiles,
  %                geometry parameters, aerodynamic coefficients, and constraints.
  %
  %   OUTPUTS:
  %       inputs - a modified structure with additional calculated fields and
  %                updated wind profile parameters.
  %
  %   The function performs the following tasks:
  %   1. Updates wind profiles if specified as strings.
  %   2. Checks if the Ostwald efficiency factor (e) is provided. If not,
  %      calculates it based on the aircraft's aspect ratio (AR) using Anderson's
  %      equation from the 6th Edition of his book, section 6.7.2.
  %   3. Calculates the wing span (b) from the aspect ratio (AR) and wing area (S).
  %   4. Sets up the maximum generator speed.
  %   5. Initializes discretized stroke length for optimization.
  %   6. Provides an initial guess for the optimization variables, including stroke
  %      length, pattern elevation, cone angle, turning radius, reel-in and reel-out
  %      speeds, and lift coefficients.
  %   7. Establishes lower and upper bounds for the optimization problem.
  %
  %   The initial guess and bounds are returned within the modified 'inputs' structure.

  % Set up selection for wind profiles
  if ischar(inputs.windProfile_h)
    inputs.windProfile_h = inputs.(inputs.windProfile_h);
  end
  if ischar(inputs.windProfile_vw)
    inputs.windProfile_vw = inputs.(inputs.windProfile_vw);
  end

  % Check if inputs.e is already present
  if inputs.givenOstwaldEff == 0
    inputs.e = 1.78 * (1 - 0.045 * inputs.AR^0.68) - 0.64; % From Anderson book (6th Edition) section 6.7.2
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

  inputs.ub = [1000, deg2rad(90), deg2rad(60), 10*inputs.b, ...
    inputs.v_d_max * inputs.nx, ...
    inputs.Cl_maxAirfoil * inputs.Cl_eff_F * inputs.nx, ...
    inputs.v_d_max * inputs.nx, 200 * inputs.nx, ...
    inputs.Cl_maxAirfoil * inputs.Cl_eff_F * inputs.nx];

end