% Input Sheet Example: AP3

% Define Input Structure
inputs      = struct();
inputs.name = 'inputFile_example_awePower';

% Model Settings
inputs.numDeltaLelems     = 5;          % Number of stroke length elements [num]
inputs.FgToggle           = 1;          % Gravity toggle: 0 = No, 1 = Yes
inputs.vertWindProfile    = 0;          % Vertical wind profile: 0 = Modelled, 1 = From dataset 

% Wind Parameters
inputs.vw_ref             = 3:1:25;     % Reference wind speeds [m/s]
inputs.h_ref              = 100;        % Reference height [m]
inputs.windShearExp       = 0.143;      % Wind shear exponent [-] (0.143 over land, 0.11 over sea)

% Sample vertical wind profile datasets (for inputs.vertWindProfile = 1)
% Profiles from: https://doi.org/10.1016/j.renene.2022.06.094
% Wind speeds normalized with value at 100m
inputs.h_windDataset          = [10, 20, 40, 60, 80, 100, 120, 140, 150, 160, 180, 200, 220, 250, 300, 500, 600]; 
inputs.v_windDataset_Cabauw   = [0.5412, 0.6074, 0.7686, 0.8685, 0.9414, 1, 1.0481, 1.0864, 1.1028, 1.1172, 1.1441, 1.1657, 1.1839, 1.2065, 1.2327, 1.2666, 1.2641]; 
inputs.v_windDataset_Ijmuiden = [0.8476, 0.8706, 0.9272, 0.9593, 0.9823, 1, 1.0138, 1.0236, 1.0277, 1.0308, 1.0366, 1.0403, 1.0428, 1.0450, 1.0446, 1.0247, 1.0108]; 
inputs.windProfile_h          = inputs.h_windDataset; 
inputs.windProfile_vw         = inputs.v_windDataset_Cabauw;

% System Parameters
inputs.S                = 12;                         % Wing surface area [m^2]
inputs.AR               = 12;                         % Aspect ratio [-]
inputs.b                = sqrt(inputs.AR * inputs.S); % Wing span [m]
inputs.P_ratedElec      = 150 * 1000;                 % Rated electrical power [W]
inputs.massOverride     = 0;                          % Mass override toggle
inputs.kiteMass         = 600;                        % (Overridden) Kite mass [kg]
inputs.peakM2E_F        = 2.5;                        % Peak mechanical to electric cycle power factor [-]

% Tether Parameters
inputs.Ft_max           = 42*1000;                  % Maximum tether force [N]
inputs.Ft_max_SF        = 0.9;                      % Safety factor for maximum tether force
inputs.maxTeLen         = 1000;                     % Maximum tether length [m]
inputs.maxHeight        = 1000;                     % Maximum operational height [m]
inputs.minGroundClear   = 100;                      % Minimum ground clearance [m]
inputs.Te_matStrength   = 7e8;                      % Tether material tensile strength [Pa]
inputs.Te_matDensity    = 980;                      % Tether material density [kg/m^3]

% Aerodynamic Parameters
inputs.Cl_maxAirfoil    = 2.5;                      % Maximum lift coefficient of airfoil [-]
inputs.Cl_eff_F         = 0.8;                      % Effective lift coefficient factor [-]
inputs.Cl0_airfoil      = 0.65;                     % Zero-lift coefficient of airfoil [-]
inputs.e                = 0.6;                      % Oswald efficiency factor [-]
inputs.Cd0              = 0.056;                    % Zero-lift drag coefficient [-]
inputs.Cd_c             = 1.2;                      % Parasitic drag coefficient [-]

% Drum operational Limits
inputs.v_d_max          = 20;                       % Maximum operational speed [m/s]
inputs.a_d_max          = 5;                        % Maximum acceleration [m/s^2]

% Drivetrain Efficiency Parameters
inputs.etaGen.param     = [0.671, -1.4141, 0.9747, 0.7233]; % Generator efficiency parameters [-]
inputs.etaGen.v_max     = inputs.v_d_max;                   % Maximum generator speed [m/s]
inputs.etaGearbox       = 0.95;                             % Gearbox efficiency [-]
inputs.etaSto           = 0.95;                             % Storage efficiency [-]
inputs.etaPE            = 0.95;                             % Power electronics efficiency [-]

% Environmental Constants
inputs.gravity          = 9.81;                     % Acceleration due to gravity [m/s^2]
inputs.airDensity       = 1.225;                    % Air density [kg/m^3]

% Optimisation Initialisation
inputs.nx               = ones(1, inputs.numDeltaLelems); % For discretized stroke length

% Initial Guess
inputs.x0 = [200, deg2rad(30), deg2rad(5), 5*inputs.b, ...                    % stroke length, pattern elevation, cone angle, initial turning radius
             inputs.v_d_max * inputs.nx, ...                          % reel-in speed
             inputs.Cl_maxAirfoil * inputs.Cl_eff_F * inputs.nx, ...  % reel-in C_L
             0.8 * inputs.nx, 90 * inputs.nx, ...                     % reel-out speed, kinematic ratio
             inputs.Cl_maxAirfoil * inputs.Cl_eff_F * inputs.nx];     % reel-out C_L

% Bounds for the initial guess
inputs.lb = [50, deg2rad(1), deg2rad(1), 5*inputs.b, ...
             1 * inputs.nx, 0.1 * inputs.nx, ...
             0.8 * inputs.nx, 1 * inputs.nx, ...
             0.1 * inputs.nx];

inputs.ub = [500, deg2rad(90), deg2rad(60), 100, ...
             inputs.v_d_max * inputs.nx, ...
             inputs.Cl_maxAirfoil * inputs.Cl_eff_F * inputs.nx, ...
             inputs.v_d_max * inputs.nx, 200 * inputs.nx, ...
             inputs.Cl_maxAirfoil * inputs.Cl_eff_F * inputs.nx];
