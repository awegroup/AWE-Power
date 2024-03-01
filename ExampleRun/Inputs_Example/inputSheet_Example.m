% Inputs sheet: AP3
inputs                = struct();

inputs.numDeltaLelems  = 5; %[num]
inputs.vertWindProfile = 0; % 0 = Modelled, 1 = From dataset 

inputs.evalPoint      = 0; % 0 = Center, 1 = Top, 2 = Side left, 3 = Bottom, 4 = Side right
inputs.FcToggle       = 0; % 0 = No, 1 = Yes
inputs.FgToggle       = 1; % 0 = No, 1 = Yes

inputs.vw_ref         = 1:1:25; %[m/s]
inputs.h_ref          = 100; %[m]
inputs.windShearExp   = 0.143; %[-] % 0.143 over land, 0.11 over sea

inputs.S              = 12; %[m^2]
inputs.AR             = 12; %[-]
inputs.b              = sqrt(inputs.AR*inputs.S);
inputs.P_ratedElec    = 150*1000; %[W]
inputs.massOverride   = 0;
inputs.kiteMass       = 600; %[kg]
inputs.peakM2E_F      = 2.5; %[-]

inputs.Ft_max            = 42; %[kN]
inputs.Ft_max_SF         = 0.8; % 0.8 for gust margin
inputs.maxTeLen          = 2000; %[m]
inputs.maxHeight         = 1000; %[m]
inputs.minGroundClear    = 100; %[m] 
inputs.Te_matStrength    = 7e8;
inputs.Te_matDensity     = 980; %[kg/m^3] 

inputs.Cl_maxAirfoil  = 2.5; %[-] % 2.7
inputs.Cl_eff_F       = 0.8; %[-] % 0.8
inputs.Cl0_airfoil    = 0.65; %[-]
inputs.e              = 0.6; %[-] % 0.6
inputs.Cd0            = 0.056; %[-]
inputs.Cd_c           = 1.2; %[-] % 1.1?

inputs.v_d_max       = 20; %[m/s]
inputs.a_d_max       = 20; %[m/s]

inputs.etaGen.param   = [0.671, -1.4141, 0.9747, 0.7233]; %[-]
inputs.etaGen.v_max   = inputs.v_d_max; %[m/s] % Or can enter a value from supplier
inputs.etaGearbox     = 0.9; %[-]
inputs.etaSto         = 0.9; %[-]
inputs.etaPE          = 0.95; %[-] % Power electronics

inputs.gravity        = 9.81; %[m/s^2]
inputs.airDensity     = 1.225; %[kg/m^3]  

% Optimisation problem data
nx = ones(1,inputs.numDeltaLelems);
% Initial guess
%               [deltaL, avgPattEle,  coneAngle,     Rp_start, v_i,               CL_i,                                    v_o,    kinematicRatio,  CL]
inputs.x0     = [200,    deg2rad(30), deg2rad(5),    50,       inputs.v_d_max*nx, inputs.Cl_maxAirfoil*inputs.Cl_eff_F*nx, 0.8*nx, 90*nx,           inputs.Cl_maxAirfoil*inputs.Cl_eff_F*nx];

% Bounds
inputs.lb     = [50,   deg2rad(1),  deg2rad(1),  50,  1*nx, 0.1*nx, 0.8*nx, 1*nx, 0.1*nx]; % 
inputs.ub     = [500,  deg2rad(90), deg2rad(60), 100, inputs.v_d_max*nx, inputs.Cl_maxAirfoil*inputs.Cl_eff_F*nx, inputs.v_d_max*nx,  200*nx, inputs.Cl_maxAirfoil*inputs.Cl_eff_F*nx]; %


