% Inputs sheet
inputs                = struct();

inputs.mainPlots      = 1; % 0 = No, 1 = Yes
inputs.numDeltaLelems = 1;

inputs.evalPoint      = 0; % 0 = Center, 1 = Top, 2 = Side left, 3 = Bottom, 4 = Side right
inputs.FcToggle       = 0; % 0 = No, 1 = Yes
inputs.FgToggle       = 1; % 0 = No, 1 = Yes

inputs.vw_ref         = 11:1:15; %[m/s]
inputs.h_ref          = 100; %[m]
inputs.windShearExp   = 0; % 0.143 over land, 0.11 over sea

inputs.S              = 12; % 12
inputs.AR             = 12;
inputs.P_ratedElec    = 150*1000; %[W]
inputs.massOverride   = 0;
inputs.kiteMass       = 600; %[kg]
inputs.peakM2E_F      = 2.5;

inputs.Ft_max            = 42; % 42 %[kN]
inputs.Ft_max_SF         = 0.8; % 0.8 for gust margin
inputs.maxTeLen          = 2000; %[m]
inputs.maxHeight         = 800; %[m]
inputs.minGroundClear    = 0; %80? % [m]
inputs.Te_matStrength    = 7e8;
inputs.Te_matDensity     = 980; %[kg/m^3] 

inputs.Cl_maxAirfoil  = 2.7; % 2.7
inputs.Cl_eff_F        = 0.8; % 0.8
inputs.Cl0_airfoil    = 0.65;
inputs.e              = 0.6; %0.6
inputs.Cd0            = 0.056;
inputs.Cd_c           = 1.2; %1.1

inputs.v_d_max       = 20;
inputs.a_d_max       = 20;

inputs.etaGen.param   = [0.671, -1.4141, 0.9747, 0.7233];
inputs.etaGen.v_max   = inputs.v_d_max; % Or can enter a value from supplier
inputs.etaGearbox     = 0.9;
inputs.etaSto         = 0.9;
inputs.etaPE          = 0.95; % Power electronics

inputs.gravity        = 9.81;
inputs.airDensity     = 1.225; %[kg/m^3]  
