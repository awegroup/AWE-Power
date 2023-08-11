% Inputs sheet: Performance model framework
inputs                = struct();

inputs.mainPlots      = 1; % 1 = yes, 0 = no
inputs.numDeltaLelems = 1;

inputs.vw_ref         = 1:1:25; %[m/s]
inputs.h_ref          = 100; %[m]
inputs.windShearExp   = 0.143; % 0.143 over land, 0.11 over sea

inputs.WA             = 80;
inputs.AR             = 12;
inputs.P_ratedElec    = 1000*1000; %[W]
inputs.massOverride   = 0;
inputs.kiteMass       = 5000; %[kg]
inputs.F_peakM2Ecyc   = 3;

inputs.Tmax              = 300; %[kN]
inputs.F_Tmax            = 0.8; % 0.8 for gust margin
inputs.maxTeLen          = 1000; %[m]
inputs.F_TeCurve         = 1; %1.05?
inputs.maxHeight         = 600; %[m]
inputs.minGroundClear    = 100; % [m]
inputs.Te_matStrength    = 7e8;
inputs.Te_matDensity     = 980; %[kg/m^3] 

inputs.CL_maxAirfoil  = 2.7; % 2.7
inputs.F_CLeff        = 0.8; % 0.8
inputs.CL0_airfoil    = 0.65;
inputs.e              = 0.6;
inputs.CD0            = 0.056;
inputs.CD_te          = 1.2; %1.1?

inputs.maxVRI         = 30;
inputs.maxAcc         = 20;

inputs.etaGen.param   = [0.671, -1.4141, 0.9747, 0.7233];
inputs.etaGen.Vmax    = max(inputs.maxVRI,25); % 25 =  Possible maximum reel-out speed
inputs.etaGearbox     = 0.9;
inputs.etaSto         = 0.9;
inputs.etaPE          = 0.98; % Power electronics

inputs.gravity        = 9.81;
inputs.airDensity     = 1.225; %[kg/m^3]  