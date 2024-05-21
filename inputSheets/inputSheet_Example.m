% Input Sheet Example (AP3)
inputs                = struct();
inputs.name           = 'inputSheet_Example';

inputs.numDeltaLelems  = 5; %[num]
inputs.FgToggle        = 1; % 0 = No, 1 = Yes
inputs.vertWindProfile = 0; % 0 = Modelled, 1 = From dataset 

inputs.vw_ref         = 1:1:25; %[m/s]
inputs.h_ref          = 100; %[m]
inputs.windShearExp   = 0.143; %[-] % 0.143 over land, 0.11 over sea

% Sample vertical wind profile datasets (i.e. for inputs.vertWindProfile = 1)
% Profiles from: https://doi.org/10.1016/j.renene.2022.06.094. Wind speeds normalized with value at 100m. 
inputs.h_windDataset          = [10,20,40,60,80,100,120,140,150,160,180,200,220,250,300,500,600]; % [m]
inputs.v_windDataset_Cabauw   = [0.541219682843206,0.607355091566827,0.768630154201962,0.868484406441142,0.941395360902529,1,1.04810058627160,1.08638854381156,1.10277338731106,1.11715868927737,1.14412258234309,1.16573551308321,1.18394938534465,1.20653423381438,1.23266397972046,1.26662287360302,1.26414483994687];
inputs.v_windDataset_Ijmuiden = [0.847612611633547,0.870603040595613,0.927240267828556,0.959346286990695,0.982291573490674,1,1.01377720773809,1.02356771954493,1.02766760602000,1.03079423355205,1.03659625208888,1.04025827758100,1.04284618416620,1.04496440015282,1.04461712713371,1.02473617783789,1.01076976884552];
inputs.windProfile_h          = inputs.h_windDataset; 
inputs.windProfile_vw         = inputs.v_windDataset_Cabauw;

inputs.S              = 12; %[m^2]
inputs.AR             = 12; %[-]
inputs.b              = sqrt(inputs.AR*inputs.S); %[m]
inputs.P_ratedElec    = 150*1000; %[W]
inputs.massOverride   = 0;
inputs.kiteMass       = 600; %[kg]
inputs.peakM2E_F      = 2.5; %[-]

inputs.Ft_max            = 42; %[kN]
inputs.Ft_max_SF         = 0.8; % 0.8 for gust margin
inputs.maxTeLen          = 1000; %[m]
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
inputs.a_d_max       = 5; %[m/s^2]

inputs.etaGen.param   = [0.671, -1.4141, 0.9747, 0.7233]; %[-]
inputs.etaGen.v_max   = inputs.v_d_max; %[m/s] % Or can enter a value from supplier
inputs.etaGearbox     = 0.95; %[-]
inputs.etaSto         = 0.95; %[-]
inputs.etaPE          = 0.95; %[-] % Power electronics

inputs.gravity        = 9.81; %[m/s^2]
inputs.airDensity     = 1.225; %[kg/m^3]  

% Optimisation problem data
inputs.nx = ones(1,inputs.numDeltaLelems);
% Initial guess
%               [deltaL, avgPattEle,  coneAngle,     Rp_start, v_i,                      CL_i,                                           v_o,           kinematicRatio,         CL]
inputs.x0     = [200,    deg2rad(30), deg2rad(5),    50,       inputs.v_d_max*inputs.nx, inputs.Cl_maxAirfoil*inputs.Cl_eff_F*inputs.nx, 0.8*inputs.nx, 90*inputs.nx,           inputs.Cl_maxAirfoil*inputs.Cl_eff_F*inputs.nx];

% Bounds
inputs.lb     = [50,   deg2rad(1),  deg2rad(1),  50,  1*inputs.nx, 0.1*inputs.nx, 0.8*inputs.nx, 1*inputs.nx, 0.1*inputs.nx]; % 
inputs.ub     = [500,  deg2rad(90), deg2rad(60), 100, inputs.v_d_max*inputs.nx, inputs.Cl_maxAirfoil*inputs.Cl_eff_F*inputs.nx,...
    inputs.v_d_max*inputs.nx,  200*inputs.nx, inputs.Cl_maxAirfoil*inputs.Cl_eff_F*inputs.nx]; %




