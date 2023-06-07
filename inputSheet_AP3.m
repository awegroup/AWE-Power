% Inputs sheet: Performance model framework
inputs                = struct();

inputs.mainPlots      = 1; % 1 = yes, 0 = no
inputs.numDeltaLelems = 5;

inputs.Vw_ref         = 1:1:25; %[m/s]
inputs.h_ref          = 100; %[m]
inputs.windShearExp   = 0.143; % 0.143 over land, 0.11 over sea
% inputs.surfRoughness  = 0.0003; %[m] 0.03 over land, 0.0002 at sea
% inputs.Vw_60          = inputs.Vw_ref.* log(60/inputs.surfRoughness)/log(inputs.h_ref/inputs.surfRoughness);

inputs.WA             = 12; % 12
inputs.AR             = 12;
inputs.P_ratedElec    = 150*1000; %[W]
inputs.massOverride   = 0;
inputs.kiteMass       = 20; %[kg]
inputs.F_peakM2Ecyc   = 3;

inputs.Tmax              = 42; % 42 %[kN]
inputs.F_Tmax            = 0.9; % 0.8 for gust margin
inputs.maxTeLen          = 1000; %[m]
inputs.F_minTeLen        = 1; % 1.2
inputs.maxHeight         = 600; %[m]
inputs.minGroundClear    = 50; %80? % [m]
inputs.Te_matStrength    = 7e8;
inputs.Te_matDensity     = 980; %[kg/m^3] 

inputs.CL_maxAirfoil  = 2.7; % 2.7
inputs.F_CLeff        = 0.8; % 0.8
inputs.CL0_airfoil    = 0.65;
inputs.e              = 0.6; %0.6
inputs.CD0            = 0.056;
inputs.CD_te          = 1.2; %1.1

inputs.maxVRI         = 30;
inputs.maxAcc         = 20;

inputs.etaGen.param   = [0.671, -1.4141, 0.9747, 0.7233];
inputs.etaGen.Vmax    = max(inputs.maxVRI,25); % 25 =  Possible maximum reel-out speed
inputs.etaGearbox     = 0.9;
inputs.etaSto         = 0.9;
inputs.etaPE          = 0.95; % Power electronics

inputs.targetPRO_mech = 0;
inputs.gravity        = 9.81;
inputs.airDensity     = 1.225; %[kg/m^3]  

%%
% z = linspace(10,500,100);
% 
% % Define the wind shear exponent for the power law
% alpha = 0.05;
% 
% % Define the wind speed at 10 m height
% u_10 = 10;
% 
% % Define the roughness length
% z_0 = 0.0003;
% 
% % Create an empty array for the wind speed
% u = zeros(size(z));
% 
% % Use the log law for heights up to 60 m
% ind = z<=60;
% u(ind) = u_10 .* log(z(ind) ./ z_0) ./ log(10/z_0);
% 
% u_60   = u_10 * log(60/ z_0) / log(10/z_0);
% 
% % Use the power law for heights above 60 m
% ind = z>60;
% u(ind) = u_60 .* (z(ind)./60).^alpha;

% % Plot the wind shear profile
% plot(u,z)
% xlabel('Wind Speed (m/s)')
% ylabel('Height (m)')

%%






