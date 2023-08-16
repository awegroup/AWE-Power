% Inputs sheet: Performance model framework
inputs                = struct();

inputs.mainPlots      = 1; % 1 = yes, 0 = no
inputs.numDeltaLelems = 5;

inputs.vw_ref         = 1:1:25; %[m/s]
inputs.h_ref          = 100; %[m]
inputs.windShearExp   = 0.143; % 0.143 over land, 0.11 over sea

inputs.S              = 12; % 12
inputs.AR             = 12;
inputs.P_ratedElec    = 150*1000; %[W]
inputs.massOverride   = 0;
inputs.kiteMass       = 50; %[kg]
inputs.peakM2E_F      = 2.5;

inputs.Ft_max            = 42; % 42 %[kN]
inputs.Ft_max_SF         = 0.8; % 0.8 for gust margin
inputs.maxTeLen          = 1000; %[m]
inputs.TeSag_F           = 1; %1.05?
inputs.maxHeight         = 800; %[m]
inputs.minGroundClear    = 80; %80? % [m]
inputs.Te_matStrength    = 7e8;
inputs.Te_matDensity     = 980; %[kg/m^3] 

inputs.CL_maxAirfoil  = 2.7; % 2.7
inputs.CLeff_F        = 0.8; % 0.8
inputs.CL0_airfoil    = 0.65;
inputs.e              = 0.6; %0.6
inputs.CD0            = 0.056;
inputs.CD_t           = 1.2; %1.1

% inputs.maxVRI         = 30;
inputs.vk_r_i_max       = 30;
inputs.winchAcc_max     = 20;

inputs.etaGen.param   = [0.671, -1.4141, 0.9747, 0.7233];
inputs.etaGen.v_max   = max(inputs.vk_r_i_max,25); % 25 =  Possible maximum reel-out speed
inputs.etaGearbox     = 0.9;
inputs.etaSto         = 0.9;
inputs.etaPE          = 0.95; % Power electronics

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






