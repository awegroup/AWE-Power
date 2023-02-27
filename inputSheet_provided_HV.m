% Inputs sheet: Performance model framework
clc;
clear;

inputs                = struct();

inputs.Vw_ref         = 1:1:25; %[m/s]
inputs.h_ref          = 100; %[m]
inputs.windShearExp   = 0.143; % 0.143 over land, 0.11 over sea
% inputs.surfRoughness  = 0.0003; %[m] 0.03 over land, 0.0002 at sea
% inputs.Vw_60          = inputs.Vw_ref.* log(60/inputs.surfRoughness)/log(inputs.h_ref/inputs.surfRoughness);

inputs.WA             = 220; %150
inputs.AR             = 12;
inputs.P_ratedElec    = 2500*1000; %[W]
inputs.massOverride   = 0;
inputs.kiteMass       = 500; %[kg]
inputs.F_peakM2Ecyc   = 2.5; %2.5

inputs.Tmax              = 900; %[kN]
inputs.F_Tmax            = 1.0; % 0.8 for gust margin
inputs.F_minTeLen        = 1; %1.2
inputs.maxTeLen          = 1200;
inputs.maxHeight         = 600;
inputs.minGroundClear    = 80; % [m]
inputs.Te_matStrength    = 7e8;
inputs.Te_matDensity     = 980; %[kg/m^3] 

inputs.CL_maxAirfoil  = 2.7; % 2.7
inputs.F_CLeff        = 0.95; % 0.8
inputs.CL0_airfoil    = 0.65;
inputs.e              = 0.6;
inputs.CD0            = 0.056;
inputs.CD_te          = 1.1; %1.2

inputs.maxVRI         = 30;
inputs.maxAcc         = 20;

inputs.etaGen.param   = [0.671, -1.4141, 0.9747, 0.7233];
inputs.etaGen.Vmax    = max(inputs.maxVRI,25); % 25 =  Possible maximum reel-out speed
inputs.etaGearbox     = 0.9;
inputs.etaSto         = 0.9;
inputs.etaPE          = 0.95; % Power electronics
inputs.targetPRO_elec = 0;
inputs.gravity        = 9.81;
inputs.airDensity     = 1.225; %[kg/m^3]  


%%
% figure()
% Vw = mean(postProRes.Vw,2);
% plot(round(Vw,2), postProRes.Pcycle_elec./10^3,'-','linewidth',1.2);
