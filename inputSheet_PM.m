% Inputs sheet: Performance model framework
clc;
clear;

inputs                = struct();
inputs.Vw_vector      = 1:1:25; %[m/s]
inputs.WA             = 80;
inputs.AR             = 12;
inputs.P_ratedElec    = 1000*1000; %[W]

inputs.Tmax           = 300; %[kN]
inputs.F_Tmax         = 0.81;
inputs.F_teLength     = 1.5;
inputs.Te_matStrength = 7e8;
inputs.Te_matDensity  = 980; %[kg/m^3] 

inputs.CL_maxAirfoil  = 2.5;
inputs.F_CLeff        = 0.84;
inputs.CL0_airfoil    = 0.65;
inputs.e              = 0.6;
inputs.CD0            = 0.06;
inputs.CD_te          = 1.2;

inputs.avgPattEle     = 22*pi()/180; % [rad] 
inputs.pattAngRadius  = 12*pi()/180; % [rad] 
inputs.maxRollAngle   = 20*pi()/180; % [rad] 
inputs.maxVRI         = 30;
inputs.maxAcc         = 20;

inputs.etaGearbox     = 0.95;
inputs.etaSto         = 0.9;
inputs.etaGrid        = 0.95;

inputs.maxNumPatterns = 1;
inputs.F_peakMechP    = 2.5;
inputs.gravity        = 9.81;
inputs.airDensity     = 1.22;   %% Make it a function of height