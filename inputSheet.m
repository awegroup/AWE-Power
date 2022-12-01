% Inputs sheet: Performance model framework
clc;
clear;

inputs                = struct();

inputs.Vw             = 1:1:25; %[m/s]
inputs.WA             = 80;
inputs.AR             = 12;
inputs.P_ratedElec    = 1000*1000; %[W]
inputs.massOverride   = 0;
inputs.kiteMass       = 20; %[kg]

inputs.Tmax              = 300; %[kN]
inputs.F_Tmax            = 0.90;
inputs.F_teLength        = 1.2;
inputs.minGroundClear    = 75; % [m]
inputs.Te_matStrength    = 7e8;
inputs.Te_matDensity     = 980; %[kg/m^3] 

inputs.CL_maxAirfoil  = 2.7; % 2.7
inputs.F_CLeff        = 0.8; % 0.8
inputs.CL0_airfoil    = 0.65;
inputs.e              = 0.6;
inputs.CD0            = 0.056;
inputs.CD_te          = 1.2;

inputs.avgPattEle     = 22*pi()/180; % [rad] 
inputs.pattAngRadius  = 12*pi()/180; % [rad] 
inputs.maxRollAngle   = 20*pi()/180; % [rad] % 20
inputs.maxVRI         = 30;
inputs.maxAcc         = 20;

inputs.etaGearbox     = 0.95;
inputs.etaSto         = 0.95;
inputs.etaPE          = 0.98; % Power electronics

inputs.F_peakElecP    = 2.5;
inputs.gravity        = 9.81;
inputs.airDensity     = 1.225; %[kg/m^3]  