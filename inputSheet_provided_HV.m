% Inputs sheet: Performance model framework
clc;
clear;

inputs                = struct();

inputs.Vw             = 1:1:25; %[m/s]
inputs.WA             = 150;
inputs.AR             = 12;
inputs.P_ratedElec    = 2500*1000; %[W]
inputs.massOverride   = 0;
inputs.kiteMass       = 500; %[kg]
inputs.F_peakM2Ecyc   = 2.5;

inputs.Tmax              = 750; %[kN]
inputs.F_Tmax            = 1.0; % 0.8 for gust margin
inputs.F_minTeLen        = 1; %1.2
inputs.maxTeLen          = 1000;
inputs.minGroundClear    = 75; % [m]
inputs.Te_matStrength    = 7e8;
inputs.Te_matDensity     = 980; %[kg/m^3] 

inputs.CL_maxAirfoil  = 2.7; % 2.7
inputs.F_CLeff        = 0.95; % 0.8
inputs.CL0_airfoil    = 0.65;
inputs.e              = 0.6;
inputs.CD0            = 0.056;
inputs.CD_te          = 1.1; %1.2

% inputs.avgPattEle     = 22*pi()/180; % [rad] 
% inputs.pattAngRadius  = 12*pi()/180; % [rad] 
% inputs.maxRollAngle   = 45*pi()/180; % [rad] % 20
inputs.maxVRI         = 30;
inputs.maxAcc         = 20;

inputs.etaGearbox     = 0.95;
inputs.etaSto         = 0.95;
inputs.etaPE          = 0.98; % Power electronics

inputs.targetPRO_elec = 0;
inputs.gravity        = 9.81;
inputs.airDensity     = 1.225; %[kg/m^3]  