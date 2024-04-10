%% WE section seminar April 2024
% Comparative Assessment of Wind Turbines and Ground-generation Airborne Wind Energy Systems
clc
clearvars
clear global

% Add folders and subfolders to path
addpath(genpath([pwd '\inputSheets']));
addpath(genpath([pwd '\outputFiles'])); 
addpath(genpath([pwd '\src']));

%% Siemens WT
% % Wind turbine characteristcis based on Gamesa SG114-2.0MW: https://en.wind-turbine-models.com/turbines/428-gamesa-g114-2.0mw
% % Model available 2010. Popular Onshore model
% % Low specific power Turbine. SP = 196 W/m^2;
% WT.S_blade   = 142.5; %[m^2] % Assuming average chord of 2.5m
% WT.SweptArea = 10207; %[m^2]
% WT.radius    = sqrt(WT.SweptArea/pi());
% % WT.Cp        = 0.46; % Optimistic value. Electrical Cp
% WT.Cp        = [0 0 0.190 0.370 0.440 0.460 0.470 0.46 0.4 0.31 0.24 0.19 0.15 0.12 0.1 0.08 0.07 0.06 0.05 0.04 0.04 0.03 0.02 0.02 0.01];
% WT.Ct        = [0 0 0.93 0.86 0.83 0.82 0.82 0.78 0.62 0.44 0.32 0.24 0.19 0.15 0.12 0.1 0.09 0.07 0.06 0.06 0.05 0.04 0.03 0.03 0.02];
% WT.P_e_rated = 2e6; %[W]
% WT.HubHt     = 100; %[m]
% WT.CutInWS   = 3; %[m/s] at 100m height
% WT.RatedWS   = 11; %[m/s] at 100m height
% WT.CutOutWS  = 25; %[m/s] at 100m height
% WT.rho_air   = 1.225; %[kg/m^3] at 100m height
% 
% vw_ref = 1:1:25; %[m/s] at 100m height
% 
% WT.Ft        = 0.5.*WT.rho_air.*WT.SweptArea.*WT.Ct.*vw_ref.^2;
% WT.P_e(vw_ref)       = 0.5.*WT.rho_air.*WT.SweptArea.*WT.Cp.*vw_ref.^3;
% WT.P_e(1:WT.CutInWS) = 0;
% % WT.P_e(WT.P_e>=WT.P_e_rated)  = WT.P_e_rated;
% % WT.P_e(WT.RatedWS:WT.CutOutWS) = WT.P_e_rated;
% WT.P_e(WT.RatedWS:21) = WT.P_e_rated;
% WT.P_e(22:25) = [1906 1681 1455 1230].*1e3;

%% Vestas WT
% Wind turbine characteristcis based on Veatas V90-2.0MW:
% https://en.wind-turbine-models.com/turbines/248-vestas-v90-gridstreamer
% https://www.vestas.com/en/pages/backup-2-mw-platform/V90-2-0-MW

WT.SweptArea = 6362; %[m^2]
WT.radius    = sqrt(WT.SweptArea/pi());
WT.S_blade   = 0.5*(3.9+1)*WT.radius; %[m^2] % Max chord is of 3.9m
% WT.Cp        = 0.46*ones(1,25); % Optimistic value. Electrical Cp
WT.Cp        = [0 0 0 0.3 0.39 0.42 0.44 0.44 0.44 0.41 0.37 0.3 0.24 0.19 0.15 0.13 0.11 0.09 0.08 0.07 0.06 0.05 0.04 0.04 0.03];
WT.Ct        = [0 0 0 0.84 0.81 0.81 0.81 0.79 0.73 0.64 0.54 0.41 0.31 0.25 0.2 0.16 0.13 0.11 0.1 0.08 0.07 0.06 0.06 0.05 0.04];
WT.P_e_rated = 2e6; %[W]
WT.HubHt     = 100; %[m]
WT.CutInWS   = 4; %[m/s] at 100m height
WT.RatedWS   = 13.5; %[m/s] at 100m height
WT.CutOutWS  = 25; %[m/s] at 100m height
WT.rho_air   = 1.225; %[kg/m^3] at 100m height

vw_ref = 1:1:25; %[m/s] at 100m height

WT.Ft          = 0.5.*WT.rho_air.*WT.SweptArea.*WT.Ct.*vw_ref.^2;
WT.P_e(vw_ref) = 0.5.*WT.rho_air.*WT.SweptArea.*WT.Cp.*vw_ref.^3;
WT.P_e(1:WT.CutInWS) = 0;
WT.P_e(WT.P_e>=WT.P_e_rated)  = WT.P_e_rated;
WT.P_e(WT.RatedWS:WT.CutOutWS) = WT.P_e_rated;

%% NREL 5MW reference
% https://www.nrel.gov/docs/fy09osti/38060.pdf
% 
% WT.radius    = 63;
% WT.S_blade   = 0.5*(3.54+1.42)*WT.radius; %[m^2] % Max chord is of 3.9m
% WT.Cp        = 0.46*ones(1,25); % Optimistic value. Electrical Cp
% WT.Ct        = 0.8*ones(1,25); % Optimistic value
% WT.P_e_rated = 5e6; %[W]
% WT.HubHt     = 100; %[m]
% WT.CutInWS   = 3; %[m/s] at 100m height
% WT.RatedWS   = 11; %[m/s] at 100m height
% WT.CutOutWS  = 25; %[m/s] at 100m height
% WT.rho_air   = 1.225; %[kg/m^3] at 100m height
% 
% vw_ref = 1:1:25; %[m/s] at 100m height
% 
% WT.Ft          = 0.5.*WT.rho_air.*pi()*WT.radius^2.*WT.Ct.*vw_ref.^2;
% WT.Ft(WT.RatedWS:WT.CutOutWS) = WT.Ft(WT.RatedWS);
% WT.P_e(vw_ref) = 0.5.*WT.rho_air.*pi()*WT.radius^2.*WT.Cp.*vw_ref.^3;
% WT.P_e(1:WT.CutInWS) = 0;
% WT.P_e(WT.P_e>=WT.P_e_rated)  = WT.P_e_rated;
% WT.P_e(WT.RatedWS:WT.CutOutWS) = WT.P_e_rated;


%% AWES
% Defined input sheet
inputSheet_MW_scale;

% Same rated power
inputs.P_ratedElec = WT.P_e_rated; %[W] 

% Equating Blade radius to wing span for a constant AR of 12 for the kite
inputs.b       = WT.radius;
inputs.AR      = 12; 
inputs.S       = inputs.b^2/inputs.AR;
inputs.Ft_max  = 3.5*inputs.S; %[kN]

% Equating 1 blade area to wing area of kite for a constant AR of 12 for the kite
% inputs.S  = WT.S_blade;
% inputs.Ft_max            = 4*inputs.S; %[kN]

% Get results
[inputs, outputs, optimDetails, processedOutputs] = main(inputs, 'inputSheet_MW_scale');

%% Comparison

% Rotor Area Specific Power (RASP)
CRAP.WT = WT.P_e./(0.5*WT.rho_air.*vw_ref.^3*3.*WT.S_blade);
CRAP.AWES = processedOutputs.P_e_avg./(0.5*inputs.airDensity.*vw_ref.^3*inputs.S);

% Weibull
% Mean wind speed at 100m height
weibull.meanWS = 9;  % in m/s

% Shape parameter (you can adjust this based on your specific scenario)
weibull.shapeParam = 2;

% Calculate the scale parameter based on mean wind speed
weibull.scaleParam = weibull.meanWS/ gamma(1 + 1/weibull.shapeParam);

% Calculate Weibull PDF values
weibull.pdf = (weibull.shapeParam / weibull.scaleParam) * ((vw_ref / weibull.scaleParam).^(weibull.shapeParam - 1)) .* exp(-(vw_ref / weibull.scaleParam).^weibull.shapeParam);

% WT AEP
AEP.WT = 8760*trapz(vw_ref, weibull.pdf .* WT.P_e)/1e6; %[MWh]
% AWES AEP
AEP.AWES = 8760*trapz(vw_ref, weibull.pdf .* processedOutputs.P_e_avg)/1e6; %[MWh]
AEP.diffPerc = (AEP.AWES-AEP.WT)*100/AEP.WT

% Moments on ground
moment.WT   =  WT.HubHt*WT.Ft/1e6;
moment.AWES = 100*outputs.d_t*mean(processedOutputs.Ft,2)/1e6;
moment.diffPerc = mean((moment.AWES-moment.WT)*100/moment.WT)
%% Plots
newcolors = [ % 0.25, 0.25, 0.25
    0 0.4470 0.7410
  0.8500 0.3250 0.0980 
  0.4660, 0.6740, 0.1880
  0.9290, 0.6940, 0.1250
  0.4940, 0.1840, 0.5560
  0.6350 0.0780 0.1840
  0.3010 0.7450 0.9330];

% % Wind profile
%   Vref = 10; % m/s
%   z = 10:10:600; % m
%   V = Vref * (z/inputs.h_ref).^inputs.windShearExp/Vref;
% %   MegAWES Onshore location Cabauw. Wind speeds normalized with value at 100m
%   z_MegAWES = [10,20,40,60,80,100,120,140,150,160,180,200,220,250,300,500,600];
%   V_MegAWES_Cabauw = [0.541219682843206,0.607355091566827,0.768630154201962,0.868484406441142,0.941395360902529,1,1.04810058627160,1.08638854381156,1.10277338731106,1.11715868927737,1.14412258234309,1.16573551308321,1.18394938534465,1.20653423381438,1.23266397972046,1.26662287360302,1.26414483994687];
%   V_MegAWES_Ijmuiden = [0.847612611633547,0.870603040595613,0.927240267828556,0.959346286990695,0.982291573490674,1,1.01377720773809,1.02356771954493,1.02766760602000,1.03079423355205,1.03659625208888,1.04025827758100,1.04284618416620,1.04496440015282,1.04461712713371,1.02473617783789,1.01076976884552];
%   figure('units','inch','Position', [15 6 2 2.2])
%   hold on
%   box on
%   grid on
%   plot(V,z,'linewidth',1.5)
%   plot(V_MegAWES_Cabauw,z_MegAWES,'--','linewidth',1.5)
%   plot(V_MegAWES_Ijmuiden,z_MegAWES,'-.','linewidth',1.5)
%   legend('α = 0.143','Cabauw,NL','Ijmuiden,NL');
%   xlim([0.5 1.5])
%   xlabel('Wind speed (-)')
%   ylabel('Height (m)')
%   hold off
% 
% % Plot the Weibull distribution
% figure('units','inch','Position', [4 4 3.5 2.2]);
% plot(vw_ref, weibull.pdf, 'LineWidth', 1.5);
% xlabel('Wind speed at 100m height (m/s)');
% ylabel('Probability (-)');
% title('k = 2, v_{w,mean} = 9m/s');
% grid on;

% Power curve
figure('units','inch','Position', [4 4 3.5 2.2])
colororder(newcolors)
hold on
grid on
box on
plot(vw_ref, processedOutputs.P_e_avg./1e6,':o','linewidth',1,'MarkerSize',3);
plot(vw_ref, WT.P_e./1e6,':^','linewidth',1,'MarkerSize',3);
ylabel('Power (MW)');
legend('AWES','WT','location','southeast');
xlabel('Wind speed at 100m height (m/s)');
hold off

% CRAP
figure('units','inch','Position', [4 4 3.5 2.2])
colororder(newcolors)
hold on
grid on
box on
plot(vw_ref, CRAP.AWES,':o','linewidth',1,'MarkerSize',3);
plot(vw_ref, CRAP.WT,':^','linewidth',1,'MarkerSize',3);
ylabel('C_{RAP} (-)');
legend('AWES','WT','location','southeast');
xlabel('Wind speed at 100m height (m/s)');
hold off

% Cp as defined for Wind turbines (with Swept area)
figure('units','inch','Position', [4 4 3.5 2.2])
colororder(newcolors)
hold on
grid on
box on
plot(vw_ref, mean(processedOutputs.Cp_e_avg,2),':o','linewidth',1,'MarkerSize',3);
plot(vw_ref, WT.Cp,':^','linewidth',1,'MarkerSize',3);
%plot(vw_ref, mean(processedOutputs.Cp_m_o,2),':^','linewidth',1,'MarkerSize',3);
ylabel('C_p (-)');
legend('AWES','WT','location','southeast');
xlabel('Wind speed at 100m height (m/s)');
hold off

% Moment on ground
figure('units','inch','Position', [4 4 3.5 2.2])
colororder(newcolors)
hold on
grid on
box on
plot(vw_ref, moment.AWES,':o','linewidth',1,'MarkerSize',3);
plot(vw_ref,moment.WT,':^','linewidth',1,'MarkerSize',3);
ylabel('Moment (MNm)');
legend('AWES','WT','location','southeast');
xlabel('Wind speed at 100m height (m/s)');
hold off



