%% Paper plots

%% Power curve comparison plot
vw = inputs.vw_ref(1):processedOutputs.vw_100m_operRange(end);
vw_loyd = vw;
% Power curve comparison plot
% Loyd
CL_loyd = inputs.Cl_maxAirfoil*inputs.Cl_eff_F;
CD_loyd = inputs.Cd0 + (CL_loyd-inputs.Cl0_airfoil)^2/(pi()*inputs.AR*inputs.e);
for i = 1:length(vw_loyd)
  P_Loyd(i) = (4/27)*(CL_loyd^3/CD_loyd^2)*(1/2)*inputs.airDensity*inputs.S*(vw_loyd(i)^3);
  %   if P_Loyd(i)>inputs.P_ratedElec
  %     P_Loyd(i) = inputs.P_ratedElec;
  %   end
end
% AP3 6DoF simulation results
AP3.PC.ws    = [7.32E+00, 8.02E+00, 9.05E+00, 1.00E+01, 1.10E+01, 1.20E+01, ...
  1.30E+01, 1.41E+01, 1.50E+01, 1.60E+01, 1.70E+01, 1.80E+01, 1.90E+01]; %[ms^{-1}]
AP3.PC.power = [0.00E+00, 7.01E+03, 2.37E+04, 4.31E+04, 6.47E+04, 8.46E+04, ...
  1.02E+05, 1.20E+05, 1.34E+05, 1.49E+05, 1.50E+05, 1.50E+05, 1.50E+05]./10^3; %[kW]

figure('units','inch','Position', [0.2 6.5 3.5 2.2])
%colororder(newcolors)
hold on
grid on
box on
plot(vw_loyd, P_Loyd./10^3,'-','linewidth',1);
plot(vw, processedOutputs.P_e_avg./10^3,':s','linewidth',1,'MarkerSize',3);
plot(AP3.PC.ws, AP3.PC.power,':^k','linewidth',1,'MarkerSize',3);
ylabel('Power (kW)');
ylim([0 155]);
legend('Loyd Reel-out','QSM (this study)','6-DOF','location','southeast');
xlabel('Wind speed at 100 m height (ms^{-1})');
xlim([1 processedOutputs.vw_100m_operRange(end)]);
hold off

%% Effect of gravity
clc; clearvars;

% Add folders and subfolders to path
addpath(genpath('C:/PhD/GitHubRepo/AWE-power/inputSheets'));
addpath(genpath('C:/PhD/GitHubRepo/AWE-power/outputFiles'));
addpath(genpath('C:/PhD/GitHubRepo/AWE-power/src'));

for i = 1:2
  inputs = loadInputs('inputFile_example_awePower.yml');
  inputs.FgToggle = i-1;
  [inputs, outputs, optimDetails, processedOutputs] = main_awePower(inputs);
  gravityEffect(i)     =  processedOutputs;
end

% Power curve plot
figure('units','inch','Position', [5 5 3.5 2.5])
hold on
grid on
box on
plot(gravityEffect(1).P_e_avg./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(gravityEffect(2).P_e_avg./1e3,'^-','linewidth',1,'MarkerSize',2);
legend('P_{no-gravity}','P_{gravity}');
ylabel('Power (kW)');
xlabel('Wind speed at 100 m height (ms^{-1})');
xlim([0 processedOutputs.vw_100m_operRange(end)]);
hold off

% Cycle power rep. plots
col_b = [0 0.4470 0.7410];
col_o = [0.8500 0.3250 0.0980];
cols = [col_b col_o];
windSpeeds = [6, 15];
for i = windSpeeds
  idx = find(inputs.vw_ref == i);
  figure('units','inch','Position', [4.2 6.5 3.5 2.2])
  hold on
  grid on
  box on
  % Plot your data
  plot(gravityEffect(1).cyclePowerRep(i).t_inst, gravityEffect(1).cyclePowerRep(i).P_e_inst,'Color',col_b,'linewidth',1.2, 'DisplayName', 'P_{no-gravity}');
  plot(gravityEffect(2).cyclePowerRep(i).t_inst, gravityEffect(2).cyclePowerRep(i).P_e_inst,'Color',col_o,'linewidth',1.2, 'DisplayName', 'P_{gravity}');
  yline(gravityEffect(1).P_e_avg(i)/10^3, '--','Color',col_b,'linewidth',1.2, 'DisplayName', 'P_{no-gravity}','HandleVisibility', 'off');
  yline(gravityEffect(2).P_e_avg(i)/10^3, '--','Color',col_o,'linewidth',1.2, 'DisplayName', 'P_{gravity}','HandleVisibility', 'off');
  % Set yline positions within the range of your data
  yline(0, '-', 'linewidth', 1, 'HandleVisibility', 'off');
  % Dummy line for average power entry in legend
  plot(NaN, NaN, '--', 'Color', 'k', 'linewidth', 1, 'DisplayName', 'P_{avg}');
  ylabel('Power (kW)');
  %ylim([1.3*min(gravityEffect(2).cyclePowerRep(i).P_e_inst) 1.3*max(gravityEffect(1).cyclePowerRep(i).P_e_inst)])
  xlabel('Time (s)');
  legend('Location', 'Best');
  title(strcat('Wind speed at 100 m:', num2str(processedOutputs.cyclePowerRep(idx).ws), ' ms^{-1}'));
  hold off
end


%% Effect of scaling on fixed wind speed of 12 ms^{-1}
% Model settings
% Zero ground clearance, no constraint on rated, max power, tether length, height, patterns, etc.
% Enforce max. tether force constraint
% No reel-in glide ratio equality constraint
% Change objective function to max. P_m_o
% Increase max. iterations in the optimisation settings
% Hard code cut-in to 7

% Fixed WA: Finding opt tether size for given WA
clc
clearvars
inputs = loadInputs('inputFile_scalingEffects_awePower.yml');
inputs.vw_ref         = 12; %[ms^{-1}]
inputs.S = 100;
for i = 1:7
  inputs.Ft_max = 4000*inputs.S + i*100000;
  fixedWA.Ft_max(i) = inputs.Ft_max;
  [inputs, outputs, optimDetails, processedOutputs] = main_awePower(inputs);
  fixedWA.m_k(i)         = processedOutputs.m_k;
  fixedWA.zetaMech(i)    = processedOutputs.zetaMech;
  fixedWA.d_te(i)        = processedOutputs.d_te;
  fixedWA.Ft_to_S(i)     = fixedWA.Ft_max(i)/inputs.S;
  clear global
end
figure('units','inch','Position', [4 4 3.5 2.2])
hold on
grid on
box on
xlabel('Tether diameter (cm)');
yyaxis left
plot(fixedWA.d_te.*100,fixedWA.zetaMech,'o-','linewidth',1,'MarkerSize',3);
ylabel('Power harvesting factor (-)');
yyaxis right
plot(fixedWA.d_te.*100,fixedWA.m_k,'^-','linewidth',1,'MarkerSize',3);
ylabel('Kite mass (kg)');
title(['S = ' num2str(inputs.S) 'm^2']);
hold off

% Fixed Ft_max: Finding opt WA for given Tether size
clc
clearvars
inputs = loadInputs('inputFile_scalingEffects_awePower.yml');
inputs.vw_ref         = 12; %[ms^{-1}]
FixedFt_max.Ft_max = 400000;
for i =1:6
  inputs.S = 10 + 10*i;
  FixedFt_max.S(i)       = inputs.S;
  inputs.Ft_max          = FixedFt_max.Ft_max;
  FixedFt_max.Ft_max     = inputs.Ft_max;
  [inputs, outputs, optimDetails, processedOutputs] = main_awePower(inputs);
  FixedFt_max.m_k(i)     = processedOutputs.m_k;
  FixedFt_max.zeta(i)    = processedOutputs.zetaMech;
  FixedFt_max.d_te       = processedOutputs.d_te;
end
figure('units','inch','Position', [12 4 3.5 2.2])
hold on
grid on
box on
xlabel('Wing area (m^2)');
yyaxis left
plot(FixedFt_max.S,FixedFt_max.zeta,'o-','linewidth',1,'MarkerSize',3);
ylabel('Power harvesting factor (-)');
yyaxis right
plot(FixedFt_max.S,FixedFt_max.m_k,'^-','linewidth',1,'MarkerSize',3);
ylabel('Kite mass (kg)');
title(['F_{t,max} = ' num2str(FixedFt_max.Ft_max) 'kN (d_{t} = ' num2str(round(FixedFt_max.d_te*100,1)) 'cm)']);
hold off


%% Effects of scaling on full operational wind speed range

% Full Power curve plots for fixed S, increasing Ft_max
clc
clearvars
inputs = loadInputs('inputFile_scalingEffects_awePower.yml');
inputs.P_ratedElec    = 1500 * 1000;
numEvals = 4;
for i = 1:numEvals
  inputs.Ft_max = 1000*inputs.S + i*100000; %inputs.S + i*100000;
  Ft_max(i)     = inputs.Ft_max;
  [inputs, outputs, optimDetails, processedOutputs] = main_awePower(inputs);
  fullPC_incrFt_max(i)  =  processedOutputs;
end
% Plotting
figSize = [5 5 3.2 2];
figure('units','inch','Position', figSize)
hold on
grid on
box on
legendLabels = cell(numEvals,1);
for i = 1:numEvals
  plot(fullPC_incrFt_max(i).P_e_avg./1e6,'o-','linewidth',1,'MarkerSize',2);
  legendLabels{i} = ['F_{t,max}= ' num2str(Ft_max(i)/1e3) 'kN'];
end
xlabel('Wind speed at 100 m height (ms^{-1})');
ylabel('Power (MW)');
title(strcat('S = ',num2str(inputs.S),'m^2'));
legend(legendLabels, 'Location', 'Best');
hold off


%% Full Power curve plots for increasing S, fixed Ft_max/S, no mechanical power limit
clc
clearvars
inputs = loadInputs('inputFile_scalingEffects_awePower.yml');
inputs.P_ratedElec    = 1000*1000; %[W]
Ft_maxByS = 3000;
numEvals = 4;
legendLabels = cell(numEvals,1);
% Computation loop
for i = 1:numEvals
  inputs.S = i*30;
  S(i)     = inputs.S;
  inputs.Ft_max = Ft_maxByS*inputs.S;
  [inputs, outputs, optimDetails, processedOutputs] = main_awePower(inputs);
  fullPC_incrS(i) =  processedOutputs;
end
% Plotting loop
figSize = [5 5 3.2 2];
figure('units','inch','Position', figSize)
hold on
grid on
box on
for i = 1:numEvals
  plot(fullPC_incrS(i).P_e_avg./1e6,'o-','linewidth',1,'MarkerSize',2);
  legendLabels{i} = ['S = ' num2str(S(i)) 'm^2'];
end
xlabel('Wind speed at 100 m height (ms^{-1})');
%xlim([0 processedOutputs.vw_100m_operRange(end)]);
ylabel('Power (MW)');
title(strcat('F_{t,max}/S = ',num2str(Ft_maxByS/1e3),'kN/m^2'));
legend(legendLabels, 'Location', 'Best');
hold off

%% Mass scaling curve plots

a1     = 0.002415;
a2     = 0.0090239;
b1     = 0.17025;
b2     = 3.2493;
k1     = 5;
c1     = 0.46608;
d1     = 0.65962;
k2     = 1.1935;
AR_ref = 12;

r  = 3.5; % Tmax/inputs.WA %[kN/m^2]
AR = 12;

a = a1*r + a2;
b = b1*r + b2;

WA = linspace(0, 310,10);
m_kite = 10*(a*WA.^2 +b*WA)*(c1*(AR/AR_ref)^2-d1*(AR/AR_ref)+k2);


w.CL = 2;
par.u0_max = 100;
q = 1/2 * 1.225 * par.u0_max^2;
W =  r*WA*1000;% q * (WA.*w.CL);% + h.A* h.CL);
ff = 0.8;
tc = 0.2;
nz = 6;

w.S = 0;
w.AR = 12;

m.w = ff * 0.036.* WA.^0.758 * (w.AR/cos(w.S)^2)^0.6 * q^0.006 * (100*tc/cos(w.S))^(-0.3) .* (nz * W).^0.49;


% Data from literature
MegAWES = [150.45, 6885];
AP2     = [3, 35];
AP3     = [12, 475];
AP4     = [35^2/12, 3500];
AP5_Luuk    = [300, 20000];
AP5_AP    = [300, 34000];
M600    = [32.9, 1730.8];
MX2     = [54, 1850];
Haas    = [138.5, 8000];
M5      = [177.3, 9900];

newcolors = [ % 0.25, 0.25, 0.25
  0 0.4470 0.7410
  0.8500 0.3250 0.0980
  0.4660, 0.6740, 0.1880
  0.4940, 0.1840, 0.5560
  0.9290, 0.6940, 0.1250

  0.6350 0.0780 0.1840
  0.3010 0.7450 0.9330];


figure()
hold on
grid on
box on
colororder(newcolors)
plot(MegAWES(1), MegAWES(2)/10^3,'p','markersize',7)
plot(AP2(1), AP2(2)/10^3,'s')
plot(AP3(1), AP3(2)/10^3,'d')
plot(AP4(1), AP4(2)/10^3,'*','markersize',7)
plot(AP5_Luuk(1), AP5_Luuk(2)/10^3,'o')
plot(AP5_AP(1), AP5_AP(2)/10^3,'x')
plot(M600(1), M600(2)/10^3,'^','markersize',5)
plot(MX2(1), MX2(2)/10^3,'v','markersize',5)
plot(M5(1), M5(2)/10^3,'<','markersize',5)
plot(Haas(1), Haas(2)/10^3,'>','markersize',5)
plot(WA, m_kite./10^3,'b','linewidth',1)
plot(WA, m.w./10^3,'k:','linewidth',1)
legend('MegAWES','AP2','AP3','AP4','AP5 low','AP5 high','M600','MX2','M5','Haas et al. 2019','Kite mass model','Aircraft wing scaling');
xlabel('Wing area (m^2)')
ylabel('Kite mass (tons)')
hold off


%% Kite mass contour plot

% Constants
a1 = 0.002415;
a2 = 0.0090239;
b1 = 0.17025;
b2 = 3.2493;
c1 = 0.46608;
d1 = 0.65962;
k2 = 1.1935;
AR_ref = 12;

% Range values
AR_range = linspace(5, 20, 50);
WA_range = linspace(5, 200, 50);
Ft_max_range = linspace(5, 1000, 50);

% Create meshgrid for AR, WA, and Ft_max
[AR_plot, WA_plot, Ft_max_plot] = ndgrid(AR_range, WA_range, Ft_max_range);

% Calculate r for each combination of Ft_max and WA
r_plot = Ft_max_plot ./ WA_plot;

% Calculate a and b for each r
a_plot = a1 * r_plot + a2;
b_plot = b1 * r_plot + b2;

% Calculate m_kite for each combination of WA, AR, and Ft_max
m_kite_plot = 10 * (a_plot .* WA_plot.^2 + b_plot .* WA_plot) .* ...
  (c1 * (AR_plot / AR_ref).^2 - d1 * (AR_plot / AR_ref) + k2);

% Convert m_kite to tons
m_kite_tons = m_kite_plot / 1000;  % Convert to tons

% Reshape variables for scatter plot
AR_vec = AR_plot(:);
WA_vec = WA_plot(:);
Ft_max_vec = Ft_max_plot(:);
m_kite_tons_vec = m_kite_tons(:);

% Plotting 3D scatter plot
figure;
scatter3( WA_vec, AR_vec, Ft_max_vec, 30, m_kite_tons_vec, 'filled');
xlabel('S [m^2]');
ylabel('F_{t,max} [kN]');
zlabel('AR [-]');
colorbar;
colormap(parula);  % Set colormap for better visualization
% Add legend to colorbar
c = colorbar;
c.Label.String = 'm_{k} [tons]';

%% Static take-off limits


AP2     = [3, 35, 1.5]; % Wing area, kite mass, CL_max
AP3     = [12, 475, 2.1]; % Wing area, kite mass, CL_max
MegAWES = [150.5, 6885, 1.9];
V3      = [19.75, 22.8, 1];
MX2     = [54, 1850, 2];
PN14    = [180, 180, 1];

AP2(4)     = calcSTOL(AP2);
AP3(4)     = calcSTOL(AP3);
MegAWES(4) = calcSTOL(MegAWES);
V3(4)      = calcSTOL(V3);
MX2(4)     = calcSTOL(MX2);
PN14(4)    = calcSTOL(PN14);

%% Comparison with Skysails

fig = openfig('skysails07.fig');
hold on;
box on;
plot(vw, processedOutputs.P_e_avg./10^3,':s','linewidth',2,'MarkerSize',5);
ylim([0 160]);
xlim([0 20]);
title('');
legend('Bin mean', 'bin edges', 'Standard deviation','150 kW system (this study)');
xlabel('Wind speeds at reference height of 100 m (ms^{-1})');
ylabel('Electrical cycle power (kW)');
xticks(1:20);
xticklabels(1:20)
yticks(10:10:160);
yticklabels(10:10:160)
hold off;

function [STOL] = calcSTOL(Kite)

  STOL = sqrt(2*Kite(2)*9.814/1.2/Kite(3)/Kite(1));

end

