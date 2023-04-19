%% Run this script to get results of the power computation model
clc
clearvars
clear global

% Defined input sheet
inputSheet_AP3;

% Get results
[optData,outputs,postProRes,timeseries] = main(inputs);

%% Main Plots
if inputs.mainPlots == 1
  Vw = inputs.Vw_ref; 

  % Wind profile
  Vref = 10; % m/s
  z = 10:10:600; % m
  V = Vref * (z/inputs.h_ref).^inputs.windShearExp/Vref;
  figure('units','inch','Position', [15 6 3.5 2.2])
  hold on
  box on
  grid on
  plot(V,z,'linewidth',1.5)
  xlim([0.5 1.5])
  xlabel('Wind Speed (-)')
  ylabel('Height (m)')

  newcolors = [ % 0.25, 0.25, 0.25
    0 0.4470 0.7410
  0.8500 0.3250 0.0980 
  0.4660, 0.6740, 0.1880
  0.9290, 0.6940, 0.1250
  0.4940, 0.1840, 0.5560
  0.6350 0.0780 0.1840
  0.3010 0.7450 0.9330];

  % Cycle timeseries plots: Pattern averages
  windSpeeds = [16];
  for i = windSpeeds
    tmax = round(max(postProRes.tCycle(windSpeeds)));
    pmax = 1.2*inputs.F_peakM2Ecyc*max(postProRes.Pcycle_elec(windSpeeds))/10^3;
  %   pmax = 1.5*max(postProRes.Pcycle_elec(windSpeeds))/10^3;
    pmin = min(-postProRes.PRIeff_elec(windSpeeds))/10^3;

    figure('units','inch','Position', [15 3 3.5 2.2])
    hold on
    grid on
    box on
    yline(0);
    yline(postProRes.Pcycle_elec(timeseries.ws(i).ws)/10^3,'--','linewidth',1);
    plot(timeseries.ws(i).t_inst, timeseries.ws(i).P_e_inst,'linewidth',1.5);
    plot(timeseries.ws(i).t_inst, timeseries.ws(i).P_m_inst,'linewidth',1.5);
    ylabel('Power (kW)');
    xlabel('Time (s)');
    %xlim([0 tmax]);
    %ylim([pmin pmax]);
    title(strcat('Wind speed at 100m:',num2str(timeseries.ws(i).ws),'m/s'));
  %   legend('Electrical','Mechanical','location','northwest');
    hold off
  end

  % Mechanical reel-out power oscillation and capping
  i = [8, 13,25];
  d.series1 = outputs.PROeff_mech_osci(i(1),:)/10^3;  
  d.series2 = outputs.PROeff_mech_osci(i(2),:)/10^3; 
  d.series3 = outputs.PROeff_mech_osci(i(3),:)/10^3; 
  d.deltaLelems = 0:1:outputs.deltaLelems-1;
  figure('units','inch','Position', [4 4 3.5 2.2])
  hold on
  grid on
  box on
  plot(d.deltaLelems,d.series1,'linewidth',1.2);
  plot(d.deltaLelems,d.series2,'linewidth',1.2);
  plot(d.deltaLelems,d.series3,'linewidth',1.2);
  ylabel('P_{m,o} (kW)');
  legend(strcat(num2str(i(1)),'m/s'),strcat(num2str(i(2)),'m/s'),strcat(num2str(i(3)),'m/s'),'location','southeast');
  xlabel('Discretized reel-out length in number of elements');
  hold off

  % VRO, VRI, tRO, tRI
  figure('units','inch','Position', [3 3 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  yyaxis left
  plot(Vw, mean(postProRes.VRO,2),'d:','markersize',3);
  plot(Vw, postProRes.VRI,'^:','markersize',3);
  ylabel('Speed (m/s)');
  yyaxis right
  plot(Vw, postProRes.tRO,'o:','markersize',3);
  plot(Vw, postProRes.tRI,'+:','markersize',3);
  ylabel('Time (s)');
  legend('v_{o,avg}','v_{i}','t_{o}','t_{i}','location','northwest');
  xlabel('Wind speed at 100m height (m/s)');
  xlim([0 25]);
  %ylim([0 160]);
  hold off

  % H_cycleAvg, Patt radius, Delta L, Avg tether length
  figure('units','inch','Position', [3 6 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  plot(Vw, postProRes.H_cycleAvg,'d:','markersize',3);
  plot(Vw, mean(postProRes.pattRad,2),'o:','markersize',3);
  plot(Vw, postProRes.deltaL,'^:','markersize',3);
  plot(Vw, postProRes.L_teMax,'+:','markersize',3);
  ylabel('Length (m)');
  legend('H_{cycleAvg}','R_{cycleAvg}','ΔL','L_{t,max}','location','northwest','Orientation','vertical');
  xlabel('Wind speed at 100m height (m/s)');
  xlim([0 25]);
  %ylim([0 160]);
  hold off

  % Roll angle, avg patt elevation
  figure('units','inch','Position', [7 3 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  plot(Vw, postProRes.avgRollAngle,'x:','markersize',4);
  plot(Vw, postProRes.avgPattEle,'^:','markersize',3);
  plot(Vw, postProRes.pattAngRadius,'o:','markersize',3);
  ylabel('Angle (deg)');
  legend('Avg. roll angle','Avg. patt ele','Patt. angular radius','location','northwest','Orientation','vertical');
  xlabel('Wind speed at 100m height (m/s)');
  xlim([0 25]);
  %ylim([0 160]);
  hold off

  % Dimension less nos: CL,CD,reel-out factor
  figure('units','inch','Position', [7 6 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  yyaxis left
  plot(Vw, postProRes.CL,'^:','markersize',3);
  plot(Vw, postProRes.CD,'+:','markersize',3);
  plot(Vw, mean(postProRes.reelOutF,2),'o:','markersize',3);
  ylabel('(-)');
  yyaxis right
  plot(Vw, mean(postProRes.TT,2)./10^3,'x:','markersize',3);
  ylabel('Tether force (kN)');
  ylim([0 1.1*max(mean(postProRes.TT,2))/10^3]);
  legend('C_{L}','C_{D}','f_{avg}','F_{t,avg}','location','northwest');
  xlabel('Wind speed at 100m height (m/s)');
  xlim([0 25]);
  %ylim([0 160]);
  hold off

  % Cycle efficiencies
  ERO_elec = sum(postProRes.PROeff_elec.*postProRes.tROeff,2)' + postProRes.PRO1_elec.*postProRes.t1;
  ERI_elec = sum(postProRes.PRIeff_elec.*postProRes.tRIeff,2)' + postProRes.PRI2_elec.*postProRes.t2;
  ERO_mech = sum(postProRes.PROeff_mech.*postProRes.tROeff,2)' + postProRes.PRO1_mech.*postProRes.t1;
  ERI_mech = sum(postProRes.PRIeff_mech.*postProRes.tRIeff,2)' + postProRes.PRI2_mech.*postProRes.t2;
  figure('units','inch','Position', [11 3 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  plot(Vw, (ERO_mech-ERI_mech)./ERO_mech*100,'o','markersize',3);
  plot(Vw, (ERO_elec-ERI_elec)./ERO_mech*100,'+','markersize',3);
  plot(Vw, postProRes.dutyCycle*100,'d','markersize',3);
  ylabel('Efficiency (%)');
  ylim([30 100]);
  legend('η_{m,avg}','η_{e,avg}','Duty cycle','location','southeast');
  xlabel('Wind speed at 100m height (m/s)');
  xlim([0 25]);
  hold off

  % Power plots
  %pastel_b = '[0 0.4470 0.7410]';
  figure('units','inch','Position', [11 6 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  plot(Vw, postProRes.PRO_mech./10^3,'-o','linewidth',1,'markersize',2.5);
  plot(Vw, postProRes.PRO_elec./10^3,':o','linewidth',1,'markersize',2.5);
  plot(Vw, postProRes.PRI_mech./10^3,'--','linewidth',1.2);
  plot(Vw, postProRes.PRI_elec./10^3,':','linewidth',1.5);
  plot(Vw, postProRes.Pcycle_elec./10^3,'-','linewidth',1.4);
  ylabel('Power (kW)');
  %title('Cycle averages');
  legend('P_{m,o,avg}','P_{e,o,avg}','P_{m,i}','P_{e,i}','P_{e,avg}','location','northwest');
  xlabel('Wind speed at 100m height (m/s)');
  xlim([0 25]);
  hold off


  % Power curve comparison plot
  % Loyd
  CL_loyd = inputs.CL_maxAirfoil*inputs.F_CLeff;
  CD_loyd = inputs.CD0 + (CL_loyd-inputs.CL0_airfoil)^2/(pi()*inputs.AR*inputs.e);
  for i = 1:length(Vw)
      P_Loyd(i) = (4/27)*(CL_loyd^3/CD_loyd^2)*(1/2)*inputs.airDensity*inputs.WA*(Vw(i)^3);
      if P_Loyd(i)>postProRes.ratedPower
          P_Loyd(i) = postProRes.ratedPower;
      end 
  end
  % AP3 6DoF simulation results
  AP3.PC.ws    = [7.32E+00, 8.02E+00, 9.05E+00, 1.00E+01, 1.10E+01, 1.20E+01, ...
    1.30E+01, 1.41E+01, 1.50E+01, 1.60E+01, 1.70E+01, 1.80E+01, 1.90E+01]; %[m/s]
  AP3.PC.power = [0.00E+00, 7.01E+03, 2.37E+04, 4.31E+04, 6.47E+04, 8.46E+04, ...
    1.02E+05, 1.20E+05, 1.34E+05, 1.49E+05, 1.50E+05, 1.50E+05, 1.50E+05]./10^3; %[kW]
  figure('units','inch','Position', [4 4 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  plot(Vw, P_Loyd./10^3,'--','linewidth',1.2);
  plot(Vw, postProRes.Pcycle_elec./10^3,'o-','linewidth',1,'MarkerSize',4);
  plot(AP3.PC.ws, AP3.PC.power,'k^--','MarkerSize',4);
  ylabel('Power (kW)');
  legend('Loyd','Model','6DOF','location','southeast');
  % legend('Loyd - Ideal','Model results','location','southeast');
  xlabel('Wind speed at 100m height (m/s)');
  xlim([0 25]);
  hold off
  
end



%% Sensitivity plots

%% Wing area
inputSheet_AP3;
inputs.massOverride   = 1;
inputs.kiteMass = 4.365695525193600e+02;
s(1).WA  = inputs.WA;
s(2).WA  = 0.7*s(1).WA;
s(3).WA  = 0.85*s(1).WA;
s(4).WA  = 1.15*s(1).WA;
s(5).WA  = 1.3*s(1).WA;
for i = 1:numel(s)
  inputs.WA = s(i).WA;
  [optData,outputs,postProRes,timeseries] = main(inputs);
  P(i).WA =  postProRes.Pcycle_elec;
end
figSize = [5 5 3.2 2];
figure('units','inch','Position', figSize)
hold on
grid on
box on
plot(P(2).WA./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(3).WA./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(1).WA./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(4).WA./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(5).WA./1e3,'o-','linewidth',1,'MarkerSize',2);
legend(strcat(num2str(s(2).WA),'m^2'),strcat(num2str(s(3).WA),'m^2'),strcat(num2str(s(1).WA),'m^2'),strcat(num2str(s(4).WA),'m^2'),strcat(num2str(s(5).WA),'m^2'),'location','southeast');
xlabel('Wind speed at pattern altitude (m/s)');
ylabel('Power (kW)');
xlim([3 20]);
hold off


%% Kite mass
inputSheet_AP3;
inputs.massOverride   = 1;
s(1).kiteMass  = 4.365695525193600e+02;
s(2).kiteMass  = 0.7*s(1).kiteMass;
s(3).kiteMass  = 0.85*s(1).kiteMass;
s(4).kiteMass  = 1.15*s(1).kiteMass;
s(5).kiteMass  = 1.3*s(1).kiteMass;
for i = 1:numel(s)
  inputs.kiteMass = s(i).kiteMass;
  [optData,outputs,postProRes,timeseries] = main(inputs);
  P(i).kiteMass =  postProRes.Pcycle_elec;
end
figSize = [5 5 3.2 2];
figure('units','inch','Position', figSize)
hold on
grid on
box on
plot(P(2).kiteMass./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(3).kiteMass./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(1).kiteMass./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(4).kiteMass./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(5).kiteMass./1e3,'o-','linewidth',1,'MarkerSize',2);
legend(strcat(num2str(s(2).kiteMass),'kg'),strcat(num2str(s(3).kiteMass),'kg'),strcat(num2str(s(1).kiteMass),'kg'),strcat(num2str(s(4).kiteMass),'kg'),strcat(num2str(s(5).kiteMass),'kg'),'location','southeast');
xlabel('Wind speed at pattern altitude (m/s)');
ylabel('Power (kW)');
xlim([3 20]);
hold off
inputs.massOverride   = 0;


%% Kite CL 
inputSheet_AP3;
s(1).CL  = inputs.CL_maxAirfoil;
s(2).CL  = 0.7*s(1).CL;
s(3).CL  = 0.85*s(1).CL;
s(4).CL  = 1.15*s(1).CL;
s(5).CL  = 1.3*s(1).CL;
for i = 1:numel(s)
  inputs.CL_maxAirfoil = s(i).CL;
  [optData,outputs,postProRes,timeseries] = main(inputs);
  P(i).CL =  postProRes.Pcycle_elec;
end
figSize = [5 5 3.2 2];
figure('units','inch','Position', figSize)
hold on
grid on
box on
plot(P(2).CL./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(3).CL./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(1).CL./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(4).CL./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(5).CL./1e3,'o-','linewidth',1,'MarkerSize',2);
legend(strcat(num2str(s(2).CL),''),strcat(num2str(s(3).CL),''),strcat(num2str(s(1).CL),''),strcat(num2str(s(4).CL),''),strcat(num2str(s(5).CL),''),'location','southeast');
xlabel('Wind speed at pattern altitude (m/s)');
ylabel('Power (kW)');
xlim([3 20]);
hold off

%% Max. Tether force
inputSheet_AP3;
inputs.massOverride   = 1;
inputs.kiteMass = 4.365695525193600e+02;
s(1).Tmax  = inputs.Tmax;
s(2).Tmax  = 0.7*s(1).Tmax;
s(3).Tmax  = 0.85*s(1).Tmax;
s(4).Tmax  = 1.15*s(1).Tmax;
s(5).Tmax  = 1.3*s(1).Tmax;
for i = 1:numel(s)
  inputs.Tmax  = s(i).Tmax;
  [optData,outputs,postProRes,timeseries] = main(inputs);
  P(i).Tmax =  postProRes.Pcycle_elec;
end
figSize = [5 5 3.2 2];
figure('units','inch','Position', figSize)
hold on
grid on
box on
plot(P(2).Tmax./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(3).Tmax./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(1).Tmax./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(4).Tmax./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(5).Tmax./1e3,'o-','linewidth',1,'MarkerSize',2);
legend(strcat(num2str(s(2).Tmax),'kN'),strcat(num2str(s(3).Tmax),'kN'),strcat(num2str(s(1).Tmax),'kN'),strcat(num2str(s(4).Tmax),'kN'),strcat(num2str(s(5).Tmax),'kN'),'location','southeast');
xlabel('Wind speed at pattern altitude (m/s)');
ylabel('Power (kW)');
xlim([3 20]);
hold off
inputs.massOverride   = 0;

%% Tether strength
inputSheet_AP3;
s(1).TeSig  = inputs.Te_matStrength;
s(2).TeSig  = 0.7*s(1).TeSig;
s(3).TeSig  = 0.85*s(1).TeSig;
s(4).TeSig  = 1.15*s(1).TeSig;
s(5).TeSig  = 1.3*s(1).TeSig;
for i = 1:numel(s)
  inputs.Te_matStrength  = s(i).TeSig;
  [optData,outputs,postProRes,timeseries] = main(inputs);
  P(i).TeSig =  postProRes.Pcycle_elec;
end
figSize = [5 5 3.2 2];
figure('units','inch','Position', figSize)
hold on
grid on
box on
plot(P(2).TeSig./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(3).TeSig./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(1).TeSig./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(4).TeSig./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(5).TeSig./1e3,'o-','linewidth',1,'MarkerSize',2);
legend(strcat(num2str(s(2).TeSig),'kg/m^3'),strcat(num2str(s(3).TeSig),'kg/m^3'),strcat(num2str(s(1).TeSig),'kg/m^3'),strcat(num2str(s(4).TeSig),'kg/m^3'),strcat(num2str(s(5).TeSig),'kg/m^3'),'location','southeast');
xlabel('Wind speed at pattern altitude (m/s)');
ylabel('Power (kW)');
xlim([3 20]);
hold off


%% Tether drag cofficient
inputSheet_AP3;
s(1).TeCd  = inputs.CD_te;
s(2).TeCd  = 0.7*s(1).TeCd;
s(3).TeCd  = 0.85*s(1).TeCd;
s(4).TeCd  = 1.15*s(1).TeCd;
s(5).TeCd  = 1.3*s(1).TeCd;
for i = 1:numel(s)
  inputs.CD_te  = s(i).TeCd;
  [optData,outputs,postProRes,timeseries] = main(inputs);
  P(i).TeCd =  postProRes.Pcycle_elec;
end
figSize = [5 5 3.2 2];
figure('units','inch','Position', figSize)
hold on
grid on
box on
plot(P(2).TeCd./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(3).TeCd./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(1).TeCd./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(4).TeCd./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(5).TeCd./1e3,'o-','linewidth',1,'MarkerSize',2);
legend(strcat(num2str(s(2).TeCd),''),strcat(num2str(s(3).TeCd),''),strcat(num2str(s(1).TeCd),''),strcat(num2str(s(4).TeCd),''),strcat(num2str(s(5).TeCd),''),'location','southeast');
xlabel('Wind speed at pattern altitude (m/s)');
ylabel('Power (kW)');
xlim([3 20]);
hold off

%% Generator size
inputSheet_AP3;
s(1).PbyA  = inputs.F_peakM2Ecyc;
s(2).PbyA  = 0.7*s(1).PbyA;
s(3).PbyA  = 0.85*s(1).PbyA;
s(4).PbyA  = 1.15*s(1).PbyA;
s(5).PbyA  = 1.3*s(1).PbyA;
for i = 1:numel(s)
  inputs.F_peakM2Ecyc  = s(i).PbyA;
  [optData,outputs,postProRes,timeseries] = main(inputs);
  P(i).PbyA =  postProRes.Pcycle_elec;
end
figSize = [5 5 3.2 2];
figure('units','inch','Position', figSize)
hold on
grid on
box on
plot(P(2).PbyA./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(3).PbyA./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(1).PbyA./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(4).PbyA./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(5).PbyA./1e3,'o-','linewidth',1,'MarkerSize',2);
legend(strcat(num2str(s(2).PbyA),''),strcat(num2str(s(3).PbyA),''),strcat(num2str(s(1).PbyA),''),strcat(num2str(s(4).PbyA),''),strcat(num2str(s(5).PbyA),''),'location','southeast');
xlabel('Wind speed at pattern altitude (m/s)');
ylabel('Power (kW)');
xlim([3 20]);
hold off


%% Maximum winch speed
inputSheet_AP3;
s(1).VRI  = inputs.maxVRI;
s(2).VRI  = 0.7*s(1).VRI;
s(3).VRI  = 0.85*s(1).VRI;
s(4).VRI  = 1.15*s(1).VRI;
s(5).VRI  = 1.3*s(1).VRI;
for i = 1:numel(s)
  inputs.maxVRI  = s(i).VRI;
  [optData,outputs,postProRes,timeseries] = main(inputs);
  P(i).VRI =  postProRes.Pcycle_elec;
end
figSize = [5 5 3.2 2];
figure('units','inch','Position', figSize)
hold on
grid on
box on
plot(P(2).VRI./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(3).VRI./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(1).VRI./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(4).VRI./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(P(5).VRI./1e3,'o-','linewidth',1,'MarkerSize',2);
legend(strcat(num2str(s(2).VRI),'m/s'),strcat(num2str(s(3).VRI),'m/s'),strcat(num2str(s(1).VRI),'m/s'),strcat(num2str(s(4).VRI),'m/s'),strcat(num2str(s(5).VRI),'m/s'),'location','southeast');
xlabel('Wind speed at pattern altitude (m/s)');
ylabel('Power (kW)');
xlim([3 20]);
hold off
 
%% Extra
% Sine wave plot for VRO oscillation
% t = linspace(270*pi/180,(270+360)*pi/180,1000);  
% x = sin(t); 
% figure() 
% hold on
% plot(t,x,'linewidth',2);
% axis off
% box on

% Oscillation theories comparison
% figure()
% hold on
% box on
% grid on
% plot(outputs.VRO_osci(7,:))
% plot(outputs.VRO_osci2(7,:))
% ylim([-3 20])
% 
% figure()
% hold on
% box on
% grid on
% plot(outputs.VRO_osci(16,:))
% plot(outputs.VRO_osci2(16,:))
% ylim([-3 20])
% 
% figure()
% hold on
% box on
% grid on
% plot(outputs.VRO_osci(25,:))
% plot(outputs.VRO_osci2(25,:))
% ylim([-3 20])
