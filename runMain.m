%% Run this script to get results of the power computation model
clc
clearvars
clear global

% Defined input sheet
inputSheet_AP3;
% inputSheet;

% Get results
[optData,outputs,processedOutputs] = main(inputs);

%% Main Plots
if inputs.mainPlots == 1
  
  vw = inputs.vw_ref;
  x_axis_limit = [0 processedOutputs.vw_100m_operRange(end)];

  newcolors = [ % 0.25, 0.25, 0.25
    0 0.4470 0.7410
  0.8500 0.3250 0.0980 
  0.4660, 0.6740, 0.1880
  0.9290, 0.6940, 0.1250
  0.4940, 0.1840, 0.5560
  0.6350 0.0780 0.1840
  0.3010 0.7450 0.9330];

%%
  % Wind profile
  Vref = 10; % m/s
  z = 10:10:600; % m
  V = Vref * (z/inputs.h_ref).^inputs.windShearExp/Vref;
%   MegAWES Onshore location Cabauw. Wind speeds normalized with value at 100m
  z_MegAWES = [10,20,40,60,80,100,120,140,150,160,180,200,220,250,300,500,600];
  V_MegAWES_Cabauw = [0.541219682843206,0.607355091566827,0.768630154201962,0.868484406441142,0.941395360902529,1,1.04810058627160,1.08638854381156,1.10277338731106,1.11715868927737,1.14412258234309,1.16573551308321,1.18394938534465,1.20653423381438,1.23266397972046,1.26662287360302,1.26414483994687];
  V_MegAWES_Ijmuiden = [0.847612611633547,0.870603040595613,0.927240267828556,0.959346286990695,0.982291573490674,1,1.01377720773809,1.02356771954493,1.02766760602000,1.03079423355205,1.03659625208888,1.04025827758100,1.04284618416620,1.04496440015282,1.04461712713371,1.02473617783789,1.01076976884552];
  figure('units','inch','Position', [15 6 2 2.2])
  hold on
  box on
  grid on
  plot(V,z,'linewidth',1.5)
  plot(V_MegAWES_Cabauw,z_MegAWES,'linewidth',1.5)
  plot(V_MegAWES_Ijmuiden,z_MegAWES,'linewidth',1.5)
  legend('α = 0.11','Cabauw,NL','Ijmuiden,NL');
  xlim([0.5 1.5])
  xlabel('Wind Speed (-)')
  ylabel('Height (m)')

  % Wind speed at 100m vs avg wind speed during reel-out
  figure('units','inch','Position', [17 6 2 2.2])
  hold on
  grid on
  box on
  plot(vw,mean(processedOutputs.vw,2),'linewidth',1);
  xlim(x_axis_limit);
  xlabel('v_w at 100m height (m/s)');
  ylabel('v_{w,avg} during operation (m/s)')

  
  % Cycle timeseries plots: Pattern averages
  windSpeeds = [processedOutputs.ratedWind];
  for i = windSpeeds
%     tmax = round(max(processedOutputs.tCycle(windSpeeds)));
%     pmax = 1.1*inputs.peakM2E_F*max(processedOutputs.Pcycle_elec(windSpeeds))/10^3;
%     pmin = 1.1*min(-processedOutputs.PRIeff_elec(windSpeeds))/10^3;
    idx = find(vw==i);
    figure('units','inch','Position', [15 3 3.5 2.2])
    hold on
    grid on
    box on
    yline(0);
    yline(processedOutputs.P_e_avg(idx)/10^3,'--','linewidth',1);
    plot(processedOutputs.cyclePowerRep(i).t_inst, processedOutputs.cyclePowerRep(i).P_e_inst,'linewidth',1.5);
    plot(processedOutputs.cyclePowerRep(i).t_inst,processedOutputs.cyclePowerRep(i).P_m_inst,'linewidth',1.5);
    ylabel('Power (kW)');
    xlabel('Time (s)');
%     xlim([0 tmax]);
%     ylim([pmin pmax]);
    title(strcat('Wind speed at 100m:',num2str(processedOutputs.cyclePowerRep(idx).ws),'m/s'));
  %   legend('Electrical','Mechanical','location','northwest');
    hold off
  end

  % Plots showing the per element results for rated wind speed


  % Reel-out, Reel-in Speeds and Time
  figure('units','inch','Position', [3 3 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  yyaxis left
  plot(vw, mean(processedOutputs.vo,2),'d:','markersize',3);
  plot(vw, mean(processedOutputs.vi,2),'^:','markersize',3);
  ylabel('Speed (m/s)');
  yyaxis right
  plot(vw, processedOutputs.to,'o:','markersize',3);
  plot(vw, processedOutputs.ti,'+:','markersize',3);
  ylabel('Time (s)');
  legend('v_{o,avg}','v_{i,avg}','t_{o}','t_{i}','location','northwest');
  xlabel('Wind speed at 100m height (m/s)');
  xlim(x_axis_limit);
  %ylim([0 160]);
  hold off

  % Lengths
  figure('units','inch','Position', [3 6 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  plot(vw, processedOutputs.h_cycleAvg,'d:','markersize',3);
  plot(vw, mean(processedOutputs.Rp,2),'o:','markersize',3);
  plot(vw, processedOutputs.deltaL,'^:','markersize',3);
  plot(vw, processedOutputs.l_t_max,'+:','markersize',3);
  ylabel('Length (m)');
  legend('h_{p,avg}','R_{p,avg}','Δl','l_{t,max}','location','northwest','Orientation','vertical');
  xlabel('Wind speed at 100m height (m/s)');
  xlim(x_axis_limit);
  %ylim([0 160]);
  hold off

  % Roll angle, avg patt elevation
  figure('units','inch','Position', [7 3 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  plot(vw, mean(processedOutputs.rollAngle,2),'x:','markersize',4);
  plot(vw, processedOutputs.beta,'^:','markersize',3);
  plot(vw, processedOutputs.gamma,'o:','markersize',3);
  ylabel('Angle (deg)');
  legend('Ψ_{avg}','β','γ','location','northwest','Orientation','vertical');
  xlabel('Wind speed at 100m height (m/s)');
  xlim(x_axis_limit);
  %ylim([0 160]);
  hold off

  % All forces 
  figure('units','inch','Position', [7 6 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  plot(vw, mean(processedOutputs.Fa,2)./10^3,'-o','linewidth',1,'markersize',2.5);
  plot(vw, mean(processedOutputs.Ft,2)./10^3,':o','linewidth',1,'markersize',2.5);
  plot(vw, mean(processedOutputs.Ft_drum,2)./10^3,'-.o','linewidth',1,'markersize',2.5);
  % plot(vw, mean(processedOutputs.W,2)./10^3,'--x','linewidth',1,'markersize',3.5);
  plot(vw, processedOutputs.W/10^3,':x','linewidth',1,'markersize',3.5);
  ylabel('Force (kN)');
  legend('F_{a}','F_{t}','F_{t,drum}','W','location','northwest');
  xlabel('Wind speed at 100m height (m/s)');
  xlim(x_axis_limit);
  %ylim([0 160]);
  hold off

  % % Non-dim factors
  % ERO_elec = sum(processedOutputs.PROeff_elec.*processedOutputs.tROeff,2)' + processedOutputs.PRO1_elec.*processedOutputs.t1;
  % ERI_elec = sum(processedOutputs.PRIeff_elec.*processedOutputs.tRIeff,2)' + processedOutputs.PRI2_elec.*processedOutputs.t2;
  % ERO_mech = sum(processedOutputs.PROeff_mech.*processedOutputs.tROeff,2)' + processedOutputs.PRO1_mech.*processedOutputs.t1;
  % ERI_mech = sum(processedOutputs.PRIeff_mech.*processedOutputs.tRIeff,2)' + processedOutputs.PRI2_mech.*processedOutputs.t2;
  % figure('units','inch','Position', [11 3 3.5 2.2])
  % colororder(newcolors)
  % hold on
  % grid on
  % box on
  % plot(vw, (ERO_mech-ERI_mech)./ERO_mech,'o','markersize',3);
  % plot(vw, (ERO_elec-ERI_elec)./ERO_mech,'+','markersize',3);
  % plot(vw, processedOutputs.dutyCycle,'d','markersize',3);
  % plot(vw, mean(processedOutputs.reelOutF,2),'x','markersize',3);
  % ylabel('(-)');
  % % ylim([30 100]);
  % legend('η_{m,avg}','η_{e,avg}','D','f','location','southeast');
  % xlabel('Wind speed at 100m height (m/s)');
  % xlim(x_axis_limit);
  % hold off

  % Power plots
  %pastel_b = '[0 0.4470 0.7410]';
  figure('units','inch','Position', [11 6 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  plot(vw, processedOutputs.P_m_o./10^3,'-o','linewidth',1,'markersize',2.5);
  plot(vw, processedOutputs.P_e_o./10^3,':o','linewidth',1,'markersize',2.5);
  plot(vw, processedOutputs.P_m_i./10^3,'--','linewidth',1.2);
  plot(vw, processedOutputs.P_e_i./10^3,':','linewidth',1.5);
  plot(vw, processedOutputs.P_e_avg./10^3,'-','linewidth',1.4);
  ylabel('Power (kW)');
  %title('Cycle averages');
  legend('P_{m,o,mean}','P_{e,o,mean}','P_{m,i,mean}','P_{e,i,mean}','P_{e,avg}','location','northwest');
  xlabel('Wind speed at 100m height (m/s)');
  xlim(x_axis_limit);
  hold off


  % Power curve comparison plot
  % Loyd
  CL_loyd = inputs.CL_maxAirfoil*inputs.CLeff_F;
  CD_loyd = inputs.CD0 + (CL_loyd-inputs.CL0_airfoil)^2/(pi()*inputs.AR*inputs.e);
  for i = 1:length(vw)
      P_Loyd(i) = (4/27)*(CL_loyd^3/CD_loyd^2)*(1/2)*inputs.airDensity*inputs.S*(vw(i)^3);
      if P_Loyd(i)>processedOutputs.ratedPower
          P_Loyd(i) = processedOutputs.ratedPower;
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
  plot(vw, P_Loyd./10^3,'--','linewidth',1.2);
  plot(vw, processedOutputs.P_e_avg./10^3,'o-','linewidth',1,'MarkerSize',4);
  plot(AP3.PC.ws, AP3.PC.power,'k^--','MarkerSize',4);
  ylabel('Power (kW)');
  legend('No effects (Loyd Reel-out)','All effects','6DOF','location','southeast');
  % legend('Loyd - Ideal','Model results','location','southeast');
  xlabel('Wind speed at 100m height (m/s)');
  xlim(x_axis_limit);
  hold off
  
end



%% Sensitivity plots

if inputs.runSensitivity == 1
  %% Wing area
  inputSheet_AP3;
  inputs.massOverride   = 1;
  inputs.kiteMass = 4.365695525193600e+02;
  s(1).S  = inputs.S;
  s(2).S  = 0.7*s(1).S;
  s(3).S  = 0.85*s(1).S;
  s(4).S  = 1.15*s(1).S;
  s(5).S  = 1.3*s(1).S;
  for i = 1:numel(s)
    inputs.S = s(i).S;
    [optData,outputs,postProRes] = main(inputs);
    P(i).S =  postProRes.Pcycle_elec;
  end
  figSize = [5 5 3.2 2];
  figure('units','inch','Position', figSize)
  hold on
  grid on
  box on
  plot(P(2).S./1e3,'o-','linewidth',1,'MarkerSize',2);
  plot(P(3).S./1e3,'o-','linewidth',1,'MarkerSize',2);
  plot(P(1).S./1e3,'o-','linewidth',1,'MarkerSize',2);
  plot(P(4).S./1e3,'o-','linewidth',1,'MarkerSize',2);
  plot(P(5).S./1e3,'o-','linewidth',1,'MarkerSize',2);
  legend(strcat(num2str(s(2).S),'m^2'),strcat(num2str(s(3).S),'m^2'),strcat(num2str(s(1).S),'m^2'),strcat(num2str(s(4).S),'m^2'),strcat(num2str(s(5).S),'m^2'),'location','southeast');
  xlabel('Wind speed at 100m height (m/s)');
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
    [optData,outputs,postProRes] = main(inputs);
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
  xlabel('Wind speed at 100m height (m/s)');
  ylabel('Power (kW)');
  xlim([3 20]);
  hold off
  
  
  %% Kite CL 
  inputSheet_AP3;
  s(1).CL  = inputs.CL_maxAirfoil;
  s(2).CL  = 0.7*s(1).CL;
  s(3).CL  = 0.85*s(1).CL;
  s(4).CL  = 1.15*s(1).CL;
  s(5).CL  = 1.3*s(1).CL;
  for i = 1:numel(s)
    inputs.CL_maxAirfoil = s(i).CL;
    [optData,outputs,postProRes] = main(inputs);
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
  xlabel('Wind speed at 100m height (m/s)');
  ylabel('Power (kW)');
  xlim([3 20]);
  hold off
  
  %% Max. Tether force
  inputSheet_AP3;
  inputs.massOverride   = 1;
  inputs.kiteMass = 4.365695525193600e+02;
  s(1).Tmax  = inputs.Ft_max;
  s(2).Tmax  = 0.7*s(1).Tmax;
  s(3).Tmax  = 0.85*s(1).Tmax;
  s(4).Tmax  = 1.15*s(1).Tmax;
  s(5).Tmax  = 1.3*s(1).Tmax;
  for i = 1:numel(s)
    inputs.Ft_max  = s(i).Tmax;
    [optData,outputs,postProRes] = main(inputs);
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
  xlabel('Wind speed at 100m height (m/s)');
  ylabel('Power (kW)');
  xlim([3 20]);
  hold off
  
  %% Tether strength
  inputSheet_AP3;
  s(1).TeSig  = inputs.Te_matStrength;
  s(2).TeSig  = 0.7*s(1).TeSig;
  s(3).TeSig  = 0.85*s(1).TeSig;
  s(4).TeSig  = 1.15*s(1).TeSig;
  s(5).TeSig  = 1.3*s(1).TeSig;
  for i = 1:numel(s)
    inputs.Te_matStrength  = s(i).TeSig;
    [optData,outputs,postProRes] = main(inputs);
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
  xlabel('Wind speed at 100m height (m/s)');
  ylabel('Power (kW)');
  xlim([3 20]);
  hold off
  
  
  %% Tether drag cofficient
  inputSheet_AP3;
  s(1).TeCd  = inputs.CD_t;
  s(2).TeCd  = 0.7*s(1).TeCd;
  s(3).TeCd  = 0.85*s(1).TeCd;
  s(4).TeCd  = 1.15*s(1).TeCd;
  s(5).TeCd  = 1.3*s(1).TeCd;
  for i = 1:numel(s)
    inputs.CD_t  = s(i).TeCd;
    [optData,outputs,postProRes] = main(inputs);
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
  xlabel('Wind speed at 100m height (m/s)');
  ylabel('Power (kW)');
  xlim([3 20]);
  hold off
  
  %% Generator size
  inputSheet_AP3;
  s(1).PbyA  = inputs.peakM2E_F;
  s(2).PbyA  = 0.7*s(1).PbyA;
  s(3).PbyA  = 0.85*s(1).PbyA;
  s(4).PbyA  = 1.15*s(1).PbyA;
  s(5).PbyA  = 1.3*s(1).PbyA;
  for i = 1:numel(s)
    inputs.peakM2E_F  = s(i).PbyA;
    [optData,outputs,postProRes] = main(inputs);
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
  xlabel('Wind speed at 100m height (m/s)');
  ylabel('Power (kW)');
  xlim([3 20]);
  hold off
  
  
  %% Maximum winch speed
  inputSheet_AP3;
  s(1).VRI  = inputs.vk_r_i_max;
  s(2).VRI  = 0.5*s(1).VRI;
  s(3).VRI  = 0.75*s(1).VRI;
  s(4).VRI  = 1.25*s(1).VRI;
  s(5).VRI  = 1.5*s(1).VRI;
  for i = 1:numel(s)
    inputs.vk_r_i_max  = s(i).VRI;
    [optData,outputs,postProRes] = main(inputs);
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
  xlabel('Wind speed at 100m height (m/s)');
  ylabel('Power (kW)');
  xlim([3 20]);
  hold off
end

%%

% lambda.NoEle    = [5.94, 5.94, 5.94, 5.94];
% lambda.NoEle_CF = [6.08, 6.08, 6.08, 6.08];
% % lambda.Ele   = [5.79,4.89, 6.18, 6.44, 5.20];
% % lambda.G     = [5.71,4.76, 6.99, 6.39, 4.06];
% % lambda.CF    = [6.68, 6.26, 4.83, 5.39];
% % lambda.All   = [6.63, 7.06, 4.67, 4.31];
% 
% zeta.NoEle    = [15.66, 15.66, 15.66, 15.66];
% zeta.NoEle_CF = [16.58, 16.59, 16.59, 16.59];
% % zeta.Ele   = [14.88,10.75, 14.30, 18.36, 14.30];
% % zeta.G     = [14.32,10, 18.4, 17.98, 9.07];
% % zeta.CF    = [19.59, 14.31, 10.4,15.31];
% % zeta.All    = [19.23, 18.39, 9.5,10.16];
% 
% figure()
% hold on
% box on
% grid on
% % yyaxis left
% plot(lambda.NoEle,'d-','markersize',5,'linewidth',1);
% plot(lambda.NoEle_CF,'^-','markersize',5,'linewidth',1);
% % plot(lambda.G,'o-','markersize',5,'linewidth',1);
% % plot(lambda.CF,'s-.','markersize',8,'linewidth',1);
% % plot(lambda.All,'<-.','markersize',8,'linewidth',1);
% legend('No Elevation', 'No Elevation + CF'); 
% ylabel('λ(-)');
% xticks([1,2,3,4]);
% xticklabels({'Top (Flying left)','Side (Flying down)','Bottom (Flying right)','Side(Flying up)'});
% title('S = 12 m^2, mk = 436 kg, vw = 20m/s, Δl = 250 m, vi = 20 m/s, β = 30 deg, γ = 10 deg, Rp_strt = 50 m, vk_r = 4m/s');
% hold off
% 
% figure()
% hold on
% box on
% grid on
% % yyaxis right
% plot(zeta.NoEle,'d-','markersize',5,'linewidth',1);
% plot(zeta.NoEle_CF,'^-','markersize',5,'linewidth',1);
% % plot(zeta.G,'o-','markersize',5,'linewidth',1);
% % plot(zeta.CF,'s-.','markersize',8,'linewidth',1);
% % plot(zeta.All,'<-.','markersize',8,'linewidth',1);
% legend('No Elevation', 'No Elevation +CF'); 
% ylabel('ζ(-)');
% xticks([1,2,3,4]);
% xticklabels({'Top (Flying left)','Side (Flying down)','Bottom (Flying right)','Side(Flying up)'});
% title('S = 12 m^2, mk = 436 kg, vw = 20m/s, Δl = 250 m, vi = 20 m/s, β = 30 deg, γ = 10 deg, Rp_strt = 50 m, vk_r = 4m/s');
% hold off
