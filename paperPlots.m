%% Paper plots


%% Effect of gravity
clc 
clearvars
clear global
for i = 1:2
  inputSheet_gravityEffect;
  inputs.FgToggle = i-1;
  [inputs, outputs, optimDetails, processedOutputs] = main(inputs, 'inputSheet_gravityEffect');
  gravityEffect(i)     =  processedOutputs;
  clear global
end

figure('units','inch','Position', [5 5 3.5 2.5])
hold on
grid on
box on
plot(gravityEffect(1).P_e_avg./1e3,'o-','linewidth',1,'MarkerSize',2);
plot(gravityEffect(2).P_e_avg./1e3,'^-','linewidth',1,'MarkerSize',2);
legend('P_{no-gravity}','P_{gravity}');
ylabel('Power (kW)');
xlabel('Wind speed at 100m height (m/s)');
xlim([0 processedOutputs.vw_100m_operRange(end)]);
hold off

figure('units','inch','Position', [5 5 3.5 2.5])
hold on
grid on
box on
plot(gravityEffect(1).dutyCycle.*100,'o-','linewidth',1,'MarkerSize',2);
plot(gravityEffect(2).dutyCycle.*100,'^-','linewidth',1,'MarkerSize',2);
legend('eta_{no-gravity}','eta_{gravity}');
ylabel('(%)');
xlabel('Wind speed at 100m height (m/s)');
xlim([0 processedOutputs.vw_100m_operRange(end)]);
hold off

%% Effect of scaling on fixed wind speed of 15 m/s
% Problem settings
% Zero ground clearance, no constraint on rated, max power, tether length, height, patterns, etc.
% Enforce max. tether force constraint
% Change objective function to max. P_m_o
% Increase max. iterations in the optimisation settings

% Fixed WA: Finding opt tether size for given WA
clc
clearvars
clear global
inputSheet_scalingEffects;
inputs.S = 60;
for i =1:6
  inputs.Ft_max = 200+200*i;
  fixedWA.Ft_max(i) = inputs.Ft_max;
  [inputs, outputs, optimDetails, processedOutputs] = main(inputs, 'inputSheet_scalingEffects');
  fixedWA.m_k(i)         = processedOutputs.m_k;
  fixedWA.zetaMech(i)    = processedOutputs.zetaMech;
  fixedWA.Ft_to_S(i)     = fixedWA.Ft_max(i)/inputs.S;
  clear global
end
figure('units','inch','Position', [4 4 3.5 2.2])
hold on
grid on
box on
xlabel('F_{t,max} (kN)');
yyaxis left
plot(fixedWA.Ft_max,fixedWA.zetaMech,'o-','linewidth',1,'MarkerSize',3);
ylabel('ζ (-)');
yyaxis right
plot(fixedWA.Ft_max,fixedWA.m_k,'^-','linewidth',1,'MarkerSize',3);
ylabel('m_k (kg)');
title('S = 100 m^2');
hold off

% Fixed Ft_max
clc
clearvars
clear global
inputSheet_scalingEffects;
FixedFt_max.Ft_max = 300;
for i =1:8
  inputs.S = 10*i;
  FixedFt_max.S(i)       = inputs.S;
  inputs.Ft_max          = FixedFt_max.Ft_max;
  FixedFt_max.Ft_max     = inputs.Ft_max;
  [inputs, outputs, optimDetails, processedOutputs] = main(inputs, 'inputSheet_scalingEffects');
  FixedFt_max.m_k(i)     = processedOutputs.m_k;
  FixedFt_max.zeta(i)    = processedOutputs.zetaMech;
  clear global
end
figure('units','inch','Position', [12 4 3.5 2.2])
hold on
grid on
box on
xlabel('S (m^2)');
yyaxis left
plot(FixedFt_max.S,FixedFt_max.zeta,'o-','linewidth',1,'MarkerSize',3);
ylabel('ζ (-)');
yyaxis right
plot(FixedFt_max.S,FixedFt_max.m_k,'^-','linewidth',1,'MarkerSize',3);
ylabel('m_k (kg)');
title('F_{t.max} = 200 kN');
hold off


%% Effects of scaling on full operational wind speed range 
% % Full Power curve plots for fixed S, increasing Ft_max
% inputSheet;
% inputs.S = 100;
% figSize = [5 5 3.2 2];
% figure('units','inch','Position', figSize)
% hold on
% grid on
% box on
% numEvals = 5;
% legendLabels = cell(numEvals,1);
% for i = 1:numEvals
%     inputs.Ft_max = 2*inputs.S + i*100;
%     [optData,outputs,processedOutputs] = main(inputs);
%     P(i).increasingFt_max =  processedOutputs.P_e_avg;
%     % Plotting
%     plot(P(i).increasingFt_max./1e6,'o-','linewidth',1,'MarkerSize',2);
%     legendLabels{i} = ['F_{t,max}= ' num2str(inputs.Ft_max) 'kN'];
% end
% xlabel('Wind speed at 100m height (m/s)');
% ylabel('Power (MW)');
% title(strcat('S = ',num2str(inputs.S),'m^2'));
% xlim([0 processedOutputs.vw_100m_operRange(end)]);
% legend(legendLabels, 'Location', 'Best');
% hold off
% 
% % Full Power curve plots for increasing S, fixed Ft_max/S, no mechanical power limit
% inputSheet;
% inputs.peakM2E_F      = 30;
% figSize = [5 5 3.2 2];
% figure('units','inch','Position', figSize)
% hold on
% grid on
% box on
% numEvals = 5;
% legendLabels = cell(numEvals,1);
% for i = 1:numEvals
%     inputs.S = i*30;
%     inputs.Ft_max = 8*inputs.S;
%     [optData,outputs,processedOutputs] = main(inputs);
%     P(i).increasingS =  processedOutputs.P_e_avg;
%     % Plotting
%     plot(P(i).increasingS./1e6,'o-','linewidth',1,'MarkerSize',2);
%     legendLabels{i} = ['S = ' num2str(inputs.S) 'm^2'];
% end
% xlabel('Wind speed at 100m height (m/s)');
% xlim([0 processedOutputs.vw_100m_operRange(end)]);
% ylabel('Power (MW)');
% title(strcat('F_{t,max}/S = ',num2str(inputs.Ft_max/inputs.S),'kN/m^2'));
% legend(legendLabels, 'Location', 'Best');
% hold off

%% Sensitivity plots example

% %% Wing area
% inputSheet_AP3;
% % inputs.massOverride   = 1;
% % inputs.kiteMass = 4.365695525193600e+02;
% s(1).S  = inputs.S;
% s(2).S  = 0.8*s(1).S;
% s(3).S  = 0.9*s(1).S;
% s(4).S  = 1.1*s(1).S;
% s(5).S  = 1.2*s(1).S;
% for i = 1:numel(s)
%   inputs.S = s(i).S;
%   [optData,outputs,postProRes] = main(inputs);
%   P(i).S =  postProRes.P_e_avg;
% end
% figSize = [5 5 3.2 2];
% figure('units','inch','Position', figSize)
% hold on
% grid on
% box on
% plot(P(2).S./1e3,'o-','linewidth',1,'MarkerSize',2);
% plot(P(3).S./1e3,'o-','linewidth',1,'MarkerSize',2);
% plot(P(1).S./1e3,'o-','linewidth',1,'MarkerSize',2);
% plot(P(4).S./1e3,'o-','linewidth',1,'MarkerSize',2);
% plot(P(5).S./1e3,'o-','linewidth',1,'MarkerSize',2);
% legend(strcat(num2str(s(2).S),'m^2'),strcat(num2str(s(3).S),'m^2'),strcat(num2str(s(1).S),'m^2'),strcat(num2str(s(4).S),'m^2'),strcat(num2str(s(5).S),'m^2'),'location','southeast');
% xlabel('Wind speed at 100m height (m/s)');
% ylabel('Power (kW)');
% xlim([3 20]);
% hold off
% 
%% Comparison for case with zero elevation (verification with base-case)

% lambda.NoEle = [5.07, 5.23, 5.36, 5.48, 5.58];
% lambda.Ele   = [4.1,4.26,4.4,4.51,4.62];
% lambda.G     = [3.45,3.83,4.08,4.26,4.41];
% 
% kByE.NoEle   = [1,1,1,1,1];
% kByE.Ele     = [1,1,1,1,1];
% kByE.G       = [0.84,0.89,0.92,0.94,0.95];
% 
% zeta.NoEle    = [15.63,15.24,14.8,14.34,13.87];
% zeta.Ele   = [10.4,10.28,10.09,9.86,9.61];
% zeta.G     = [6.59,7.7,8.21,8.43,8.49];
% ws = 11:15;
% 
% figure()
% subplot(1,2,1)
% hold on
% box on
% grid on
% xlabel('Wind speed (m/s)');
% plot(ws,kByE.NoEle,'-d','markersize',3,'linewidth',1);
% plot(ws,kByE.Ele,'-o','markersize',3,'linewidth',1);
% plot(ws,kByE.G,'-s','markersize',3,'linewidth',1);
% legend('No Elevation', 'Elevation',' Elevation + Gravity'); 
% ylabel('k/E(-)');
% % xticks([1,2,3,4,5]);
% % xticklabels({'Top (Flying left)','Side (Flying down)','Bottom (Flying right)','Side(Flying up)','Center (flying left)'});
% hold off
% subplot(1,2,2)
% % figure()
% hold on
% box on
% grid on
% xlabel('Wind speed (m/s)');
% plot(ws,zeta.NoEle,'-d','markersize',3,'linewidth',1);
% plot(ws,zeta.Ele,'-o','markersize',3,'linewidth',1);
% plot(ws,zeta.G,'-s','markersize',3,'linewidth',1);
% legend('No Elevation', 'Elevation',' Elevation + Gravity'); 
% ylabel('ζ(-)');
% % xticks([1,2,3,4,5]);
% % xticklabels({'Top (Flying left)','Side (Flying down)','Bottom (Flying right)','Side(Flying up)','Center (flying left)'});
% hold off