%% AWEC2024: Performance model parametric study
clc
clearvars

%% Baseline case
clear global
inputSheet_AWEC2024_ParamStudy;
[inputs, outputs, optimDetails, processedOutputs] = main(inputs, 'inputSheet_AWEC2024_ParamStudy');
baseline     =  processedOutputs;
[baseline_AEP, pdf] = calcAEP(9,2,inputs.vw_ref,baseline.P_e_avg);
% Plotting
% Define pastel colors
pastelBlue = [0 0.4470 0.7410];  % RGB values for pastel blue
pastelOrange = [0.8500 0.3250 0.0980];  % RGB values for pastel orange
% Plotting
figSize = [5 5 5 4]; % Increased width for the first subplot
figure('units','inch','Position', figSize)
hold on
grid on
box on
% Plot Power on left y-axis with pastel blue color
yyaxis left
powerPlot = plot(baseline.P_e_avg/1e6, 'o-', 'linewidth', 1, 'MarkerSize', 2, 'Color', pastelBlue);
ylabel('Power (MW)', 'Color', pastelBlue);
set(gca, 'YColor', pastelBlue);  % Set y-axis color to pastel blue
% Plot Probability on right y-axis with pastel orange color
yyaxis right
probPlot = plot(pdf*100, 'o-', 'linewidth', 1, 'MarkerSize', 2, 'Color', pastelOrange);
ylabel('Probability (%)', 'Color', pastelOrange);
set(gca, 'YColor', pastelOrange);  % Set y-axis color to pastel orange
xlabel('Wind speed at 100m height (m/s)', 'Color', 'black'); % Change x-axis label color to black
legend('Power', 'Wind pdf', 'Location', 'Best');
% Increase font sizes
set(gca,'FontSize',24); % Increase font size of axes labels
set(get(gca,'XLabel'),'FontSize',24); % Increase font size of x-axis label
set(get(gca,'YLabel'),'FontSize',24); % Increase font size of y-axis label
set(get(gca,'Title'),'FontSize',24); % Increase font size of title
legendFontSize = 16; % Adjust the legend font size as needed
set(legend,'FontSize',legendFontSize); % Increase font size of legend
hold off

%% Wing area sensitivity, fixed Ft_max/S
% clc
% clearvars
clear global
inputSheet_AWEC2024_ParamStudy;
inputs.peakM2E_F = 2.5;
Ft_maxByS = 3;
numEvals = 4;
% Computation loop
for i = 1:numEvals
    inputs.S = (i-1)*50 + 50;
    S(i)     = inputs.S;
    inputs.Ft_max = Ft_maxByS*inputs.S;
    [inputs, outputs, optimDetails, processedOutputs] = main(inputs, 'inputSheet_AWEC2024_ParamStudy;');
    wingArea_sens(i) =  processedOutputs;
    wingArea_sens_AEP(i) = calcAEP(9,2,inputs.vw_ref,wingArea_sens(i).P_e_avg);
    clear global
end
% Normalize AEP by max value
wingArea_sens_AEP = wingArea_sens_AEP/baseline_AEP*100;

% Plotting
figSize = [5 5 10 4]; % Increased width for the first subplot
figure('units','inch','Position', figSize)
% Plot first subplot with increased width
subplot(1,2,1)
hold on
grid on
box on
legendLabels = cell(numEvals,1);
for i = 1:numEvals
    plot(wingArea_sens(i).P_e_avg./1e6,'o-','linewidth',1,'MarkerSize',2);
    legendLabels{i} = ['S = ' num2str(S(i)) 'm^2'];
end
title(strcat('F_{t,max}/S = ',num2str(Ft_maxByS),' kN/m^2'));
subplot1Formatting(processedOutputs.vw_100m_operRange(end), 2.5, legendLabels);
hold off
% Plot second subplot
subplot(1,2,2)
hold on
grid on
box on
for i = 1:numEvals
  scatter(S(i),wingArea_sens_AEP(i),100, "filled");
end
subplot2Formatting(S,wingArea_sens_AEP,'Wing area (m^2)','(%)');
hold off


%% Tether diameter variation
% clc
% clearvars
clear global
inputSheet_AWEC2024_ParamStudy;
numEvals = 4;
for i = 1:numEvals
    inputs.Ft_max = 200 + i*100;
    Ft_max(i)     = inputs.Ft_max;
    [inputs, outputs, optimDetails, processedOutputs] = main(inputs, 'inputSheet_AWEC2024_ParamStudy');
    Ft_max_sens(i)  =  processedOutputs;
    Ft_max_sens_AEP(i) = calcAEP(9,2,inputs.vw_ref,Ft_max_sens(i).P_e_avg);
    clear global
end
% Normalize AEP by max value
Ft_max_sens_AEP = Ft_max_sens_AEP/baseline_AEP*100;

% Plotting
figSize = [5 5 10 4]; % Increased width for the first subplot
figure('units','inch','Position', figSize)
% Plot first subplot with increased width
subplot(1,2,1)
hold on
grid on
box on
legendLabels = cell(numEvals,1);
for i = 1:numEvals
    plot(Ft_max_sens(i).P_e_avg./1e6,'o-','linewidth',1,'MarkerSize',2);
    legendLabels{i} = ['d_{t}= ' num2str(round(Ft_max_sens(i).Dia_te*100,1)) 'cm'];
end
title(strcat('S = ',num2str(inputs.S),'m^2, ','C_{d,t} = ', num2str(1.2)));
subplot1Formatting(processedOutputs.vw_100m_operRange(end), 2.5, legendLabels);
hold off
% Plot second subplot
subplot(1,2,2)
hold on
grid on
box on
for i = 1:numEvals
  d_t(i) = round(Ft_max_sens(i).Dia_te*100,1);
end
for i = 1:numEvals
  scatter(d_t(i),Ft_max_sens_AEP(i),100,"filled");
end
subplot2Formatting(d_t,Ft_max_sens_AEP,'Tether diameter (cm)','(%)');
hold off

%% Tether drag coefficient sensitivity
% clc
% clearvars
clear global
inputSheet_AWEC2024_ParamStudy;
inputs.Cd_c = 0.8;
numEvals = 4;
for i = 1:numEvals
    inputs.Cd_c = 0.8 + (i-1)*0.1;
    Cd_t(i)     = inputs.Cd_c;
    [inputs, outputs, optimDetails, processedOutputs] = main(inputs, 'inputSheet_AWEC2024_ParamStudy');
    Cd_t_sens(i)  =  processedOutputs;
    Cd_t_sens_AEP(i) = calcAEP(9,2,inputs.vw_ref,Cd_t_sens(i).P_e_avg);
    clear global
end
% Normalize AEP by max value
Cd_t_sens_AEP = Cd_t_sens_AEP/baseline_AEP*100;

% Plotting
figSize = [5 5 10 4]; % Increased width for the first subplot
figure('units','inch','Position', figSize)
% Plot first subplot with increased width
subplot(1,2,1)
hold on
grid on
box on
legendLabels = cell(numEvals,1);
for i = 1:numEvals
    plot(Cd_t_sens(i).P_e_avg./1e6,'o-','linewidth',1,'MarkerSize',2);
    legendLabels{i} = ['C_{d,t}= ' num2str(Cd_t(i))];
end
title(strcat('S = ',num2str(inputs.S),'m^2'));
subplot1Formatting(processedOutputs.vw_100m_operRange(end), 2.5, legendLabels);
hold off
% Plot second subplot
subplot(1,2,2)
hold on
grid on
box on
for i = 1:numEvals
  scatter(Cd_t(i),Cd_t_sens_AEP(i),100, "filled");
end
subplot2Formatting(Cd_t,Cd_t_sens_AEP,'Cross-section drag coeff.','(%)');
hold off




%% CL,max sensitivity
clear global
inputSheet_AWEC2024_ParamStudy;
numEvals = 4;
% Computation loop
for i = 1:numEvals
    inputs.Cl_maxAirfoil = 2 + (i-1)*0.3;
    CL(i) = inputs.Cl_eff_F*inputs.Cl_maxAirfoil;
    % Optimisation problem data
    inputs.nx = ones(1,inputs.numDeltaLelems);
    % Initial guess
    %               [deltaL, avgPattEle,  coneAngle,     Rp_start,         v_i,                      CL_i,                                           v_o,           kinematicRatio, CL]
    inputs.x0     = [200,    deg2rad(30), deg2rad(5),   3*inputs.b,       inputs.v_d_max*inputs.nx, inputs.Cl_maxAirfoil*inputs.Cl_eff_F*inputs.nx, 0.8*inputs.nx, 90*inputs.nx,   inputs.Cl_maxAirfoil*inputs.Cl_eff_F*inputs.nx];
    % Bounds
    inputs.lb     = [50,   deg2rad(1),  deg2rad(1),  3*inputs.b,  1*inputs.nx, 0.1*inputs.nx, 0.8*inputs.nx, 1*inputs.nx, 0.1*inputs.nx]; % 
    inputs.ub     = [500,  deg2rad(90), deg2rad(60), 500, inputs.v_d_max*inputs.nx, inputs.Cl_maxAirfoil*inputs.Cl_eff_F*inputs.nx,...
      inputs.v_d_max*inputs.nx,  200*inputs.nx, inputs.Cl_maxAirfoil*inputs.Cl_eff_F*inputs.nx]; %
    [inputs, outputs, optimDetails, processedOutputs] = main(inputs, 'inputSheet_AWEC2024_ParamStudy;');
    CL_sens(i) =  processedOutputs;
    CL_sens_AEP(i) = calcAEP(9,2,inputs.vw_ref,CL_sens(i).P_e_avg);
    clear global
end
% Normalize AEP by max value
CL_sens_AEP = CL_sens_AEP/baseline_AEP*100;

% Plotting
figSize = [5 5 10 4]; % Increased width for the first subplot
figure('units','inch','Position', figSize)
% Plot first subplot with increased width
subplot(1,2,1)
hold on
grid on
box on
legendLabels = cell(numEvals,1);
for i = 1:numEvals
    plot(CL_sens(i).P_e_avg./1e6,'o-','linewidth',1,'MarkerSize',2);
    legendLabels{i} = ['C_{L,max} = ' num2str(CL(i))];
end
title(strcat('S = ',num2str(inputs.S),'m^2'));
subplot1Formatting(processedOutputs.vw_100m_operRange(end), 2.5, legendLabels);
hold off
% Plot second subplot
subplot(1,2,2)
hold on
grid on
box on
for i = 1:numEvals
  scatter(CL(i),CL_sens_AEP(i),100, "filled");
end
subplot2Formatting(CL,CL_sens_AEP,'C_{L,max}','(%)');
hold off

%% Generator size sensitivity
clear global
inputSheet_AWEC2024_ParamStudy;
numEvals = 4;
% Computation loop
for i = 1:numEvals
    inputs.peakM2E_F = 0.5 + i*0.5;
    Gen(i) = inputs.peakM2E_F*inputs.P_ratedElec/1e6; %[MW]
    [inputs, outputs, optimDetails, processedOutputs] = main(inputs, 'inputSheet_AWEC2024_ParamStudy;');
    Gen_sens(i) =  processedOutputs;
    Gen_sens_AEP(i) = calcAEP(9,2,inputs.vw_ref,Gen_sens(i).P_e_avg);
    clear global
end
% Normalize AEP by max value
Gen_sens_AEP = Gen_sens_AEP/baseline_AEP*100;

% Plotting
figSize = [5 5 10 4]; % Increased width for the first subplot
figure('units','inch','Position', figSize)
% Plot first subplot with increased width
subplot(1,2,1)
hold on
grid on
box on
legendLabels = cell(numEvals,1);
for i = 1:numEvals
    plot(Gen_sens(i).P_e_avg./1e6,'o-','linewidth',1,'MarkerSize',2);
    legendLabels{i} = ['Gen size = ' num2str(Gen(i)) 'MW'];
end
title(strcat('S = ',num2str(inputs.S),'m^2'));
subplot1Formatting(processedOutputs.vw_100m_operRange(end), 2.5, legendLabels);
hold off
% Plot second subplot
subplot(1,2,2)
hold on
grid on
box on
for i = 1:numEvals
  scatter(Gen(i),Gen_sens_AEP(i),100, "filled");
end
subplot2Formatting(Gen,Gen_sens_AEP,'Gen size (MW)','(%)');
hold off


%% Mass sensitivity
% % clc
% % clearvars
% clear global
% inputSheet_AWEC2024_ParamStudy;
% inputs.massOverride   = 1;
% numEvals = 4;
% for i = 1:numEvals
%     inputs.kiteMass   = 6000 + (i-1)*1000; %[kg]
%     mk(i)     = inputs.kiteMass;
%     [inputs, outputs, optimDetails, processedOutputs] = main(inputs, 'inputSheet_AWEC2024_ParamStudy');
%     mk_sens(i)     =  processedOutputs;
%     mk_sens_AEP(i) = calcAEP(9,2,inputs.vw_ref,mk_sens(i).P_e_avg);
%     clear global
% end
% % Normalize AEP by max value
% mk_sens_AEP = mk_sens_AEP/baseline_AEP*100;
% 
% % Plotting
% figSize = [5 5 10 4]; % Increased width for the first subplot
% figure('units','inch','Position', figSize)
% % Plot first subplot with increased width
% subplot(1,2,1)
% hold on
% grid on
% box on
% legendLabels = cell(numEvals,1);
% for i = 1:numEvals
%     plot(mk_sens(i).P_e_avg./1e6,'o-','linewidth',1,'MarkerSize',2);
%     legendLabels{i} = ['Kite mass= ' num2str(mk(i)) 'kg'];
% end
% title(strcat('S = ',num2str(inputs.S),'m^2'));
% subplot1Formatting(processedOutputs.vw_100m_operRange(end), 2.5, legendLabels);
% hold off
% % Plot second subplot
% subplot(1,2,2)
% hold on
% grid on
% box on
% for i = 1:numEvals
%   scatter(mk(i),mk_sens_AEP(i),100, "filled");
% end
% subplot2Formatting(mk,mk_sens_AEP,'Kite mass (kg)','(%)');
% hold off


%% Mass scaling curve plots

% a1     = 0.002415;
% a2     = 0.0090239;
% b1     = 0.17025;
% b2     = 3.2493;
% k1     = 5;
% c1     = 0.46608;
% d1     = 0.65962;
% k2     = 1.1935;
% AR_ref = 12;
% 
% r  = 3.5; % Tmax/inputs.WA %[kN/m^2]
% AR = 12;
% 
% a = a1*r + a2;
% b = b1*r + b2;
% 
% WA = linspace(0, 310,10);
% m_kite = 10*(a*WA.^2 +b*WA)*(c1*(AR/AR_ref)^2-d1*(AR/AR_ref)+k2); 
% 
% 
% w.CL = 2;
% par.u0_max = 100;
% q = 1/2 * 1.225 * par.u0_max^2;
% W =  r*WA*1000;% q * (WA.*w.CL);% + h.A* h.CL); 
% ff = 0.8;
% tc = 0.2;
% nz = 6;
% 
% w.S = 0;
% w.AR = 12;
% 
% m.w = ff * 0.036.* WA.^0.758 * (w.AR/cos(w.S)^2)^0.6 * q^0.006 * (100*tc/cos(w.S))^(-0.3) .* (nz * W).^0.49;
% 
% 
% % Data from literature
% MegAWES = [150.45, 6885];
% AP2     = [3, 35];
% AP3     = [12, 475];
% AP4     = [35^2/12, 3500];
% AP5_Luuk    = [300, 20000];
% AP5_AP    = [300, 34000];
% M600    = [32.9, 1730.8];
% MX2     = [54, 1850];
% Haas    = [138.5, 8000];
% M5      = [177.3, 9900];
% 
% newcolors = [ % 0.25, 0.25, 0.25
%   0 0.4470 0.7410
% 0.8500 0.3250 0.0980 
% 0.4660, 0.6740, 0.1880
% 0.4940, 0.1840, 0.5560
% 0.9290, 0.6940, 0.1250
% 
% 0.6350 0.0780 0.1840
% 0.3010 0.7450 0.9330];
% 
% 
% figSize = [5 5 6 4];
% figure('units','inch','Position', figSize)
% hold on
% grid on
% box on
% colororder(newcolors)
% plot(MegAWES(1), MegAWES(2)/10^3,'p','markersize',7) % Set marker size to 7
% plot(AP2(1), AP2(2)/10^3,'s','markersize',7) % Set marker size to 7
% plot(AP3(1), AP3(2)/10^3,'d','markersize',7) % Set marker size to 7
% plot(AP4(1), AP4(2)/10^3,'*','markersize',7) % Set marker size to 7
% plot(AP5_Luuk(1), AP5_Luuk(2)/10^3,'o','markersize',7) % Set marker size to 7
% plot(AP5_AP(1), AP5_AP(2)/10^3,'x','markersize',7) % Set marker size to 7
% plot(M600(1), M600(2)/10^3,'^','markersize',7) % Set marker size to 7
% plot(MX2(1), MX2(2)/10^3,'v','markersize',7) % Set marker size to 7
% plot(M5(1), M5(2)/10^3,'<','markersize',7) % Set marker size to 7
% plot(Haas(1), Haas(2)/10^3,'>','markersize',7) % Set marker size to 7
% plot(WA, m_kite./10^3,'k','linewidth',1)
% plot(WA, m.w./10^3,'k--','linewidth',1)
% legend('MegAWES','AP2','AP3','AP4','AP5 low','AP5 high','M600','MX2','M5','Haas et al. 2019','AP scaling','Aircraft wing scaling','NumColumns',2); 
% xlabel('Wing area (m^2)')
% ylabel('Kite mass (tons)')
% % % Increase font sizes
% set(gca,'FontSize',14); % Increase font size of axes labels
% set(get(gca,'XLabel'),'FontSize',16); % Increase font size of x-axis label
% set(get(gca,'YLabel'),'FontSize',16); % Increase font size of y-axis label
% set(get(gca,'Title'),'FontSize',16); % Increase font size of title
% legendFontSize = 12; % Adjust the legend font size as needed
% set(legend,'FontSize',legendFontSize); % Increase font size of legend
% hold off

function [AEP, pdf] = calcAEP(meanWS, shapeParam, vw_ref, power)

  % Weibull
  
  % Calculate the scale parameter based on mean wind speed
  scaleParam = meanWS/ gamma(1 + 1/shapeParam);
  
  % Calculate Weibull PDF values
  pdf = (shapeParam / scaleParam) * ((vw_ref / scaleParam).^(shapeParam - 1)) .* exp(-(vw_ref / scaleParam).^shapeParam);
  
  % AWES AEP
  AEP = 8760*trapz(vw_ref, pdf .* power)/1e6; %[MWh]

end

function subplot1Formatting(XLIM, YLIM, legendLabels)
  xlabel('Wind speed at 100m height (m/s)');
  ylabel('Power (MW)');
  xlim([0 XLIM]);
  ylim([0 YLIM]);
  legend(legendLabels, 'Location', 'Best');
  % Increase font sizes
  set(gca,'FontSize',18); % Increase font size of axes labels
  set(get(gca,'XLabel'),'FontSize',24); % Increase font size of x-axis label
  set(get(gca,'YLabel'),'FontSize',24); % Increase font size of y-axis label
  set(get(gca,'Title'),'FontSize',24); % Increase font size of title
  legendFontSize = 16; % Adjust the legend font size as needed
  set(legend,'FontSize',legendFontSize); % Increase font size of legend
end

function subplot2Formatting(param,paramAEP, XLABEL,YLABEL)
  title('Normalized AEP');
  xlim([0.9*min(param) 1.1*max(param)])
  ylim([0.95*min(paramAEP) 100])
  xticks(param)
  xticklabels(cellstr(num2str(param(:))))
  xlabel(XLABEL)
  ylabel(YLABEL)
  % Increase font sizes
  set(gca,'FontSize',18); % Increase font size of axes labels
  set(get(gca,'XLabel'),'FontSize',24); % Increase font size of x-axis label
  set(get(gca,'YLabel'),'FontSize',24); % Increase font size of y-axis label
  set(get(gca,'Title'),'FontSize',24); % Increase font size of title
end
