%% AWEC2024: Performance model parametric study
clc
clearvars

% Add folders and subfolders to path
addpath(genpath([pwd '\inputSheets']));
addpath(genpath([pwd '\outputFiles'])); 
addpath(genpath([pwd '\src']));

%% Baseline case
clear global
inputSheet_AWEC2024_ParamStudy;
[inputs, outputs, optimDetails, processedOutputs] = main(inputs);
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
powerPlot = plot(inputs.vw_ref,baseline.P_e_avg/1e6, 'o-', 'linewidth', 1, 'MarkerSize', 2, 'Color', pastelBlue);
ylabel('Power (MW)', 'Color', pastelBlue);
ylim([0 2.1])
set(gca, 'YColor', pastelBlue);  % Set y-axis color to pastel blue
% Plot Probability on right y-axis with pastel orange color
yyaxis right
probPlot = plot(pdf*100, 'o-', 'linewidth', 1, 'MarkerSize', 2, 'Color', pastelOrange);
ylabel('Probability (%)', 'Color', pastelOrange);
set(gca, 'YColor', pastelOrange);  % Set y-axis color to pastel orange
xlabel('Wind speed at 100m height (m/s)', 'Color', 'black'); % Change x-axis label color to black
legend('Power', 'Wind pdf', 'Location', 'Best');
% Increase font sizes
set(gca,'FontSize',20); % Increase font size of axes labels
set(get(gca,'XLabel'),'FontSize',20); % Increase font size of x-axis label
set(get(gca,'YLabel'),'FontSize',20); % Increase font size of y-axis label
set(get(gca,'Title'),'FontSize',20); % Increase font size of title
legendFontSize = 18; % Adjust the legend font size as needed
set(legend,'FontSize',legendFontSize); % Increase font size of legend
hold off

%% Wing area sensitivity, fixed Ft_max/S
clear global
inputSheet_AWEC2024_ParamStudy;
inputs.peakM2E_F = 2.5;
Ft_maxByS = 3;
numEvals = 4;
legendLabels = cell(numEvals,1);
% Computation loop
for i = 1:numEvals
    inputs.S = (i-1)*50 + 50;
    S(i)     = inputs.S;
    inputs.Ft_max = Ft_maxByS*inputs.S;
    [inputs, outputs, optimDetails, processedOutputs] = main(inputs);
    wingArea_sens(i) =  processedOutputs;
    wingArea_sens_AEP(i) = calcAEP(9,2,inputs.vw_ref,wingArea_sens(i).P_e_avg);
    legendLabels{i} = ['S = ' num2str(S(i)) 'm^2'];
    clear global
end
% Normalize AEP by max value
wingArea_sens_AEP = wingArea_sens_AEP/baseline_AEP*100;

% Plotting
Title = strcat('F_{t,max}/S = ',num2str(Ft_maxByS),' kN/m^2');
XLABELsubPlot2 = 'Wing area (m^2)';
XLIMsubPlot1 = processedOutputs.vw_100m_operRange(end);
plotSubplots(numEvals, inputs.vw_ref, wingArea_sens, S, Title, XLIMsubPlot1, wingArea_sens_AEP, legendLabels, XLABELsubPlot2);

%% Tether diameter variation
clear global
inputSheet_AWEC2024_ParamStudy;
numEvals = 4;
legendLabels = cell(numEvals,1);
for i = 1:numEvals
    inputs.Ft_max = 200 + i*100;
    Ft_max(i)     = inputs.Ft_max;
    [inputs, outputs, optimDetails, processedOutputs] = main(inputs);
    Ft_max_sens(i)  =  processedOutputs;
    Ft_max_sens_AEP(i) = calcAEP(9,2,inputs.vw_ref,Ft_max_sens(i).P_e_avg);
    legendLabels{i} = ['d_{t}= ' num2str(round(Ft_max_sens(i).Dia_te*100,1)) 'cm'];
    clear global
end
% Normalize AEP by max value
Ft_max_sens_AEP = Ft_max_sens_AEP/baseline_AEP*100;

% Plotting
Title = strcat('S = ',num2str(inputs.S),'m^2');
XLABELsubPlot2 = 'Tether diameter (cm)';
XLIMsubPlot1 = processedOutputs.vw_100m_operRange(end);
plotSubplots(numEvals, inputs.vw_ref, Ft_max_sens, Ft_max, Title, XLIMsubPlot1, Ft_max_sens_AEP, legendLabels, XLABELsubPlot2);



%% CL,max sensitivity
clear global
inputSheet_AWEC2024_ParamStudy;
numEvals = 4;
legendLabels = cell(numEvals,1);
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
    [inputs, outputs, optimDetails, processedOutputs] = main(inputs);
    CL_sens(i) =  processedOutputs;
    CL_sens_AEP(i) = calcAEP(9,2,inputs.vw_ref,CL_sens(i).P_e_avg);
    legendLabels{i} = ['C_{L,max} = ' num2str(CL(i))];
    clear global
end
% Normalize AEP by max value
CL_sens_AEP = CL_sens_AEP/baseline_AEP*100;

% Plotting
Title = strcat('S = ',num2str(inputs.S),'m^2');
XLABELsubPlot2 = 'C_{L,max}';
XLIMsubPlot1 = processedOutputs.vw_100m_operRange(end);
plotSubplots(numEvals, inputs.vw_ref, CL_sens, CL, Title, XLIMsubPlot1, CL_sens_AEP, legendLabels, XLABELsubPlot2);

%% Generator size sensitivity
clear global
inputSheet_AWEC2024_ParamStudy;
numEvals = 4;
legendLabels = cell(numEvals,1);
% Computation loop
for i = 1:numEvals
    inputs.peakM2E_F = 0.5 + i*0.5;
    Gen(i) = inputs.peakM2E_F*inputs.P_ratedElec/1e6; %[MW]
    [inputs, outputs, optimDetails, processedOutputs] = main(inputs);
    Gen_sens(i) =  processedOutputs;
    Gen_sens_AEP(i) = calcAEP(9,2,inputs.vw_ref,Gen_sens(i).P_e_avg);
    legendLabels{i} = ['Gen size = ' num2str(Gen(i)) 'MW'];
    clear global
end
% Normalize AEP by max value
Gen_sens_AEP = Gen_sens_AEP/baseline_AEP*100;

% Plotting
Title = strcat('S = ',num2str(inputs.S),'m^2');
XLIMsubPlot1 = processedOutputs.vw_100m_operRange(end);
XLABELsubPlot2 = 'Gen size (MW)';
plotSubplots(numEvals, inputs.vw_ref, Gen_sens, Gen, Title, XLIMsubPlot1, Gen_sens_AEP, legendLabels, XLABELsubPlot2);

%% Wind shear sensitivity
clear global
inputSheet_AWEC2024_ParamStudy;
inputs.vw_ref         = 2:1:25; %[m/s]
inputs.windShearExp   = 0; %[-] % 0.143 over land, 0.11 over sea
numEvals = 4;
legendLabels = cell(numEvals,1);
for i = 1:numEvals
    inputs.windShearExp   = 0 + (i-1)*0.05; %[-] % 0.143 over land, 0.11 over sea
    windShear(i) = inputs.windShearExp;
    [inputs, outputs, optimDetails, processedOutputs] = main(inputs);
    windShear_sens(i)     =  processedOutputs;
    windShear_sens_AEP(i) = calcAEP(9,2,inputs.vw_ref,windShear_sens(i).P_e_avg);
    legendLabels{i} = ['α = ' num2str(windShear(i))];
    clear global
end
% Normalize AEP by max value
windShear_sens_AEP =windShear_sens_AEP/baseline_AEP*100;

% Plotting
Title = strcat('S = ',num2str(inputs.S),'m^2');
XLIMsubPlot1 = processedOutputs.vw_100m_operRange(end);
XLABELsubPlot2 = 'α (-)';
plotSubplots(numEvals, inputs.vw_ref,windShear_sens, windShear, Title, XLIMsubPlot1, windShear_sens_AEP, legendLabels, XLABELsubPlot2);

% %% Mass sensitivity
% clear global
% inputSheet_AWEC2024_ParamStudy;
% inputs.massOverride   = 1;
% numEvals = 4;
% legendLabels = cell(numEvals,1);
% for i = 1:numEvals
%     inputs.kiteMass   = 6000 + (i-1)*1000; %[kg]
%     mk(i)     = inputs.kiteMass/1e3;
%     [inputs, outputs, optimDetails, processedOutputs] = main(inputs);
%     mk_sens(i)     =  processedOutputs;
%     mk_sens_AEP(i) = calcAEP(9,2,inputs.vw_ref,mk_sens(i).P_e_avg);
%     legendLabels{i} = ['Kite mass = ' num2str(mk(i)) ' tons'];
%     clear global
% end
% % Normalize AEP by max value
% mk_sens_AEP = mk_sens_AEP/baseline_AEP*100;
% 
% % Plotting
% Title = strcat('S = ',num2str(inputs.S),'m^2');
% XLIMsubPlot1 = processedOutputs.vw_100m_operRange(end);
% XLABELsubPlot2 = 'Kite mass (tons)';
% plotSubplots(numEvals, inputs.vw_ref, mk_sens, mk, Title, XLIMsubPlot1, mk_sens_AEP, legendLabels, XLABELsubPlot2);
% 
% %% Tether drag coefficient sensitivity
% clear global
% inputSheet_AWEC2024_ParamStudy;
% inputs.Cd_c = 0.8;
% numEvals = 4;
% legendLabels = cell(numEvals,1);
% for i = 1:numEvals
%     inputs.Cd_c = 0.8 + (i-1)*0.1;
%     Cd_t(i)     = inputs.Cd_c;
%     [inputs, outputs, optimDetails, processedOutputs] = main(inputs);
%     Cd_t_sens(i)  =  processedOutputs;
%     Cd_t_sens_AEP(i) = calcAEP(9,2,inputs.vw_ref,Cd_t_sens(i).P_e_avg);
%     legendLabels{i} = ['C_{d,t}= ' num2str(Cd_t(i))];
%     clear global
% end
% % Normalize AEP by max value
% Cd_t_sens_AEP = Cd_t_sens_AEP/baseline_AEP*100;
% 
% % Plotting
% Title = strcat('S = ',num2str(inputs.S),'m^2');
% XLABELsubPlot2 = 'Cross-section drag coeff.';
% XLIMsubPlot1 = processedOutputs.vw_100m_operRange(end);
% plotSubplots(numEvals, inputs.vw_ref, Cd_t_sens, Cd_t, Title, XLIMsubPlot1, Cd_t_sens_AEP, legendLabels, XLABELsubPlot2);

%% Functions

function [AEP, pdf] = calcAEP(meanWS, shapeParam, vw_ref, power)

  % Weibull
  
  % Calculate the scale parameter based on mean wind speed
  scaleParam = meanWS/ gamma(1 + 1/shapeParam);
  
  % Calculate Weibull PDF values
  pdf = (shapeParam / scaleParam) * ((vw_ref / scaleParam).^(shapeParam - 1)) .* exp(-(vw_ref / scaleParam).^shapeParam);
  
  % AWES AEP
  AEP = 8760*trapz(vw_ref, pdf .* power)/1e6; %[MWh]

end

function plotSubplots(numEvals, vw_ref, param_sens, param, Title, XLIMsubPlot1, param_sens_AEP, legendLabels, XLABELsubplot2)
    % Define the total plot width
    totalWidth = 13; % inches
    % Define the widths for each subplot
    subplot1Width = 0.5 * totalWidth; % 55% of total width
    subplot2Width = 0.22 * totalWidth; % 20% of total width

    % Define the figure size
    figSize = [1 1 totalWidth 8]; % Use totalWidth for the figure width

    % Create the figure with the specified size
    figure('units','inch','Position', figSize)

    % Subplot 1
    subplot('Position', [0.1 0.16 subplot1Width/totalWidth 0.7]) % Adjust the position to take 55% width
    hold on
    grid on
    box on
    for i = 1:numEvals
        plot(vw_ref, param_sens(i).P_e_avg./1e6, 'o-', 'linewidth', 1.5, 'MarkerSize', 2);
    end
    title(Title);
    xlabel('Wind speed at 100m height (m/s)');
    ylabel('Power (MW)');
    xlim([0 XLIMsubPlot1]);
    ylim([0 2.1]);
    legend(legendLabels, 'Location', 'Best');
    % Increase font sizes
    set(gca,'FontSize',26); % Increase font size of axes labels
    set(get(gca,'XLabel'),'FontSize',34); % Increase font size of x-axis label
    set(get(gca,'YLabel'),'FontSize',34); % Increase font size of y-axis label
    set(get(gca,'Title'),'FontSize',34); % Increase font size of title
    legendFontSize = 24; % Adjust the legend font size as needed
    set(legend,'FontSize',legendFontSize); % Increase font size of legend
    hold off

    % Subplot 2
    subplot('Position', [0.72 0.16 subplot2Width/totalWidth 0.7]) % Adjust the position to take 20% width
    hold on
    grid on
    box on
    for i = 1:numEvals
        scatter(param(i),param_sens_AEP(i), 250, "filled");
    end
    % Fit a linear curve
    p = polyfit(param, param_sens_AEP, 2);
    fitline = polyval(p, param);
    plot(param, fitline, 'k--', 'LineWidth', 2);  % Plot the fit line in red dashed line
    title('Normalized AEP');
    xlim([0.9*min(param) 1.1*max(param)])
    ylim([20 140])
    xticks(param)
    xticklabels(cellstr(num2str(param(:))))
    xlabel(XLABELsubplot2)
    ylabel('(%)')
    % Increase font sizes
    set(gca,'FontSize',22); % Increase font size of axes labels
    set(get(gca,'XLabel'),'FontSize',34); % Increase font size of x-axis label
    set(get(gca,'YLabel'),'FontSize',34); % Increase font size of y-axis label
    set(get(gca,'Title'),'FontSize',34); % Increase font size of title
    hold off
end

