function plotResults_awePower(inputs)        

  % Load processed outputs
  load(['outputFiles/' inputs.name '_' 'processedOutputs' '.mat']);

  %% Plot settings

  vw           = processedOutputs.vw_100m_operRange;
  cutIn        = processedOutputs.cutIn;
  x_axis_limit = [1 processedOutputs.vw_100m_operRange(end)];
  
  newcolors = [ % 0.25, 0.25, 0.25
    0 0.4470 0.7410
  0.8500 0.3250 0.0980 
  0.4660, 0.6740, 0.1880
  0.9290, 0.6940, 0.1250
  0.4940, 0.1840, 0.5560
  0.6350 0.0780 0.1840
  0.3010 0.7450 0.9330];

  %% Plots
  % Plots showing the per element results for rated wind speed
  ws = processedOutputs.ratedWind;
  fig = figure();
  colororder(newcolors)
  % Lengths
  subplot(2,2,1)
  hold on
  grid on
  box on
  plot(processedOutputs.Rp(ws,:),':o','linewidth',1,'markersize',3);
  plot(processedOutputs.l_t_inCycle(ws,:),':^','linewidth',1,'markersize',3);
  plot(processedOutputs.h_inCycle(ws,:),':s','linewidth',1,'markersize',3);
  ylabel('Length (m)');
  legend('R_{p}','l_{t}','h_{p}','location','northwest');
  hold off
  % Speeds
  subplot(2,2,2)
  hold on
  grid on
  box on
  yyaxis left
  plot(processedOutputs.lambda(ws,:),':o','linewidth',1,'markersize',3);
  plot(processedOutputs.f(ws,:),':s','linewidth',1,'markersize',3);
  plot(processedOutputs.f_i(ws,:),':v','linewidth',1,'markersize',3);
  ylabel('Velocity factors (-)');
  yyaxis right
  plot(processedOutputs.vw(ws,:),':s','linewidth',1,'markersize',3);
  plot(processedOutputs.vw_i(ws,:),':v','linewidth',1,'markersize',3);
  ylabel('Wind speed (m/s)');
  legend('λ','f_o','f_i','v_{w,o}','v_{w,i}','location','northwest');
  hold off
  % Glide ratios
  subplot(2,2,3)
  hold on
  grid on
  box on
  plot(processedOutputs.E(ws,:),':s','linewidth',1,'markersize',3);
  plot(processedOutputs.E_i(ws,:),':v','linewidth',1,'markersize',3);
  % ylim([0 2.5]);
  ylabel('Glide ratio (-)');
  legend('E_o','E_i','location','northwest');
  hold off
  % Forces
  subplot(2,2,4)
  hold on
  grid on
  box on
  plot(processedOutputs.Fa(ws,:)/1e3,':^','linewidth',1,'markersize',3);
  plot(processedOutputs.Fa_i(ws,:)/1e3,':o','linewidth',1,'markersize',3);
  plot(processedOutputs.Ft(ws,:)/1e3,':s','linewidth',1,'markersize',3);
  plot(processedOutputs.Ft_i(ws,:)/1e3,':v','linewidth',1,'markersize',3);
  plot(processedOutputs.W(ws,:)/1e3,':x','linewidth',1,'markersize',3);
  legend('F_{a,o}','F_{a,i}','F_{t,o}','F_{t,i}','F_{g}','location','northwest');
  ylabel('Force (kN)');
  hold off
  han=axes(fig,'visible','off'); 
  han.Title.Visible='on';
  han.XLabel.Visible='on';
  han.YLabel.Visible='off';
  ylabel(han,'yourYLabel');
  xlabel(han,'Discretized reeling distance over cycle');
  title(han,['Wind speed at 100 m = ' num2str(ws) 'm/s']);

  % Speeds
  figure('units','inch','Position', [0.2 0.5 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  plot(vw, calcRowMean(processedOutputs.lambda, cutIn),':o','linewidth',1,'markersize',3);
  plot(vw, calcRowMean(processedOutputs.f, cutIn),':s','linewidth',1,'markersize',3);
  plot(vw, calcRowMean(processedOutputs.f_i, cutIn),':v','linewidth',1,'markersize',3);
  ylabel('Velocity factors (-)');
  legend('λ','f_{o}','f_{i}','location','northwest');
  xlabel('Wind speed at 100 m height (m/s)');
  xlim(x_axis_limit)
  hold off

  % CL
  figure('units','inch','Position', [4.2 0.5 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  plot(vw, calcRowMean(processedOutputs.CL, cutIn),'o:','linewidth',1,'markersize',3);
  plot(vw, calcRowMean(processedOutputs.CL_i, cutIn),'d:','linewidth',1,'markersize',3);
  ylabel('Lift coefficient (-)');
  legend('C_{L,o}','C_{L,i}','location','northwest','Orientation','vertical');
  xlabel('Wind speed at 100 m height (m/s)');
  xlim(x_axis_limit)
  hold off

  % CD
  figure('units','inch','Position', [8.2 0.5 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  plot(vw, calcRowMean(processedOutputs.CD, cutIn),'o:','linewidth',1,'markersize',3);
  plot(vw, calcRowMean(processedOutputs.CD_i, cutIn),'d:','linewidth',1,'markersize',3);
  plot(vw, calcRowMean(processedOutputs.CD_k, cutIn),'^:','linewidth',1,'markersize',3);
  plot(vw, calcRowMean(processedOutputs.CD_k_i, cutIn),'v:','linewidth',1,'markersize',3);
  plot(vw, calcRowMean(processedOutputs.CD_t, cutIn),'s:','linewidth',1,'markersize',3);
  ylabel('Drag coefficient (-)');
  legend('C_{D,o}','C_{D,i}','C_{D,k,o}','C_{D,k,i}','C_{D,t}','location','northwest','Orientation','vertical');
  xlabel('Wind speed at 100 m height (m/s)');
  xlim(x_axis_limit)
  hold off


  % Lengths
  figure('units','inch','Position', [12.2 0.5 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  plot(vw, processedOutputs.h_cycleAvg(cutIn:end),'d:','linewidth',1,'markersize',3);
  plot(vw, calcRowMean(processedOutputs.Rp, cutIn),'o:','linewidth',1,'markersize',3);
  plot(vw, processedOutputs.deltaL(cutIn:end),'^:','linewidth',1,'markersize',3);
  plot(vw, processedOutputs.l_t_max(cutIn:end),'s:','linewidth',1,'markersize',3);
  plot(vw, processedOutputs.l_t_min(cutIn:end),'x:','linewidth',1,'markersize',3);
  ylabel('Lengths (m)');
  legend('h_{p,avg}','R_{p,avg}','Δl','l_{t,max}','l_{t,0}','location','northwest','Orientation','vertical');
  xlabel('Wind speed at 100 m height (m/s)');
  xlim(x_axis_limit)
  hold off

  % Roll angle, avg patt elevation
  figure('units','inch','Position', [16.2 0.5 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  plot(vw, calcRowMean(processedOutputs.rollAngle, cutIn),'s:','linewidth',1,'markersize',3);
  plot(vw, processedOutputs.beta(cutIn:end),'o:','linewidth',1,'markersize',3);
  plot(vw, processedOutputs.gamma(cutIn:end),'d:','linewidth',1,'markersize',3);
  ylabel('Angles (deg.)');
  legend('$$\Psi_\mathrm{p}$$','$$\beta_\mathrm{p}$$','$$\gamma_\mathrm{p}$$','location','northwest','Orientation','vertical','Interpreter','latex');
  xlabel('Wind speed at 100 m height (m/s)');
  xlim(x_axis_limit)
  hold off

  % Reel-out forces
  figure('units','inch','Position', [0.2 3.5 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  plot(vw, calcRowMean(processedOutputs.W, cutIn)./10^3,':o','markersize',3)
  plot(vw, calcRowMean(processedOutputs.Fa, cutIn)./10^3,':s','linewidth',1,'markersize',3);
  plot(vw, calcRowMean(processedOutputs.Ft, cutIn)./10^3,':^','linewidth',1,'markersize',3);
  ylabel('Force (kN)');
  legend('F_{g}','F_{a,o}','F_{t,o}','location','northwest');
  xlabel('Wind speed at 100 m height (m/s)');
  xlim(x_axis_limit)
  hold off

  % Reel-in forces
  figure('units','inch','Position', [4.2 3.5 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  plot(vw, calcRowMean(processedOutputs.W, cutIn)./10^3,':o','markersize',3)
  plot(vw, calcRowMean(processedOutputs.Fa_i, cutIn)./10^3,':s','linewidth',1,'markersize',3);
  plot(vw, calcRowMean(processedOutputs.Ft_i, cutIn)./10^3,':^','linewidth',1,'markersize',3);
  ylabel('Force (kN)');
  legend('F_{g}','F_{a,i}','F_{t,i}','location','northwest');
  xlabel('Wind speed at 100 m height (m/s)');
  xlim(x_axis_limit)
  hold off


  % Time, num of patterns
  figure('units','inch','Position', [8.2 3.5 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  yyaxis left
  plot(vw, processedOutputs.to(cutIn:end),':s','linewidth',1,'markersize',3);
  plot(vw, processedOutputs.ti(cutIn:end),':d','linewidth',1,'markersize',3);
  plot(vw, calcRowMean(processedOutputs.tPatt, cutIn),':o','linewidth',1,'markersize',3);
  ylabel('Time durations (s)');
  yyaxis right
  plot(vw, calcRowMean(processedOutputs.numOfPatt,cutIn),':^','linewidth',1,'markersize',3);
  ylabel('Number of patterns(-)');
  xlabel('Wind speed at 100 m height (m/s)');
  legend('t_{o}','t_{i}','t_{patt,avg}','N_{p}','location','northwest');
  xlim(x_axis_limit)
  hold off

  % Reel-out power
  figure('units','inch','Position', [12.2 3.5 3.5 2.2])
  colororder(newcolors)
  hold on
%   grid on
  box on
  plot(vw, processedOutputs.P_m_o(cutIn:end)./10^3,':o','linewidth',1,'markersize',3);
  plot(vw, processedOutputs.P_e_o(cutIn:end)./10^3,':s','linewidth',1,'markersize',3);
  plot(vw, processedOutputs.P_e_avg(cutIn:end)./10^3,':x','linewidth',1,'markersize',5);
  [~, maxIndex] = max(round(mean(processedOutputs.Ft,2)));
  xline(maxIndex(1));
  xline(processedOutputs.ratedWind);
  ylabel('Power (kW)');
  legend('P_{m,o}','P_{e,o}','P_{e,avg}','location','northwest');
  xlabel('Wind speed at 100 m height (m/s)');
  xlim(x_axis_limit)
  hold off

  % Reel-in power
  figure('units','inch','Position', [16.2 3.5 3.5 2.2])
  colororder(newcolors)
  hold on
  grid on
  box on
  plot(vw, processedOutputs.P_m_i(cutIn:end)./10^3,':o','linewidth',1,'markersize',3);
  plot(vw, processedOutputs.P_e_i(cutIn:end)./10^3,':s','linewidth',1,'markersize',3);
  ylabel('Power (kW)');
  %title('Cycle averages');
  legend('P_{m,i}','P_{e,i}','location','northwest');
  xlabel('Wind speed at 100 m height (m/s)');
  xlim(x_axis_limit)
  hold off
 
  % Cycle timeseries plots: Pattern averages
  windSpeeds = [processedOutputs.cutIn, processedOutputs.ratedWind, processedOutputs.cutOut];
  for i = windSpeeds
      figure('units','inch','Position', [4.2 6.5 3.5 2.2])
      hold on
      grid on
      box on
      % Plot your data
      plot(processedOutputs.cyclePowerRep(i).t_inst, processedOutputs.cyclePowerRep(i).P_e_inst,'linewidth',1.5, 'DisplayName', 'P_{e}');
      plot(processedOutputs.cyclePowerRep(i).t_inst,processedOutputs.cyclePowerRep(i).P_m_inst,':','linewidth',1.5, 'DisplayName', 'P_{m}');
      yline(processedOutputs.P_e_avg(i)/10^3, '-.', 'linewidth',1.2, 'DisplayName', 'P_{e,avg}');
      yline(processedOutputs.P_m_avg(i)/10^3, ':', 'linewidth',1.2, 'DisplayName', 'P_{m,avg}');
      % Set yline positions within the range of your data
      yline(0, '-', 'linewidth', 1, 'HandleVisibility', 'off');      
      ylabel('Power (kW)');
      xlabel('Time (s)');      
      % Create a custom legend with 'Data1' and 'Data2' only
      legend('Location', 'Best');      
      % Customize legend colors
      legend('TextColor', 'k', 'Color', 'w'); % 'TextColor' sets the text color, 'Color' sets the background color      
      xlim([0 processedOutputs.cyclePowerRep(i).t_inst(end)]);
      ylim([1.5*min(processedOutputs.cyclePowerRep(i).P_e_inst) 1.05*max(processedOutputs.cyclePowerRep(i).P_m_inst)]);
      title(strcat('Wind speed at 100 m:', num2str(i), ' m/s'));      
      hold off
  end

  % Wind profile
  if inputs.vertWindProfile == 0
    Vref = 10; % m/s
    z = 10:10:600; % m
    V = Vref * (z/inputs.h_ref).^inputs.windShearExp/Vref;
    %   MegAWES Onshore location Cabauw. Wind speeds normalized with value at 100 m. Profiles from: https://doi.org/10.1016/j.renene.2022.06.094
    z_MegAWES = [10,20,40,60,80,100,120,140,150,160,180,200,220,250,300,500,600];
    V_MegAWES_Cabauw = [0.541219682843206,0.607355091566827,0.768630154201962,0.868484406441142,0.941395360902529,1,1.04810058627160,1.08638854381156,1.10277338731106,1.11715868927737,1.14412258234309,1.16573551308321,1.18394938534465,1.20653423381438,1.23266397972046,1.26662287360302,1.26414483994687];
    V_MegAWES_Ijmuiden = [0.847612611633547,0.870603040595613,0.927240267828556,0.959346286990695,0.982291573490674,1,1.01377720773809,1.02356771954493,1.02766760602000,1.03079423355205,1.03659625208888,1.04025827758100,1.04284618416620,1.04496440015282,1.04461712713371,1.02473617783789,1.01076976884552];
    figure('units','inch','Position', [8.2 6.5 2 2.2])
    hold on
    box on
    grid on
    plot(V,z,'linewidth',1.2)
    plot(V_MegAWES_Cabauw,z_MegAWES,'--','linewidth',1.2)
    plot(V_MegAWES_Ijmuiden,z_MegAWES,'-.','linewidth',1.2)
    legend('α = 0.143','Cabauw,NL','Ijmuiden,NL');
    xlim([0.5 1.5])
    xlabel('Wind Speed (-)')
    ylabel('Height (m)')
    hold off
  else

    % Continuous function from vertical wind profile
    % Define the heights at which you want to interpolate
    heights = [0:10:1000]; % Adjust the desired height values as needed
    % Interpolate at the specified heights
    V_interpolated = interp1(inputs.windProfile_h, inputs.windProfile_vw, heights, 'linear', 'extrap');
    % Plot the interpolated function
    figure('units','inch','Position', [15 6 2 2.2])
    hold on 
    box on
    grid on
    plot(inputs.windProfile_vw, inputs.windProfile_h, 'o', V_interpolated, heights, '-');
    xlabel('Wind speed (m/s)');
    ylabel('Height (m)');
    xlim([0.5 1.5])
    title('Wind profile');
    legend('Data', 'Interpolated', 'Location', 'Best');
    hold off
  end

   
end

function rowMean = calcRowMean(array, cutIn)
  rowMean = mean(array,2);
  rowMean = rowMean(cutIn:end);
end