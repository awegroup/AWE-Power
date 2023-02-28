% Performance model
clc
clearvars
clear global

inputSheet_AP3;
% inputSheet;
% inputSheet_provided_HV;
% inputSheet_softKite;

%% Optimise operation for every wind speed: Uncapped electrical power
%      [deltaL, VRI, CL, avgPattEle, pattAngRadius, startPattRadius]
x0      = [100, 4, 1.5, deg2rad(20),deg2rad(12),200]; % 
for i=1:length(inputs.Vw_ref)
  % Output of previous wind speed as input to next wind speed
  x_init = x0;  
%  x_init    = [100, 4, 1.5, deg2rad(20),deg2rad(12),200]; % 
  x0     = x_init./x_init;
  lb     = [20, 2, inputs.CL0_airfoil, deg2rad(5), deg2rad(5),sqrt(inputs.AR*inputs.WA)/2]./x_init; % 
  ub     = [500, inputs.maxVRI, inputs.CL_maxAirfoil*inputs.F_CLeff, deg2rad(80),deg2rad(80),500]./x_init; % 
  options                           = optimoptions('fmincon');
  options.Display                   = 'iter-detailed';
  options.Algorithm                 = 'sqp';
  options.FiniteDifferenceType      = 'central';
  %options.FiniteDifferenceStepSize  = [1e-3 1e-3 1e-6 1e-7 1e-7 1e-3];
%   options.FiniteDifferenceStepSize  = 1e-6;
%   options.OptimalityTolerance       = 1e-12;
%   options.StepTolerance             = 1e-6;
  options.MaxFunctionEvaluations    = 5000*numel(x_init);
  options.MaxIterations             = 1000*numel(x_init);
%   options.ConstraintTolerance      = 1e-3;
%   options.FunctionTolerance        = 1e-9;
%   options.DiffMaxChange            = 1e-1;
%   options.DiffMinChange            = 0;
%   options.OutputFcn                = @O_outfun;
%    options.PlotFcns                 = {@optimplotfval, @optimplotx, @optimplotfirstorderopt,...
%                                @optimplotconstrviolation, @optimplotfunccount, @optimplotstepsize};
  con = @(x) constraints(i,inputs);
  
  [x,fval,exitflag(i),optHist(i),lambda(i)] = fmincon(@(x) objective(x,x_init,i,inputs),x0,[],[],[],[],lb,ub,con,options);
  
  % Storing final results
   [~,inputs,outputs] = objective(x,x_init,i,inputs);
   x0 = x.*x_init;
     
   
   % Changing initial guess if previous wind speed evaluation is infeasible
   if outputs.P_cycleElec(i) <= 0
       x0 = [100, 5, 1.5, deg2rad(20),deg2rad(12),200]; % 
   end  
end

%% Second optimisation iteration for following mean of capped electrical power from the first iteration
outputs1 = outputs;
clear outputs
inputs.targetPRO_elec = outputs1.PROeff_elec_cap;
%      [deltaL, VRI, CL, avgPattEle, pattAngRadius, startPattRadius]
x0      = [100, 4, 1.5, deg2rad(20),deg2rad(12),200]; % 
for i=1:length(inputs.Vw_ref)
  % Output of previous wind speed as input to next wind speed
  x_init = x0;  
%  x_init    = [100, 4, 1.5, deg2rad(20),deg2rad(12),200]; % 
  x0     = x_init./x_init;
  lb     = [20, 2, inputs.CL0_airfoil, deg2rad(5), deg2rad(5),sqrt(inputs.AR*inputs.WA)/2]./x_init; % 
  ub     = [500, inputs.maxVRI, inputs.CL_maxAirfoil*inputs.F_CLeff, deg2rad(80),deg2rad(80),500]./x_init; % 
  options                           = optimoptions('fmincon');
  options.Display                   = 'iter-detailed';
  options.Algorithm                 = 'sqp';
  options.FiniteDifferenceType      = 'central';
  options.MaxFunctionEvaluations    = 5000*numel(x_init);
  options.MaxIterations             = 1000*numel(x_init);
  con = @(x) constraints(i,inputs);
  
  [x,fval,exitflag(i),optHist(i),lambda(i)] = fmincon(@(x) objective(x,x_init,i,inputs),x0,[],[],[],[],lb,ub,con,options);
  
  % Storing final results
   [~,inputs,outputs] = objective(x,x_init,i,inputs);
   x0 = x.*x_init;
     
   
   % Changing initial guess if previous wind speed evaluation is infeasible
   if outputs.P_cycleElec(i) <= 0
       x0 = [100, 5, 1.5, deg2rad(20),deg2rad(12),200]; % 
   end  
end
% Storing back the capped powerresults from first optimisation
outputs.PROeff_elec_osci = outputs1.PROeff_elec_osci_cap;

%% Post processing
Vw = inputs.Vw_ref; 

%% Cut-in wind speed
% temp1         = Vw((outputs.VRO-outputs.VRO_osciAmp)>0);
temp2         = Vw(outputs.P_cycleElec > 0);
temp5         = Vw(mean(outputs.VRO(:,:),2) >=1);
% postProRes.cutIn = max(temp1(1), temp2(1),temp(5));
postProRes.cutIn = max(temp2(1),temp5(1));

%% Rated wind and power
temp3 = round(outputs.P_cycleElec./max(outputs.P_cycleElec),2);
temp4 = Vw(temp3==1);
postProRes.ratedWind  = temp4(1);
postProRes.ratedPower = outputs.P_cycleElec(postProRes.ratedWind);

%% System data
% system.VRO_osci         = zeros(length(Vw),length(outputs.VRO_osci(1,:)));
% system.PROeff_mech_osci = zeros(length(Vw),length(outputs.VRO_osci(1,:)));
% system.PROeff_elec_osci = zeros(length(Vw),length(outputs.VRO_osci(1,:)));
% system.numPattParts     = outputs.numPattParts;
postProRes.D_te         = outputs.D_te;
for i=1:length(Vw)
    if Vw(i)>=postProRes.cutIn
        postProRes.Vw(i,:)           = outputs.Vw(i,:);
        postProRes.PRO1_mech(i)    = outputs.PRO1_mech(i);
        postProRes.PRI2_mech(i)    = outputs.PRI2_mech(i); 
        postProRes.PROeff_mech(i,:)  = outputs.PROeff_mech(i,:);
        postProRes.PRIeff_mech(i,:)  = outputs.PRIeff_mech(i,:);  
        postProRes.PRO1_elec(i)    = outputs.PRO1_elec(i);
        postProRes.PRI2_elec(i)    = outputs.PRI2_elec(i); 
        postProRes.PROeff_elec(i,:)  = outputs.PROeff_elec(i,:);
        postProRes.PRIeff_elec(i,:)  = outputs.PRIeff_elec(i,:);
        postProRes.PRO_mech(i)    = outputs.PRO_mech(i);
        postProRes.PRI_mech(i)    = outputs.PRI_mech(i);
        postProRes.PRO_elec(i)    = outputs.PRO_elec(i);
        postProRes.PRI_elec(i)    = outputs.PRI_elec(i);
        postProRes.deltaL(i)      = outputs.deltaL(i);
        postProRes.t1(i)          = outputs.t1(i);
        postProRes.t2(i)          = outputs.t2(i);
        postProRes.tROeff(i,:)      = outputs.tROeff(i,:);
        postProRes.tRIeff(i,:)      = outputs.tRIeff(i,:);
        postProRes.tRO(i)         = outputs.tRO(i);
        postProRes.tRI(i)         = outputs.tRI(i);
        postProRes.tCycle(i)      = outputs.tCycle(i);
        postProRes.TT(i,:)          = outputs.T(i,:);
        postProRes.VRO(i,:)         = outputs.VRO(i,:);
        postProRes.VRI(i)         = outputs.VRI(i);
        postProRes.VRI_app(i,:)     = outputs.VA_RI(i,:);
        postProRes.Pcycle_elec(i) = outputs.P_cycleElec(i);
        postProRes.Pcycle_mech(i) = outputs.P_cycleMech(i);
%         system.VRO_osci(i,:)  = outputs.VRO_osci(i,:);
%         system.PROeff_mech_osci(i,:)  = outputs.PROeff_mech_osci(i,:);
%         system.PROeff_elec_osci(i,:)  = outputs.PROeff_elec_osci(i,:);
        postProRes.pattRad(i,:)       = outputs.pattRadius(i,:);
        postProRes.H_cycleStart(i)    = outputs.H_cycleStart(i);
        postProRes.H_cycleAvg(i)      = outputs.H_cycleAvg(i);
        postProRes.L_teMax(i)         = outputs.L_teMax(i);
        postProRes.avgRollAngle(i)    = rad2deg(outputs.avgRollAngle(i));
        postProRes.avgPattEle(i)      = rad2deg(outputs.avgPattEle(i));
        postProRes.pattAngRadius(i)   = rad2deg(outputs.pattAngRadius(i));
        postProRes.CL(i)              = outputs.CL(i);
        postProRes.CD(i)              = outputs.CD(i);
        postProRes.numOfPatt(i,:)       = outputs.numOfPatt(i,:);
        postProRes.reelOutF(i,:)      =  postProRes.VRO(i,:)/postProRes.Vw(i,:);  
        postProRes.dutyCycle(i)       = postProRes.tRO(i)/postProRes.tCycle(i);
        postProRes.tPatt(i,:)         = 2*pi()*postProRes.pattRad(i,:)/outputs.VC(i);
    end
%     for j = 1:outputs.numPattParts
%       system.tPatt(i)  = system.tPatt(i) + 2*pi()*system.pattRadius(i)/outputs.numPattParts/outputs.VC_osci(i,j);
%     end 
%     system.osciFactor(i,:) = outputs.osciFactor(i,:);

%     %% Coefficient of power (Cp): For comparing with wind turbines

%     outputs.sweptArea(i)       = pi()*((outputs.pattRadius(i)+outputs.wingSpan/2)^2 - (outputs.pattRadius(i)-outputs.wingSpan/2)^2);

%     outputs.Cp_reelOut(i) = outputs.PRO_mech(i)/(0.5*outputs.rho_air(i)*outputs.sweptArea(i)*inputs.Vw(i)^3);
%     % Axial induction factor maximising Cp is a = 0.5 and is independent of reel-out factor
%     % Ref paper: The Betz limit applied to Airborne Wind Energy, https://doi.org/10.1016/j.renene.2018.04.034
%     f(i)    = outputs.reelOutF(i);
%     a_ideal = 0.5;
%     outputs.Cp_max_ifaIsHalf(i)  = (1-f(i))^2*4*a_ideal*(1-a_ideal)*f(i);
%     a                            = abs(roots([1 -1 outputs.Cp_reelOut(i)/(4*f(i)*(1-f(i))^2)]));
%     outputs.axialInductionF(i,1) = a(1);
%     outputs.axialInductionF(i,2) = a(2);
end

%% Representative instantaneous cycle data
timeseries = struct();
for i = postProRes.cutIn:length(Vw)
  [timeseries.ws(i)] = createTimeseries(i,postProRes);
end

%% Plots

% Cycle timeseries plots, Check reel-in representation: Time and power in each regime should add to total reel-in energy
windSpeeds = [postProRes.ratedWind];
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
  ylabel('Electrical power (kW)');
  xlabel('Time (s)');
  %xlim([0 tmax]);
  %ylim([pmin pmax]);
  title(strcat('Wind speed at 100m:',num2str(timeseries.ws(i).ws),'m/s'));
%   legend('Electrical','Mechanical','location','northwest');
  hold off
end

% Parameter plots

newcolors = [ % 0.25, 0.25, 0.25
  0 0.4470 0.7410
0.8500 0.3250 0.0980 
0.4660, 0.6740, 0.1880
0.9290, 0.6940, 0.1250
0.4940, 0.1840, 0.5560
0.6350 0.0780 0.1840
0.3010 0.7450 0.9330];

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
legend('V_{RO,avg}','V_{RI}','t_{RO}','t_{RI}','location','northwest');
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
legend('η_{m,cycle}','η_{e,cycle}','Duty cycle','location','southeast');
xlabel('Wind speed at 100m height (m/s)');
xlim([0 25]);
hold off

% Power plots
pastel_b = '[0 0.4470 0.7410]';
figure('units','inch','Position', [11 6 3.5 2.2])
colororder(newcolors)
hold on
grid on
box on
plot(Vw, postProRes.PRI_mech./10^3,'--','linewidth',1.2);
plot(Vw, postProRes.PRO_mech./10^3,'-.','linewidth',1.2);
plot(Vw, postProRes.PRI_elec./10^3,'--','linewidth',1.2);
plot(Vw, postProRes.PRO_elec./10^3,'-.','linewidth',1.2);
plot(Vw, postProRes.Pcycle_elec./10^3,'-','linewidth',1.2);
ylabel('Power (kW)');
%title('Cycle averages');
legend('P_{m,RI}','P_{m,RO}','P_{e,RI}','P_{e,RO}','P_{e,avg}','location','northwest');
xlabel('Wind speed at 100m height (m/s)');
xlim([0 25]);
hold off


%% Power curve comparison plot
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
plot(Vw, P_Loyd./10^3,'--','linewidth',1.5);
plot(Vw, postProRes.Pcycle_elec./10^3,'-','linewidth',1.5);
plot(AP3.PC.ws, AP3.PC.power,'k^--','MarkerSize',4);
ylabel('Power (kW)');
legend('Loyd','Model','6DOF','location','southeast');
% legend('Loyd - Ideal','Model results','location','southeast');
xlabel('Wind speed at 100m height (m/s)');
xlim([0 25]);
hold off

%% Plot oscillating VRO and PRO_mech
% i = [7,12,18, 24];
% d.series1 = system.PROeff_elec_osci(i(1),:)/10^3;  
% d.series2 = system.PROeff_elec_osci(i(2),:)/10^3; 
% d.series3 = system.PROeff_elec_osci(i(3),:)/10^3; 
% d.series4 = system.PROeff_elec_osci(i(4),:)/10^3; 
% 
% figure('units','inch','Position', [4 4 3.5 2.2])
% hold on
% grid on
% box on
% plot(d.series1,'linewidth',1.2);
% plot(d.series2,'linewidth',1.2);
% plot(d.series3,'linewidth',1.2);
% plot(d.series4,'linewidth',1.2);
% ylabel('Electrical capped power (kW)');
% %ylim([30 100]);
% legend(strcat(num2str(i(1)),'m/s'),strcat(num2str(i(2)),'m/s'),strcat(num2str(i(3)),'m/s'),strcat(num2str(i(4)),'m/s'),'location','southeast');
% xlabel('Positions in a single pattern');
% %xlim([0 25]);
% hold off



%% GA
%    outputs = struct();
%     lb     = [50, 5, inputs.CL0_airfoil, deg2rad(2), deg2rad(2),deg2rad(2)];
%     ub     = [250, inputs.maxVRI, inputs.CL_maxAirfoil*inputs.F_CLeff, deg2rad(90),deg2rad(60),deg2rad(90)];
%     options                           = optimoptions('ga');
%     con = @(x) constraints_GA(i,inputs,outputs);
%   obj = @(x) objective_GA(i,inputs,outputs);
%   [x, fval,exitflag(i),optHist(i)] = ga(obj ,numel(lb),[],[],[],[],lb,ub,con,options);
%      %Storing final results
%    [~,inputs,outputs] = objective_GA(x,i,inputs,outputs);



 