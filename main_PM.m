% Performance model
clc;
clear;

inputSheet_AP3_PM;
%inputSheet_PM;

outputs = struct();

for i=1:length(inputs.Vw)
  
  % Optimisation
  lb = [10,5,inputs.Tmax];
  ub = [500,inputs.maxVRI, inputs.Tmax];
  %x0 = (lb + ub)/2; 
  x0 = [50,5, inputs.Tmax];
  [x,fval] = fmincon(@(x) handleFn_PM(x,i,inputs,outputs),x0,[],[],[],[],lb,ub);
  
  % Storing final results
  [~, inputs,outputs] = handleFn_PM(x,i,inputs,outputs);
  
end
%% 
hold on
for i =1:25

        plot(outputs.VRO_osc(i,:));

end

%% Post processing

Vw = inputs.Vw;

% Cut-in wind speed
temp = Vw(outputs.P_cycleElec>0);
system.cutIn = temp(1);

% Rated electrical power
system.ratedPower = min(inputs.P_ratedElec,max(outputs.P_cycleElec));

% Rated wind speed
temp = Vw(outputs.P_cycleElec>=system.ratedPower);
system.ratedWind  = temp(1);

% tRO,tRI, TT, VRO, PRO,PRI,P_cycle
for i=1:length(Vw)
    if Vw(i)<system.cutIn
        system.PRO_mech(i)    = 0;
        system.PRI_mech(i)    = 0;
        system.PRO_elec(i)    = 0;
        system.PRI_elec(i)    = 0;
        system.deltaL(i)      = 0;
        system.tRO(i)         = 0;
        system.tRI(i)         = 0;
        system.TT(i)          = 0;
        system.VRO(i)         = 0;
        system.VRI(i)         = 0;
        system.Pcycle_elec(i) = 0;
    elseif Vw(i)>=system.ratedWind
        system.PRO_mech(i)    = outputs.PRO_mech(system.ratedWind);
        system.PRI_mech(i)    = outputs.PRI_mech(system.ratedWind);
        system.PRO_elec(i)    = outputs.PRO_elec(system.ratedWind);
        system.PRI_elec(i)    = outputs.PRI_elec(system.ratedWind);
        system.deltaL(i)      = outputs.deltaL(system.ratedWind);
        system.tRO(i)         = outputs.tRO(system.ratedWind);
        system.tRI(i)         = outputs.tRI(system.ratedWind);
        system.TT(i)          = outputs.T(i);
        system.VRO(i)         = outputs.VRO(system.ratedWind);
        system.VRI(i)         = outputs.VRI(system.ratedWind);
        system.Pcycle_elec(i) = system.ratedPower;
    else
        system.PRO_mech(i)    = outputs.PRO_mech(i);
        system.PRI_mech(i)    = outputs.PRI_mech(i);
        system.PRO_elec(i)    = outputs.PRO_elec(i);
        system.PRI_elec(i)    = outputs.PRI_elec(i);
        system.deltaL(i)      = outputs.deltaL(i);
        system.tRO(i)         = outputs.tRO(i);
        system.tRI(i)         = outputs.tRI(i);
        system.TT(i)          = outputs.T(i);
        system.VRO(i)         = outputs.VRO(i);
        system.VRI(i)         = outputs.VRI(i);
        system.Pcycle_elec(i) = outputs.P_cycleElec(i);
    end
    system.tCycle(i) = system.tRO(i)+system.tRI(i);    
end

%% Plots

% TT,VRO,VRI plot
figure('units','inch','Position', [5 5 3.5 2.2])
hold on
grid on
box on
yyaxis left
plot(Vw, system.TT./10^3,'x','markersize',4);
ylabel('Tether tension (kN)');
ylim([0 1.1*max(system.TT)/10^3]);
yyaxis right
plot(Vw, system.VRO,'+','markersize',3);
plot(Vw, system.VRI,'o','markersize',3);
legend('T','V_{RO}','V_{RI}','location','southeast');
xlabel('Wind speed at avg. pattern altitude (m/s)');
ylabel('Speed (m/s)');
xlim([0 25]);
ylim([0 1.05*inputs.maxVRI]);
hold off


% Cycle times, reel-out length
figure('units','inch','Position', [5 5 3.5 2.2])
hold on
grid on
box on
yyaxis left
plot(Vw, system.deltaL,'x','markersize',4);
ylabel('Reel-out length(m)');
yyaxis right
plot(Vw, system.tRO,'+','markersize',3);
plot(Vw, system.tRI,'o','markersize',3);
ylabel('Time(s)');
legend('ΔL','t_{RO}','t_{RI}','location','northwest');
xlabel('Wind speed at avg. pattern altitude (m/s)');
xlim([0 25]);
%ylim([0 160]);
hold off


% Cycle efficiency due to reel-in and drivetrain losses
figure('units','inch','Position', [5 5 3.5 2.2])
hold on
grid on
box on
plot(Vw, system.Pcycle_elec./system.PRO_mech*100,'x','markersize',4);
ylabel('Efficiency (%)');
ylim([0 100]);
legend('η_{Cycle}','location','southeast');
xlabel('Wind speed at avg. pattern altitude (m/s)');
xlim([0 25]);
hold off


%Power curve comparison plot
% Loyd
for i = 1:length(Vw)
    P_Loyd(i) = (4/27)*(outputs.CL^3/outputs.CD(i)^2)*(1/2)*outputs.rho_air(i)*inputs.WA*(Vw(i)^3);
    if P_Loyd(i)>system.ratedPower
        P_Loyd(i) = system.ratedPower;
    end 
end
% AP3 6DoF simulation results
AP3.PC.ws    = [7.32E+00, 8.02E+00, 9.05E+00, 1.00E+01, 1.10E+01, 1.20E+01, ...
  1.30E+01, 1.41E+01, 1.50E+01, 1.60E+01, 1.70E+01, 1.80E+01, 1.90E+01]; %[m/s]
AP3.PC.power = [0.00E+00, 7.01E+03, 2.37E+04, 4.31E+04, 6.47E+04, 8.46E+04, ...
  1.02E+05, 1.20E+05, 1.34E+05, 1.49E+05, 1.50E+05, 1.50E+05, 1.50E+05]./10^3; %[kW]

figure('units','inch','Position', [5 5 3.5 2.2])
hold on
grid on
box on
plot(Vw, P_Loyd./10^3,'--','linewidth',1.2);
plot(Vw, system.Pcycle_elec./10^3,'linewidth',1.5);
plot(AP3.PC.ws, AP3.PC.power,'k^--','MarkerSize',4);
legend('Loyd','Model results','6DOF simulation','location','southeast');
xlabel('Wind speed at avg. pattern altitude (m/s)');
ylabel('Power (kW)');
xlim([0 25]);
%ylim([0 160]);
hold off






