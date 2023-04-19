%% Results plot
% Power curve comparison and sensitivities
clc
clearvars

% modelSpreadsheet    = "Ampyx LCOE model (with KPIs) 2021-05-28 FN MK2 AP3 trial_copy.xlsx";
sensitivityStudyAP4 = "Sensitivity study AP4_copy.xlsx";

%% Model results

%% AP3 6DoF simulation results
AP3.PC.ws    = [7.32E+00, 8.02E+00, 9.05E+00, 1.00E+01, 1.10E+01, 1.20E+01, ...
  1.30E+01, 1.41E+01, 1.50E+01, 1.60E+01, 1.70E+01, 1.80E+01, 1.90E+01]; %[m/s]

AP3.PC.power = [0.00E+00, 7.01E+03, 2.37E+04, 4.31E+04, 6.47E+04, 8.46E+04, ...
  1.02E+05, 1.20E+05, 1.34E+05, 1.49E+05, 1.50E+05, 1.50E+05, 1.50E+05]./10^3; %[kW]

%% Sensitivity plots
ws_sensitivity = xlsread(sensitivityStudyAP4, 'Sheet1','B5:B24');

genEff.p1 = xlsread(sensitivityStudyAP4, 'Sheet1','C5:C24')./10^6;
genEff.p2 = xlsread(sensitivityStudyAP4, 'Sheet1','D5:D24')./10^6;

genPA.p1 = xlsread(sensitivityStudyAP4, 'Sheet1','Z5:Z24')./10^6;
genPA.p2 = xlsread(sensitivityStudyAP4, 'Sheet1','AA5:AA24')./10^6;
genPA.p3 = xlsread(sensitivityStudyAP4, 'Sheet1','AB5:AB24')./10^6;
genPA.p4 = xlsread(sensitivityStudyAP4, 'Sheet1','AC5:AC24')./10^6;
genPA.p5 = xlsread(sensitivityStudyAP4, 'Sheet1','AD5:AD24')./10^6;

RI.p1 = xlsread(sensitivityStudyAP4, 'Sheet1','AG5:AG24')./10^6;
RI.p2 = xlsread(sensitivityStudyAP4, 'Sheet1','AH5:AH24')./10^6;
RI.p3 = xlsread(sensitivityStudyAP4, 'Sheet1','AI5:AI24')./10^6;
RI.p4 = xlsread(sensitivityStudyAP4, 'Sheet1','AJ5:AJ24')./10^6;
RI.p5 = xlsread(sensitivityStudyAP4, 'Sheet1','AK5:AK24')./10^6;

Te.p1 = xlsread(sensitivityStudyAP4, 'Sheet1','G5:G24')./10^6;
Te.p2 = xlsread(sensitivityStudyAP4, 'Sheet1','H5:H24')./10^6;
Te.p3 = xlsread(sensitivityStudyAP4, 'Sheet1','I5:I24')./10^6;
Te.p4 = xlsread(sensitivityStudyAP4, 'Sheet1','J5:J24')./10^6;

CD.p1 = xlsread(sensitivityStudyAP4, 'Sheet1','L5:L24')./10^6;
CD.p2 = xlsread(sensitivityStudyAP4, 'Sheet1','M5:M24')./10^6;
CD.p3 = xlsread(sensitivityStudyAP4, 'Sheet1','N5:N24')./10^6;

mass.p1 = xlsread(sensitivityStudyAP4, 'Sheet1','Q5:Q24')./10^6;
mass.p2 = xlsread(sensitivityStudyAP4, 'Sheet1','R5:R24')./10^6;
mass.p3 = xlsread(sensitivityStudyAP4, 'Sheet1','S5:S24')./10^6;
mass.p4 = xlsread(sensitivityStudyAP4, 'Sheet1','T5:T24')./10^6;
mass.p5 = xlsread(sensitivityStudyAP4, 'Sheet1','U5:U24')./10^6;
mass.p6 = xlsread(sensitivityStudyAP4, 'Sheet1','V5:V24')./10^6;

CLmax.p1 = xlsread(sensitivityStudyAP4, 'Sheet1','AM5:AM24')./10^6;
CLmax.p2 = xlsread(sensitivityStudyAP4, 'Sheet1','AN5:AN24')./10^6;
CLmax.p3 = xlsread(sensitivityStudyAP4, 'Sheet1','AO5:AO24')./10^6;
CLmax.p4 = xlsread(sensitivityStudyAP4, 'Sheet1','AP5:AP24')./10^6;
CLmax.p5 = xlsread(sensitivityStudyAP4, 'Sheet1','AQ5:AQ24')./10^6;

CLimprov.p1 = xlsread(sensitivityStudyAP4, 'Sheet1','AS5:AS24')./10^6;
CLimprov.p2 = xlsread(sensitivityStudyAP4, 'Sheet1','AT5:AT24')./10^6;
CLimprov.p3 = xlsread(sensitivityStudyAP4, 'Sheet1','AU5:AU24')./10^6;
CLimprov.p4 = xlsread(sensitivityStudyAP4, 'Sheet1','AV5:AV24')./10^6;


%%

figSize = [5 5 3.2 2];

% Generator eff
figure('units','inch','Position', figSize)
hold on
grid on
box on
plot(ws_sensitivity, genEff.p1,'linewidth',1);
plot(ws_sensitivity, genEff.p2,'--','linewidth',1);
legend('Variable efficiency','Constant efficiency', 'location','southeast');
xlabel('Wind speed at pattern altitude (m/s)');
ylabel('Power (MW)');
xlim([3 19]);
ylim([0 1.05]);
hold off

% Generator oversizing Peak/rated
figure('units','inch','Position', figSize)
hold on
grid on
box on
plot(ws_sensitivity, genPA.p1,'linewidth',1);
plot(ws_sensitivity, genPA.p2','linewidth',1);
plot(ws_sensitivity, genPA.p3,'linewidth',1);
plot(ws_sensitivity, genPA.p4,'linewidth',1);
plot(ws_sensitivity, genPA.p5,'linewidth',1);
legend('P/A=3','P/A=2.75','P/A=2.5','P/A=2.25','P/A=2', 'location','southeast');
xlabel('Wind speed at pattern altitude (m/s)');
ylabel('Power (MW)');
xlim([3 19]);
ylim([0 1.05]);
hold off

% Reel-in speed effect
figure('units','inch','Position', figSize)
hold on
grid on
box on
plot(ws_sensitivity, RI.p1,'linewidth',1);
plot(ws_sensitivity, RI.p2','linewidth',1);
plot(ws_sensitivity, RI.p3,'linewidth',1);
plot(ws_sensitivity, RI.p4,'linewidth',1);
plot(ws_sensitivity, RI.p5,'linewidth',1);
legend('15m/s','20m/s','25m/s','30m/s','35m/s', 'location','southeast');
xlabel('Wind speed at pattern altitude (m/s)');
ylabel('Power (MW)');
xlim([3 19]);
ylim([0 1.05]);
hold off

% Tether material strength
figure('units','inch','Position', figSize)
hold on
grid on
box on
plot(ws_sensitivity, Te.p1,'linewidth',1);
plot(ws_sensitivity, Te.p2','linewidth',1);
plot(ws_sensitivity, Te.p3,'linewidth',1);
plot(ws_sensitivity, Te.p4,'linewidth',1);
legend('900 MPa','700 MPa','500 MPa','300 MPa', 'location','southeast');
xlabel('Wind speed at pattern altitude (m/s)');
ylabel('Power (MW)');
xlim([3 19]);
ylim([0 1.05]);
hold off

% Effect of drag coefficients
figure('units','inch','Position', figSize)
hold on
grid on
box on
plot(ws_sensitivity, CD.p1,'linewidth',1);
plot(ws_sensitivity, CD.p2','linewidth',1);
plot(ws_sensitivity, CD.p3,'linewidth',1);
legend('C_D=1.2','C_D=1','C_D=0.8', 'location','southeast');
xlabel('Wind speed at pattern altitude (m/s)');
ylabel('Power (MW)');
xlim([3 19]);
ylim([0 1.05]);
hold off

% CL max effect
figure('units','inch','Position', figSize)
hold on
grid on
box on
plot(ws_sensitivity, CLmax.p1,'linewidth',1);
plot(ws_sensitivity, CLmax.p2','linewidth',1);
plot(ws_sensitivity, CLmax.p3,'linewidth',1);
plot(ws_sensitivity, CLmax.p4,'linewidth',1);
plot(ws_sensitivity, CLmax.p5,'linewidth',1);
legend('C_L=3.2','C_L=3.0','C_L=2.7','C_L=2.4','C_L=2.1', 'location','southeast');
xlabel('Wind speed at pattern altitude (m/s)');
ylabel('Power (MW)');
xlim([3 19]);
ylim([0 1.05]);
hold off

% CL improvement potential
figure('units','inch','Position', figSize)
hold on
grid on
box on
plot(ws_sensitivity, CLmax.p1,'linewidth',1);
plot(ws_sensitivity, CLmax.p2','linewidth',1);
plot(ws_sensitivity, CLmax.p3,'linewidth',1);
plot(ws_sensitivity, CLmax.p4,'linewidth',1);
legend('C_L=2.1','C_L=2.7','C_L=3.3','C_L=3.9', 'location','southeast');
xlabel('Wind speed at pattern altitude (m/s)');
ylabel('Power (MW)');
xlim([3 19]);
ylim([0 1.05]);
hold off

% Mass sensitivity
figure('units','inch','Position', figSize)
hold on
grid on
box on
plot(ws_sensitivity, mass.p1,'linewidth',1);
plot(ws_sensitivity, mass.p2','linewidth',1);
plot(ws_sensitivity, mass.p3,'linewidth',1);
plot(ws_sensitivity, mass.p4,'linewidth',1);
plot(ws_sensitivity, mass.p5,'linewidth',1);
plot(ws_sensitivity, mass.p6,'linewidth',1);
legend('m=1000kg','m=2000kg','m=3000kg','m=4000kg','m=5000kg','m=6000kg', 'location','southeast');
xlabel('Wind speed at pattern altitude (m/s)');
ylabel('Power (MW)');
xlim([3 19]);
ylim([0 1.05]);
hold off


% %% Wind speed vector
% ws                       = xlsread(modelSpreadsheet, 'Performance', 'A83:A108'); % [m/s]
% 
% %% TT plot
% TT   = xlsread(modelSpreadsheet, 'Performance', 'B83:B108')./10^3; %[kN]
% V_RO = xlsread(modelSpreadsheet, 'Performance', 'G83:G108'); %[m/s]
% 
% figure('units','inch','Position', [5 5 3.5 2.2])
% hold on
% grid on
% box on
% yyaxis left
% plot(ws, TT,'--','linewidth',1.5);
% ylabel('Tether tension (kN)');
% yyaxis right
% plot(ws, V_RO,'o','linewidth',1.5,'markersize',2);
% legend('Tether tension','Reel-out speed');
% xlabel('Wind speed at pattern altitude (m/s)');
% ylabel('Reel-out speed (m/s)');
% xlim([0 26]);
% %ylim([0 160]);
% hold off
% 
% %% RO and RI time plot
% V_RO = xlsread(modelSpreadsheet, 'Performance', 'G112:G137');   %[m/s]
% a    = xlsread(modelSpreadsheet, 'Inputs', 'P525'); %[m/s]
% V_RI = xlsread(modelSpreadsheet, 'Inputs', 'P521');  %[m/s]
% deltaL = xlsread(modelSpreadsheet, 'Performance', 'C39') - xlsread(modelSpreadsheet, 'Performance', 'C38');  %[m]
% 
% t_RO = deltaL./V_RO + V_RO./a;
% t_RI = deltaL/V_RI + V_RI/a;
% t_cycle = t_RO+t_RI;
% 
% t_RO(t_RO==Inf) = 0;
% t_RI = t_RI.*ones(length(t_RO),1);
% 
% figure('units','inch','Position', [5 5 3.5 2.2])
% hold on
% grid on
% box on
% plot(ws(6:end), t_RO(5:end),'--','linewidth',1.5);
% plot(ws(6:end), t_RI(5:end),':','linewidth',1.5);
% plot(ws(6:end), t_cycle(5:end),'-.','linewidth',1.5);
% legend('Reel-out time','Reel-in time','Cycle time');
% xlabel('Wind speed at pattern altitude (m/s)');
% ylabel('Time(s)');
% xlim([0 26]);
% %ylim([0 160]);
% hold off
% 
% 
% %% Variable efficiencies plot
% eta_DT    = xlsread(modelSpreadsheet, 'Performance', 'Z83:Z108').*100; %[%]
% eta_RI    = xlsread(modelSpreadsheet, 'Performance', 'AC83:AC108').*100; %[%]
% eta_Cycle = xlsread(modelSpreadsheet, 'Performance', 'W83:W108').*100; %[%]
% eta_tot   = xlsread(modelSpreadsheet, 'Performance', 'AE83:AE108').*100; %[%]
% 
% figure('units','inch','Position', [5 5 3.5 2.2])
% hold on
% grid on
% box on
% plot(ws, eta_DT,'--','linewidth',1.5);
% plot(ws, eta_RI,':','linewidth',1.5);
% plot(ws, eta_Cycle,'-.','linewidth',1.5);
% plot(ws, eta_tot,'-','linewidth',1.5);
% legend('Drivetrain eff.','Reel-in eff.','Cycle eff.','Total eff.');
% xlabel('Wind speed at pattern altitude (m/s)');
% ylabel('Efficiency (%)');
% xlim([0 26]);
% %ylim([0 160]);
% hold off
% 
% 
% %% Power curve comparison plot
% analytical_PC_simple     = xlsread(modelSpreadsheet, 'Simple models', 'K5:K29')./10^3;   %[kW]
% analytical_PC_allEffects = xlsread(modelSpreadsheet, 'Performance', 'AJ83:AJ108')./10^3; %[kW]
% 
% figure('units','inch','Position', [5 5 3.5 2.2])
% hold on
% grid on
% box on
% plot(ws(2:end), analytical_PC_simple,'--','linewidth',1.5);
% plot(ws, analytical_PC_allEffects,'linewidth',1.5);
% plot(AP3.PC.ws, AP3.PC.power,'k^--','MarkerSize',4);
% legend('Analytical ideal','Analytical all effects','6DOF simulation');
% xlabel('Wind speed at pattern altitude (m/s)');
% ylabel('Power (kW)');
% xlim([0 26]);
% ylim([0 160]);
% hold off

%% 
