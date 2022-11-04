%% Analytical model paper plots script
clc
clearvars

% To access PM model results 
path.datasets  = addpath(genpath('C:\Users\rjoshi1\OneDrive - Delft University of Technology\Rishikesh\PhD\Ampyx\Datasets'));

% To access the analytical model
ExcelModelName = 'C:\Users\rjoshi1\OneDrive - Delft University of Technology\Rishikesh\PhD\Ampyx\GitLabRepo\ap4_sensitivity\Excel\Ampyx LCOE model.xlsx';

%% Point mass model
% AP3res_PMmodel     = load('AP3results_v3.mat');
% AP4res_PMmodel     = load('AP4results.mat');
AP3res_PMmodel     = load('C:\PhD\Ampyx\Datasets\PM-model-results\AP3_S12_S202\results\results.mat');
AP4res_PMmodel     = load('C:\PhD\Ampyx\Datasets\PM-model-results\AP3_Scaled_S85_SEP2013\results\results.mat');

AP3CycleRes_PMmodel = PM_postProcess(AP3res_PMmodel);
AP4CycleRes_PMmodel = PM_postProcess(AP4res_PMmodel);

%% Analytical model

% To match rated power 
inputParamAP3.gen  = 172; %[kW]
inputParamAP4.gen  = 1072; %[kW]

inputParamAP3.name = 'new';
inputParamAP4.name = 'new';

% Running the analytical model
[AP3mechPower_AnModel, inputParamAP3] = run_AnalyticalModel(inputParamAP3, ExcelModelName, AP3res_PMmodel, AP3CycleRes_PMmodel);
[AP4mechPower_AnModel, inputParamAP4] = run_AnalyticalModel(inputParamAP4, ExcelModelName, AP4res_PMmodel, AP4CycleRes_PMmodel);

%% Plots: Power curve comparisons

% AP4
figure()
hold on
grid on
plot(AP4CycleRes_PMmodel.windSpeedRange, AP4CycleRes_PMmodel.mechPowerAvg,'o--','MarkerSize',4)
plot(AP4CycleRes_PMmodel.windSpeedRange, AP4mechPower_AnModel,'linewidth',1.5)
xlabel('Wind speeds at pattern altitude (m/s)');
ylabel('Power (MW)');
legend('Dynamic (point-mass) model', 'Analytical model','location','North west');
title('Mechanical power (capped): AP4');
hold off  

% AP3
figure()
hold on
grid on
plot(AP3CycleRes_PMmodel.windSpeedRange, AP3CycleRes_PMmodel.mechPowerAvg,'o--','MarkerSize',4)
plot(AP3CycleRes_PMmodel.windSpeedRange, AP3mechPower_AnModel,'linewidth',1.5)
xlabel('Wind speeds at pattern altitude (m/s)');
ylabel('Power (MW)');
legend('Dynamic (point-mass) model', 'Analytical model','location','North west');
title('Mechanical power (capped): AP3');
hold off  

%% Running analytical model for chosen input parameters
function [mechPower, inputParam] = run_AnalyticalModel(inputParam, ExcelModelName, PMmodelResults, PMres_cycle)
  
  % Aliging inputs which vary per wind speed
  for i =1:numel(PMmodelResults.concept.modelInterfaces.performanceModel.powerCurve.data)
    allWS.CLwing(i)     = mean(PMmodelResults.concept.modelInterfaces.performanceModel.powerCurve.data(i).results.CL);
    allWS.CDwing(i)     = mean(PMmodelResults.concept.modelInterfaces.performanceModel.powerCurve.data(i).results.CD);
    allWS.pattEle(i)    = rad2deg(mean(PMmodelResults.concept.modelInterfaces.performanceModel.powerCurve.data(i).results.theta)); 
    allWS.pattRad(i)    = rad2deg(mean(abs(PMmodelResults.concept.modelInterfaces.performanceModel.powerCurve.data(i).results.phi)));
    allWS.maxRollAng(i) = rad2deg(max(PMmodelResults.concept.modelInterfaces.performanceModel.powerCurve.data(i).results.roll)); 
  end
  inputParam.pattEle    = 90 - allWS.pattEle(PMres_cycle.ratedWindSpeed-2); %[deg] % In Matlab model it is measured from the vertical % At rated
  inputParam.pattRad    = allWS.pattRad(PMres_cycle.ratedWindSpeed-2); %[deg] % At rated
  inputParam.maxRollAng = allWS.maxRollAng(PMres_cycle.ratedWindSpeed-2); %[deg] % At rated
  inputParam.CLwing     = mean(allWS.CLwing);
  inputParam.CDwing     = mean(allWS.CDwing);
  
 if inputParam.name == 'old'
  inputParam.mass       = PMmodelResults.concept.Parameters.value.RPA.mass.value; %[kg]
 else
  inputParam.mass       = PMmodelResults.concept.Results.value.RPA.mass.value; %[kg]
 end
  inputParam.TT         = PMmodelResults.concept.modelInputs.Facility.value.RPA.nominalTension.value*0.96; % [kN] % Max TT is never reached
  inputParam.WA         = PMmodelResults.concept.modelInputs.Facility.value.wing.area.value; % [m^2]
  inputParam.AR         = PMmodelResults.concept.modelInputs.Facility.value.wing.aspectRatio.value;  %[-]
  inputParam.numPattCyc = PMmodelResults.concept.modelInterfaces.performanceModel.powerCurve.data(PMres_cycle.ratedWindSpeed-2).results.numCycles  ;  %[-] % At rated
   
  % Writing input param to the Excel
  PowerRatio = inputParam.gen/inputParam.WA;
  TmaxpRatio = inputParam.TT/inputParam.gen;
  writematrix(10,ExcelModelName,'Sheet','Lists','Range','E15');
  writematrix(inputParam.WA,ExcelModelName,'Sheet','Scenarios','Range','R20');
  writematrix(PowerRatio,ExcelModelName,'Sheet','Scenarios','Range','R22');
  writematrix(TmaxpRatio,ExcelModelName,'Sheet','Scenarios','Range','R24');
  writematrix(inputParam.mass,ExcelModelName,'Sheet','Scenarios','Range','R64');
  writematrix(inputParam.CLwing,ExcelModelName,'Sheet','Scenarios','Range','R69');
  writematrix(inputParam.CDwing,ExcelModelName,'Sheet','Performance','Range','C46');
  writematrix(inputParam.pattEle,ExcelModelName,'Sheet','Scenarios','Range','R77');
  writematrix(inputParam.pattRad,ExcelModelName,'Sheet','Scenarios','Range','R86');
  writematrix(inputParam.AR,ExcelModelName,'Sheet','Scenarios','Range','R94');
  writematrix(inputParam.maxRollAng,ExcelModelName,'Sheet','Scenarios','Range','R95');
  writematrix(inputParam.numPattCyc,ExcelModelName,'Sheet','Scenarios','Range','R98');
  
  % Reading Mechanical power curve
  mechPower = xlsread(ExcelModelName,'Performance','AF81:AF106')./10e5;   % [MW]

end

%% Post processing PM model results file
function [cycleData] = PM_postProcess(DataFile)
  
  % User defined values for simplifying calculations
  paramForCalc                 = struct(); %
  paramForCalc.windSpeedsCount = 26; % [0m/s to 25m/s]
  
  temp1 = DataFile.concept.modelInterfaces.performanceModel.powerCurve.data;
  
  for i=1:numel(temp1)
     tf = isreal(temp1(i).results.winch.mechanicalPower);
     if tf == 1 
      paramForCalc.cutIn = i+2;
      break
     end
  end

  % Instantaneous data 
  dataInput               = struct();
  for i=1:paramForCalc.windSpeedsCount  % Wind speed index: 0m/s to 25m/s
    if i<=paramForCalc.cutIn % P-M model calculates only from 3m/s wind speed and Complex values till cut-in are ignored
      dataInput(i).simTimeSteps  = length(DataFile.concept.modelInterfaces.performanceModel.powerCurve.data(4).results.winch.mechanicalPower); % [s]
      dataInput(i).instMechPower = zeros(dataInput(i).simTimeSteps,1); % [W]
      dataInput(i).instTime      = DataFile.concept.modelInterfaces.performanceModel.powerCurve.data(4).results.winch.t; % [s]
    else
      dataInput(i).simTimeSteps  = length(DataFile.concept.modelInterfaces.performanceModel.powerCurve.data(i-3).results.winch.mechanicalPower); % [s]
      dataInput(i).instMechPower = DataFile.concept.modelInterfaces.performanceModel.powerCurve.data(i-3).results.winch.mechanicalPower/10^6; % [MW]
      dataInput(i).instTime      = DataFile.concept.modelInterfaces.performanceModel.powerCurve.data(i-3).results.winch.t; % [s]
    end
  end

  % Cycle average mechanical energy and power per wind speed
  cycleData               = struct(); 
  cycleData.mechPowerAvg  = zeros(paramForCalc.windSpeedsCount,1);
  cycleData.mechEnergyAvg = zeros(paramForCalc.windSpeedsCount,1);
  cycleData.time          = zeros(paramForCalc.windSpeedsCount,1);

  for i=1:paramForCalc.windSpeedsCount 
    if i<=paramForCalc.cutIn
      cycleData.mechEnergyAvg(i) = 0;
    else
      for j=1:dataInput(i).simTimeSteps-1
        cycleData.mechEnergyAvg(i)   = cycleData.mechEnergyAvg(i) + dataInput(i).instMechPower(j)*(dataInput(i).instTime(j+1)-dataInput(i).instTime(j));
      end
      cycleData.time(i)          = dataInput(i).instTime(end);
      cycleData.mechPowerAvg(i)  = cycleData.mechEnergyAvg(i)/cycleData.time(i);
    end
  end
  
    % Extracting rated wind speed
    temp(1)   = cycleData.mechPowerAvg(1);
  for i = 2:numel(cycleData.mechPowerAvg)-1
    temp(i) = cycleData.mechPowerAvg(i+1) - cycleData.mechPowerAvg(i);
  end
  temp = temp./max(temp);
  for i = paramForCalc.cutIn:numel(cycleData.mechPowerAvg)-1
    if temp(i) < 0.04
      paramForCalc.ratedWindSpeed = i-1;
      break
    else
      paramForCalc.ratedWindSpeed = i;
    end
  end
  
  cycleData.windSpeedsCount = paramForCalc.windSpeedsCount;
  cycleData.cutIn           = paramForCalc.cutIn;
  cycleData.ratedWindSpeed  = paramForCalc.ratedWindSpeed;
  cycleData.windSpeedRange = linspace(0,25,paramForCalc.windSpeedsCount);

end
