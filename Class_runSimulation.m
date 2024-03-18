% Run the below function to generate results
% Make sure to clean the variables from the workspace
clc
clearvars
%% Required
inputs.windSpeedReference  = 1:1:25; %[m/s]
inputs.heightWindReference = 100; %[m]
inputs.areaWing            = 12; %[m^2]
inputs.span                = 12; %[m]
inputs.forceTether_max     = 42*1000; %[N]
inputs.P_ratedElec         = 150*1000; %[W]
inputs.Cl_maxAirfoil       = 2.5; %[-]
inputs.Cl_eff_F            = 0.8; %[-]
inputs.Cl0_airfoil         = 0.65; %[-]
inputs.Cd0                 = 0.056; %[-]

%% Optional (has default value built in)
inputs.doMassOverride   = true;
inputs.kiteMass         = 600; %[kg]
inputs.numDeltaLelems   = 2; %[num]

%% Create class
simulationQSM = KiteQSMsimulation2(inputs);

options                           = optimoptions('fmincon');
options.Display                   = 'none'; %'iter-detailed', 'notify-detailed'
options.Algorithm                 = 'sqp';
options.FiniteDifferenceType      = 'central'; % 'forward'
options.ScaleProblem              = true;
% options.FiniteDifferenceStepSize  = [1e-12 1e-12 1e-12 1e-12 1e-6*nx 1e-6*nx 1e-6*nx 1e-6*nx];% 1e-6*nx];
% options.FiniteDifferenceStepSize  = 1e-9;
% options.ConstraintTolerance       = 1e-4;
% options.OptimalityTolerance       = 1e-9;
% options.StepTolerance             = 1e-4;
% options.MaxFunctionEvaluations    = 5000*numel(x0);
% options.MaxIterations             = 500*numel(x0);
% options.FunctionTolerance        = 1e-9;
% options.DiffMaxChange            = 1e-1;
% options.DiffMinChange            = 0;

%% Run QSM
simulationQSM = simulationQSM.runsimulation(simulationQSM.inputs, options);

%% Post processing
simulationQSM = simulationQSM.processoutputs(simulationQSM.inputs, simulationQSM.outputs);

%% Plotting
simulationQSM.plotresults(simulationQSM.inputs, simulationQSM.processedOutputs);