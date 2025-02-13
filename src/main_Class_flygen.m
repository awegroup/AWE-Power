function [inputs, outputs, optimDetails, processedOutputs, simulationQSM] = main_Class_flygen(inputs)
  % main_Class_flygen: Optimizes the power generation for a fixed-wing airborne wind energy system
  %
  % This function optimizes and computes the power output of a fly-gen fixed-wing airborne wind energy
  % system across a range of wind speeds. It involves mass estimation, system optimization using
  % the `fmincon` optimizer, and post-processing of outputs. The function uses input parameters
  % and constants provided in the `inputs` structure, which are validated and appended before
  % the optimization process.
  %
  % Inputs:
  %   inputs (struct): Structure containing necessary input parameters and constants for the
  %                    airborne wind energy system, including wing characteristics, tether force,
  %                    and aerodynamic parameters.
  %
  % Outputs:
  %   outputs          (struct): Structure containing the results of the optimization process
  %                              for various wind speeds.
  %   optimDetails     (struct): Structure containing details and intermediate values from the
  %                              optimization process, including convergence information.
  %   processedOutputs (struct): Structure containing post-processed results such as power
  %                              output, efficiency, and additional operational details.
  %   simulationQSM    (object): Class object of the QSM simulation.
  %
  % The function also handles the kite mass estimation, depending on whether a mass override is
  % specified. After performing the optimization, it post-processes the outputs to represent cycle
  % power for each wind speed within the operational range and saves the results.

  %% Set additional input parameters
  loggerID = 'QSMlogger';
  logger = mlog.Logger(loggerID,'QSMlogger.log');
  logger.clearLog();
  logger.CommandWindowThreshold = "INFO";
  
  options                           = optimoptions('fmincon');
  options.Display                   = 'none'; %'iter-detailed', 'notify-detailed'
  options.Algorithm                 = 'sqp';
  options.FiniteDifferenceType      = 'central'; % 'forward'
  options.ScaleProblem              = true;
  options.FiniteDifferenceStepSize  = 1e-3;
  options.MaxFunctionEvaluations    = 15000;
  options.MaxIterations             = 10000;
  options.FunctionTolerance         = 1e-3;
  
  % Optimisation Initialisation
  elevationAngle_deg = inputs.elevationAngle_deg;
  coneAngle_deg = inputs.coneAngle_deg;
  tetherLength_m = inputs.tetherLength_m;
  kinematicRatio_lambda = inputs.kinematicRatio_lambda;
  thrustToKiteDragRatio = inputs.thrustToKiteDragRatio;

  % Initial Guess
    inputs.x0 = [deg2rad(elevationAngle_deg(2)), ...        % pattern elevation
                deg2rad(coneAngle_deg(2)), ...              % cone angle
                tetherLength_m(2), ...                      % tether length
                kinematicRatio_lambda(2), ...               % kinematicRatio
                inputs.Cl_maxAirfoil*inputs.Cl_eff_F, ...   % CL
                thrustToKiteDragRatio(2)];                  % CT/CDkite

    inputs.lb = [deg2rad(elevationAngle_deg(1)), ...
                deg2rad(coneAngle_deg(1)), ...
                tetherLength_m(1), ...
                kinematicRatio_lambda(1), ...
                0.1,...
                thrustToKiteDragRatio(1)];
            
    inputs.ub = [deg2rad(elevationAngle_deg(3)), ...
                deg2rad(coneAngle_deg(3)), ...
                tetherLength_m(3), ...
                kinematicRatio_lambda(3), ...
                inputs.Cl_maxAirfoil*inputs.Cl_eff_F,...
                thrustToKiteDragRatio(3)];

  %% Instantiate class
  simulationQSM = KiteQSMsimulation_flygen(inputs, loggerID);

  %% Run simulation
  % Accessing class properties here once and let matlab make a deep-copy is
  % faster than accessing each variable separately inside the class method.
  simulationQSM = simulationQSM.runsimulation(simulationQSM.inputs, options);

  % Post processing
  simulationQSM = simulationQSM.processoutputs();
    
  % Extract results from class object
  inputs = simulationQSM.inputs;
  outputs = simulationQSM.outputs;
  optimDetails = simulationQSM.optimDetails;
  processedOutputs = simulationQSM.processedOutputs;

  %% Save outputs  
  saveResults(inputs, optimDetails, outputs, processedOutputs);

end