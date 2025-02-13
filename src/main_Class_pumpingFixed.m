function [inputs, outputs, optimDetails, processedOutputs, simulationQSM] = main_Class_pumpingFixed(inputs)
  % main_Class_pumpingFixed: Optimizes the power generation for a fixed-wing airborne wind energy system
  %
  % This function optimizes and computes the power output of a fixed-wing airborne wind energy
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
  options.FiniteDifferenceType      = 'forward'; % 'forward'
  options.ScaleProblem              = true;
  options.ConstraintTolerance       = 1e-3;
  options.FiniteDifferenceStepSize  = 1e-6;
  options.MaxIterations             = 200;
  
  % Optimisation Initialisation
  nx = ones(1, inputs.numDeltaLelems); % For discretized stroke length

  stroke_length = inputs.stroke_length_m;
  elevationAngle_deg = inputs.elevationAngle_deg;
  coneAngle_deg = inputs.coneAngle_deg;
  turningRadius_over_span_ratio = inputs.turningRadius_over_span_ratio;
  reel_in_speed = inputs.reel_in_speed;
  reel_out_speed = inputs.reel_out_speed;
  kinematicRatio_lambda = inputs.kinematicRatio_lambda;

  % Initial Guess
  inputs.x0 = [stroke_length(2), ...                     % stroke length
      deg2rad(elevationAngle_deg(2)), ...                % pattern elevation
      deg2rad(coneAngle_deg(2)), ...                     % cone angle
      turningRadius_over_span_ratio(2)*inputs.span, ...  % initial turning radius
      reel_in_speed(2) .* nx, ...                        % reel-in speed
      inputs.Cl_maxAirfoil * inputs.Cl_eff_F .* nx, ...  % reel-in C_L
      reel_out_speed(2) .* nx, ...                       % reel-out speed
      kinematicRatio_lambda(2) .* nx, ...                % kinematic ratio
      inputs.Cl_maxAirfoil * inputs.Cl_eff_F .* nx];     % reel-out C_L

  % Bounds for the initial guess
  inputs.lb = [stroke_length(1), ...                     % stroke length
      deg2rad(elevationAngle_deg(1)), ...                % pattern elevation
      deg2rad(coneAngle_deg(1)), ...                     % cone angle
      turningRadius_over_span_ratio(1)*inputs.span, ...  % initial turning radius
      reel_in_speed(1) .* nx, ...                        % reel-in speed
      0.1 .* nx, ...                                     % reel-in C_L
      reel_out_speed(1) .* nx, ...                       % reel-out speed
      kinematicRatio_lambda(1) .* nx, ...                % kinematic ratio
      0.1 .* nx];                                        % reel-out C_L
  
  inputs.ub = [stroke_length(3), ...                     % stroke length
      deg2rad(elevationAngle_deg(3)), ...                % pattern elevation
      deg2rad(coneAngle_deg(3)), ...                     % cone angle
      turningRadius_over_span_ratio(3)*inputs.span, ...  % initial turning radius
      reel_in_speed(3) .* nx, ...                        % reel-in speed
      inputs.Cl_maxAirfoil * inputs.Cl_eff_F .* nx, ...  % reel-in C_L
      reel_out_speed(3) .* nx, ...                       % reel-out speed
      kinematicRatio_lambda(3) .* nx, ...                % kinematic ratio
      inputs.Cl_maxAirfoil * inputs.Cl_eff_F .* nx];     % reel-out C_L

  %% Instantiate class
  simulationQSM = KiteQSMsimulation_pumpingFixed(inputs, loggerID);

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

  % Cycle power representation for wind speeds in the operational range
  for i = processedOutputs.cutIn:processedOutputs.cutOut
    [processedOutputs.cyclePowerRep(i)] = createCyclePowerRep(i, processedOutputs, processedOutputs.vw_reference);
  end

  %% Save outputs  
  saveResults(inputs, optimDetails, outputs, processedOutputs);

end