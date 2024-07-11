function [inputs, outputs, optimDetails, processedOutputs] = main_awePower(inputs)
  % main_awePower Optimizes the power generation for a fixed-wing airborne wind energy system
  %
  % This function performs mass estimation, optimization, and post-processing
  % to compute the power output of a fixed-wing airborne wind energy system across
  % a range of wind speeds. It utilizes the fmincon optimizer to find optimal
  % operation parameters for each wind speed reference value.
  %
  % Inputs:
  %   inputs - Structure containing necessary input parameters and constants
  %
  % Outputs:
  %   outputs          - Structure containing computed values from the optimization
  %   optimDetails     - Structure containing details of the optimization process
  %   processedOutputs - Structure containing post-processed outputs

  clear global outputs

  % Validate input data
  validateInput(inputs, inputValidators());

  % Calculate additional input parameters
  inputs = appendInputs(inputs);

  % Kite mass estimate
  if inputs.massOverride == 1
    kiteMass = inputs.kiteMass;
  else
    kiteMass = estimateKiteMass(inputs.Ft_max, inputs.S, inputs.AR);
  end

  % Optimisation
  [optimDetails, outputs] = optProblemFormulation(inputs, kiteMass);

  % Post processing
  [processedOutputs] = postProcessOutputs(inputs, outputs);

  % Cycle power representation for wind speeds in the operational range
  for i = processedOutputs.cutIn:processedOutputs.cutOut
    [processedOutputs.cyclePowerRep(i)] = createCyclePowerRep(i, processedOutputs, inputs.vw_ref);
  end

  % Save outputs
  saveResults(inputs, optimDetails, outputs, processedOutputs);

  clear global outputs

end