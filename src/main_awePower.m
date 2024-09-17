function [inputs, outputs, optimDetails, processedOutputs] = main_awePower(inputs)
  % main_awePower: Optimizes the power generation for a fixed-wing airborne wind energy system
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
  %
  % The function also handles the kite mass estimation, depending on whether a mass override is
  % specified. After performing the optimization, it post-processes the outputs to represent cycle
  % power for each wind speed within the operational range and saves the results.


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