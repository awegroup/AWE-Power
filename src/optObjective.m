function [fval,outputs] = optObjective(x,i,inputs,kiteMass)
  % optObjective_awePower Objective function for optimizing power generation in AWE system
  %
  % This function computes the objective value for optimizing the power generation
  % of an Airborne Wind Energy (AWE) system at a specific wind speed index.
  %
  % Inputs:
  %   x        - Vector of decision variables for optimization
  %   i        - Index indicating the current wind speed iteration
  %   inputs   - Structure containing input parameters and constants
  %   kiteMass - Mass of the kite used in the system
  %
  % Outputs:
  %   fval     - Objective function value (negative of average electrical power)
  %   outputs  - Structure containing updated output variables after computation

  global outputs

  % Assign Kite mass
  outputs.m_k = kiteMass;

  % Assign variable values
  outputs.deltaL(i)                               = x(1);
  outputs.beta(i)                                 = x(2);
  outputs.gamma(i)                                = x(3);
  outputs.Rp_start(i)                             = x(4);
  outputs.vk_r_i(i,1:inputs.numDeltaLelems)       = x(0*inputs.numDeltaLelems+5:1*inputs.numDeltaLelems+4);
  outputs.CL_i(i,1:inputs.numDeltaLelems)         = x(1*inputs.numDeltaLelems+5:2*inputs.numDeltaLelems+4);
  outputs.vk_r(i,1:inputs.numDeltaLelems)         = x(2*inputs.numDeltaLelems+5:3*inputs.numDeltaLelems+4);
  outputs.kRatio(i,1:inputs.numDeltaLelems)       = x(3*inputs.numDeltaLelems+5:4*inputs.numDeltaLelems+4);
  outputs.CL(i,1:inputs.numDeltaLelems)           = x(4*inputs.numDeltaLelems+5:5*inputs.numDeltaLelems+4);

  % Main Analysis
  computePower(i,inputs);

  % Objective: Maximise electrical cycle average power
  fval = - outputs.P_e_avg(i);

end                    