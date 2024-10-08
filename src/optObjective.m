function [fval,outputs] = optObjective(x,i,inputs,kiteMass)
  % optObjective_awePower: Objective function for optimizing power generation in AWE system
  %
  % This function computes the objective value for maximizing the power generation of a
  % fixed-wing airborne wind energy (AWE) system at a specific wind speed index by optimizing
  % decision variables such as tether length, angle of attack, and glide ratio. The goal is
  % to maximize the electrical power output by minimizing the negative of the cycle average
  % electrical power.
  %
  % Inputs:
  %   x        (vector): Vector of decision variables for optimization, which include tether
  %                      length ratios, velocity components, and lift coefficients.
  %   i        (int): Index indicating the current wind speed iteration.
  %   inputs   (struct): Structure containing input parameters and system constants such as
  %                      the number of decision variables and wind profiles.
  %   kiteMass (double): Mass of the kite used in the AWE system, either provided or calculated.
  %
  % Outputs:
  %   fval     (double): Objective function value, which is the negative of the average electrical
  %                      power for a given wind speed and set of decision variables.
  %   outputs  (struct): Updated structure containing output variables such as tether forces,
  %                      power values, and velocities after the power computation step.
  %
  % The objective function aims to find the set of decision variables that maximize the
  % system's power output while ensuring operational constraints such as tether forces
  % and velocities are met.


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