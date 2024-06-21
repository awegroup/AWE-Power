function [c, ceq] = optConstraints_awePower(i,inputs)
  % optConstraints_awePower Nonlinear constraints for optimizing power generation in AWE system
  %
  % This function computes the nonlinear inequality and equality constraints
  % for optimizing the power generation of an Airborne Wind Energy (AWE) system
  % at a specific wind speed index.
  %
  % Inputs:
  %   i       - Index indicating the current wind speed iteration
  %   inputs  - Structure containing input parameters and constants
  %
  % Outputs:
  %   c       - Vector of nonlinear inequality constraint values
  %   ceq     - Vector of nonlinear equality constraint values

  global outputs

  %% Inequality constraints

  % Min clearance between ground and bottom point of pattern
  c(1)   = (inputs.minGroundClear - outputs.pattStartGrClr(i));
  c(2)   = (inputs.minGroundClear - outputs.pattEndGrClr(i));

  % Capping for requested electrical rated power
  c(3)   = (outputs.P_e_avg(i) - inputs.P_ratedElec);

  % Tether length limit
  c(4) = (outputs.l_t_max(i) - inputs.maxTeLen);

  % Min number of patterns to get into transition
  c(5) = (1 - mean(outputs.numOfPatt(i,:)));

  % Max. cycle avg height
  c(6) = (outputs.h_cycleEnd(i) - inputs.maxHeight);

  % Peak mechanical power limit
  c(1,7:inputs.numDeltaLelems+6) = (outputs.P_m_o_eff(i,:) - inputs.peakM2E_F*inputs.P_ratedElec);

  % % Maximum tether force
  c(1,7+inputs.numDeltaLelems:2*inputs.numDeltaLelems+6) = (outputs.Ft(i,:) - inputs.Ft_max*inputs.Ft_max_SF);

  % Apparent speed cannot be negative
  c(1,7+2*inputs.numDeltaLelems:3*inputs.numDeltaLelems+6) = (0 - outputs.va(i,:));

  % Tangential velocity factor cannot be negative
  c(1,7+3*inputs.numDeltaLelems:4*inputs.numDeltaLelems+6) = (0 - outputs.lambda(i,:));

  % Tether force during reel-in cannot be negative
  c(1,7+4*inputs.numDeltaLelems:5*inputs.numDeltaLelems+6) = (0 - outputs.Ft_i(i,:));


  %% Equality constraints

  % Glide ratio during reel-out
  ceq(1,0*inputs.numDeltaLelems+1:1*inputs.numDeltaLelems) = (outputs.E_result(i,:) - outputs.E(i,:));

  % Glide ratio during reel-in
  ceq(1,1*inputs.numDeltaLelems+1:2*inputs.numDeltaLelems) = (outputs.E_result_i(i,:) - outputs.E_i(i,:));

end