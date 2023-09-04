function [c, ceq] = constraints(i,inputs)
  
  global outputs
   
  %% Inequality constraints
  % c = 0;
  % Min clearance between ground and bottom point of pattern
  c(1)   = (inputs.minGroundClear - outputs.pattStartGrClr(i))/1000;
  c(2)   = (inputs.minGroundClear - outputs.pattEndGrClr(i))/1000;

  % Capping for requested electrical rated power
  c(3)   = (outputs.P_cycleElec(i) - inputs.P_ratedElec)/inputs.P_ratedElec/1000; 

  % Tether length limit
  c(4) = (outputs.l_t_max(i) - inputs.maxTeLen)/inputs.maxTeLen/100; 

  % Min number of patterns to get into transition 
  c(5) = (1 - mean(outputs.numOfPatt(i,:)))/10;

  % Max. cycle avg height
  c(6) = (outputs.h_cycleEnd(i) - inputs.maxHeight)/inputs.maxHeight/1000;

  % Peak mechanical power limit
  c(1,7:inputs.numDeltaLelems+6) = (outputs.PROeff_mech(i,:) - inputs.peakM2E_F*inputs.P_ratedElec)/(inputs.peakM2E_F*inputs.P_ratedElec*1000);

  % Maximum tether force
  c(1,7+inputs.numDeltaLelems:2*inputs.numDeltaLelems+6) = (outputs.Ft(i,:) - inputs.Ft_max*inputs.Ft_max_SF*1000)/inputs.Ft_max*inputs.Ft_max_SF*100000;

  % Apparent speed cannot be negative
  c(1,7+2*inputs.numDeltaLelems:3*inputs.numDeltaLelems+6) = (0 - outputs.va(i,:))/1000;

  % Tangential velocity factor cannot be negative
  c(1,7+3*inputs.numDeltaLelems:4*inputs.numDeltaLelems+6) = (0 - outputs.lambda(i,:))/1000;

  % Tether force during reel-in 
  c(1,7+4*inputs.numDeltaLelems:5*inputs.numDeltaLelems+6) = (0 - outputs.Ft_drum_i(i,:))/10000;
  
% %   Kite speed limit
%   c(1,7+3*inputs.numDeltaLelems:4*inputs.numDeltaLelems+6) = (outputs.vk_omega(i,:) - 40)/1000;

    
  %% Equality constraints
  
  % Kinematic ratio
  ceq(1,0*inputs.numDeltaLelems+1:1*inputs.numDeltaLelems) = (outputs.G_result(i,:) - (outputs.CL(i,:)/outputs.CD(i,:)))/100;

  ceq(1,1*inputs.numDeltaLelems+1:2*inputs.numDeltaLelems) = (outputs.G_result_i(i,:) - (outputs.CL_i(i,:)/outputs.CD_i(i,:)))/100;
 

 
end