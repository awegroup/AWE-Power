function [c, ceq] = constraints(i,inputs)
  
  global outputs
   
  %% Inequality constraints
  
  % Min clearance between ground and bottom point of pattern
  c(1)   = (inputs.minGroundClear - outputs.pattStartGrClr(i))/inputs.minGroundClear/100;
  c(2)   = (inputs.minGroundClear - outputs.pattEndGrClr(i))/inputs.minGroundClear/100;
  
  % Capping for requested electrical rated power
  c(3)   = (outputs.P_cycleElec(i) - inputs.P_ratedElec)/inputs.P_ratedElec/1000; 
  
  % Tether length limit
  c(4) = (outputs.L_teMax(i) - inputs.maxTeLen)/inputs.maxTeLen/100; 
  
  % Min number of patterns to get into transition 
  c(5) = (1 - mean(outputs.numOfPatt(i,:)))/10;
  
  % Max. cycle avg height
  c(6) = (outputs.H_cycleEnd(i) - inputs.maxHeight)/inputs.maxHeight/1000;
  
  % Peak mechanical power limit
  c(1,7:inputs.numDeltaLelems+6) = (outputs.PROeff_mech(i,:) - inputs.peakM2E_F*inputs.P_ratedElec)/(inputs.peakM2E_F*inputs.P_ratedElec*1000);
  
  % Maximum tether force
  c(1,7+inputs.numDeltaLelems:2*inputs.numDeltaLelems+6) = (outputs.Ft(i,:) - inputs.Ft_max*inputs.Ft_max_SF*1000)/inputs.Ft_max*inputs.Ft_max_SF*100000;
  
  % Kite speed limit
%   c(1,7+2*inputs.numDeltaLelems:3*inputs.numDeltaLelems+6) = (outputs.vk_omega(i,:) - 40)/1000;
    
  %% Equality constraints
  
  % Kinematic ratio
  ceq(1:inputs.numDeltaLelems) = (outputs.k_result(i,:) - (outputs.CL(i,:)/outputs.CD(i,:)))/100;

 
end