function [c, ceq] = constraints(i,inputs)
  
  global outputs
   
  %% Inequality constraints
  
  % Min clearance between ground and bottom point of pattern
  c(1)   = (inputs.minGroundClear - outputs.pattStartGrClr(i))/inputs.minGroundClear;
  
  c(2)   = (inputs.minGroundClear - outputs.pattEndGrClr(i))/inputs.minGroundClear;
  
  % Capping for requested electrical rated power
  c(3)   = (outputs.P_cycleElec(i) - inputs.P_ratedElec)/inputs.P_ratedElec; 
  
  % Tether length limit
  c(4) = (outputs.L_teMax(i) - inputs.maxTeLen)/inputs.maxTeLen; 
  
  % Min number of patterns to get into transition 
  c(5) = (1 - mean(outputs.numOfPatt(i,:)))/3;
  
  % Max. cycle avg height
  c(6) = (outputs.H_cycleEnd(i) - inputs.maxHeight)/inputs.maxHeight;
  
  % Peak mech to cycle elec ratio
 % c(1,7:outputs.deltaLelems+6) = (outputs.PROeff_mech(i,:) - inputs.F_peakM2Ecyc*inputs.P_ratedElec)/(inputs.F_peakM2Ecyc*inputs.P_ratedElec);
  
  %% Equality constraints
  
  if inputs.targetPRO_mech ~= 0 % Only run in Second Optimisation to follow capped mean mech. power
      ceq(1)    = (mean(outputs.PROeff_mech(i,:)) - inputs.targetPRO_mech(i))/inputs.targetPRO_mech(i)/10;
  else
      ceq(1)    = 0; % For First Optimisation
  end
 
end