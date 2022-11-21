function [c, ceq] = constraints(i,inputs)
  
  global outputs
   
  c(1)   = - outputs.H_avg(i) + outputs.pattRadius(i) + inputs.minGroundClear; % Min clearance between ground and bottom point of pattern
  c(2)   = outputs.P_cycleElec(i) - inputs.P_ratedElec; % Capping for requested electrical rated power
  ceq    = 0;  
 
end