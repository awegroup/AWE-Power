function [c, ceq] = constraints(i,inputs)
  
  global outputs
   
  c(1)   = inputs.minGroundClear - outputs.H_avg(i) + outputs.pattRadius(i); % Min clearance between ground and bottom point of pattern
  c(3)   = max(outputs.PROeff_elec_osci(i,:)) - inputs.F_peakElecP*inputs.P_ratedElec; % Capping for max. electrical power
  c(2)   = outputs.P_cycleElec(i) - inputs.P_ratedElec; % Capping for requested electrical rated power
  ceq    = 0;  
 
end