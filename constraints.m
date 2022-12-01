function [c, ceq] = constraints(i,inputs)
  
  global outputs
   
  %% Min clearance between ground and bottom point of pattern
  c(1)   = inputs.minGroundClear - outputs.H_avg(i) + outputs.pattRadius(i); 
  
  %% Capping for max. reel-out electrical power
  c(3)   = max(outputs.PROeff_elec_osci(i,:)) - inputs.F_peakElecP*inputs.P_ratedElec; 
  
  %% Capping for requested electrical rated power
  c(2)   = outputs.P_cycleElec(i) - inputs.P_ratedElec; 
  
  %% Required argument
  ceq    = 0;  
end