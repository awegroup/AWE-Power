function [c, ceq] = constraints(i,inputs)
  
  global outputs
   
  %% Inequality constraints
  
  % Min clearance between ground and bottom point of pattern
  c(1)   = inputs.minGroundClear - outputs.L_teMin(i)*sin(outputs.avgPattEle(i)-outputs.pattAngRadius(i)); 
  
  % Capping for requested electrical rated power
  c(2)   = outputs.P_cycleElec(i) - inputs.P_ratedElec; 
  
  % Unphysical turn radius
  c(3)   = outputs.wingSpan/2 - outputs.pattRadius(i);
  
  % Avg. patt. elevation > patt. angular radius
  c(4)   = outputs.pattAngRadius(i) - outputs.avgPattEle(i);
  
  % Tether length limit
  c(5) = outputs.L_teMax(i) - inputs.maxTeLen; 
  
  %% Equality constraints
  
  % Consistency in roll angles due to centrifugal and gravity effects
  ceq(1)    = outputs.rollAngleC(i) - outputs.rollAngleG(i);
  
%   % Only used in second iter of optimisation when capping max. electrical power
%   if inputs.targetPRO_elec ~=0
%       % Capping for max. reel-out electrical power
%       ceq(2)    = outputs.PROeff_elec(i) - inputs.targetPRO_elec(i);
%   end
 
end