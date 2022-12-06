function [c, ceq] = constraints_GA(i,inputs,outputs)
  
%   global outputs
   
%   %% Min clearance between ground and bottom point of pattern
%   c(1)   = inputs.minGroundClear - outputs.H_avg(i) + outputs.pattRadius(i); 
%   
%   %% Capping for max. reel-out electrical power
%   c(2)   = max(outputs.PROeff_elec_osci(i,:)) - inputs.F_peakElecP*inputs.P_ratedElec; 
%   
%   %% Capping for requested electrical rated power
%   c(3)   = outputs.P_cycleElec(i) - inputs.P_ratedElec; 
%   
%   %% Min required radius to turn
%   c(4)   = outputs.wingSpan/2 - outputs.pattRadius(i);
  
%%  
% Assign variable values
%     outputs.deltaL(i)        = x(1);
%     outputs.VRI(i)           = x(2);
%     outputs.CL(i)            = x(3);
%     outputs.avgPattEle(i)    = x(4);
%     outputs.pattAngRadius(i) = x(5);
%     outputs.maxRollAngle(i)  = x(6);

    % Main computation
%     compute(i,inputs, outputs);
  %%
  c = [inputs.minGroundClear - outputs.H_avg(i) + outputs.pattRadius(i); max(outputs.PROeff_elec_osci(i,:)) - inputs.F_peakElecP*inputs.P_ratedElec;...
        outputs.P_cycleElec(i) - inputs.P_ratedElec; outputs.wingSpan/2 - outputs.pattRadius(i)];
  
  %% Required argument
  ceq    = outputs.rollAngleC(i) - outputs.rollAngleG(i);
 
end