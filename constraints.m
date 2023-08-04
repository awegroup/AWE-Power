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
  % If not running second optimisation
  c(1,7:inputs.numDeltaLelems+6) = (outputs.PROeff_mech(i,:) - inputs.F_peakM2Ecyc*inputs.P_ratedElec)/(inputs.F_peakM2Ecyc*inputs.P_ratedElec);
  
  % Reel-out speed should be less than tether direction component of wind
  % Top pint
%  c(1,inputs.numDeltaLelems+7:2*inputs.numDeltaLelems+6) = outputs.VRO_top(i,:) - outputs.Vw(i,:)*cos(outputs.avgPattEle(i)+outputs.pattAngRadius(i));
  
  % Bottom pint
%   c(1,inputs.numDeltaLelems+7:2*inputs.numDeltaLelems+6) = outputs.VRO_top(i,:) - outputs.Vw(i,:)*cos(outputs.avgPattEle(i)-outputs.pattAngRadius(i));
  
  % Tangential speed constraint
%    c(1,2*inputs.numDeltaLelems+7:3*inputs.numDeltaLelems+6) = 0 - outputs.Vc_top(i,:);
   
   
  %% Equality constraints
  
  ceq(1) = 0;
  
  % Force balance condition 
  
  % Top point of the pattern
  %ceq(1:inputs.numDeltaLelems) = outputs.Fa_top(i,:).*sin(outputs.rollAngleTop(i,:))-...
  %                    (outputs.W(i).*cos(outputs.avgPattEle(i)+outputs.pattAngRadius(i))+outputs.Fc_top(i,:).*cos(outputs.pattAngRadius(i))); 
  
  %Bottom point of the pattern
 % ceq(1:inputs.numDeltaLelems) = outputs.Fa_top(i,:).*sin(outputs.rollAngleTop(i,:))-...
%                     (outputs.W(i).*cos(outputs.avgPattEle(i)-outputs.pattAngRadius(i))-outputs.Fc_top(i,:).*cos(outputs.pattAngRadius(i))); 
  
                  
                  
                  
%   if inputs.targetPRO_mech ~= 0 % Only run in Second Optimisation to follow capped mean mech. power
%       ceq(1)    = (mean(outputs.PROeff_mech(i,:)) - inputs.targetPRO_mech(i))/inputs.targetPRO_mech(i)/10;
%   else
%       ceq(1)    = 0; % For First Optimisation
%   end
  
  % Kite tangential speed at top pattern point
   ceq(1,2:inputs.numDeltaLelems+1) = (outputs.Va_top(i,:).^2 - outputs.Vc_top(i,:).^2 - (outputs.VSR_top(i,:)).^2)./50.^2; %./cos(outputs.rollAngleTop(i,:))
  % Kite tangential speed at bottom pattern point
%   ceq(1,inputs.numDeltaLelems+2:2*inputs.numDeltaLelems+1) = (outputs.Va_top(i,:).^2 - outputs.Vc_top(i,:).^2 - (outputs.VSR_top(i,:)/outputs.rollAngleTop(i,:)).^2)./50.^2;

end