function [fval,inputs,outputs] = objective_GA(i,inputs, outputs)

%     global outputs
  
    

     % Assign variable values
    outputs.deltaL(i)        = x(1);
    outputs.VRI(i)           = x(2);
    outputs.CL(i)            = x(3);
    outputs.avgPattEle(i)    = x(4);
    outputs.pattAngRadius(i) = x(5);
    outputs.maxRollAngle(i)  = x(6);

    % Main computation
    [inputs, outputs]  = compute(i,inputs, outputs);
    
    % Objective
    fval = -outputs.P_cycleElec(i);
    
end