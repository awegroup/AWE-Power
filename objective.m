function [fval,inputs,outputs] = objective(x,x_init,i,inputs)

    global outputs
  
    % Denormalize x
    x = x.*x_init;
    
    % Assign variable values
    outputs.deltaL(i)                               = x(1);
    outputs.vk_r_i(i)                               = x(2);
    outputs.beta(i)                                 = x(3);
    outputs.gamma(i)                                = x(4);
    outputs.Rp_start(i)                             = x(5);
    outputs.CL(i,1:inputs.numDeltaLelems)           = x(6:inputs.numDeltaLelems+6-1);
    outputs.f(i,1:inputs.numDeltaLelems)            = x(1*inputs.numDeltaLelems+6:2*inputs.numDeltaLelems+6-1); 
    outputs.kRatio(i,1:inputs.numDeltaLelems)       = x(2*inputs.numDeltaLelems+6:3*inputs.numDeltaLelems+6-1); 
    
    % Main computation
    [inputs]  = compute(i,inputs);
    
    % Objective
    fval = - outputs.P_cycleElec(i)/inputs.P_ratedElec/1000;
        
end                    