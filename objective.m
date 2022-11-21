function [fval,inputs,outputs] = objective(x,x_init,i,inputs)

    global outputs
  
    % Denormalize
    x_input = x.*x_init; 

    [inputs]  = updateVars(x_input,i,inputs);
    [inputs]  = compute(i,inputs);
    
    fval = -outputs.P_cycleElec(i);
    
end