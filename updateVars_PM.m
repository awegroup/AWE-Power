function [inputs,outputs] = updateVars_PM(x,i,inputs,outputs)
    
    outputs.deltaL(i) = x(1);
    outputs.VRI(i)    = x(2);
    inputs.Tmax       = x(3);

end