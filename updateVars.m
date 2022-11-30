function [inputs] = updateVars(x,i,inputs)
  
    global outputs
  
    outputs.deltaL(i)     = x(1);
    outputs.VRI(i)        = x(2);
    inputs.Tmax           = x(3);
    outputs.CL_airfoil(i) = x(4);

end