function [inputs] = updateVars(x,i,inputs)
  
    global outputs
  
    outputs.deltaL(i)     = x(1);
    outputs.VRI(i)        = x(2);
    outputs.CL(i)         = x(3);

end