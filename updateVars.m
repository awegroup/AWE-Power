function [inputs] = updateVars(x,i,inputs)
  
    global outputs
  
    outputs.deltaL(i)        = x(1);
    outputs.VRI(i)           = x(2);
    outputs.CL(i)            = x(3);
    outputs.avgPattEle(i)    = x(4);
    outputs.pattAngRadius(i) = x(5);
    outputs.maxRollAngle(i)  = x(6);

end