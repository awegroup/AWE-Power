function [fval,inputs,outputs] = handleFn_PM(x,i,inputs,outputs)
       
    [inputs,outputs] = updateVars_PM(x,i,inputs,outputs);
    [inputs,outputs] = evaluation_PM(i,inputs,outputs);
    %[inputs,outputs] = Copy_of_evaluation_PM(i,inputs,outputs);

    fval = -outputs.P_cycleElec(i);
    
end