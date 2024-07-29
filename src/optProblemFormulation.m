function [optimDetails, outputs] = optProblemFormulation(inputs, kiteMass)

  % Initial guess for optimisation
  x0     = inputs.x0;

  % Bounds
  lb = inputs.lb;
  ub = inputs.ub;

  % Optimise operation for every wind speed
  for i=1:length(inputs.vw_ref)

    options                           = optimoptions('fmincon');
    options.Display                   = 'none';
    options.Algorithm                 = 'sqp';
    options.FiniteDifferenceType      = 'central';
    %   options.FiniteDifferenceType      = 'forward';
    options.ScaleProblem              = true;
    %           options.MaxFunctionEvaluations    = 5000*numel(x0);
    %           options.MaxIterations             = 500*numel(x0);


    con = @(x) optConstraints(i,inputs);
    obj = @(x) optObjective(x,i,inputs,kiteMass);


    [x,~,exitflag(i),optHist(i),lambda(i)] = fmincon(obj,x0,[],[],[],[],lb,ub,con,options);

    % Storing final results
    [~,outputs] = optObjective(x,i,inputs,kiteMass);

    % Changing initial guess if previous wind speed evaluation is considered infeasible
    if abs(mean((outputs.E_result - outputs.E),2)) > 1 %0.01
      x0 = inputs.x0;
    else
      % Solution of current wind speed as the initial guess for the next wind speed
      x0 = x;
    end

  end
%  disp(exitflag)
  % Store optimisation results data
  optimDetails.optHist  = optHist;
  optimDetails.exitflag = exitflag;
  optimDetails.lambda   = lambda;


end