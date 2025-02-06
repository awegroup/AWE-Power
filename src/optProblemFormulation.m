function [optimDetails, outputs] = optProblemFormulation(inputs, kiteMass)

  % optProblemFormulation Optimizes the power generation for an AWE system across different wind speeds
  %
  % This function formulates and solves an optimization problem for a fixed-wing airborne wind
  % energy (AWE) system at multiple wind speed reference values. The optimization aims to
  % maximize the electrical power output by adjusting operational parameters within specified
  % bounds. The function uses the `fmincon` solver with the Sequential Quadratic Programming (SQP)
  % algorithm to find optimal parameters for each wind speed.
  %
  % Inputs:
  %   inputs    (struct): Structure containing input parameters including initial guesses,
  %                       bounds, and wind speed references.
  %   kiteMass  (double): Mass of the kite used in the AWE system, either provided or calculated.
  %
  % Outputs:
  %   optimDetails (struct): Structure containing details of the optimization process such as
  %                          optimization history, exit flags, and Lagrange multipliers.
  %   outputs      (struct): Structure containing computed values from the optimization, including
  %                          updated parameters and results for each wind speed.
  %
  % The function iterates over each wind speed reference value, performing optimization with
  % constraints and an objective function tailored for that specific wind speed. The results from
  % each optimization are used to update the initial guess for the next wind speed iteration.
  % The function captures the optimization details and final results for further analysis.


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
    %       options.FiniteDifferenceType      = 'forward';
    options.ScaleProblem              = true;
    %           options.MaxFunctionEvaluations    = 5000*numel(x0);
    %           options.MaxIterations             = 500*numel(x0);

    %     options.ConstraintTolerance       = 1e-3;


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