function [optData,outputs,processedOutputs] = main(inputs)
  nx = ones(1,inputs.numDeltaLelems);
  
  %% Optimise operation for every wind speed
  %%        [deltaL, VRI,               avgPattEle,  pattAngRadius, startPattRadius, reelOutSpeed, kinematicRatio]
  % All free
    x0     = [200,    inputs.vk_r_i_max, deg2rad(30), deg2rad(5),    50,              0.2*nx,     80*nx];
  
  % Fixed
  % x0     = [250,    20,   deg2rad(30), deg2rad(10),   50, 4*nx, 80*nx];
  % x_init = [250,    20,   deg2rad(30), deg2rad(10),   80, 4*nx, 80*nx];

  
  for i=1:length(inputs.vw_ref)
   
%    Output of previous wind speed as input to next wind speed
    x_init = x0;
    x0     = x_init./x_init;
    % Bounds
    % All free
    lb     = [50,   1,                 deg2rad(1),  deg2rad(1),  50,  0.1*nx,   1*nx]./x_init; % 
    ub     = [600,  inputs.vk_r_i_max, deg2rad(90), deg2rad(60), 100,                          25*nx,    300*nx]./x_init; % 

    % x0 = x0./x_init;
    % Fixed
    % lb     = [250,    20,   deg2rad(30), deg2rad(10),   50, 4*nx, 1*nx]./x_init;   % 
    % ub     = [250,    20,   deg2rad(30), deg2rad(10),   50, 4*nx, 250*nx]./x_init; % 


    options                           = optimoptions('fmincon');
    options.Display                   = 'iter-detailed';
    options.Algorithm                 = 'sqp';
    options.FiniteDifferenceType      = 'central';
%    options.FiniteDifferenceStepSize  = [1e-12 1e-12 1e-8 1e-8 1e-12 1e-6*nx 1e-6*nx 1e-6*nx];% 1e-6*nx];
%     options.FiniteDifferenceStepSize  = 1e-5.*[1 1 1 1 1 nx nx nx nx];
%     options.ConstraintTolerance      = 1e-9;
  %  options.OptimalityTolerance       = 1e-9;
  %     options.StepTolerance             = 1e-6;
    options.MaxFunctionEvaluations    = 3000*numel(x_init);
    options.MaxIterations             = 300*numel(x_init);
  %   options.FunctionTolerance        = 1e-9;
  %   options.DiffMaxChange            = 1e-1;
  %   options.DiffMinChange            = 0;
  %   options.OutputFcn                = @O_outfun;
  %    options.PlotFcns                 = {@optimplotfval, @optimplotx, @optimplotfirstorderopt,...
  %                                @optimplotconstrviolation, @optimplotfunccount, @optimplotstepsize};
    con = @(x) constraints(i,inputs);

    Objective = @(x) objective(x,x_init,i,inputs);

  % % Multistart options
  % ms_opts = MultiStart('UseParallel',true);
  % problem = createOptimProblem('fmincon', 'objective', Objective, 'x0', x0, 'lb', lb, 'ub', ub, 'nonlcon', con, 'options', options);
  % ms = GlobalSearch(ms_opts);
  % [x,~] = run(ms, problem);

    [x,~,exitflag(i),optHist(i),lambda(i)] = fmincon(Objective,x0,[],[],[],[],lb,ub,con,options);

    % Storing final results
    [~,inputs,outputs] = objective(x,x_init,i,inputs);
    x0 = x.*x_init;

    % Changing initial guess if previous wind speed evaluation is infeasible
    if outputs.P_cycleElec(i) <= 0    
        % All free
        x0 = [200, inputs.vk_r_i_max, deg2rad(30), deg2rad(5),    50,             0.2*nx,       80*nx];        
        % Fixed 
        % x0      = [250,    20,   deg2rad(30), deg2rad(10),   50,               4*nx,      80*nx];
    end  
  end
  disp(exitflag)
  % Store optimisation results data
  optData(1).optHist  = optHist;
  optData(1).exitflag = exitflag;
  optData(1).lambda   = lambda;

  %% Post processing
  vw = inputs.vw_ref; 

  %% Cut-in wind speed
  temp1                  = vw(outputs.P_cycleElec > 0);
  temp2                  = vw(outputs.vk_r(:,1) > 0.3); % min speed control limit of DT
  processedOutputs.cutIn = max(temp1(1),temp2(1));
%  processedOutputs.cutIn = max(temp1(1));

  %% Rated wind and power
  temp3                       = round(outputs.P_cycleElec./max(outputs.P_cycleElec),2);
  temp4                       = vw(temp3==1);
  processedOutputs.ratedWind  = temp4(1);
  processedOutputs.ratedPower = outputs.P_cycleElec(vw==temp4(1));

  %% System data
  processedOutputs.Dia_te = outputs.d_t;
   for i=1:length(vw)
      if vw(i)>=processedOutputs.cutIn
          processedOutputs.Vw(i,:)           = outputs.vw(i,:);
          processedOutputs.PRO1_mech(i)      = outputs.PRO1_mech(i);
          processedOutputs.PRI2_mech(i)      = outputs.PRI2_mech(i); 
          processedOutputs.PROeff_mech(i,:)  = outputs.PROeff_mech(i,:);
          processedOutputs.PRIeff_mech(i,:)  = outputs.PRIeff_mech(i,:);  
          processedOutputs.PRO1_elec(i)      = outputs.PRO1_elec(i);
          processedOutputs.PRI2_elec(i)      = outputs.PRI2_elec(i); 
          processedOutputs.PROeff_elec(i,:)  = outputs.PROeff_elec(i,:);
          processedOutputs.PRIeff_elec(i,:)  = outputs.PRIeff_elec(i,:);
          processedOutputs.PRO_mech(i)       = outputs.PRO_mech(i);
          processedOutputs.PRI_mech(i)       = outputs.PRI_mech(i);
          processedOutputs.PRO_elec(i)       = outputs.PRO_elec(i);
          processedOutputs.PRI_elec(i)       = outputs.PRI_elec(i);
          processedOutputs.deltaL(i)         = outputs.deltaL(i);
          processedOutputs.t1(i)             = outputs.t1(i);
          processedOutputs.t2(i)             = outputs.t2(i);
          processedOutputs.tROeff(i,:)       = outputs.tROeff(i,:);
          processedOutputs.tRIeff(i,:)       = outputs.tRIeff(i,:);
          processedOutputs.tRO(i)            = outputs.tRO(i);
          processedOutputs.tRI(i)            = outputs.tRI(i);
          processedOutputs.tCycle(i)         = outputs.tCycle(i);
          processedOutputs.Ft(i,:)           = outputs.Ft(i,:);
          processedOutputs.Fa(i,:)           = outputs.Fa(i,:);
          processedOutputs.Fc(i,:)           = outputs.Fc(i,:);
          processedOutputs.W(i)              = outputs.W(i);
          processedOutputs.VRO(i,:)          = outputs.vk_r(i,:);
          processedOutputs.VRI(i)            = outputs.vk_r_i(i);
          processedOutputs.VRI_app(i,:)      = outputs.va_i(i,:);
          processedOutputs.Pcycle_elec(i)    = outputs.P_cycleElec(i);
          processedOutputs.Pcycle_mech(i)    = outputs.P_cycleMech(i);
          processedOutputs.pattRad(i,:)      = outputs.Rp(i,:);
          processedOutputs.H_cycleStart(i)   = outputs.h_cycleStart(i);
          processedOutputs.H_cycleAvg(i)     = outputs.h_cycleAvg(i);
          processedOutputs.L_teMax(i)        = outputs.l_t_max(i);
          processedOutputs.Vc(i,:)           = outputs.vk_omega(i,:); 
          processedOutputs.Va(i,:)           = outputs.va(i,:);
          processedOutputs.avgPattEle(i)     = rad2deg(outputs.beta(i));
          processedOutputs.rollAngle(i,:)    = rad2deg(-outputs.rollAngle(i,:));
          processedOutputs.pattAngRadius(i)  = rad2deg(outputs.gamma(i));
          processedOutputs.CL(i,:)           = outputs.CL(i,:);
          processedOutputs.CD(i,:)           = outputs.CD(i,:);
          processedOutputs.numOfPatt(i,:)    = outputs.numOfPatt(i,:);
          processedOutputs.reelOutF(i,:)     = outputs.f(i,:);  
          processedOutputs.dutyCycle(i)      = processedOutputs.tRO(i)/processedOutputs.tCycle(i);
          processedOutputs.tPatt(i,:)        = 2*pi()*processedOutputs.pattRad(i,:)/outputs.vk_omega(i);
          
          % Coefficient of power (Cp) as defined for HAWTs
%           outputs.sweptArea(i)         = pi()*((outputs.Rp(i)+outputs.b/2)^2 - (outputs.Rp(i)-outputs.b/2)^2);
%           outputs.Cp_reelOut(i)        = outputs.PRO_mech(i)/(0.5*outputs.rho_air(i)*outputs.sweptArea(i)*vw(i)^3);
%           % Axial induction factor maximising Cp is a = 0.5 and is independent of reel-out factor
%           % Ref paper: The Betz limit applied to Airborne Wind Energy, https://doi.org/10.1016/j.renene.2018.04.034
%           f(i)                         = outputs.f(i);
%           a_ideal                      = 0.5;
%           outputs.Cp_max_ifaIsHalf(i)  = (1-f(i))^2*4*a_ideal*(1-a_ideal)*f(i);
%           a                            = abs(roots([1 -1 outputs.Cp_reelOut(i)/(4*f(i)*(1-f(i))^2)]));
%           outputs.axialInductionF(i,1) = a(1);
%           outputs.axialInductionF(i,2) = a(2);      
      end      
  end

  %% Cycle power representation for wind speeds in the operational range
  for i = processedOutputs.cutIn:vw(end)
    [processedOutputs.cyclePowerRep(i)] = createCyclePowerRep(i,processedOutputs); 
  end
  function [cpr] = createCyclePowerRep(ws, system)
              
    cpr.ws     = ws;
    cpr.idx    = find(vw==ws);
    cpr.t1     = round(system.t1(cpr.idx),2);
    cpr.tROeff = round(sum(system.tROeff(cpr.idx,:),2),2);
    cpr.t2     = round(system.t2(cpr.idx),2);
    cpr.tRIeff = round(sum(system.tRIeff(cpr.idx,:),2),2);
    cpr.tRO    = cpr.t1 + cpr.tROeff; %[s]
    cpr.tRI    = cpr.t2 + cpr.tRIeff; %[s]
    cpr.tCycle = cpr.tRO+cpr.tRI; %[s]
    
    cpr.t_inst = cumsum([0 cpr.t1 system.tROeff(cpr.idx,:) cpr.t1 cpr.t2 system.tRIeff(cpr.idx,:) cpr.t2 0]);
    
    cpr.P_e_inst = [0 system.PROeff_elec(cpr.idx,1) system.PROeff_elec(cpr.idx,:) 0 ...
                -flip(system.PRIeff_elec(cpr.idx,1)) -flip(system.PRIeff_elec(cpr.idx,:)) 0 0]./10^3;
              
    cpr.P_m_inst = [0 system.PROeff_mech(cpr.idx,1) system.PROeff_mech(cpr.idx,:) 0 ...
                -flip(system.PRIeff_mech(cpr.idx,1)) -flip(system.PRIeff_mech(cpr.idx,:)) 0 0]./10^3;
    
  end
  
end