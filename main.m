function [optData,outputs,processedOutputs] = main(inputs)
  nx = ones(1,inputs.numDeltaLelems);
  
  %% Optimise operation for every wind speed
  %%        [deltaL,  avgPattEle,  pattAngRadius, startPattRadius, VRI,   CL_ri,  reelOutSpeed, CL_ro,                                       kinematicRatio]
  % All free
    x0     = [200,    deg2rad(30), deg2rad(5),    50,              20*nx, 1.5*nx, 0.5*nx,      inputs.CL_maxAirfoil*inputs.CLeff_F*nx,      120*nx];
    x_init = [500,  deg2rad(90), deg2rad(60), 100, inputs.vk_r_i_max*nx, inputs.CL_maxAirfoil*inputs.CLeff_F*nx, 20*nx,   inputs.CL_maxAirfoil*inputs.CLeff_F*nx, 500*nx];
 
  % Fixed
  % x0     = [250,    deg2rad(30), deg2rad(5),    50,              20*nx, 1.5*nx, 2*nx,      inputs.CL_maxAirfoil*inputs.CLeff_F*nx,      300*nx];
  % x_init = [250,    deg2rad(30), deg2rad(5),    50,              20*nx, 1.5*nx, 2*nx,      inputs.CL_maxAirfoil*inputs.CLeff_F*nx,     120*nx];

  
  for i=1:length(inputs.vw_ref)
   
%    Output of previous wind speed as input to next wind speed
    x0     = x0./x_init;

    % Bounds
    % All free
    lb     = [50,   deg2rad(1),  deg2rad(1),  50,  0*nx,                0*nx,                                   0*nx, 0*nx,                                   1*nx]./x_init; % 
    ub     = [500,  deg2rad(90), deg2rad(60), 100, inputs.vk_r_i_max*nx, inputs.CL_maxAirfoil*inputs.CLeff_F*nx, 20*nx,   inputs.CL_maxAirfoil*inputs.CLeff_F*nx, 500*nx]./x_init; % 

    % Fixed
    % lb     = [250,    deg2rad(30), deg2rad(5),    50,              20*nx, 0*nx, 2*nx,      0*nx,      1*nx]./x_init;   % 
    % ub     = [250,    deg2rad(30), deg2rad(5),    50,              20*nx, inputs.CL_maxAirfoil*inputs.CLeff_F*nx, 2*nx, inputs.CL_maxAirfoil*inputs.CLeff_F*nx,      600*nx]./x_init; % 


    options                           = optimoptions('fmincon');
    options.Display                   = 'iter-detailed';
    options.Algorithm                 = 'sqp';
    options.FiniteDifferenceType      = 'forward';
   % options.FiniteDifferenceStepSize  = [1e-12 1e-12 1e-12 1e-12 1e-6*nx 1e-6*nx 1e-6*nx 1e-6*nx];% 1e-6*nx];
    % options.FiniteDifferenceStepSize  = 1e-12.*[1 1 1 1 nx nx nx nx];
    options.ConstraintTolerance      = 1e-5;
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
  % ms_opts = MultiStart('UseParallel',true,'StartPointsToRun', 'random');
  % ms = GlobalSearch(ms_opts);
  % % ms = MultiStart('UseParallel',true,'StartPoints', 5);
  % problem = createOptimProblem('fmincon', 'objective', Objective, 'x0', x0, 'lb', lb, 'ub', ub, 'nonlcon', con, 'options', options);
  % [x,fval,exitflag(i),Output(i)] = run(ms, problem);

    [x,~,exitflag(i),optHist(i),lambda(i)] = fmincon(Objective,x0,[],[],[],[],lb,ub,con,options);

    % Storing final results
    [~,inputs,outputs] = objective(x,x_init,i,inputs);
    x0 = x.*x_init;

    % Changing initial guess if previous wind speed evaluation is infeasible
    if outputs.P_e_avg(i) <= 0    
        % All free
        x0 = [200,    deg2rad(30), deg2rad(5),    50, 20*nx, 1.5*nx, 0.5*nx, inputs.CL_maxAirfoil*inputs.CLeff_F*nx, 120*nx];        
        % Fixed 
        % x0      = [250,    deg2rad(30), deg2rad(5),    50,              20*nx, 1.5*nx, 2*nx,      inputs.CL_maxAirfoil*inputs.CLeff_F*nx,      300*nx];
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
  temp1                  = vw(outputs.P_e_avg > 0);
  temp2                  = vw(mean((outputs.E_result - outputs.E),2) < 0.001); % Glide ratio constraint violation acceptance
  processedOutputs.cutIn = max(temp1(1),temp2(1)); % At Reference height
 % processedOutputs.cutIn = max(temp1(1));
 

  %% Rated wind and power
  temp3                       = round(outputs.P_e_avg./max(outputs.P_e_avg),2);
  temp4                       = vw(temp3==1);
  processedOutputs.ratedWind  = temp4(1); % At Reference height 
  processedOutputs.ratedPower = outputs.P_e_avg(vw==temp4(1));

  %% Cut-out wind speed
  processedOutputs.cutOut   = 25; % At operational height
  processedOutputs.vw_100m_operRange = vw(mean(outputs.vw,2)<=processedOutputs.cutOut);

  %% System data
  processedOutputs.Dia_te = outputs.d_t;
   for i=1:length(vw)
      if vw(i)>=processedOutputs.cutIn
          processedOutputs.vw(i,:)           = outputs.vw(i,:);
          processedOutputs.P1_m_o(i)         = outputs.P1_m_o(i);
          processedOutputs.P2_m_i(i)         = outputs.P2_m_i(i); 
          processedOutputs.P_m_o_eff(i,:)    = outputs.P_m_o_eff(i,:);
          processedOutputs.P_m_i_eff(i,:)    = outputs.P_m_i_eff(i,:);  
          processedOutputs.P1_e_o(i)         = outputs.P1_e_o(i);
          processedOutputs.P2_e_i(i)         = outputs.P2_e_i(i); 
          processedOutputs.P_e_o_eff(i,:)    = outputs.P_e_o_eff(i,:);
          processedOutputs.P_e_i_eff(i,:)    = outputs.P_e_i_eff(i,:);
          processedOutputs.P_m_o(i)          = outputs.P_m_o(i);
          processedOutputs.P_m_i(i)          = outputs.P_m_i(i);
          processedOutputs.P_e_o(i)          = outputs.P_e_o(i);
          processedOutputs.P_e_i(i)          = outputs.P_e_i(i);
          processedOutputs.deltaL(i)         = outputs.deltaL(i);
          processedOutputs.t1(i)             = outputs.t1(i);
          processedOutputs.t2(i)             = outputs.t2(i);
          processedOutputs.to_eff(i,:)       = outputs.to_eff(i,:);
          processedOutputs.ti_eff(i,:)       = outputs.ti_eff(i,:);
          processedOutputs.to(i)             = outputs.to(i);
          processedOutputs.ti(i)             = outputs.ti(i);
          processedOutputs.tCycle(i)         = outputs.tCycle(i);
          processedOutputs.Ft(i,:)           = outputs.Ft(i,:);
          processedOutputs.Ft_drum(i,:)      = outputs.Ft_drum(i,:);
          processedOutputs.Ft_drum_i(i,:)    = outputs.Ft_drum_i(i,:);
          processedOutputs.Fa(i,:)           = outputs.Fa(i,:);
          processedOutputs.Fa_i(i,:)         = outputs.Fa_i(i,:);
          processedOutputs.Fc(i,:)           = outputs.Fc(i,:);
          processedOutputs.W(i)              = outputs.W(i);
          processedOutputs.vo(i,:)           = outputs.vk_r(i,:);
          processedOutputs.vk_tau(i,:)       = outputs.vk_tau(i,:);
          processedOutputs.vi(i,:)           = outputs.vk_r_i(i,:);
          processedOutputs.P_e_avg(i)        = outputs.P_e_avg(i);
          processedOutputs.P_m_avg(i)        = outputs.P_m_avg(i);
          processedOutputs.Rp(i,:)           = outputs.Rp(i,:);
          processedOutputs.h_cycleStart(i)   = outputs.h_cycleStart(i);
          processedOutputs.h_cycleAvg(i)     = outputs.h_cycleAvg(i);
          processedOutputs.l_t_max(i)        = outputs.l_t_max(i);
          processedOutputs.l_t_inCycle(i,:)  = outputs.l_t_inCycle(i,:);
          processedOutputs.h_inCycle(i,:)    = outputs.h_inCycle(i,:);
          processedOutputs.vk_omega(i,:)     = outputs.vk_omega(i,:); 
          processedOutputs.va(i,:)           = outputs.va(i,:);
          processedOutputs.beta(i)           = rad2deg(outputs.beta(i));
          processedOutputs.rollAngle(i,:)    = rad2deg(-outputs.rollAngle(i,:));
          processedOutputs.gamma(i)          = rad2deg(outputs.gamma(i));
          processedOutputs.CL(i,:)           = outputs.CL(i,:);
          processedOutputs.CD(i,:)           = outputs.CD(i,:);
          processedOutputs.CL_i(i,:)         = outputs.CL_i(i,:);
          processedOutputs.CD_i(i,:)         = outputs.CD_i(i,:);
          processedOutputs.E(i,:)            = outputs.E(i,:);
          processedOutputs.E_i(i,:)          = outputs.E_i(i,:);
          processedOutputs.numOfPatt(i,:)    = outputs.numOfPatt(i,:);
          processedOutputs.lambda(i,:)       = outputs.lambda(i,:); 
          processedOutputs.f(i,:)            = outputs.f(i,:); 
          processedOutputs.f_i(i,:)          = outputs.f_i(i,:); 
          processedOutputs.tPatt(i,:)        = outputs.tPatt(i,:);
          processedOutputs.dutyCycle(i)      = processedOutputs.to(i)/processedOutputs.tCycle(i);
          
          % Coefficient of power (Cp) as defined for HAWTs
%           outputs.sweptArea(i)         = pi()*((outputs.Rp(i)+outputs.b/2)^2 - (outputs.Rp(i)-outputs.b/2)^2);
%           outputs.Cp_reelOut(i)        = outputs.P_m_o(i)/(0.5*outputs.rho_air(i)*outputs.sweptArea(i)*vw(i)^3);
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
  for i = processedOutputs.cutIn:processedOutputs.cutOut
    [processedOutputs.cyclePowerRep(i)] = createCyclePowerRep(i,processedOutputs); 
  end
  function [cyclePowerRep] = createCyclePowerRep(ws, system)
              
    cyclePowerRep.ws     = ws;
    cyclePowerRep.idx    = find(vw==ws);
    cyclePowerRep.t1     = round(system.t1(cyclePowerRep.idx),2);
    cyclePowerRep.to_eff = round(sum(system.to_eff(cyclePowerRep.idx,:),2),2);
    cyclePowerRep.t2     = round(system.t2(cyclePowerRep.idx),2);
    cyclePowerRep.ti_eff = round(sum(system.ti_eff(cyclePowerRep.idx,:),2),2);
    cyclePowerRep.to     = cyclePowerRep.t1 + cyclePowerRep.to_eff; %[s]
    cyclePowerRep.ti     = cyclePowerRep.t2 + cyclePowerRep.ti_eff; %[s]
    cyclePowerRep.tCycle = cyclePowerRep.to + cyclePowerRep.to; %[s]
    
    % cyclePowerRep.t_inst = cumsum([0 cyclePowerRep.t1 system.to_eff(cyclePowerRep.idx,:) cyclePowerRep.t1 cyclePowerRep.t2 system.ti_eff(cyclePowerRep.idx,:) cyclePowerRep.t2 0]);
    % 
    % cyclePowerRep.P_e_inst = [0 system.P_e_o_eff(cyclePowerRep.idx,1) system.P_e_o_eff(cyclePowerRep.idx,:) 0 ...
    %             -flip(system.P_e_i_eff(cyclePowerRep.idx,1)) -flip(system.P_e_i_eff(cyclePowerRep.idx,:)) 0 0]./10^3;
    % 
    % cyclePowerRep.P_m_inst = [0 system.P_m_o_eff(cyclePowerRep.idx,1) system.P_m_o_eff(cyclePowerRep.idx,:) 0 ...
    %             -flip(system.P_m_i_eff(cyclePowerRep.idx,1)) -flip(system.P_m_i_eff(cyclePowerRep.idx,:)) 0 0]./10^3;

    cyclePowerRep.t_inst = cumsum([0 cyclePowerRep.t1 (cyclePowerRep.to_eff-cyclePowerRep.t1) cyclePowerRep.t1 cyclePowerRep.t2 (cyclePowerRep.ti_eff-cyclePowerRep.t2) cyclePowerRep.t2]);

    cyclePowerRep.P_e_inst = [0 system.P_e_o_eff(cyclePowerRep.idx,1) system.P_e_o_eff(cyclePowerRep.idx,end) 0 ...
                -system.P_e_i_eff(cyclePowerRep.idx,end) -system.P_e_i_eff(cyclePowerRep.idx,1) 0]./10^3;

    cyclePowerRep.P_m_inst = [0 system.P_m_o_eff(cyclePowerRep.idx,1) system.P_m_o_eff(cyclePowerRep.idx,end) 0 ...
                -system.P_m_i_eff(cyclePowerRep.idx,end) -system.P_m_i_eff(cyclePowerRep.idx,1) 0]./10^3;
    
  end
  
end