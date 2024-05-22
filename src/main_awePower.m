function [inputs, outputs, optimDetails, processedOutputs] = main_awePower(inputs)


  %% Kite mass 
    if inputs.massOverride == 1
      kiteMass = inputs.kiteMass;
    else
      kiteMass = estimateKiteMass_awePower(inputs.Ft_max, inputs.S, inputs.AR); 
    end

  %% Optimisation

  % Initial guess for optimisation
  x0     = inputs.x0;

  % Bounds
  lb = inputs.lb;
  ub = inputs.ub;

  % Optimise operation for every wind speed
  for i=1:length(inputs.vw_ref)
   
    options                           = optimoptions('fmincon');
    %     options.Display                   = 'iter-detailed';
    options.Display                   = 'final-detailed';
    % options.Display                   = 'notify-detailed';
    options.Algorithm                 = 'sqp';
    options.FiniteDifferenceType      = 'central';
    %   options.FiniteDifferenceType      = 'forward';
    options.ScaleProblem              = true; 
    %   options.MaxFunctionEvaluations    = 5000*numel(x0);
    %   options.MaxIterations             = 500*numel(x0);
 
    
    con = @(x) optConstraints_awePower(i,inputs);
    obj = @(x) optObjective_awePower(x,i,inputs,kiteMass);

    
    [x,~,exitflag(i),optHist(i),lambda(i)] = fmincon(obj,x0,[],[],[],[],lb,ub,con,options);

    % Storing final results
    [~,inputs,outputs] = optObjective_awePower(x,i,inputs,kiteMass);

    % Changing initial guess if previous wind speed evaluation is considered infeasible
    if abs(mean((outputs.E_result - outputs.E),2)) > 1 %0.01 
        x0 = inputs.x0;     
    else
        % Solution of current wind speed as the initial guess for the next wind speed
        x0 = x;
    end 

  end
  disp(exitflag)
  % Store optimisation results data
  optimDetails.optHist  = optHist;
  optimDetails.exitflag = exitflag;
  optimDetails.lambda   = lambda;

  %% Post processing
  vw = inputs.vw_ref; % Wind speed at ref. height

  %% Cut-in wind speed defined at ref. height of 100m
  % Glide ratio constraint violation acceptance for feasible solution
  temp1                  = vw(abs(mean((outputs.E_result - outputs.E),2)) < 0.01); %0.01
  processedOutputs.cutIn = max(temp1(1));

  %% Rated wind and power defined at ref. height of 100m
  temp3                       = round(outputs.P_e_avg./max(outputs.P_e_avg),2);
  temp4                       = vw(temp3>=0.99);
  processedOutputs.ratedWind  = temp4(1); % At Reference height 
  processedOutputs.ratedPower = outputs.P_e_avg(vw==temp4(1));

  %% Cut-out wind speed defined at ref. height of 100m
  % Operating range traslated to wind speed at ref. height
  processedOutputs.vw_100m_operRange = vw(mean(outputs.vw,2)<=25); % Cut-out at 25m/s at operational height (Assumption)
  processedOutputs.cutOut   = processedOutputs.vw_100m_operRange(end);

  %% Extract feasible results in the operational range
  processedOutputs.Dia_te = outputs.d_t;
   for i=1:length(vw)
      if vw(i)>=processedOutputs.cutIn && vw(i)<=processedOutputs.cutOut
          processedOutputs.vw(i,:)           = outputs.vw(i,:);
          processedOutputs.vw_i(i,:)         = outputs.vw_i(i,:);
          processedOutputs.m_k               = outputs.m_k;  
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
          processedOutputs.Ft_i(i,:)         = outputs.Ft_i(i,:);
          processedOutputs.Fa(i,:)           = outputs.Fa(i,:);
          processedOutputs.Fa_i(i,:)         = outputs.Fa_i(i,:);
          processedOutputs.W(i,:)            = outputs.W(i,:);
          processedOutputs.vo(i,:)           = outputs.vk_r(i,:);
          processedOutputs.vk_tau(i,:)       = outputs.vk_tau(i,:);
          processedOutputs.vi(i,:)           = outputs.vk_r_i(i,:);
          processedOutputs.P_e_avg(i)        = outputs.P_e_avg(i);
          processedOutputs.P_m_avg(i)        = outputs.P_m_avg(i);
          processedOutputs.Rp(i,:)           = outputs.Rp(i,:);
          processedOutputs.h_cycleStart(i)   = outputs.h_cycleStart(i);
          processedOutputs.h_cycleAvg(i)     = outputs.h_cycleAvg(i);
          processedOutputs.l_t_max(i)        = outputs.l_t_max(i);
          processedOutputs.l_t_min(i)        = outputs.l_t_min(i);
          processedOutputs.l_t_inCycle(i,:)  = outputs.l_t_inCycle(i,:);
          processedOutputs.h_inCycle(i,:)    = outputs.h_inCycle(i,:); 
          processedOutputs.va(i,:)           = outputs.va(i,:);
          processedOutputs.beta(i)           = rad2deg(outputs.beta(i));
          processedOutputs.rollAngle(i,:)    = rad2deg(-outputs.rollAngle(i,:));
          processedOutputs.gamma(i)          = rad2deg(outputs.gamma(i));
          processedOutputs.CL(i,:)           = outputs.CL(i,:);
          processedOutputs.CD(i,:)           = outputs.CD(i,:);
          processedOutputs.CD_k(i,:)         = outputs.CD_k(i,:);
          processedOutputs.CD_k_i(i,:)       = outputs.CD_k_i(i,:);
          processedOutputs.CD_t(i,:)         = outputs.CD_t(i,:);
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
          processedOutputs.zetaMech(i,:)     = outputs.zetaMech(i,:);  
          processedOutputs.d_te              = outputs.d_t;  

          % Cycle efficiency
          processedOutputs.cycleEff_elec(i)  = (processedOutputs.P_e_avg(i)*processedOutputs.tCycle(i))/...
                                                    (processedOutputs.P_e_o(i)*processedOutputs.to(i));
          processedOutputs.cycleEff_mech(i)  = (processedOutputs.P_m_avg(i)*processedOutputs.tCycle(i))/...
                                                    (processedOutputs.P_m_o(i)*processedOutputs.to(i));

          % Intermediate storage (Ultracapacitor) sizing for power smoothing
          processedOutputs.ultraCapSize(i)   = processedOutputs.P_m_avg(i)*processedOutputs.ti(i)/3.6e6;  %[kWh]

          
          % Coefficient of power (Cp) as defined for HAWTs
          processedOutputs.sweptArea(i,:)     = pi().*((outputs.Rp(i,:)+inputs.b/2).^2 - (outputs.Rp(i,:)-inputs.b/2).^2);
          processedOutputs.Cp_m_o(i,:)        = outputs.P_m_o(i)./(0.5.*inputs.airDensity.*processedOutputs.sweptArea(i,:).*outputs.vw(i,:).^3);
          processedOutputs.Cp_e_avg(i,:)      = outputs.P_e_avg(i)./(0.5.*inputs.airDensity.*processedOutputs.sweptArea(i,:).*outputs.vw(i,:).^3);
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
    
    cyclePowerRep.t_inst   = cumsum([0 cyclePowerRep.t1 (cyclePowerRep.to_eff-cyclePowerRep.t1) cyclePowerRep.t1 ...
                              cyclePowerRep.t2 (cyclePowerRep.ti_eff-cyclePowerRep.t2) cyclePowerRep.t2]);

    cyclePowerRep.P_e_inst = [0 system.P_e_o_eff(cyclePowerRep.idx,1) system.P_e_o_eff(cyclePowerRep.idx,end) 0 ...
                              -system.P_e_i_eff(cyclePowerRep.idx,end) -system.P_e_i_eff(cyclePowerRep.idx,1) 0]./10^3;

    cyclePowerRep.P_m_inst = [0 system.P_m_o_eff(cyclePowerRep.idx,1) system.P_m_o_eff(cyclePowerRep.idx,end) 0 ...
                              -system.P_m_i_eff(cyclePowerRep.idx,end) -system.P_m_i_eff(cyclePowerRep.idx,1) 0]./10^3;
    
  end

  

  %% Save outputs
    
  mkdir 'outputFiles';

  % Change names to associate with specific input file
  save(['outputFiles\' inputs.name '_' 'optimDetails' '.mat'], 'optimDetails');
  save(['outputFiles\' inputs.name '_' 'outputs' '.mat'], 'outputs');
  save(['outputFiles\' inputs.name '_' 'processedOutputs' '.mat'], 'processedOutputs');

end