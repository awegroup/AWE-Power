
 %% Optimise operation for every wind speed

function [optData,outputs,postProRes,timeseries] = main(inputs)
  
  nx = ones(1,inputs.numDeltaLelems);
  %         [deltaL, VRI, avgPattEle,  pattAngRadius, startPattRadius, CL,     rollAngleTop,   reelOutFactor kinematicRatio]
  x0      = [200,    4,   deg2rad(20), deg2rad(12),   40,              inputs.CL_maxAirfoil*inputs.F_CLeff*nx, deg2rad(20)*nx, 0.33*nx,100*nx];

  for i=1:length(inputs.vw_ref)
    
    % Output of previous wind speed as input to next wind speed
    x_init = x0;
    x0     = x_init./x_init;
     
    % Bounds
    % inputs.CL0_airfoil
    lb     = [100, 2,             deg2rad(3),  deg2rad(3),  30,  inputs.CL_maxAirfoil*inputs.F_CLeff*nx,  deg2rad(2)*nx,   0.01*nx, 0*nx]./x_init; % 
    ub     = [600, inputs.maxVRI, deg2rad(80), deg2rad(80), 500, inputs.CL_maxAirfoil*inputs.F_CLeff*nx,   deg2rad(80)*nx, 1*nx, 150*nx]./x_init; % 
    
    options                           = optimoptions('fmincon');
    options.Display                   = 'iter-detailed';
    options.Algorithm                 = 'sqp';
    options.FiniteDifferenceType      = 'central';
%    options.FiniteDifferenceStepSize  = [1e-12 1e-12 1e-8 1e-8 1e-12 1e-6*nx 1e-6*nx 1e-6*nx];% 1e-6*nx];
%     options.FiniteDifferenceStepSize  = 1e-5.*[1 1 1 1 1 nx nx nx nx];
%     options.ConstraintTolerance      = 1e-9;
  %  options.OptimalityTolerance       = 1e-9;
  %     options.StepTolerance             = 1e-6;
    options.MaxFunctionEvaluations    = 800*numel(x_init);
    options.MaxIterations             = 300*numel(x_init);
  %   options.FunctionTolerance        = 1e-9;
  %   options.DiffMaxChange            = 1e-1;
  %   options.DiffMinChange            = 0;
  %   options.OutputFcn                = @O_outfun;
  %    options.PlotFcns                 = {@optimplotfval, @optimplotx, @optimplotfirstorderopt,...
  %                                @optimplotconstrviolation, @optimplotfunccount, @optimplotstepsize};
    con = @(x) constraints(i,inputs);

    [x,~,exitflag(i),optHist(i),lambda(i)] = fmincon(@(x) objective(x,x_init,i,inputs),x0,[],[],[],[],lb,ub,con,options);

    % Storing final results
    [~,inputs,outputs] = objective(x,x_init,i,inputs);
    x0 = x.*x_init;

    % Changing initial guess if previous wind speed evaluation is infeasible
    if outputs.P_cycleElec(i) <= 0
        x0 = [200, 4, deg2rad(20), deg2rad(12), 40,  inputs.CL_maxAirfoil*inputs.F_CLeff*nx, deg2rad(20)*nx, 0.33*nx, 100*nx]; % 
    end  
  end
  disp(exitflag)
  % Store optimisation results data
  optData(1).optHist  = optHist;
  optData(1).exitflag = exitflag;
  optData(1).lambda   = lambda;

  %% Post processing
  Vw = inputs.vw_ref; 

  %% Cut-in wind speed
  temp1            = Vw(outputs.P_cycleElec > 0);
  temp2            = Vw(mean(outputs.VRO(:,:),2) >= 0);
  postProRes.cutIn = max(temp1(1),temp2(1));

  %% Rated wind and power
  temp3 = round(outputs.P_cycleElec./max(outputs.P_cycleElec),2);
  temp4 = Vw(temp3==1);
  postProRes.ratedWind  = temp4(1);
  postProRes.ratedPower = outputs.P_cycleElec(postProRes.ratedWind);

  %% System data
  postProRes.D_te         = outputs.D_te;
  for i=1:length(Vw)
      if Vw(i)>=postProRes.cutIn
          postProRes.Vw(i,:)           = outputs.Vw_top(i,:);
          postProRes.PRO1_mech(i)    = outputs.PRO1_mech(i);
          postProRes.PRI2_mech(i)    = outputs.PRI2_mech(i); 
          postProRes.PROeff_mech(i,:)  = outputs.PROeff_mech(i,:);
          postProRes.PRIeff_mech(i,:)  = outputs.PRIeff_mech(i,:);  
          postProRes.PRO1_elec(i)    = outputs.PRO1_elec(i);
          postProRes.PRI2_elec(i)    = outputs.PRI2_elec(i); 
          postProRes.PROeff_elec(i,:)  = outputs.PROeff_elec(i,:);
          postProRes.PRIeff_elec(i,:)  = outputs.PRIeff_elec(i,:);
          postProRes.PRO_mech(i)    = outputs.PRO_mech(i);
          postProRes.PRI_mech(i)    = outputs.PRI_mech(i);
          postProRes.PRO_elec(i)    = outputs.PRO_elec(i);
          postProRes.PRI_elec(i)    = outputs.PRI_elec(i);
          postProRes.deltaL(i)      = outputs.deltaL(i);
          postProRes.t1(i)          = outputs.t1(i);
          postProRes.t2(i)          = outputs.t2(i);
          postProRes.tROeff(i,:)      = outputs.tROeff(i,:);
          postProRes.tRIeff(i,:)      = outputs.tRIeff(i,:);
          postProRes.tRO(i)         = outputs.tRO(i);
          postProRes.tRI(i)         = outputs.tRI(i);
          postProRes.tCycle(i)      = outputs.tCycle(i);
          postProRes.TT(i,:)          = outputs.T(i,:);
          postProRes.VRO(i,:)         = outputs.VRO(i,:);
          postProRes.VRI(i)         = outputs.VRI(i);
          postProRes.VRI_app(i,:)     = outputs.VA_RI(i,:);
          postProRes.Pcycle_elec(i) = outputs.P_cycleElec(i);
          postProRes.Pcycle_mech(i) = outputs.P_cycleMech(i);
          postProRes.pattRad(i,:)       = outputs.pattRadius(i,:);
          postProRes.H_cycleStart(i)    = outputs.H_cycleStart(i);
          postProRes.H_cycleAvg(i)      = outputs.H_cycleAvg(i);
          postProRes.L_teMax(i)         = outputs.L_teMax(i);
          postProRes.Vc(i,:)            = outputs.VC(i,:); 
          postProRes.Va(i,:)            = outputs.Va_top(i,:);
          postProRes.avgRollAngle(i,:)    = rad2deg(outputs.rollAngleTop(i,:));
          postProRes.avgPattEle(i)      = rad2deg(outputs.avgPattEle(i));
          postProRes.pattAngRadius(i)   = rad2deg(outputs.pattAngRadius(i));
          postProRes.CL(i,:)            = outputs.CL(i,:);
          postProRes.CD(i,:)            = outputs.CD(i,:);
          postProRes.numOfPatt(i,:)     = outputs.numOfPatt(i,:);
          postProRes.reelOutF(i,:)      = outputs.reelOutF(i,:);  
          postProRes.dutyCycle(i)       = postProRes.tRO(i)/postProRes.tCycle(i);
          postProRes.tPatt(i,:)         = 2*pi()*postProRes.pattRad(i,:)/outputs.VC(i);
         
          
          % Coefficient of power (Cp): For comparing with wind turbines

          outputs.sweptArea(i)       = pi()*((outputs.pattRadius(i)+outputs.wingSpan/2)^2 - (outputs.pattRadius(i)-outputs.wingSpan/2)^2);

          outputs.Cp_reelOut(i) = outputs.PRO_mech(i)/(0.5*outputs.rho_air(i)*outputs.sweptArea(i)*Vw(i)^3);
          % Axial induction factor maximising Cp is a = 0.5 and is independent of reel-out factor
          % Ref paper: The Betz limit applied to Airborne Wind Energy, https://doi.org/10.1016/j.renene.2018.04.034
          f(i)    = outputs.reelOutF(i);
          a_ideal = 0.5;
          outputs.Cp_max_ifaIsHalf(i)  = (1-f(i))^2*4*a_ideal*(1-a_ideal)*f(i);
          a                            = abs(roots([1 -1 outputs.Cp_reelOut(i)/(4*f(i)*(1-f(i))^2)]));
          outputs.axialInductionF(i,1) = a(1);
          outputs.axialInductionF(i,2) = a(2);
          
      end      
  end

  %% Representative pattern average cycle timeseries data
  timeseries = struct();
  for i = postProRes.cutIn:length(Vw)
    [timeseries.ws(i)] = createTimeseries(i,postProRes);
  end
  
end