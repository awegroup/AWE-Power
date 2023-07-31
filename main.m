function [optData,outputs,postProRes,timeseries] = main(inputs)
  
  %% Optimise operation for every wind speed: Uncapped electrical power
  
  nx = ones(1,inputs.numDeltaLelems);
  %         [deltaL, VRI, avgPattEle,  pattAngRadius, startPattRadius, CL,     rollAngleTop,   kiteSpeedTangTop, rollAngleBot]
  x0      = [200,    4,   deg2rad(20), deg2rad(12),   40,              1.5*nx, deg2rad(50)*nx, 20*nx,            deg2rad(70)*nx]; 
  
  for i=1:length(inputs.Vw_ref)
    % Output of previous wind speed as input to next wind speed
    x_init = x0;  
    x0     = x_init./x_init;
    % sqrt(inputs.AR*inputs.WA)/2
    lb     = [100, 2,             deg2rad(5),  deg2rad(3),  30,  inputs.CL0_airfoil*nx,                  deg2rad(5)*nx,  20*nx, deg2rad(5)*nx]./x_init; % 
    ub     = [600, inputs.maxVRI, deg2rad(80), deg2rad(80), 500, inputs.CL_maxAirfoil*inputs.F_CLeff*nx, deg2rad(80)*nx, 50*nx, deg2rad(80)*nx]./x_init; % 
    options                           = optimoptions('fmincon');
    options.Display                   = 'iter-detailed';
    options.Algorithm                 = 'sqp';
    options.FiniteDifferenceType      = 'central';
%    options.FiniteDifferenceStepSize  = [1e-12 1e-12 1e-8 1e-8 1e-12 1e-6*nx 1e-6*nx 1e-6*nx];
    options.FiniteDifferenceStepSize  = 1e-5.*[1 1 1 1 1 nx nx nx nx];
  %  options.OptimalityTolerance       = 1e-9;
  %     options.StepTolerance             = 1e-6;
    options.MaxFunctionEvaluations    = 800*numel(x_init);
    options.MaxIterations             = 200*numel(x_init);
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
        x0 = [200, 4, deg2rad(20), deg2rad(12), 60, 1.5*nx, deg2rad(50)*nx, 20*nx, deg2rad(70)*nx]; % 
    end  
  end
  disp(exitflag)
  % Store optimisation results data
  optData(1).optHist  = optHist;
  optData(1).exitflag = exitflag;
  optData(1).lambda   = lambda;
  
  %% Second optimisation iteration for following mean of capped mech. power from the first iteration
%   outputs1 = outputs;
%   clear outputs
%   inputs.targetPRO_mech = outputs1.PROeff_mech_cap;
%   %      [deltaL, VRI, CL, avgPattEle, pattAngRadius, startPattRadius]
%   x0      = [200, 4, 1.5, deg2rad(20),deg2rad(12),200]; % 
%   for i=1:length(inputs.Vw_ref)
%     % Output of previous wind speed as input to next wind speed
%     x_init = x0;  
%   %  x_init    = [200, 4, 1.5, deg2rad(20),deg2rad(12),200]; % 
%     x0     = x_init./x_init;
%     lb     = [100, 2, inputs.CL0_airfoil, deg2rad(5), deg2rad(5),sqrt(inputs.AR*inputs.WA)/2]./x_init; % 
%     ub     = [600, inputs.maxVRI, inputs.CL_maxAirfoil*inputs.F_CLeff, deg2rad(80),deg2rad(80),500]./x_init; % 
%     options                           = optimoptions('fmincon');
%     options.Display                   = 'iter-detailed';
%     options.Algorithm                 = 'sqp';
%     options.FiniteDifferenceType      = 'central';
%   %   options.FiniteDifferenceStepSize  = [1e-12 1e-12 eps^(1/3) eps^(1/3) eps^(1/3) eps^(1/3)];
%   %options.FiniteDifferenceStepSize  = [eps^(1/3) eps^(1/3) eps^(1/3) 1e-12 eps^(1/3) eps^(1/3)];
%   %   options.OptimalityTolerance       = 1e-9;
%   %   options.StepTolerance             = 1e-6;
%     options.MaxFunctionEvaluations    = 1000*numel(x_init);
%     options.MaxIterations             = 1000*numel(x_init);
%     options.ConstraintTolerance       = 1e-3;
% 
%     con = @(x) constraints(i,inputs);
% 
%     [x,fval,exitflag(i),optHist(i),lambda(i)] = fmincon(@(x) objective(x,x_init,i,inputs),x0,[],[],[],[],lb,ub,con,options);
% 
%     % Storing final results
%      [~,inputs,outputs] = objective(x,x_init,i,inputs);
%      x0 = x.*x_init;
% 
% 
%      % Changing initial guess if previous wind speed evaluation is infeasible
%      if outputs.P_cycleElec(i) <= 0
%          x0 = [200, 5, 1.5, deg2rad(20),deg2rad(12),200]; % 
%      end  
%   end
% 
%   % Storing back the capped oscillating power results from first optimisation
%   outputs.PROeff_mech_osci = outputs1.PROeff_mech_osci_cap;
%   outputs.PROeff_elec_osci = outputs1.PROeff_elec_osci_cap;
%   % Store optimisation results data
%   optData(2).optHist  = optHist;
%   optData(2).exitflag = exitflag;

  %% Post processing
  Vw = inputs.Vw_ref; 

  %% Cut-in wind speed
  % temp1         = Vw((outputs.VRO-outputs.VRO_osciAmp)>0);
  temp2         = Vw(outputs.P_cycleElec > 0);
  temp5         = Vw(mean(outputs.VRO(:,:),2) >=1);
  % postProRes.cutIn = max(temp1(1), temp2(1),temp(5));
  postProRes.cutIn = max(temp2(1),temp5(1));

  %% Rated wind and power
  temp3 = round(outputs.P_cycleElec./max(outputs.P_cycleElec),2);
  temp4 = Vw(temp3==1);
  postProRes.ratedWind  = temp4(1);
  postProRes.ratedPower = outputs.P_cycleElec(postProRes.ratedWind);

  %% System data
  % system.VRO_osci         = zeros(length(Vw),length(outputs.VRO_osci(1,:)));
  % system.PROeff_mech_osci = zeros(length(Vw),length(outputs.VRO_osci(1,:)));
  % system.PROeff_elec_osci = zeros(length(Vw),length(outputs.VRO_osci(1,:)));
  % system.numPattParts     = outputs.numPattParts;
  postProRes.D_te         = outputs.D_te;
  for i=1:length(Vw)
      if Vw(i)>=postProRes.cutIn
          postProRes.Vw(i,:)           = outputs.Vw(i,:);
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
  %         system.VRO_osci(i,:)  = outputs.VRO_osci(i,:);
  %         system.PROeff_mech_osci(i,:)  = outputs.PROeff_mech_osci(i,:);
  %         system.PROeff_elec_osci(i,:)  = outputs.PROeff_elec_osci(i,:);
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
          postProRes.reelOutF(i,:)      = postProRes.VRO(i,:)/postProRes.Vw(i,:);  
          postProRes.dutyCycle(i)       = postProRes.tRO(i)/postProRes.tCycle(i);
          postProRes.tPatt(i,:)         = 2*pi()*postProRes.pattRad(i,:)/outputs.VC(i);
      end
  %     for j = 1:outputs.numPattParts
  %       system.tPatt(i)  = system.tPatt(i) + 2*pi()*system.pattRadius(i)/outputs.numPattParts/outputs.VC_osci(i,j);
  %     end 
  %     system.osciFactor(i,:) = outputs.osciFactor(i,:);

  %     %% Coefficient of power (Cp): For comparing with wind turbines

  %     outputs.sweptArea(i)       = pi()*((outputs.pattRadius(i)+outputs.wingSpan/2)^2 - (outputs.pattRadius(i)-outputs.wingSpan/2)^2);

  %     outputs.Cp_reelOut(i) = outputs.PRO_mech(i)/(0.5*outputs.rho_air(i)*outputs.sweptArea(i)*inputs.Vw(i)^3);
  %     % Axial induction factor maximising Cp is a = 0.5 and is independent of reel-out factor
  %     % Ref paper: The Betz limit applied to Airborne Wind Energy, https://doi.org/10.1016/j.renene.2018.04.034
  %     f(i)    = outputs.reelOutF(i);
  %     a_ideal = 0.5;
  %     outputs.Cp_max_ifaIsHalf(i)  = (1-f(i))^2*4*a_ideal*(1-a_ideal)*f(i);
  %     a                            = abs(roots([1 -1 outputs.Cp_reelOut(i)/(4*f(i)*(1-f(i))^2)]));
  %     outputs.axialInductionF(i,1) = a(1);
  %     outputs.axialInductionF(i,2) = a(2);
  end

  %% Test
%   figure()
%   plot(mean(postProRes.PROeff_mech,2)' - inputs.targetPRO_mech);
% 
%   figure()
%   plot(inputs.targetPRO_mech)
%   hold on
%   plot(mean(postProRes.PROeff_mech,2)')

  %% Representative pattern average cycle timeseries data
  timeseries = struct();
  for i = postProRes.cutIn:length(Vw)
    [timeseries.ws(i)] = createTimeseries(i,postProRes);
  end
  
end