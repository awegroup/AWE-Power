function [inputs] = compute(i,inputs)
 global outputs
    
    %% Kite mass 
    if inputs.massOverride == 1
      outputs.m_kite = inputs.kiteMass;
    else
      %Vincent's simple mass model      
      a1     = 0.002415;       a2     = 0.0090239;       b1     = 0.17025;       b2     = 3.2493;
      k1     = 5;              c1     = 0.46608;         d1     = 0.65962;       k2     = 1.1935;
      AR_ref = 12;
      a = a1*(inputs.Tmax/inputs.WA) + a2;
      b = b1*(inputs.Tmax/inputs.WA) + b2;
      outputs.m_kite = 10*(a*inputs.WA^2 +b*inputs.WA-k1)*(c1*(inputs.AR/AR_ref)^2-d1*(inputs.AR/AR_ref)+k2); 
    end
    
    %% Minimum Tether length and minimum patt. radius - Centrifugal force balance - no gravity
    outputs.wingSpan           = sqrt(inputs.AR*inputs.WA);     
    outputs.L_teMin         = outputs.startPattRadius/sin(outputs.pattAngRadius)*inputs.F_minTeLen;
    outputs.pattStartGrClr  = outputs.L_teMin*sin(outputs.avgPattEle-outputs.pattAngRadius);
    outputs.H_cycleStart    = outputs.L_teMin*cos(outputs.pattAngRadius)*sin(outputs.avgPattEle);
    outputs.L_teMax         = outputs.L_teMin+outputs.deltaL/cos(outputs.pattAngRadius); 
    outputs.L_teAvg         = (outputs.L_teMax+outputs.L_teMin)/2; %[m]
    outputs.H_cycleAvg      = outputs.L_teAvg*cos(outputs.pattAngRadius)*sin(outputs.avgPattEle);
    outputs.H_cycleEnd      = outputs.L_teMax*cos(outputs.pattAngRadius)*sin(outputs.avgPattEle);
    outputs.D_te               = sqrt(inputs.Tmax*1000/inputs.Te_matStrength*4/pi()); %[m] safety factor could be added (say *1.1)
    
    %% Effective mass (Kite + tether)
    outputs.m_te  = inputs.Te_matDensity*pi()/4*outputs.D_te^2*outputs.L_teAvg; % Could add (say 0.85) as safety factor on material density
    outputs.m_eff = outputs.m_kite+sin(outputs.avgPattEle)*outputs.m_te;

    %% Effective CD
    outputs.CD_kite   = inputs.CD0 + (outputs.CL-inputs.CL0_airfoil)^2/(pi()*inputs.AR*inputs.e);
    outputs.CD_tether = (1/4)*inputs.CD_te*outputs.D_te*outputs.L_teAvg/inputs.WA;
    outputs.CD        = outputs.CD_kite + outputs.CD_tether;

    %% Tmax and W
    outputs.Tmax_act = inputs.Tmax*inputs.F_Tmax*1000; %[N]
    outputs.W     = outputs.m_eff*inputs.gravity;
    
    %% Cycle avg. mech reel-out power considering vertical wind shear
    
    intermRes = deltaLDiscret(i, inputs);
    
%     merged_structure = outputs;
    fields = fieldnames(intermRes);
    for k = 1:numel(fields)
        field = fields{k};
        outputs.(field) = intermRes.(field);
    end

  function [intermRes] = deltaLDiscret(i, inputs)
 
     
    % Divide the reel-out length in equal number of elements of around 5m lengths
     intermRes.deltaLelems = 10; %round(outputs.deltaL/5);
     intermRes.elemDeltaL = outputs.deltaL/intermRes.deltaLelems;
  
     % For each element of deltaL
     for j = 1:intermRes.deltaLelems
       if j == 1
         intermRes.h_inCycle(j) = outputs.H_cycleStart + intermRes.elemDeltaL/2*sin(outputs.avgPattEle);
       else
         intermRes.h_inCycle(j) = intermRes.h_inCycle(j-1) + intermRes.elemDeltaL*sin(outputs.avgPattEle);
       end
       
       intermRes.Vw(j) = inputs.Vw_ref(i)*(intermRes.h_inCycle(j)/inputs.h_ref)^inputs.windShearExp;
        % Air density as a function of height   % Ref: https://en.wikipedia.org/wiki/Density_of_air
        M = 0.0289644; % [kg/mol]
        R = 8.3144598; % [N·m/(mol·K)]
        T = 288.15; %[Kelvin]
        L = 0.0065; %[Kelvin/m] 
        intermRes.rho_air(j) = inputs.airDensity*(1-L*(5+ intermRes.h_inCycle(j))/T)^(inputs.gravity*M/R/L-1); % +5 is for approx. platform height

        % Pattern radius at point of interest on deltaL
        if j == 1
          intermRes.pattRadius(j) = (outputs.startPattRadius + (intermRes.elemDeltaL*tan(outputs.pattAngRadius) + outputs.startPattRadius))/2;
        else
          intermRes.pattRadius(j) = (intermRes.pattRadius(j-1) + (j*intermRes.elemDeltaL*tan(outputs.pattAngRadius) + outputs.startPattRadius))/2;
        end

        % Roll angle calculation for the particular element
        intermRes.avgRollAngle(j) = asin(outputs.m_eff*cos(outputs.pattAngRadius)/...
                                (0.5*intermRes.rho_air(j)*inputs.WA*outputs.CL*intermRes.pattRadius(j)));

        % Tether tension: Reduction in lift due to roll
        intermRes.T(j)   = min(outputs.Tmax_act, (4/9)*(outputs.CL*cos(intermRes.avgRollAngle(j)))^3/outputs.CD^2*(1/2).*intermRes.rho_air(j)*...
                          inputs.WA*(intermRes.Vw(j).*cos(outputs.avgPattEle))^2);

        % Sink rate: Effect of gravity
        intermRes.J(j)   = sqrt(intermRes.T(j)^2 + outputs.W^2 + 2*intermRes.T(j)*outputs.W*sin(outputs.avgPattEle)); % Equilibrium force for Aerodynamic resultant
        intermRes.VSR(j) = outputs.CD/(outputs.CL*cos(intermRes.avgRollAngle(j)))^(3/2)*...
                             (intermRes.T(j)+outputs.W*sin(outputs.avgPattEle))/sqrt(0.5*intermRes.rho_air(j)*inputs.WA*intermRes.J(j));

        % VRO
        intermRes.VRO(j) = intermRes.Vw(j)*cos(outputs.avgPattEle)-intermRes.VSR(j);

        % PRO 
        intermRes.PROeff_mech(j) = intermRes.T(j)*intermRes.VRO(j); %[W]

        % PRO elec
        % Generator efficiency. As a function of RPM/RPM_max, where RPM_max is driven by winch i.e Max VRI
        intermRes.genEff_RO(j)  = (inputs.etaGen.param(1)*(intermRes.VRO(j)/inputs.etaGen.Vmax)^3 + ...
                                      inputs.etaGen.param(2)*(intermRes.VRO(j)/inputs.etaGen.Vmax)^2 + ...
                                        inputs.etaGen.param(3)*(intermRes.VRO(j)/inputs.etaGen.Vmax)+inputs.etaGen.param(4))^sign(1);
        intermRes.PROeff_elec(j) = intermRes.PROeff_mech(j)*inputs.etaGearbox*intermRes.genEff_RO(j)*inputs.etaPE;

        % Airspeed and tangential speed
         intermRes.lambda(j) = 0.5*intermRes.rho_air(j)*inputs.WA*outputs.CL;
         intermRes.delta(j)  = 0.5*intermRes.rho_air(j)*inputs.WA*outputs.CD;
         intermRes.a(j)      = sqrt((intermRes.lambda(j)^2+intermRes.delta(j)^2))/intermRes.J(j);
         intermRes.VA(j)     = 1/sqrt(intermRes.a(j)); % Apparent speed [m/s];
         intermRes.VC(j)     = sqrt(intermRes.VA(j)^2 - intermRes.VSR(j)^2);

        % PRI effective mech
        intermRes.VA_RI(j)       = sqrt(intermRes.Vw(j)^2 +outputs.VRI^2 +2*intermRes.Vw(j)*outputs.VRI*cos(outputs.avgPattEle));
        intermRes.CL_RI(j)       = 2*outputs.W/(intermRes.rho_air(j)*intermRes.VA_RI(j)^2*inputs.WA);
        intermRes.CD_RI(j)       = inputs.CD0+(intermRes.CL_RI(j)- inputs.CL0_airfoil)^2/(pi()*inputs.AR*inputs.e) + outputs.CD_tether;

        intermRes.PRIeff_mech(j) = 0.5*intermRes.rho_air(j)*intermRes.CD_RI(j)*inputs.WA*intermRes.VA_RI(j)^3;

        % Generator efficiency during RI: As a function of RPM/RPM_max, where RPM_max is driven by winch i.e Max VRI
        intermRes.genEff_RI = (inputs.etaGen.param(1)*(outputs.VRI/inputs.etaGen.Vmax)^3 + ...
                                 inputs.etaGen.param(2)*(outputs.VRI/inputs.etaGen.Vmax)^2 + ...
                                  inputs.etaGen.param(3)*(outputs.VRI/inputs.etaGen.Vmax)+inputs.etaGen.param(4))^sign(1);
        % PRI effective elec
        intermRes.PRIeff_elec(j) = intermRes.PRIeff_mech(j)/inputs.etaGearbox/inputs.etaSto/intermRes.genEff_RI/inputs.etaPE;

     end
     
   end
    
      %% Cycle simulation

      % tRO
       if outputs.VRO(:)<0
        outputs.t1             = 0;
        outputs.tROeff(:)       = 0./outputs.VRO(:);
        outputs.tRO            = 0;
       else
        outputs.t1       = outputs.VRO(1)/inputs.maxAcc;
        outputs.tROeff(:) = outputs.elemDeltaL./outputs.VRO(:);
        outputs.tRO      = outputs.t1 + sum(outputs.tROeff(:));
       end
     

      % Reel-out power during transition
      outputs.PRO1_mech = outputs.PROeff_mech(1)/2;
      outputs.PRO1_elec = outputs.PROeff_elec(1)/2;

      % PRO 
      if outputs.VRO(:)<0
        outputs.PRO_mech = 1e-9;
      else
        outputs.PRO_mech = (sum(outputs.PROeff_mech(:).*outputs.tROeff(:)) + outputs.PRO1_mech*outputs.t1)/outputs.tRO;
        outputs.PRO_elec = (sum(outputs.PROeff_elec(:).*outputs.tROeff(:)) + outputs.PRO1_elec*outputs.t1)/outputs.tRO;
      end

      
      % tRI
      if outputs.VRO(:)<0
        outputs.t2             = 0;
        outputs.tRIeff(:)       = ones(1,outputs.deltaLelems).*0;
        outputs.tRI            = 0;
      else
        outputs.t2       = outputs.VRI/inputs.maxAcc;
        outputs.tRIeff(:) = ones(1,outputs.deltaLelems).*outputs.elemDeltaL/outputs.VRI;
        outputs.tRI      = outputs.t2 + sum(outputs.tRIeff(:));
      end
      
      % Transition reel-in mech power
      outputs.PRI2_mech     = outputs.PRIeff_mech(1)/2;
      outputs.PRI2_elec     = outputs.PRIeff_elec(1)/2;

      % PRI 
      outputs.PRI_mech = (sum(outputs.PRIeff_mech(:).*outputs.tRIeff(:)) + outputs.PRI2_mech*outputs.t2)/outputs.tRI;
      outputs.PRI_elec = (sum(outputs.PRIeff_elec(:).*outputs.tRIeff(:)) + outputs.PRI2_elec*outputs.t2)/outputs.tRI;

      % tCycle
      outputs.tCycle = outputs.tRO+outputs.tRI;

      outputs.tPatt(:)     = 2*pi()*outputs.pattRadius(:)./outputs.VC(:);
      outputs.numOfPatt(:) = outputs.tRO./outputs.tPatt(:);

      % P_cycleElec 
      if outputs.VRO(:)<0
        outputs.P_cycleElec = 0;
      else
       outputs.P_cycleElec = (sum(outputs.tROeff(:).*outputs.PROeff_elec(:)) + outputs.t1*outputs.PRO1_elec - ...
                                   sum(outputs.tRIeff(:).*outputs.PRIeff_elec(:)) -outputs.t2*outputs.PRI2_elec)/outputs.tCycle;    
      end

      % Without drivetrain eff
      outputs.P_cycleMech = (sum(outputs.tROeff(:).*outputs.PROeff_mech(:)) + outputs.t1*outputs.PRO1_mech - ...
                                   sum(outputs.tRIeff(:).*outputs.PRIeff_mech(:)) - outputs.t2*outputs.PRI2_mech)/outputs.tCycle;    
 
                   
  
%    %% Reel-out speed osci due to gravity. 
%    % MK theory of SR oscillation
%      outputs.T_simple   = min(outputs.Tmax_act, (4/9)*outputs.CL^3/outputs.CD^2*(1/2)*outputs.rho_air*...
%                             inputs.WA*inputs.Vw^2);
%     outputs.VSR_simple  = outputs.CD/outputs.CL^(3/2)*sqrt(outputs.T_simple/(0.5*outputs.rho_air*inputs.WA));
%      outputs.VSR_simple2 = outputs.CD/(outputs.CL*cos(outputs.avgRollAngle))^(3/2)*sqrt(outputs.T/(0.5*outputs.rho_air*inputs.WA));
%      outputs.a_simple = sqrt((outputs.lambda^2+outputs.delta^2)/(outputs.T_simple^2+outputs.W^2));
%      outputs.VA_simple = 1/sqrt(outputs.a_simple);
%      outputs.s = outputs.a_simple*(outputs.delta*outputs.T_simple-outputs.lambda*outputs.W)/(outputs.lambda^2+outputs.delta^2);
%     outputs.VSR_down = outputs.VA_simple*outputs.s;
%      outputs.VRO_osciAmp1 = abs(outputs.VSR_down - outputs.VSR_simple);
%      outputs.VRO_osciAmp2 = abs(outputs.VSR - outputs.VSR_simple2);
%      outputs.VRO_osciAmp  = 0;
%      outputs.numPattParts           = 31;
%     for i=1:outputs.numPattParts
%       outputs.VRO_osci(i,i)          = outputs.VRO + outputs.VRO_osciAmp*sin((i-1)*2*pi()/(outputs.numPattParts-1));      
%       outputs.PROeff_mech_osci(i,i) = outputs.T*outputs.VRO_osci(i,i); %[W]
%       % Generator efficiency. As a function of RPM/RPM_max, where RPM_max is driven by winch i.e Max VRI
%       a = 0.671;
%       b = -1.4141;
%       c = 0.9747;
%       d = 0.7233;
%       V_maxGen = max(inputs.maxVRI,25); % 25 =  Possible maximum reel-out speed
%       outputs.genEff_RO(i,i)        = (a*(outputs.VRO_osci(i,i)/V_maxGen)^3+b*(outputs.VRO_osci(i,i)/V_maxGen)^2+c*(outputs.VRO_osci(i,i)/V_maxGen)+d)^sign(1);
%       outputs.PROeff_elec_osci(i,i) = outputs.PROeff_mech_osci(i,i)*inputs.etaGearbox*outputs.genEff_RO(i,i)*inputs.etaPE;
%       outputs.PROeff_elec_osci_cap(i,i) = min(inputs.F_peakM2Ecyc*inputs.P_ratedElec,...
%                                      outputs.PROeff_mech_osci(i,i)*inputs.etaGearbox*outputs.genEff_RO(i,i)*inputs.etaPE);
%     end
%     
%     outputs.PROeff_mech     = mean(outputs.PROeff_mech_osci(:));
%     outputs.PROeff_elec     = mean(outputs.PROeff_elec_osci(:));
%     outputs.PROeff_elec_cap = mean(outputs.PROeff_elec_osci_cap(:));       

end 

% % Air density as a function of height
%     % Ref: https://en.wikipedia.org/wiki/Density_of_air
%     M = 0.0289644; % [kg/mol]
%     R = 8.3144598; % [N·m/(mol·K)]
%     T = 288.15; %[Kelvin]
%     L = 0.0065; %[Kelvin/m] 
%     outputs.rho_air = inputs.airDensity*(1-L*(5+outputs.H_avgPatt)/T)^(inputs.gravity*M/R/L-1); % +5 is for approx. platform height
%     
%     % Roll angle calculation base on new mass, airdensity and pattern radius
%     outputs.avgRollAngle = asin(outputs.m_eff*cos(outputs.pattAngRadius)/...
%                             (0.5*outputs.rho_air*inputs.WA*outputs.CL*outputs.pattRadius));
%     
%     % Tether tension, sinkrate, reel-out speed, PRO_mech 
%     
%     outputs.Tmax_act = inputs.Tmax*inputs.F_Tmax*1000; %[N]
%     outputs.W     = outputs.m_eff*inputs.gravity;
% 
%     % Tether tension: Reduction in lift due to roll
%     outputs.T   = min(outputs.Tmax_act, (4/9)*(outputs.CL*cos(outputs.avgRollAngle))^3/outputs.CD^2*(1/2)*outputs.rho_air*...
%                       inputs.WA*(inputs.Vw*cos(outputs.avgPattEle))^2);
%     
%     % Sink rate: Effect of gravity
%     outputs.i   = sqrt(outputs.T^2 + outputs.W^2 + 2*outputs.T*outputs.W*sin(outputs.avgPattEle)); % Equilibrium force for Aerodynamic resultant
%     outputs.VSR = outputs.CD/(outputs.CL*cos(outputs.avgRollAngle))^(3/2)*...
%                          (outputs.T+outputs.W*sin(outputs.avgPattEle))/sqrt(0.5*outputs.rho_air*inputs.WA*outputs.i);
%          
%     % VRO
%     outputs.VRO = inputs.Vw*cos(outputs.avgPattEle)-outputs.VSR;
%      
%     % PRO without oscillation theory
%     outputs.PROeff_mech = outputs.T*outputs.VRO; %[W]
%     % Generator efficiency. As a function of RPM/RPM_max, where RPM_max is driven by winch i.e Max VRI
%     a = 0.671;
%     b = -1.4141;
%     c = 0.9747;
%     d = 0.7233;
%     V_maxGen = max(inputs.maxVRI,25); % 25 =  Possible maximum reel-out speed
%     outputs.genEff_RO        = (a*(outputs.VRO/V_maxGen)^3+b*(outputs.VRO/V_maxGen)^2+c*(outputs.VRO/V_maxGen)+d)^sign(1);
