function [inputs] = compute(i,inputs)
 global outputs
    
    %% Kite mass 
    if inputs.massOverride == 1
      outputs.m_kite = inputs.kiteMass;
    else
      %Vincent's simple mass model
%       a = 1.747e-2;
%       b = 3.84518;
%       k1 = 5;
%       c = 0.46608;
%       d = 0.65962;
%       k2 =1.1935;
%       outputs.m_kite = 10*(a*inputs.WA^2 +b*inputs.WA-k1)*(c*(inputs.AR/12)^2-d*(inputs.AR/12)+k2); 
      
      a1     = 0.002415;
      a2     = 0.0090239;
      b1     = 0.17025;
      b2     = 3.2493;
      k1     = 5;
      c1     = 0.46608;
      d1     = 0.65962;
      k2     = 1.1935;
      AR_ref = 12;
      
      a = a1*(inputs.Tmax/inputs.WA) + a2;
      b = b1*(inputs.Tmax/inputs.WA) + b2;
      
      outputs.m_kite = 10*(a*inputs.WA^2 +b*inputs.WA-k1)*(c1*(inputs.AR/AR_ref)^2-d1*(inputs.AR/AR_ref)+k2); 
      
    end
    
    %% Minimum Tether length and minimum patt. radius - Centrifugal force balance - no gravity
    outputs.wingSpan           = sqrt(inputs.AR*inputs.WA);
    %outputs.pattRadius(i)      = outputs.m_kite*cos(outputs.pattAngRadius(i))/...
                            %     (0.5*inputs.airDensity*inputs.WA*outputs.CL(i)*sin(outputs.maxRollAngle(i)));
    outputs.sweptArea(i)       = pi()*((outputs.pattRadius(i)+outputs.wingSpan/2)^2 - (outputs.pattRadius(i)-outputs.wingSpan/2)^2); 
    outputs.L_teMin(i)         = outputs.pattRadius(i)/sin(outputs.pattAngRadius(i))*inputs.F_minTeLen;
    outputs.pattStartGrClr(i)  = outputs.L_teMin(i)*sin(outputs.avgPattEle(i)-outputs.pattAngRadius(i));
    outputs.H_minPatt(i)       = outputs.L_teMin(i)*cos(outputs.pattAngRadius(i))*sin(outputs.avgPattEle(i));
    outputs.L_teMax(i)         = outputs.L_teMin(i)+outputs.deltaL(i); 
    outputs.L_teAvg(i)         = (outputs.L_teMax(i)+outputs.L_teMin(i))/2; %[m]
    outputs.H_avgPatt(i)       = outputs.L_teAvg(i)*cos(outputs.pattAngRadius(i))*sin(outputs.avgPattEle(i));
    outputs.D_te               = sqrt(inputs.Tmax*1000/inputs.Te_matStrength*4/pi()); %[m] safety factor could be added (say *1.1)
    
    %% Effective mass (Kite + tether)
    outputs.m_te(i)  = inputs.Te_matDensity*pi()/4*outputs.D_te^2*outputs.L_teAvg(i); % Could add (say 0.85) as safety factor on material density
    outputs.m_eff(i) = outputs.m_kite+sin(outputs.avgPattEle(i))*outputs.m_te(i);

    %% Effective CD
    outputs.CD_kite(i)   = inputs.CD0 + (outputs.CL(i)-inputs.CL0_airfoil)^2/(pi()*inputs.AR*inputs.e);
    outputs.CD_tether(i) = (1/4)*inputs.CD_te*outputs.D_te*outputs.L_teAvg(i)/inputs.WA;
    outputs.CD(i)        = outputs.CD_kite(i) + outputs.CD_tether(i);

    %% Pattern avg. mech reel-out power considering vertical wind shear
    
%   function [intermRes] = pattAvgPROmech(i, inputs)
%     
%     
%   end
    
   % u(ind) = u_60 .* (z(ind)./60).^alpha;
    
    % Air density as a function of height
    % Ref: https://en.wikipedia.org/wiki/Density_of_air
    M = 0.0289644; % [kg/mol]
    R = 8.3144598; % [N·m/(mol·K)]
    T = 288.15; %[Kelvin]
    L = 0.0065; %[Kelvin/m] 
    outputs.rho_air(i) = inputs.airDensity*(1-L*(5+outputs.H_avgPatt(i))/T)^(inputs.gravity*M/R/L-1); % +5 is for approx. platform height
    
    % Roll angle calculation base on new mass, airdensity and pattern radius
    outputs.avgRollAngle(i) = asin(outputs.m_eff(i)*cos(outputs.pattAngRadius(i))/...
                            (0.5*outputs.rho_air(i)*inputs.WA*outputs.CL(i)*outputs.pattRadius(i)));
    
    % Tether tension, sinkrate, reel-out speed, PRO_mech 
    
    outputs.Tmax_act = inputs.Tmax*inputs.F_Tmax*1000; %[N]
    outputs.W(i)     = outputs.m_eff(i)*inputs.gravity;

    % Tether tension: Reduction in lift due to roll
    outputs.T(i)   = min(outputs.Tmax_act, (4/9)*(outputs.CL(i)*cos(outputs.avgRollAngle(i)))^3/outputs.CD(i)^2*(1/2)*outputs.rho_air(i)*...
                      inputs.WA*(inputs.Vw(i)*cos(outputs.avgPattEle(i)))^2);
    
    % Sink rate: Effect of gravity
    outputs.J(i)   = sqrt(outputs.T(i)^2 + outputs.W(i)^2 + 2*outputs.T(i)*outputs.W(i)*sin(outputs.avgPattEle(i))); % Equilibrium force for Aerodynamic resultant
    outputs.VSR(i) = outputs.CD(i)/(outputs.CL(i)*cos(outputs.avgRollAngle(i)))^(3/2)*...
                         (outputs.T(i)+outputs.W(i)*sin(outputs.avgPattEle(i)))/sqrt(0.5*outputs.rho_air(i)*inputs.WA*outputs.J(i));
         
    % VRO
    outputs.VRO(i) = inputs.Vw(i)*cos(outputs.avgPattEle(i))-outputs.VSR(i);
     
    % PRO without oscillation theory
    outputs.PROeff_mech(i) = outputs.T(i)*outputs.VRO(i); %[W]
    % Generator efficiency. As a function of RPM/RPM_max, where RPM_max is driven by winch i.e Max VRI
    a = 0.671;
    b = -1.4141;
    c = 0.9747;
    d = 0.7233;
    V_maxGen = max(inputs.maxVRI,25); % 25 =  Possible maximum reel-out speed
    outputs.genEff_RO(i)        = (a*(outputs.VRO(i)/V_maxGen)^3+b*(outputs.VRO(i)/V_maxGen)^2+c*(outputs.VRO(i)/V_maxGen)+d)^sign(1);
    
    
    %% PRO elec
    outputs.PROeff_elec(i) = outputs.PROeff_mech(i)*inputs.etaGearbox*outputs.genEff_RO(i)*inputs.etaPE;
       
    % Airspeed and tangential speed
     outputs.lambda(i) = 0.5*outputs.rho_air(i)*inputs.WA*outputs.CL(i);
     outputs.delta(i)  = 0.5*outputs.rho_air(i)*inputs.WA*outputs.CD(i);
     outputs.a(i)      = sqrt((outputs.lambda(i)^2+outputs.delta(i)^2))/outputs.J(i);
     outputs.VA(i) = 1/sqrt(outputs.a(i)); % Apparent speed [m/s];
     outputs.VC(i) = sqrt(outputs.VA(i)^2 - outputs.VSR(i)^2);
        
    %% Cycle simulation
       
    % tRO
    outputs.t1(i)     = outputs.VRO(i)/inputs.maxAcc;
    outputs.tROeff(i) = outputs.deltaL(i)/outputs.VRO(i);
    outputs.tRO(i)    = outputs.t1(i)+outputs.tROeff(i);
    if outputs.VRO(i)<0
      outputs.t1(i)           = 0;
      outputs.tROeff(i)       = 0;
      outputs.tRO(i)          = 0;
    end
    
    % Reel-out power during transition
    outputs.PRO1_mech(i) = outputs.PROeff_mech(i)/2;
    outputs.PRO1_elec(i) = outputs.PROeff_elec(i)/2;
    
    % PRO 
    outputs.PRO_mech(i) = (outputs.PROeff_mech(i)* outputs.tROeff(i) + outputs.PRO1_mech(i)*outputs.t1(i))/outputs.tRO(i);
    outputs.PRO_elec(i) = (outputs.PROeff_elec(i)* outputs.tROeff(i) + outputs.PRO1_elec(i)*outputs.t1(i))/outputs.tRO(i);
    if outputs.VRO(i)<0
      outputs.PRO_mech(i) = 1e-9;
    end
    
    % PRI effective mech
    
    % UPDATE Vw
    
    outputs.VA_RI(i)      = sqrt(inputs.Vw(i)^2 +outputs.VRI(i)^2 +2*inputs.Vw(i)*outputs.VRI(i)*cos(outputs.avgPattEle(i)));
    outputs.CL_RI(i)       = 2*outputs.W(i)/(outputs.rho_air(i)*outputs.VA_RI(i)^2*inputs.WA);
    outputs.CD_RI(i)       = inputs.CD0+(outputs.CL_RI(i)- inputs.CL0_airfoil)^2/(pi()*inputs.AR*inputs.e);
    outputs.PRIeff_mech(i) = 0.5*outputs.rho_air(i)*outputs.CD_RI(i)*inputs.WA*outputs.VA_RI(i)^3;
    
    % Generator efficiency during RI: As a function of RPM/RPM_max, where RPM_max is driven by winch i.e Max VRI
    a = 0.671;
    b = -1.4141;
    c = 0.9747;
    d = 0.7233;
    %RPM_max = max(inputs.maxVRI,25); % 25 =  Possible maximum reel-out speed
    outputs.genEff_RI(i)        = (a*(outputs.VRI(i)/V_maxGen)^3+b*(outputs.VRI(i)/V_maxGen)^2+c*(outputs.VRI(i)/V_maxGen)+d)^sign(1);
      
    % PRI effective elec
    outputs.PRIeff_elec(i) = outputs.PRIeff_mech(i)/inputs.etaGearbox/inputs.etaSto/outputs.genEff_RI(i)/inputs.etaPE;
    
    % tRI
    outputs.t2(i)     = outputs.VRI(i)/inputs.maxAcc;
    outputs.tRIeff(i) = outputs.deltaL(i)/outputs.VRI(i);
    outputs.tRI(i)    = outputs.t2(i)+outputs.tRIeff(i);
   
    % Transition reel-in mech power
    outputs.PRI2_mech(i)     = outputs.PRIeff_mech(i)/2;
    outputs.PRI2_elec(i)     = outputs.PRIeff_elec(i)/2;
    
    % PRI 
    outputs.PRI_mech(i) = (outputs.PRIeff_mech(i)* outputs.tRIeff(i) + outputs.PRI2_mech(i)*outputs.t2(i))/outputs.tRI(i);
    outputs.PRI_elec(i) = (outputs.PRIeff_elec(i)* outputs.tRIeff(i) + outputs.PRI2_elec(i)*outputs.t2(i))/outputs.tRI(i);
    
    % tCycle
    outputs.tCycle(i) = outputs.tRO(i)+outputs.tRI(i);
    
    outputs.tPatt(i)     = 2*pi()*outputs.pattRadius(i)/outputs.VC(i);
    outputs.numOfPatt(i) = outputs.tRO(i)/outputs.tPatt(i);
    
    % P_cycleElec 
     outputs.P_cycleElec(i) = (outputs.tROeff(i)*outputs.PROeff_elec(i)+outputs.t1(i)*outputs.PRO1_elec(i) - ...
                     outputs.tRIeff(i)*outputs.PRIeff_elec(i)-outputs.t2(i)*outputs.PRI2_elec(i))/outputs.tCycle(i);    
    if outputs.VRO(i)<0
      outputs.P_cycleElec(i) = 0;
    end
     
    % Without drivetrain eff
    outputs.P_cycleMech(i) = (outputs.tROeff(i)*outputs.PROeff_mech(i)+outputs.t1(i)*outputs.PRO1_mech(i) - ...
                     outputs.tRIeff(i)*outputs.PRIeff_mech(i)-outputs.t2(i)*outputs.PRI2_mech(i))/outputs.tCycle(i);    
 
%    %% Reel-out speed osci due to gravity. 
%    % MK theory of SR oscillation
%      outputs.T_simple(i)   = min(outputs.Tmax_act, (4/9)*outputs.CL(i)^3/outputs.CD(i)^2*(1/2)*outputs.rho_air(i)*...
%                             inputs.WA*inputs.Vw(i)^2);
%     outputs.VSR_simple(i)  = outputs.CD(i)/outputs.CL(i)^(3/2)*sqrt(outputs.T_simple(i)/(0.5*outputs.rho_air(i)*inputs.WA));
%      outputs.VSR_simple2(i) = outputs.CD(i)/(outputs.CL(i)*cos(outputs.avgRollAngle(i)))^(3/2)*sqrt(outputs.T(i)/(0.5*outputs.rho_air(i)*inputs.WA));
%      outputs.a_simple(i) = sqrt((outputs.lambda(i)^2+outputs.delta(i)^2)/(outputs.T_simple(i)^2+outputs.W(i)^2));
%      outputs.VA_simple(i) = 1/sqrt(outputs.a_simple(i));
%      outputs.s(i) = outputs.a_simple(i)*(outputs.delta(i)*outputs.T_simple(i)-outputs.lambda(i)*outputs.W(i))/(outputs.lambda(i)^2+outputs.delta(i)^2);
%     outputs.VSR_down(i) = outputs.VA_simple(i)*outputs.s(i);
%      outputs.VRO_osciAmp1(i) = abs(outputs.VSR_down(i) - outputs.VSR_simple(i));
%      outputs.VRO_osciAmp2(i) = abs(outputs.VSR(i) - outputs.VSR_simple2(i));
%      outputs.VRO_osciAmp(i)  = 0;
%      outputs.numPattParts           = 31;
%     for j=1:outputs.numPattParts
%       outputs.VRO_osci(i,j)          = outputs.VRO(i) + outputs.VRO_osciAmp(i)*sin((j-1)*2*pi()/(outputs.numPattParts-1));      
%       outputs.PROeff_mech_osci(i,j) = outputs.T(i)*outputs.VRO_osci(i,j); %[W]
%       % Generator efficiency. As a function of RPM/RPM_max, where RPM_max is driven by winch i.e Max VRI
%       a = 0.671;
%       b = -1.4141;
%       c = 0.9747;
%       d = 0.7233;
%       V_maxGen = max(inputs.maxVRI,25); % 25 =  Possible maximum reel-out speed
%       outputs.genEff_RO(i,j)        = (a*(outputs.VRO_osci(i,j)/V_maxGen)^3+b*(outputs.VRO_osci(i,j)/V_maxGen)^2+c*(outputs.VRO_osci(i,j)/V_maxGen)+d)^sign(1);
%       outputs.PROeff_elec_osci(i,j) = outputs.PROeff_mech_osci(i,j)*inputs.etaGearbox*outputs.genEff_RO(i,j)*inputs.etaPE;
%       outputs.PROeff_elec_osci_cap(i,j) = min(inputs.F_peakM2Ecyc*inputs.P_ratedElec,...
%                                      outputs.PROeff_mech_osci(i,j)*inputs.etaGearbox*outputs.genEff_RO(i,j)*inputs.etaPE);
%     end
%     
%     outputs.PROeff_mech(i)     = mean(outputs.PROeff_mech_osci(i,:));
%     outputs.PROeff_elec(i)     = mean(outputs.PROeff_elec_osci(i,:));
%     outputs.PROeff_elec_cap(i) = mean(outputs.PROeff_elec_osci_cap(i,:));       

end 