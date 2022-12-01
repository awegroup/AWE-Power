function [inputs] = compute(i,inputs)
 global outputs
    
    %% Kite mass 
    if inputs.massOverride == 1
      outputs.m_kite = inputs.kiteMass;
    else
      %Vincent's simple mass model
      a = 1.747e-2;
      b = 3.84518;
      k1 = 5;
      c = 0.46608;
      d = 0.65962;
      k2 =1.1935;
      outputs.m_kite = 10*(a*inputs.WA^2 +b*inputs.WA-k1)*(c*(inputs.AR/12)^2-d*(inputs.AR/12)+k2); 
    end
    
    %% Minimum Tether length and minimum patt. radius - Centrifugal force balance - no gravity
    outputs.L_teMin_req(i) = outputs.m_kite*cos(inputs.pattAngRadius)/...
        (0.5*inputs.airDensity*inputs.WA*outputs.CL(i)*sin(inputs.maxRollAngle)*sin(inputs.pattAngRadius));
    outputs.minPattRad(i)  = outputs.L_teMin_req(i)*sin(inputs.pattAngRadius);
    outputs.L_teMin(i)     = outputs.L_teMin_req(i)*inputs.F_teLength; %[m] could add 100 additional safety margin
    outputs.L_teMax(i)     = outputs.L_teMin(i)+outputs.deltaL(i); 
    outputs.L_teAvg(i)     = (outputs.L_teMax(i)+outputs.L_teMin(i))/2; %[m]
    outputs.pattRadius(i)  = outputs.L_teAvg(i)*sin(inputs.pattAngRadius); %[m]
    outputs.H_avg(i)       = outputs.L_teAvg(i)*sin(inputs.avgPattEle);
    outputs.D_te           = sqrt(inputs.Tmax*1000/inputs.Te_matStrength*4/pi())*1.1; %[m] 1.1 is safety factor
    
    %% Effective mass (Kite + tether)
    outputs.m_te(i)  = inputs.Te_matDensity*0.85*pi()/4*outputs.D_te^2*outputs.L_teAvg(i); % 0.85 is safety factor on material density
    outputs.m_eff(i) = outputs.m_kite+sin(inputs.avgPattEle)*outputs.m_te(i);

    %% Effective CD
    outputs.CD_kite(i)   = inputs.CD0 + (outputs.CL(i)-inputs.CL0_airfoil)^2/(pi()*inputs.AR*inputs.e);
    outputs.CD_tether(i) = (1/4)*inputs.CD_te*outputs.D_te*outputs.L_teAvg(i)/inputs.WA;
    outputs.CD(i)        = outputs.CD_kite(i) + outputs.CD_tether(i);

    %% Air density as a function of average pattern altitude
    % Ref: https://en.wikipedia.org/wiki/Density_of_air
    M = 0.0289644; % [kg/mol]
    R = 8.3144598; % [N·m/(mol·K)]
    T = 288.15; %[Kelvin]
    L = 0.0065; %[Kelvin/m] 
    outputs.rho_air(i) = inputs.airDensity*(1-L*(5+outputs.H_avg(i))/T)^(inputs.gravity*M/R/L-1); % +5 is for approx. platform height
    
    %% Roll angle calculation base on new mass, airdensity and pattern radius
    outputs.rollAngle(i) = asin(outputs.m_kite*cos(inputs.pattAngRadius)/...
        (0.5*inputs.airDensity*inputs.WA*outputs.CL(i)*outputs.pattRadius(i)));
    
    %% Tether tension, sinkrate, reel-out speed

    outputs.Tmax_act = inputs.Tmax*inputs.F_Tmax*1000; %[N]
    outputs.W(i)     = outputs.m_eff(i)*inputs.gravity;

    % Tether tension
    outputs.T(i)   = min(outputs.Tmax_act, (4/9)*outputs.CL(i)^3*cos(outputs.rollAngle(i))/outputs.CD(i)^2*(1/2)*outputs.rho_air(i)*...
                      inputs.WA*(inputs.Vw(i)*cos(outputs.rollAngle(i)-inputs.pattAngRadius))^2);

    % Sink rate
    outputs.J(i)   = sqrt(outputs.T(i)^2+ outputs.W(i)^2+2*outputs.T(i)*outputs.W(i)*sin(inputs.avgPattEle));
    outputs.VSR(i) = outputs.CD(i)/outputs.CL(i)^(3/2)*cos(outputs.rollAngle(i))*outputs.T(i)/sqrt(0.5*outputs.rho_air(i)*inputs.WA*outputs.J(i));
     
    % VRO
    outputs.VRO(i) = inputs.Vw(i)*cos(outputs.rollAngle(i)-inputs.avgPattEle)-outputs.VSR(i);
    
    % Reel-out factor
    outputs.reelOutF(i) = outputs.VRO(i)/inputs.Vw(i);
    
    %% Airspeed and tangential speed
     lambda        = 0.5*outputs.rho_air(i)*inputs.WA*outputs.CL(i);
     delta         = 0.5*outputs.rho_air(i)*inputs.WA*outputs.CD(i);
     a             = sqrt((lambda^2+delta^2)/(outputs.T(i)^2+outputs.W(i)^2+2*outputs.T(i)*outputs.W(i)*sin(inputs.avgPattEle)));
     outputs.VA(i) = 1/sqrt(a); % Apparent wind speed [m/s];
     outputs.VC(i) = sqrt(outputs.VA(i)^2-outputs.VSR(i)^2); 
   
    %% Reel-out speed osci due to gravity. Since Tether has a tension limit the osci cannot be translated to TT 
     % Theory of centripetal force with gravity. V_C osci -> V_A osci -> VRO osci 
     % VC oscillation amplitude
     outputs.VC_osciAmp(i)          = sqrt(outputs.VC(i)^2+outputs.pattRadius(i)*inputs.gravity) - outputs.VC(i);
     outputs.numPattParts           = 31;
    for j=1:outputs.numPattParts
      outputs.VC_osci(i,j)          = outputs.VC(i) + outputs.VC_osciAmp(i)*sin((j-1)*2*pi()/(outputs.numPattParts-1));
      outputs.VA_osci(i,j)          = sqrt(outputs.VC_osci(i,j)^2 + outputs.VSR(i)^2); 
      outputs.osciFactor(i,j)       = (outputs.VA_osci(i,j)/outputs.VA(i))^2*cos(outputs.rollAngle(i));
      outputs.VRO_osci(i,j)         = outputs.osciFactor(i,j)*outputs.VRO(i);
      outputs.PROeff_mech_osci(i,j) = outputs.T(i)*outputs.VRO_osci(i,j); %[W]
      % Generator efficiency. As a function of RPM/RPM_max, where RPM_max is driven by winch i.e Max VRI
      a = 0.671;
      b = -1.4141;
      c = 0.9747;
      d = 0.7233;
      RPM_max = max(inputs.maxVRI,25); % 25 =  Possible maximum reel-out speed
      outputs.genEff_RO(i,j)        = (a*(outputs.VRO_osci(i,j)/RPM_max)^3+b*(outputs.VRO_osci(i,j)/RPM_max)^2+c*(outputs.VRO_osci(i,j)/RPM_max)+d)^sign(1);
      outputs.PROeff_elec_osci(i,j) = outputs.PROeff_mech_osci(i,j)*inputs.etaGearbox*outputs.genEff_RO(i,j)*inputs.etaPE;
    end
    outputs.PROeff_mech(i) = mean(outputs.PROeff_mech_osci(i,:));
    outputs.PROeff_elec(i) = mean(outputs.PROeff_elec_osci(i,:));
    
    %% Cycle simulation
    
    % tRO
    outputs.t1(i)     = outputs.VRO(i)/inputs.maxAcc;
    outputs.tROeff(i) = outputs.deltaL(i)/outputs.VRO(i);
    outputs.tRO(i)    = outputs.t1(i)+outputs.tROeff(i);
    
    % Reel-out power during transition
    outputs.PRO1_mech(i) = outputs.PROeff_mech(i)/2;
    outputs.PRO1_elec(i) = outputs.PROeff_elec(i)/2;
    
    % PRO 
    outputs.PRO_mech(i) = (outputs.PROeff_mech(i)* outputs.tROeff(i) + outputs.PRO1_mech(i)*outputs.t1(i))/outputs.tRO(i);
    outputs.PRO_elec(i) = (outputs.PROeff_elec(i)* outputs.tROeff(i) + outputs.PRO1_elec(i)*outputs.t1(i))/outputs.tRO(i);
    
    % PRI effective mech
    outputs.VAS_RI(i)      = inputs.Vw(i)*cos(inputs.avgPattEle)+outputs.VRI(i);
    outputs.CL_RI(i)       = 2*outputs.W(i)/(outputs.rho_air(i)*outputs.VAS_RI(i)^2*inputs.WA);
    outputs.CD_RI(i)       = inputs.CD0+(outputs.CL_RI(i)- inputs.CL0_airfoil)^2/(pi()*inputs.AR*inputs.e);
    outputs.PRIeff_mech(i) = 0.5*outputs.rho_air(i)*outputs.CD_RI(i)*inputs.WA*outputs.VAS_RI(i)^3;
    
    % Generator efficiency. As a function of RPM/RPM_max, where RPM_max is driven by winch i.e Max VRI
    outputs.genEff_RI(i)   = a*(outputs.VRI(i)/inputs.maxVRI)^3+b*(outputs.VRI(i)/inputs.maxVRI)^2+c*(outputs.VRI(i)/inputs.maxVRI)+d;
      
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
    
    % P_cycleElec 
     outputs.P_cycleElec(i) = (outputs.tROeff(i)*outputs.PROeff_elec(i)+outputs.t1(i)*outputs.PRO1_elec(i) - ...
                     outputs.tRIeff(i)*outputs.PRIeff_elec(i)-outputs.t2(i)*outputs.PRI2_elec(i))/outputs.tCycle(i); 
     
    % Without drivetrain eff
    outputs.P_cycleMech(i) = (outputs.tROeff(i)*outputs.PROeff_mech(i)+outputs.t1(i)*outputs.PRO1_mech(i) - ...
                     outputs.tRIeff(i)*outputs.PRIeff_mech(i)-outputs.t2(i)*outputs.PRI2_mech(i))/outputs.tCycle(i);    

    %% Coefficient of power (Cp): For comparing with wind turbines
    % Swept area 
    outputs.wingSpan      = sqrt(inputs.AR*inputs.WA);
    outputs.sweptArea(i)  = pi()*((outputs.pattRadius(i)+outputs.wingSpan/2)^2 - (outputs.pattRadius(i)-outputs.wingSpan/2)^2);
    outputs.Cp_reelOut(i) = outputs.PRO_mech(i)/(0.5*outputs.rho_air(i)*outputs.sweptArea(i)*inputs.Vw(i)^3);
    % Axial induction factor maximising Cp is a = 0.5 and is independent of reel-out factor
    % Ref paper: The Betz limit applied to Airborne Wind Energy, https://doi.org/10.1016/j.renene.2018.04.034
    f(i)    = outputs.reelOutF(i);
    a_ideal = 0.5;
    outputs.Cp_max_ifaIsHalf(i)  = (1-f(i))^2*4*a_ideal*(1-a_ideal)*f(i);
    a                            = abs(roots([1 -1 outputs.Cp_reelOut(i)/(4*f(i)*(1-f(i))^2)]));
    outputs.axialInductionF(i,1) = a(1);
    outputs.axialInductionF(i,2) = a(2);
    
    %% Without considering oscillations
%     %outputs.PRO_mech(i) = outputs.T(i)*outputs.VRO(i); %[W]
%     % Generator efficiency. As a function of RPM/RPM_max, where RPM_max is driven by winch i.e Max VRI
%     a = 0.671;
%     b = -1.4141;
%     c = 0.9747;
%     d = 0.7233;
%     outputs.genEff_RO(i) = a*(outputs.VRO(i)/inputs.maxVRI)^3+b*(outputs.VRO(i)/inputs.maxVRI)^2+c*(outputs.VRO(i)/inputs.maxVRI)+d;
%     outputs.genEff_RI(i) = a*(outputs.VRI(i)/inputs.maxVRI)^3+b*(outputs.VRI(i)/inputs.maxVRI)^2+c*(outputs.VRI(i)/inputs.maxVRI)+d;
%     % PRO elec
%     outputs.PRO_elec(i) = outputs.PRO_mech(i)*inputs.etaGearbox*outputs.genEff_RO(i);

end