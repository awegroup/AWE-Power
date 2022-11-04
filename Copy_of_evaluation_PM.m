function [inputs,outputs] = Copy_of_evaluation_PM(i,inputs,outputs)
  
    %% Wing CL
    outputs.CL             = inputs.CL_maxAirfoil*inputs.F_CLeff;
    
    %% Kite mass estimate: Vincent's simple mass model
    a = 1.747e-2;
    b = 3.84518;
    k1 = 5;
    c = 0.46608;
    d = 0.65962;
    k2 =1.1935;
    outputs.m_kite   = 10*(a*inputs.WA^2 +b*inputs.WA-k1)*(c*(inputs.AR/12)^2-d*(inputs.AR/12)+k2);
    
    %% Minimum Tether length- Centrifugal force balance
    outputs.L_teMin_req    = outputs.m_kite*cos(inputs.pattAngRadius)/...
        (0.5*inputs.airDensity*inputs.WA*outputs.CL*sin(inputs.maxRollAngle)*sin(inputs.pattAngRadius));
    outputs.minPattRad   = outputs.L_teMin_req*sin(inputs.pattAngRadius);
    outputs.L_teMin        = outputs.L_teMin_req*inputs.F_teLength +100; %[m] + 100 additional safety margin
    outputs.L_teMax(i)     = outputs.L_teMin+outputs.deltaL(i);
    outputs.L_teAvg(i)     = (outputs.L_teMax(i)+outputs.L_teMin)/2; %[m]
    outputs.D_te           = sqrt(inputs.Tmax*1000/inputs.Te_matStrength*4/pi())*1.1; %[m] 1.1 is safety factor
    
    %% Swept area 
    outputs.wingSpan = sqrt(inputs.AR*inputs.WA);
    outputs.sweptArea = pi()*((outputs.minPattRad+outputs.wingSpan/2)^2 - (outputs.minPattRad-outputs.wingSpan/2)^2);
    
    %% Effective mass (Kite + tether)
    outputs.m_te(i)  = inputs.Te_matDensity*0.85*pi()/4*outputs.D_te^2*outputs.L_teAvg(i); % 0.85 is safety factor on material density
    outputs.m_eff(i) = outputs.m_kite+sin(inputs.avgPattEle)*outputs.m_te(i);

    %% Wing CD
    outputs.CD(i) = inputs.CD0 + (outputs.CL-inputs.CL0_airfoil)^2/(pi()*inputs.AR*inputs.e) + 0.25*inputs.CD_te...
        *outputs.D_te*outputs.L_teAvg(i)/inputs.WA;

    %% Air density as a function of average pattern altitude
    
    % Ref: https://en.wikipedia.org/wiki/Density_of_air
    M = 0.0289644; % [kg/mol]
    R = 8.3144598; % [N·m/(mol·K)]
    T = 288.15; %[Kelvin]
    L = 0.0065; %[Kelvin/m] 
    outputs.rho_air(i) = inputs.airDensity*(1-L*(5+outputs.L_teAvg(i)*sin(inputs.avgPattEle))/T)^(inputs.gravity*M/R/L-1); % +5 is for approx. platform height

    
    %% Roll angle
    outputs.rollAngle(i) = asin(outputs.m_eff(i)/(tan(inputs.pattAngRadius)*0.5*outputs.rho_air(i)*inputs.WA*outputs.CL*outputs.L_teAvg(i))); %[rad]

    %% Power evaluation

    outputs.Tmax_act = inputs.Tmax*inputs.F_Tmax*1000; %[N]
    outputs.W(i)        = outputs.m_eff(i)*inputs.gravity;

    % Tether tension
    outputs.T(i)   = min(outputs.Tmax_act, (4/9)*(outputs.CL*cos(outputs.rollAngle(i)))^3/outputs.CD(i)^2*(1/2)*outputs.rho_air(i)*...
    inputs.WA*(inputs.Vw(i)*cos(inputs.pattAngRadius))^2);

    % Sink rate
    outputs.J(i)   = sqrt(outputs.T(i)^2+ outputs.W(i)^2+2*outputs.T(i)*outputs.W(i)*sin(inputs.avgPattEle));
    outputs.VSR(i) = outputs.CD(i)/(outputs.CL*cos(outputs.rollAngle(i)))^(3/2)*outputs.T(i)/sqrt(0.5*outputs.rho_air(i)*inputs.WA*outputs.J(i));

    % VRO
    outputs.VRO(i) = inputs.Vw(i)*cos(inputs.avgPattEle)-outputs.VSR(i);
    
    % Reel-out factor
    outputs.reelOutF(i) = outputs.VRO(i)/inputs.Vw(i);
    
    % VRO oscillation due to gravity: Fluctuations in reel-out power
    amp = sqrt(outputs.VRO(i)^2+outputs.minPattRad*inputs.gravity) - outputs.VRO(i);
    for j=1:21
      VRO_osc = outputs.VRO(i) + amp*sin((j-1)*2*pi()/20);
      outputs.VRO_osc(i,j) = VRO_osc;
    end
    
    % PRO mech
    outputs.PRO_mech(i) = outputs.T(i)*outputs.VRO(i); %[W]
    
    % Coefficient of power (Cp): For comparing with wind turbines
    outputs.Cp_reelOut(i) = outputs.PRO_mech(i)/(0.5*outputs.rho_air(i)*outputs.sweptArea*inputs.Vw(i)^3);
    
    % Axial induction factor maximising Cp is a = 0.5 and is independent of reel-out factor
    % Ref paper: The Betz limit applied to Airborne Wind Energy, https://doi.org/10.1016/j.renene.2018.04.034
    f(i) = outputs.reelOutF(i);
    a_ideal = 0.5;
    outputs.Cp_max_ifaIsHalf(i) = (1-f(i))^2*4*a_ideal*(1-a_ideal)*f(i);
    a = abs(roots([1 -1 outputs.Cp_reelOut(i)/(4*f(i)*(1-f(i))^2)]));
    outputs.axialInductionF(i,1) = a(1);
    outputs.axialInductionF(i,2) = a(2);
    
    % Generator efficiency. As a function of RPM/RPM_max, where RPM_max is driven by winch i.e Max VRI
    a = 0.671;
    b = -1.4141;
    c = 0.9747;
    d = 0.7233;
    outputs.genEff_RO(i) = a*(outputs.VRO(i)/inputs.maxVRI)^3+b*(outputs.VRO(i)/inputs.maxVRI)^2+c*(outputs.VRO(i)/inputs.maxVRI)+d;
    outputs.genEff_RI(i) = a*(outputs.VRI(i)/inputs.maxVRI)^3+b*(outputs.VRI(i)/inputs.maxVRI)^2+c*(outputs.VRI(i)/inputs.maxVRI)+d;
    
    % PRO elec
    outputs.PRO_elec(i) = outputs.PRO_mech(i)*inputs.etaGearbox*outputs.genEff_RO(i);
    
    % tRO
    outputs.tRO(i) = (outputs.deltaL(i)/outputs.VRO(i)+outputs.VRO(i)/inputs.maxAcc);
    
    % PRI mech
    outputs.PRI_mech(i) = 0.5*outputs.rho_air(i)*(inputs.CD0+((2*outputs.W(i)/outputs.rho_air(i)/(inputs.Vw(i)+outputs.VRI(i))^2/inputs.WA)...
        -inputs.CL0_airfoil)^2/pi()/inputs.AR/inputs.e)*inputs.WA*(inputs.Vw(i)*cos(inputs.avgPattEle)+outputs.VRI(i))^3;
    
    % PRI elec
    outputs.PRI_elec(i) = outputs.PRI_mech(i)/inputs.etaGearbox/inputs.etaSto/outputs.genEff_RI(i);
    
    % tRI
    outputs.tRI(i) = outputs.deltaL(i)/outputs.VRI(i)+outputs.VRI(i)/inputs.maxAcc;
    
    % tCycle
    outputs.tCycle(i) = outputs.tRO(i)+outputs.tRI(i);
    
    %P_cycleElec 
     outputs.P_cycleElec(i) = (outputs.tRO(i)*outputs.PRO_elec(i) - outputs.tRI(i)*outputs.PRI_elec(i))/outputs.tCycle(i)+inputs.etaGrid;    
    
    % Without drivetrain eff
    outputs.P_cycleMech(i) = (outputs.tRO(i)*outputs.PRO_mech(i) - outputs.tRI(i)*outputs.PRI_mech(i))/(outputs.tRO(i)+outputs.tRI(i));    




end