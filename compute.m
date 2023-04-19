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
    outputs.L_teMin(i)         = outputs.startPattRadius(i)/sin(outputs.pattAngRadius(i))*inputs.F_minTeLen;
    outputs.pattStartGrClr(i)  = outputs.L_teMin(i)*sin(outputs.avgPattEle(i)-outputs.pattAngRadius(i));
    outputs.H_cycleStart(i)    = outputs.L_teMin(i)*cos(outputs.pattAngRadius(i))*sin(outputs.avgPattEle(i));
    outputs.L_teMax(i)         = outputs.L_teMin(i)+outputs.deltaL(i)/cos(outputs.pattAngRadius(i)); 
    outputs.pattEndGrClr(i)    = outputs.L_teMax(i)*sin(outputs.avgPattEle(i)-outputs.pattAngRadius(i));
    outputs.L_teAvg(i)         = (outputs.L_teMax(i)+outputs.L_teMin(i))/2; %[m]
    outputs.H_cycleAvg(i)      = outputs.L_teAvg(i)*cos(outputs.pattAngRadius(i))*sin(outputs.avgPattEle(i));
    outputs.H_cycleEnd(i)      = outputs.L_teMax(i)*cos(outputs.pattAngRadius(i))*sin(outputs.avgPattEle(i));
    outputs.D_te               = sqrt(inputs.Tmax*1000/inputs.Te_matStrength*4/pi()); %[m] safety factor could be added (say *1.1)
    
    %% Effective mass (Kite + tether)
    outputs.m_te(i)  = inputs.Te_matDensity*pi()/4*outputs.D_te^2*outputs.L_teAvg(i); % Could add (say 0.85) as safety factor on material density
    outputs.m_eff(i) = outputs.m_kite+sin(outputs.avgPattEle(i))*outputs.m_te(i);

    %% Effective CD
    outputs.CD_kite(i)   = inputs.CD0 + (outputs.CL(i)-inputs.CL0_airfoil)^2/(pi()*inputs.AR*inputs.e);
    outputs.CD_tether(i) = (1/4)*inputs.CD_te*outputs.D_te*outputs.L_teAvg(i)/inputs.WA;
    outputs.CD(i)        = outputs.CD_kite(i) + outputs.CD_tether(i);

    %% Tmax and W
    outputs.Tmax_act = inputs.Tmax*inputs.F_Tmax*1000; %[N]
    outputs.W(i)     = outputs.m_eff(i)*inputs.gravity;
    
    %% Cycle avg. mech reel-out power considering vertical wind shear
     
    % Divide the reel-out length in equal number of elements of around 5m lengths
     outputs.deltaLelems = 30; % Found to be not sensitive to the number of elements once element size ~<=1.5m
     outputs.elemDeltaL(i) = outputs.deltaL(i)/outputs.deltaLelems;
  
     % For each element of deltaL
     for j = 1:outputs.deltaLelems
       if j == 1
         outputs.h_inCycle(i,j) = outputs.H_cycleStart(i) + outputs.elemDeltaL(i)/2*sin(outputs.avgPattEle(i));
       else
         outputs.h_inCycle(i,j) = outputs.h_inCycle(i,j-1) + outputs.elemDeltaL(i)*sin(outputs.avgPattEle(i));
       end
       
       outputs.Vw(i,j) = inputs.Vw_ref(i)*(outputs.h_inCycle(i,j)/inputs.h_ref)^inputs.windShearExp;
        % Air density as a function of height   % Ref: https://en.wikipedia.org/wiki/Density_of_air
        M = 0.0289644; % [kg/mol]
        R = 8.3144598; % [N·m/(mol·K)]
        T = 288.15; %[Kelvin]
        L = 0.0065; %[Kelvin/m] 
        outputs.rho_air(i,j) = inputs.airDensity*(1-L*(5+ outputs.h_inCycle(i,j))/T)^(inputs.gravity*M/R/L-1); % +5 is for approx. platform height

        % Pattern radius at point of interest on deltaL
        if j == 1
          outputs.pattRadius(i,j) = (outputs.startPattRadius(i) + (outputs.elemDeltaL(i)*tan(outputs.pattAngRadius(i)) + outputs.startPattRadius(i)))/2;
        else
          outputs.pattRadius(i,j) = (outputs.pattRadius(i,j-1) + (j*outputs.elemDeltaL(i)*tan(outputs.pattAngRadius(i)) + outputs.startPattRadius(i)))/2;
        end

        % Roll angle calculation for the particular element
        outputs.avgRollAngle(i,j) = asin(outputs.m_eff(i)*cos(outputs.pattAngRadius(i))/...
                                (0.5*outputs.rho_air(i,j)*inputs.WA*outputs.CL(i)*outputs.pattRadius(i,j)));

        % Tether tension: Reduction in lift due to roll
        outputs.T(i,j)   = min(outputs.Tmax_act, (4/9)*(outputs.CL(i)*cos(outputs.avgRollAngle(i,j)))^3/outputs.CD(i)^2*(1/2).*outputs.rho_air(i,j)*...
                          inputs.WA*(outputs.Vw(i,j).*cos(outputs.avgPattEle(i)))^2);

        % Sink rate: Effect of gravity
        outputs.J(i,j)   = sqrt(outputs.T(i,j)^2 + outputs.W(i)^2 + 2*outputs.T(i,j)*outputs.W(i)*sin(outputs.avgPattEle(i))); % Equilibrium force for Aerodynamic resultant
        outputs.VSR(i,j) = outputs.CD(i)/(outputs.CL(i)*cos(outputs.avgRollAngle(i,j)))^(3/2)*...
                             (outputs.T(i,j)+outputs.W(i)*sin(outputs.avgPattEle(i)))/sqrt(0.5*outputs.rho_air(i,j)*inputs.WA*outputs.J(i,j));

        % Airspeed and tangential speed
         outputs.lambda(i,j) = 0.5*outputs.rho_air(i,j)*inputs.WA*outputs.CL(i);
         outputs.delta(i,j)  = 0.5*outputs.rho_air(i,j)*inputs.WA*outputs.CD(i);
         outputs.a(i,j)      = sqrt((outputs.lambda(i,j)^2+outputs.delta(i,j)^2))/outputs.J(i,j);
         outputs.VA(i,j)     = 1/sqrt(outputs.a(i,j)); % Apparent speed [m/s];
         outputs.VC(i,j)     = sqrt(outputs.VA(i,j)^2 - outputs.VSR(i,j)^2);
                           
        % VRO
        outputs.VRO(i,j) = outputs.Vw(i,j)*cos(outputs.avgPattEle(i))-outputs.VSR(i,j);

        % PRO 
        outputs.PROeff_mech(i,j) = outputs.T(i,j)*outputs.VRO(i,j); %[W]

        % PRO elec
        % Generator efficiency. As a function of RPM/RPM_max, where RPM_max is driven by winch i.e Max VRI
        outputs.genEff_RO(i,j)  = (inputs.etaGen.param(1)*(outputs.VRO(i,j)/inputs.etaGen.Vmax)^3 + ...
                                      inputs.etaGen.param(2)*(outputs.VRO(i,j)/inputs.etaGen.Vmax)^2 + ...
                                        inputs.etaGen.param(3)*(outputs.VRO(i,j)/inputs.etaGen.Vmax)+inputs.etaGen.param(4))^sign(1);
        outputs.PROeff_elec(i,j) = outputs.PROeff_mech(i,j)*inputs.etaGearbox*outputs.genEff_RO(i,j)*inputs.etaPE;

        % PRI effective mech
        outputs.VA_RI(i,j)       = sqrt(outputs.Vw(i,j)^2 +outputs.VRI(i)^2 +2*outputs.Vw(i,j)*outputs.VRI(i)*cos(outputs.avgPattEle(i)));
        outputs.CL_RI(i,j)       = 2*outputs.W(i)/(outputs.rho_air(i,j)*outputs.VA_RI(i,j)^2*inputs.WA);
        outputs.CD_RI(i,j)       = inputs.CD0+(outputs.CL_RI(i,j)- inputs.CL0_airfoil)^2/(pi()*inputs.AR*inputs.e) + outputs.CD_tether(i);

        outputs.PRIeff_mech(i,j) = 0.5*outputs.rho_air(i,j)*outputs.CD_RI(i,j)*inputs.WA*outputs.VA_RI(i,j)^3;

        % Generator efficiency during RI: As a function of RPM/RPM_max, where RPM_max is driven by winch i.e Max VRI
        outputs.genEff_RI(i) = (inputs.etaGen.param(1)*(outputs.VRI(i)/inputs.etaGen.Vmax)^3 + ...
                                 inputs.etaGen.param(2)*(outputs.VRI(i)/inputs.etaGen.Vmax)^2 + ...
                                  inputs.etaGen.param(3)*(outputs.VRI(i)/inputs.etaGen.Vmax)+inputs.etaGen.param(4))^sign(1);
        % PRI effective elec
        outputs.PRIeff_elec(i,j) = outputs.PRIeff_mech(i,j)/inputs.etaGearbox/inputs.etaSto/outputs.genEff_RI(i)/inputs.etaPE;

     end
         
      %% Cycle simulation
     
      % tRO
       if outputs.VRO(i,:)<0
        outputs.t1(i)             = 0;
        outputs.tROeff(i,:)       = 0./outputs.VRO(i,:);
        outputs.tRO(i)            = 0;
       else
        outputs.t1(i)       = outputs.VRO(i,1)/inputs.maxAcc;
        outputs.tROeff(i,:) = outputs.elemDeltaL(i)./outputs.VRO(i,:);
        outputs.tRO(i)      = outputs.t1(i) + sum(outputs.tROeff(i,:));
       end
     
      % Reel-out power during transition
      outputs.PRO1_mech(i) = outputs.PROeff_mech(i,1)/2;
      outputs.PRO1_elec(i) = outputs.PROeff_elec(i,1)/2;

      % PRO 
      if outputs.VRO(i,:)<0
        outputs.PRO_mech(i) = 1e-9;
      else
        outputs.PRO_mech(i) = (sum(outputs.PROeff_mech(i,:).*outputs.tROeff(i,:)) + outputs.PRO1_mech(i)*outputs.t1(i))/outputs.tRO(i);
        outputs.PRO_elec(i) = (sum(outputs.PROeff_elec(i,:).*outputs.tROeff(i,:)) + outputs.PRO1_elec(i)*outputs.t1(i))/outputs.tRO(i);
      end

      % tRI
      if outputs.VRO(i,:)<0
        outputs.t2(i)             = 0;
        outputs.tRIeff(i,:)       = ones(1,outputs.deltaLelems).*0;
        outputs.tRI(i)            = 0;
      else
        outputs.t2(i)       = outputs.VRI(i)/inputs.maxAcc;
        outputs.tRIeff(i,:) = ones(1,outputs.deltaLelems).*outputs.elemDeltaL(i)/outputs.VRI(i);
        outputs.tRI(i)      = outputs.t2(i) + sum(outputs.tRIeff(i,:));
      end
      
      % Transition reel-in mech power
      outputs.PRI2_mech(i)     = outputs.PRIeff_mech(i,1)/2;
      outputs.PRI2_elec(i)     = outputs.PRIeff_elec(i,1)/2;

      % PRI 
      outputs.PRI_mech(i) = (sum(outputs.PRIeff_mech(i,:).*outputs.tRIeff(i,:)) + outputs.PRI2_mech(i)*outputs.t2(i))/outputs.tRI(i);
      outputs.PRI_elec(i) = (sum(outputs.PRIeff_elec(i,:).*outputs.tRIeff(i,:)) + outputs.PRI2_elec(i)*outputs.t2(i))/outputs.tRI(i);

      % tCycle
      outputs.tCycle(i) = outputs.tRO(i)+outputs.tRI(i);

      outputs.tPatt(i,:)     = 2*pi()*outputs.pattRadius(i,:)./outputs.VC(i,:);
      outputs.numOfPatt(i,:) = outputs.tRO(i)./outputs.tPatt(i,:);
      
      %% Reel-out speed oscillation due to gravity
      
      if inputs.targetPRO_mech == 0 % Only run in First Optimisation-run
      
        % Vertical and Horizontal Sink rate difference theory
        outputs.T_simple(i)    = min(outputs.Tmax_act, (4/9)*outputs.CL(i)^3/outputs.CD(i)^2*(1/2)*mean(outputs.rho_air(i,:))*...
                                  inputs.WA*mean(outputs.Vw(i,:))^2);
        outputs.VSR_simple(i)  = outputs.CD(i)/outputs.CL(i)^(3/2)*sqrt(outputs.T_simple(i)/(0.5*mean(outputs.rho_air(i,:))*inputs.WA));
        outputs.a_simple(i)    = sqrt((mean(outputs.lambda(i,:))^2+mean(outputs.delta(i,:))^2)/(outputs.T_simple(i)^2+outputs.W(i)^2));
        outputs.VA_simple(i)   = 1/sqrt(outputs.a_simple(i));
        outputs.s(i)           = outputs.a_simple(i)*(mean(outputs.delta(i,:))*outputs.T_simple(i)-mean(outputs.lambda(i,:))*outputs.W(i))/...
                                    (mean(outputs.lambda(i,:))^2+mean(outputs.delta(i,:))^2);
        outputs.VSR_down(i)    = outputs.VA_simple(i)*outputs.s(i);
        outputs.VRO_osciAmp(i) = abs(abs(outputs.VSR_down(i)) - outputs.VSR_simple(i))*cos(outputs.avgPattEle(i));
%         outputs.VRO_osciAmp(i) = 0;
        % To match the number of deltaL elements to the number of elements considering pattern 
        if mean(outputs.numOfPatt(i,:)) == 0 
          outputs.numOfPatt(i,:) = 1;
        end
        outputs.numPattParts(i)   = outputs.deltaLelems/mean(outputs.numOfPatt(i,:));  
        for j=1:round(outputs.numPattParts(i)*mean(outputs.numOfPatt(i,:)))
          
          % Vertical and Horizontal Sink rate difference theory
          outputs.VRO_osci(i,j)             = outputs.VRO(i,j) + outputs.VRO_osciAmp(i)*sin((j-1)*2*pi()/(outputs.numPattParts(i))+270/180*pi()); 
          
          % KE + PE exchange theory
          outputs.VA_osciAmp(i,j) = sqrt(outputs.VA(i,j)^2+2*inputs.gravity*outputs.pattRadius(i,j)) - outputs.VA(i,j);
          outputs.VA_osci(i,j)    = outputs.VA(i,j) + outputs.VA_osciAmp(i)*sin((j-1)*2*pi()/(outputs.numPattParts(i))+270/180*pi());
          outputs.osciFactor(i,j) = (outputs.VA_osci(i,j)/outputs.VA(i,j))^2*cos(outputs.avgPattEle(i));
          outputs.VRO_osci2(i,j)   = outputs.VRO(i,j)*outputs.osciFactor(i,j);
            
          outputs.PROeff_mech_osci(i,j)     = outputs.T(i)*outputs.VRO_osci(i,j); %[W]
          outputs.PROeff_mech_osci_cap(i,j) = min(inputs.F_peakM2Ecyc*inputs.P_ratedElec, outputs.PROeff_mech_osci(i,j));
          outputs.genEff_RO_osci(i,j)  = (inputs.etaGen.param(1)*(outputs.VRO_osci(i,j)/inputs.etaGen.Vmax)^3 + ...
                                           inputs.etaGen.param(2)*(outputs.VRO_osci(i,j)/inputs.etaGen.Vmax)^2 + ...
                                            inputs.etaGen.param(3)*(outputs.VRO_osci(i,j)/inputs.etaGen.Vmax)+inputs.etaGen.param(4))^sign(1);
          outputs.PROeff_elec_osci_cap(i,j) = outputs.PROeff_mech_osci_cap(i,j)*inputs.etaGearbox*outputs.genEff_RO_osci(i,j)*inputs.etaPE;
        end
        outputs.PROeff_mech_cap(i) = mean(outputs.PROeff_mech_osci_cap(i,:)); 
        outputs.PROeff_elec_cap(i) = mean(outputs.PROeff_elec_osci_cap(i,:)); 
      end

      %% P_cycleElec 
      if outputs.VRO(i,:)<0
        outputs.P_cycleElec(i) = 0;
      else
       outputs.P_cycleElec(i) = (sum(outputs.tROeff(i,:).*outputs.PROeff_elec(i,:)) + outputs.t1(i)*outputs.PRO1_elec(i) - ...
                                   sum(outputs.tRIeff(i,:).*outputs.PRIeff_elec(i,:)) -outputs.t2(i)*outputs.PRI2_elec(i))/outputs.tCycle(i);    
      end

      % Without drivetrain eff
      outputs.P_cycleMech(i) = (sum(outputs.tROeff(i,:).*outputs.PROeff_mech(i,:)) + outputs.t1(i)*outputs.PRO1_mech(i) - ...
                                   sum(outputs.tRIeff(i,:).*outputs.PRIeff_mech(i,:)) - outputs.t2(i)*outputs.PRI2_mech(i))/outputs.tCycle(i);
end 
