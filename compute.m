function [inputs] = compute(i,inputs)
 global outputs
    
    %% Kite mass 
    if inputs.massOverride == 1
      outputs.m_k(i) = inputs.kiteMass;
    else
      %Vincent Bonnin's simple mass model developed at Ampyx Power. Based on AP3 data and projected data for larger systems (AP4-AP5)      
      a1     = 0.002415;       a2     = 0.0090239;       b1     = 0.17025;       b2     = 3.2493;
      k1     = 5;              c1     = 0.46608;         d1     = 0.65962;       k2     = 1.1935;
      AR_ref = 12;
      a = a1*(inputs.Ft_max/inputs.S) + a2;
      b = b1*(inputs.Ft_max/inputs.S) + b2;
      outputs.m_k(i) = 10*(a*inputs.S^2 +b*inputs.S-k1)*(c1*(inputs.AR/AR_ref)^2-d1*(inputs.AR/AR_ref)+k2); 
    end
    
    %% Constants
    outputs.wingSpan           = sqrt(inputs.AR*inputs.S);     
    outputs.l_t_min(i)         = outputs.Rp_start(i)/sin(outputs.gamma(i));
    outputs.pattStartGrClr(i)  = outputs.l_t_min(i)*sin(outputs.beta(i)-outputs.gamma(i));
    outputs.h_cycleStart(i)    = outputs.l_t_min(i)*cos(outputs.gamma(i))*sin(outputs.beta(i));
    outputs.l_t_max(i)         = outputs.l_t_min(i)+outputs.deltaL(i)/cos(outputs.gamma(i)); 
    outputs.pattEndGrClr(i)    = outputs.l_t_max(i)*sin(outputs.beta(i)-outputs.gamma(i));
    outputs.l_t_avg(i)         = (outputs.l_t_max(i)+outputs.l_t_min(i))/2; %[m]
    outputs.h_cycleAvg(i)      = outputs.l_t_avg(i)*cos(outputs.gamma(i))*sin(outputs.beta(i));
    outputs.h_cycleEnd(i)      = outputs.l_t_max(i)*cos(outputs.gamma(i))*sin(outputs.beta(i));
    outputs.d_t                = sqrt(inputs.Ft_max*1000/inputs.Te_matStrength*4/pi()); %[m] safety factor could be added (say *1.1)
    
    %% Cycle avg. mech reel-out power considering vertical wind shear
     
    % Discretizing the reel-out length in chosen number of elements
    % Found to be not highly sensitive to the number of elements
     outputs.deltaLelems   = inputs.numDeltaLelems; 
     outputs.elemDeltaL(i) = outputs.deltaL(i)/outputs.deltaLelems;
  
     % Assigning and evaluating a single flight state equilibrium for each length element 
     for j = 1:outputs.deltaLelems
        
        % Defining kite positioning parameters 

        if j == 1
          outputs.l_t_inCycle(i,j) = outputs.l_t_min(i) + outputs.elemDeltaL(i)/2/cos(outputs.gamma(i));
        else
          outputs.l_t_inCycle(i,j) = outputs.l_t_inCycle(i,j-1) + outputs.elemDeltaL(i)/cos(outputs.gamma(i));
        end

        % Pattern ground clearance
        outputs.pattGrClr(i,j) = outputs.l_t_inCycle(i,j)*sin(outputs.beta(i)-outputs.gamma(i));

        % Effective mass lumped at kite (Kite + tether)
        if inputs.FgToggle == 0
            outputs.m_t(i,j) = 0;
            outputs.m_k(i)   = 0;
        else
            outputs.m_t(i,j)  = inputs.Te_matDensity*pi()/4*outputs.d_t^2*outputs.l_t_inCycle(i,j); % Could add (say 0.85) as safety factor on material density
        end
        outputs.m_eff(i,j) = outputs.m_k(i)+sin(outputs.beta(i))*outputs.m_t(i,j);

        if inputs.evalPoint      == 0
            % Bottom point
            outputs.theta(i,j) = pi()/2 - (outputs.beta(i)-outputs.gamma(i));
            outputs.phi(i,j)   = 0;
            outputs.chi(i,j)   = deg2rad(270);
        elseif inputs.evalPoint      == 1
            % Top point
            outputs.theta(i,j) = pi()/2 - (outputs.beta(i)+outputs.gamma(i));
            outputs.phi(i,j)   = 0;
            outputs.chi(i,j)   = deg2rad(90);
        else 
            % Representative flight state
             outputs.theta(i,j) = pi()/2 - (outputs.beta(i));
             outputs.phi(i,j)   = 0;
             outputs.chi(i,j)   = deg2rad(90);
        end
        
        % Effective CD
        outputs.CL(i,j)     = inputs.CL_maxAirfoil*inputs.CLeff_F;
        outputs.CD_k(i,j)   = inputs.CD0 + (outputs.CL(i,j)-inputs.CL0_airfoil)^2/(pi()*inputs.AR*inputs.e);
        outputs.CD_t(i,j)   = (1/4)*inputs.CD_t*outputs.d_t*outputs.l_t_inCycle(i,j)/inputs.S;
        outputs.CD(i,j)     = outputs.CD_k(i,j) + outputs.CD_t(i,j);

        %% Increasing radius
        % Pattern average height
        if j == 1
          outputs.h_inCycle(i,j) = outputs.h_cycleStart(i) + outputs.elemDeltaL(i)/2*sin(outputs.beta(i));
        else
          outputs.h_inCycle(i,j) = outputs.h_inCycle(i,j-1) + outputs.elemDeltaL(i)*sin(outputs.beta(i));
        end

        % Pattern radius at point of interest on deltaL
        if j == 1
          outputs.Rp(i,j) = outputs.Rp_start(i) + outputs.elemDeltaL(i)/2*tan(outputs.gamma(i));
        else
          outputs.Rp(i,j) = outputs.Rp(i,j-1) + outputs.elemDeltaL(i)*tan(outputs.gamma(i));
        end
        
        % Wind speed at evaluation point                       
        outputs.vw(i,j) = inputs.vw_ref(i)*((outputs.h_inCycle(i,j)-outputs.Rp(i,j)*cos(outputs.beta(i)))...
                            /inputs.h_ref)^inputs.windShearExp;
        % outputs.vw(i,j) = inputs.vw_ref(i)*((outputs.h_inCycle(i,j))/inputs.h_ref)^inputs.windShearExp;
        
        % Air density as a function of height   % Ref: https://en.wikipedia.org/wiki/Density_of_air
        % M = 0.0289644; % [kg/mol]
        % R = 8.3144598; % [N·m/(mol·K)]
        % T = 288.15;    % [Kelvin]
        % L = 0.0065;    % [Kelvin/m] 
        % outputs.rho_air(i,j) = inputs.airDensity*(1-L*(outputs.h_inCycle(i,j)-outputs.Rp(i,j)*cos(outputs.beta(i)))/T)^(inputs.gravity*M/R/L-1);
        % outputs.rho_air(i,j) = inputs.airDensity*(1-L*(outputs.h_inCycle(i,j))/T)^(inputs.gravity*M/R/L-1);
        outputs.rho_air(i,j)   = inputs.airDensity;

        % Intermediate calculation for brevity
        outputs.CR(i,j)       = sqrt(outputs.CL(i,j)^2+outputs.CD(i,j)^2);
        outputs.halfRhoS(i,j) = 0.5*outputs.rho_air(i,j)*inputs.S;

        % apparent wind velocity magnitude
        outputs.va(i,j) = (sin(outputs.theta(i,j))*cos(outputs.phi(i,j))-outputs.f(i,j))*outputs.vw(i,j)*sqrt(1+outputs.kRatio(i,j)^2);
        
        % aerodynamic force magnitude
        outputs.Fa(i,j) = outputs.halfRhoS(i,j)*outputs.CR(i,j)*outputs.va(i,j)^2;
        
        % tangential kite velocity factor
        a = cos(outputs.theta(i,j))*cos(outputs.phi(i,j))*cos(outputs.chi(i,j))-sin(outputs.phi(i,j))*sin(outputs.chi(i,j));
        b = sin(outputs.theta(i,j))*cos(outputs.phi(i,j));
        outputs.lambda(i,j) = a + sqrt(a^2 + b^2 - 1 + outputs.kRatio(i,j)^2 *(b - outputs.f(i,j))^2);
                          
        % Circumferential velocity
        outputs.vk_omega(i,j)   = outputs.lambda(i,j)*outputs.vw(i,j);  
        
        % magnitude of centripetal force
        if inputs.FcToggle == 0
            outputs.Fc(i,j) = 0;
        else
            outputs.Fc(i,j) = outputs.m_eff(i)*outputs.vk_omega(i,j)^2/outputs.Rp(i,j);
        end

        if inputs.evalPoint == 0
          % centrifugal force vector: Bottom point
          outputs.Fc_r(i,j)     = +outputs.Fc(i,j)*sin(outputs.gamma(i));
          outputs.Fc_theta(i,j) = +outputs.Fc(i,j)*cos(outputs.gamma(i));
          outputs.Fc_phi(i,j)   = 0;
        else
        % centrifugal force vector: Top point
           outputs.Fc_r(i,j)     = +outputs.Fc(i,j)*sin(outputs.gamma(i));
           outputs.Fc_theta(i,j) = -outputs.Fc(i,j)*cos(outputs.gamma(i));
           outputs.Fc_phi(i,j)   = 0;
        end
     
        % gravitational force vector (kite + tether)
        outputs.W(i,j)        = outputs.m_eff(i,j)*inputs.gravity;
        outputs.Fg_r(i,j)     = -outputs.W(i)*cos(outputs.theta(i,j));
        outputs.Fg_theta(i,j) = outputs.W(i)*sin(outputs.theta(i,j));
        outputs.Fg_phi(i,j)   = 0;
        
        
        % aerodynamic force vector
        outputs.Fa_theta(i,j) = -outputs.Fg_theta(i,j) -outputs.Fc_theta(i,j);
        outputs.Fa_phi(i,j) = 0;
        outputs.Fa_r(i,j) = sqrt(outputs.Fa(i,j)^2-outputs.Fa_theta(i,j)^2-outputs.Fa_phi(i,j)^2);
        
        % Roll angle
        outputs.rollAngle(i,j) = atan(outputs.Fa_theta(i,j)/outputs.Fa_r(i,j));
        
        
        % wind velocity vector
        outputs.vw_r(i,j)     = outputs.vw(i,j)*sin(outputs.theta(i,j))*cos(outputs.phi(i,j));
        outputs.vw_theta(i,j) = outputs.vw(i,j)*cos(outputs.theta(i,j))*cos(outputs.phi(i,j));
        outputs.vw_phi(i,j)   = outputs.vw(i,j)*(-sin(outputs.phi(i,j)));
        
        % apparent wind velocity vector
        outputs.va_r(i,j)     = outputs.vw(i,j)*(sin(outputs.theta(i,j))*cos(outputs.phi(i,j))-outputs.f(i,j));
        outputs.va_theta(i,j) = outputs.vw(i,j)*(cos(outputs.theta(i,j))*cos(outputs.phi(i,j))-outputs.lambda(i,j)*cos(outputs.chi(i,j)));
        outputs.va_phi(i,j)   = outputs.vw(i,j)*(-sin(outputs.phi(i,j))-outputs.lambda(i,j)*sin(outputs.chi(i,j)));
        
        % dot product F_a*v_a;
        outputs.F_dot_v(i,j) = (outputs.Fa_r(i,j)*outputs.va_r(i,j) + outputs.Fa_theta(i,j)*outputs.va_theta(i,j) + outputs.Fa_phi(i,j)*outputs.va_phi(i,j));
        
        % drag vector
        outputs.D_r(i,j)     = (outputs.F_dot_v(i,j)/outputs.va(i,j)^2)*outputs.va_r(i,j);
        outputs.D_theta(i,j) = (outputs.F_dot_v(i,j)/outputs.va(i,j)^2)*outputs.va_theta(i,j);
        outputs.D_phi(i,j)   = (outputs.F_dot_v(i,j)/outputs.va(i,j)^2)*outputs.va_phi(i,j);
        
        % lift vector
        outputs.L_r(i,j) = outputs.Fa_r(i,j) - outputs.D_r(i,j);
        outputs.L_theta(i,j) = outputs.Fa_theta(i,j) - outputs.D_theta(i,j);
        outputs.L_phi(i,j) = outputs.Fa_phi(i,j) - outputs.D_phi(i,j);
        
        % drag magnitude
        outputs.D(i,j) = sqrt(outputs.D_r(i,j)^2 + outputs.D_theta(i,j)^2 + outputs.D_phi(i,j)^2);
        
        % lift magnitude
        outputs.L(i,j) = sqrt(outputs.L_r(i,j)^2 + outputs.L_theta(i,j)^2 + outputs.L_phi(i,j)^2);
        
        % Straight-tether force 
        outputs.Ft(i,j) = outputs.Fa_r(i,j) + outputs.Fg_r(i,j) + outputs.Fc_r(i,j);

        % Loss in Ft due to tether sag
        outputs.Ft_k_theta(i,j) = -(1/2)*sin(outputs.theta(i,j))*outputs.m_t(i,j)*inputs.gravity;
        outputs.Ft_k_r(i,j)     = sqrt(outputs.Ft(i,j)^2 - outputs.Ft_k_theta(i,j)^2);
        outputs.Ft_k_phi(i,j)   = 0; 

        % Tether force at the drum
        outputs.Ft_drum(i,j) = sqrt((sqrt(outputs.Ft_k_r(i,j)^2 - outputs.Ft_k_theta(i,j)^2) - ...
                                cos(outputs.theta(i,j))*outputs.m_t(i)*inputs.gravity)^2 + outputs.Ft_k_theta(i,j)^2);

        % reel-out speed
        outputs.vk_r(i,j) = outputs.f(i,j)*outputs.vw(i,j);
        
        % lift-to-drag ratio that follows from the chosen kinematic ratio
        outputs.k_result(i,j) = sqrt(((outputs.Fa(i,j)*outputs.va(i,j))/outputs.F_dot_v(i,j))^2-1);
          
        % Effective mechanical reel-out power
        outputs.PROeff_mech(i,j) = outputs.Ft_drum(i,j)*outputs.vk_r(i,j); %[W]        
       
        % Effective electrical reel-out power
        % Generator efficiency. As a function of RPM/RPM_max, where RPM_max is driven by winch i.e Max VRI
        outputs.genEff_RO(i,j)  = (inputs.etaGen.param(1)*(outputs.vk_r(i,j)/inputs.etaGen.v_max)^3 + ...
                                      inputs.etaGen.param(2)*(outputs.vk_r(i,j)/inputs.etaGen.v_max)^2 + ...
                                        inputs.etaGen.param(3)*(outputs.vk_r(i,j)/inputs.etaGen.v_max)+inputs.etaGen.param(4))^sign(1);
        outputs.PROeff_elec(i,j) = outputs.PROeff_mech(i,j)*inputs.etaGearbox*outputs.genEff_RO(i,j)*inputs.etaPE;

        % Effective mechanical reel-in power
        outputs.va_i(i,j)       = sqrt(outputs.vw(i,j)^2 +outputs.vk_r_i(i)^2 + 2*outputs.vw(i,j)*outputs.vk_r_i(i)*cos(outputs.beta(i)));
        outputs.CL_i(i,j)       = 2*outputs.W(i)/(outputs.rho_air(i,j)*outputs.va_i(i,j)^2*inputs.S);
        outputs.CD_i(i,j)       = inputs.CD0+(outputs.CL_i(i,j)- inputs.CL0_airfoil)^2/(pi()*inputs.AR*inputs.e) + outputs.CD_t(i);
        outputs.PRIeff_mech(i,j) = 0.5*outputs.rho_air(i,j)*outputs.CD_i(i,j)*inputs.S*outputs.va_i(i,j)^3;

        % Generator efficiency during RI: As a function of RPM/RPM_max, where RPM_max is driven by winch i.e Max VRI
        outputs.genEff_RI(i) = (inputs.etaGen.param(1)*(outputs.vk_r_i(i)/inputs.etaGen.v_max)^3 + ...
                                 inputs.etaGen.param(2)*(outputs.vk_r_i(i)/inputs.etaGen.v_max)^2 + ...
                                  inputs.etaGen.param(3)*(outputs.vk_r_i(i)/inputs.etaGen.v_max)+inputs.etaGen.param(4))^sign(1);
        
        % Effective electrical reel-in power
        outputs.PRIeff_elec(i,j) = outputs.PRIeff_mech(i,j)/inputs.etaGearbox/inputs.etaSto/outputs.genEff_RI(i)/inputs.etaPE;

     end
         
    %% Cycle simulation
     
    % Reel-out time
     if outputs.vk_r(i,:)<0
      outputs.t1(i)             = 0;
      outputs.tROeff(i,:)       = 0./outputs.vk_r(i,:);
      outputs.tRO(i)            = 0;
     else
      outputs.t1(i)       = outputs.vk_r(i,1)/inputs.winchAcc_max;
      outputs.tROeff(i,:) = outputs.elemDeltaL(i)./outputs.vk_r(i,:);
      outputs.tRO(i)      = outputs.t1(i) + sum(outputs.tROeff(i,:));
     end
     
      % Reel-out power during transition
      outputs.PRO1_mech(i) = outputs.PROeff_mech(i,1)/2;
      outputs.PRO1_elec(i) = outputs.PROeff_elec(i,1)/2;

      % Reel-out power 
      if outputs.vk_r(i,:)<0
        outputs.PRO_mech(i) = 1e-9;
      else
        outputs.PRO_mech(i) = (sum(outputs.PROeff_mech(i,:).*outputs.tROeff(i,:)) + outputs.PRO1_mech(i)*outputs.t1(i))/outputs.tRO(i);
        outputs.PRO_elec(i) = (sum(outputs.PROeff_elec(i,:).*outputs.tROeff(i,:)) + outputs.PRO1_elec(i)*outputs.t1(i))/outputs.tRO(i);
      end

      % Reel-in time
      if outputs.vk_r(i,:)<0
        outputs.t2(i)             = 0;
        outputs.tRIeff(i,:)       = ones(1,outputs.deltaLelems).*0;
        outputs.tRI(i)            = 0;
      else
        outputs.t2(i)       = outputs.vk_r_i(i)/inputs.winchAcc_max;
        outputs.tRIeff(i,:) = ones(1,outputs.deltaLelems).*outputs.elemDeltaL(i)/outputs.vk_r_i(i);
        outputs.tRI(i)      = outputs.t2(i) + sum(outputs.tRIeff(i,:));
      end
      
      % Reel-in power duing transition
      outputs.PRI2_mech(i)     = outputs.PRIeff_mech(i,1)/2;
      outputs.PRI2_elec(i)     = outputs.PRIeff_elec(i,1)/2;

      % Reel-in power 
      outputs.PRI_mech(i) = (sum(outputs.PRIeff_mech(i,:).*outputs.tRIeff(i,:)) + outputs.PRI2_mech(i)*outputs.t2(i))/outputs.tRI(i);
      outputs.PRI_elec(i) = (sum(outputs.PRIeff_elec(i,:).*outputs.tRIeff(i,:)) + outputs.PRI2_elec(i)*outputs.t2(i))/outputs.tRI(i);

      % Cycle time
      outputs.tCycle(i) = outputs.tRO(i)+outputs.tRI(i);

      % Time for one pattern revolution and number of patterns in the cycle
      outputs.tPatt(i,:)     = 2*pi()*outputs.Rp(i,:)./outputs.vk_omega(i,:);
      outputs.numOfPatt(i,:) = outputs.tRO(i)./outputs.tPatt(i,:);

      % Electrical cycle power
      if outputs.vk_r(i,:)<=0
        outputs.P_cycleElec(i) = 0;
      else
       outputs.P_cycleElec(i) = (sum(outputs.tROeff(i,:).*outputs.PROeff_elec(i,:)) + outputs.t1(i)*outputs.PRO1_elec(i) - ...
                                   sum(outputs.tRIeff(i,:).*outputs.PRIeff_elec(i,:)) -outputs.t2(i)*outputs.PRI2_elec(i))/outputs.tCycle(i);    
      end

      % Mechanical cycle power - without drivetrain eff
      outputs.P_cycleMech(i) = (sum(outputs.tROeff(i,:).*outputs.PROeff_mech(i,:)) + outputs.t1(i)*outputs.PRO1_mech(i) - ...
                                   sum(outputs.tRIeff(i,:).*outputs.PRIeff_mech(i,:)) - outputs.t2(i)*outputs.PRI2_mech(i))/outputs.tCycle(i);

end      