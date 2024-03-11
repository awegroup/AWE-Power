function [inputs] = compute(i,inputs)
 global outputs
    
    %% Kite mass 
    if inputs.massOverride == 1
      outputs.m_k(i) = inputs.kiteMass;
    else
      % Vincent Bonnin's simple mass model developed at Ampyx Power. Based on AP3 data and projected data for larger systems (AP4-AP5)      
      a1     = 0.002415;       a2     = 0.0090239;       b1     = 0.17025;       b2     = 3.2493;
      k1     = 5;              c1     = 0.46608;         d1     = 0.65962;       k2     = 1.1935;
      AR_ref = 12;
      a = a1*(inputs.Ft_max/inputs.S) + a2;
      b = b1*(inputs.Ft_max/inputs.S) + b2;
      outputs.m_k(i) = 10*(a*inputs.S^2 +b*inputs.S-k1)*(c1*(inputs.AR/AR_ref)^2-d1*(inputs.AR/AR_ref)+k2); 
    end
    
    %% Constant for the cycle     
    outputs.l_t_min(i)         = outputs.Rp_start(i)/sin(outputs.gamma(i));
    outputs.pattStartGrClr(i)  = outputs.l_t_min(i)*sin(outputs.beta(i)-outputs.gamma(i));
    outputs.h_cycleStart(i)    = outputs.l_t_min(i)*cos(outputs.gamma(i))*sin(outputs.beta(i));
    outputs.l_t_max(i)         = outputs.l_t_min(i)+outputs.deltaL(i)/cos(outputs.gamma(i)); 
    outputs.pattEndGrClr(i)    = outputs.l_t_max(i)*sin(outputs.beta(i)-outputs.gamma(i));
    outputs.l_t_avg(i)         = (outputs.l_t_max(i)+outputs.l_t_min(i))/2; %[m]
    outputs.h_cycleAvg(i)      = outputs.l_t_avg(i)*cos(outputs.gamma(i))*sin(outputs.beta(i));
    outputs.h_cycleEnd(i)      = outputs.l_t_max(i)*cos(outputs.gamma(i))*sin(outputs.beta(i));
    outputs.d_t                = sqrt(inputs.Ft_max*1000/inputs.Te_matStrength*4/pi()); %[m] safety factor could be added (say *1.1)
    
    %% Discretizing the reel-out length in chosen number of elements
    % Found to be not highly sensitive to the number of elements
     outputs.deltaLelems   = inputs.numDeltaLelems; 
     outputs.elemDeltaL(i) = outputs.deltaL(i)/outputs.deltaLelems;
  
     %% Assigning and evaluating a single flight state equilibrium for each length element for reel-out and reel-in phase
     for j = 1:outputs.deltaLelems
        %% Reel-out phase:
        % Tether length at jth element
        if j == 1
          outputs.l_t_inCycle(i,j) = outputs.l_t_min(i) + outputs.elemDeltaL(i)/2/cos(outputs.gamma(i));
        else
          outputs.l_t_inCycle(i,j) = outputs.l_t_inCycle(i,j-1) + outputs.elemDeltaL(i)/cos(outputs.gamma(i));
        end

        % Pattern ground clearance
        outputs.pattGrClr(i,j) = outputs.l_t_inCycle(i,j)*sin(outputs.beta(i)-outputs.gamma(i));

        % Effective mass lumped at kite point (Kite + tether)
        outputs.m_t(i,j)   = inputs.Te_matDensity*pi()/4*outputs.d_t^2*outputs.l_t_inCycle(i,j);
        outputs.m_eff(i,j) = outputs.m_k(i) + 0.5*outputs.m_t(i,j);

        % Coordinates of the Kite's position and orientation in the Spherical ref. frame
        % Center point (representative point)
        outputs.theta(i,j) = pi()/2 - (outputs.beta(i));
        outputs.phi(i,j)   = 0;
        outputs.chi(i,j)   = deg2rad(90);
        
        % Effective CD
        outputs.CD_k(i,j)   = inputs.Cd0 + (outputs.CL(i,j)-inputs.Cl0_airfoil)^2/(pi()*inputs.AR*inputs.e);
        outputs.CD_t(i,j)   = (1/4)*inputs.Cd_c*outputs.d_t*outputs.l_t_inCycle(i,j)/inputs.S;
        outputs.CD(i,j)     = outputs.CD_k(i,j) + outputs.CD_t(i,j);

        % Effective Glide ratio
        outputs.E(i,j)      = outputs.CL(i,j)/outputs.CD(i,j);

        % Average pattern height at point of interest on deltaL
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
        
        if inputs.vertWindProfile == 0
          % Modelled 
          outputs.vw(i,j) = inputs.vw_ref(i)*((outputs.h_inCycle(i,j))/inputs.h_ref)^inputs.windShearExp;
        else
          % Extrapolated from dataset
          outputs.vw(i,j) = inputs.vw_ref(i)*interp1(inputs.windProfile_h, inputs.windProfile_vw, (outputs.h_inCycle(i,j)), 'linear', 'extrap');
        end

        % Intermediate calculation for brevity
        outputs.halfRhoS = 0.5*inputs.airDensity*inputs.S;
        
        % Wind velocity vector
        outputs.vw_r(i,j)     = outputs.vw(i,j)*sin(outputs.theta(i,j))*cos(outputs.phi(i,j));
        outputs.vw_theta(i,j) = outputs.vw(i,j)*cos(outputs.theta(i,j))*cos(outputs.phi(i,j));
        outputs.vw_phi(i,j)   = -outputs.vw(i,j)*sin(outputs.phi(i,j));

        % Reel-out factor
        outputs.f(i,j) = outputs.vk_r(i,j)/outputs.vw(i,j);

        % Apparent wind velocity magnitude
        outputs.va(i,j) = (outputs.vw_r(i,j)-outputs.vk_r(i,j))*sqrt(1+outputs.kRatio(i,j)^2);
        
        % Aerodynamic force magnitude
        outputs.Fa(i,j) = outputs.halfRhoS*sqrt(outputs.CL(i,j)^2+outputs.CD(i,j)^2)*outputs.va(i,j)^2;
        
        % Tangential kite velocity factor
        a = cos(outputs.theta(i,j))*cos(outputs.phi(i,j))*cos(outputs.chi(i,j))-sin(outputs.phi(i,j))*sin(outputs.chi(i,j));
        b = sin(outputs.theta(i,j))*cos(outputs.phi(i,j));
        outputs.lambda(i,j) = a + sqrt(a^2 + b^2 - 1 + outputs.kRatio(i,j)^2*(b - outputs.f(i,j))^2);
        % Tangential kite velocity
        outputs.vk_tau(i,j)     = outputs.lambda(i,j)*outputs.vw(i,j); 
        % Tangential kite velocity components in theta and phi directions
        outputs.vk_theta(i,j) = outputs.vk_tau(i,j)*cos(outputs.chi(i,j));
        outputs.vk_phi(i,j)   = outputs.vk_tau(i,j)*sin(outputs.chi(i,j));
     
        % Gravity toggle
        if inputs.FgToggle == 0
            outputs.W(i,j) = 0;
        else
            outputs.W(i,j)        = outputs.m_eff(i,j)*inputs.gravity;
        end
        
        % Gravitational force vector (kite + tether)
        outputs.Fg_r(i,j)     = -outputs.W(i)*cos(outputs.theta(i,j));
        outputs.Fg_theta(i,j) = outputs.W(i)*sin(outputs.theta(i,j));
        outputs.Fg_phi(i,j)   = 0;
        
        % Aerodynamic force vector
        outputs.Fa_theta(i,j) = -outputs.Fg_theta(i,j);
        outputs.Fa_phi(i,j)   = -outputs.Fg_phi(i,j); 
        outputs.Fa_r(i,j)     = sqrt(outputs.Fa(i,j)^2-outputs.Fa_theta(i,j)^2-outputs.Fa_phi(i,j)^2);
        
        % Roll angle in case of Gravity (N.A when CF is included)
        outputs.rollAngle(i,j) = atan(outputs.Fa_theta(i,j)/outputs.Fa_r(i,j));
        
        % Apparent wind velocity vector
        outputs.va_r(i,j)     = outputs.vw(i,j)*(sin(outputs.theta(i,j))*cos(outputs.phi(i,j))-outputs.f(i,j));
        outputs.va_theta(i,j) = outputs.vw(i,j)*(cos(outputs.theta(i,j))*cos(outputs.phi(i,j))-outputs.lambda(i,j)*cos(outputs.chi(i,j)));
        outputs.va_phi(i,j)   = outputs.vw(i,j)*(-sin(outputs.phi(i,j))-outputs.lambda(i,j)*sin(outputs.chi(i,j)));

        % Dot product F_a*v_a;
        outputs.F_dot_v(i,j) = outputs.Fa_r(i,j)*outputs.va_r(i,j) + outputs.Fa_theta(i,j)*outputs.va_theta(i,j) + outputs.Fa_phi(i,j)*outputs.va_phi(i,j);
        
        % Drag vector from dot product
        outputs.D_r(i,j)     = (outputs.F_dot_v(i,j)/outputs.va(i,j)^2)*outputs.va_r(i,j);
        outputs.D_theta(i,j) = (outputs.F_dot_v(i,j)/outputs.va(i,j)^2)*outputs.va_theta(i,j);
        outputs.D_phi(i,j)   = (outputs.F_dot_v(i,j)/outputs.va(i,j)^2)*outputs.va_phi(i,j);
        
        % Lift vector
        outputs.L_r(i,j)     = outputs.Fa_r(i,j) - outputs.D_r(i,j);
        outputs.L_theta(i,j) = outputs.Fa_theta(i,j) - outputs.D_theta(i,j);
        outputs.L_phi(i,j)   = outputs.Fa_phi(i,j) - outputs.D_phi(i,j);
        
        % Drag magnitude
        outputs.D(i,j) = sqrt(outputs.D_r(i,j)^2 + outputs.D_theta(i,j)^2 + outputs.D_phi(i,j)^2);
        
        % Lift magnitude
        outputs.L(i,j) = sqrt(outputs.L_r(i,j)^2 + outputs.L_theta(i,j)^2 + outputs.L_phi(i,j)^2);
        
        % Straight-tether force 
        outputs.Ft(i,j) = outputs.Fa_r(i,j) + outputs.Fg_r(i,j);
        
        % Lift-to-drag ratio that follows from the chosen kinematic ratio
        outputs.E_result(i,j) = sqrt(((outputs.Fa(i,j)*outputs.va(i,j))/outputs.F_dot_v(i,j))^2-1);

        % Ratio of Kinemactic Ratio and Glide Ratio
        outputs.kByE(i,j) = outputs.kRatio(i,j)/outputs.E_result(i,j);
          
        % Effective mechanical reel-out power
        outputs.P_m_o_eff(i,j) = outputs.Ft(i,j)*outputs.vk_r(i,j); %[W]   
        
        outputs.zetaMech(i,j)    = outputs.P_m_o_eff(i,j)/(outputs.halfRhoS*outputs.vw(i,j)^3);
       
        % Effective electrical reel-out power
        % Generator efficiency. As a function of RPM/RPM_max, where RPM_max is driven by winch 
        outputs.etaGen_o(i,j)  = (inputs.etaGen.param(1)*(outputs.vk_r(i,j)/inputs.etaGen.v_max)^3 + ...
                                      inputs.etaGen.param(2)*(outputs.vk_r(i,j)/inputs.etaGen.v_max)^2 + ...
                                        inputs.etaGen.param(3)*(outputs.vk_r(i,j)/inputs.etaGen.v_max)+inputs.etaGen.param(4))^sign(1);
        outputs.P_e_o_eff(i,j) = outputs.P_m_o_eff(i,j)*inputs.etaGearbox*outputs.etaGen_o(i,j)*inputs.etaPE;

        %% Retraction Phase
        % Kite position in spherical coordinates
        % Reel-in is assumed to start from the top of the pattern
        outputs.theta_i(i,j) = pi()/2 - (outputs.beta(i)+outputs.gamma(i));
        outputs.phi_i(i,j)   = 0;
        
        if inputs.vertWindProfile == 0
          % Modelled 
          outputs.vw_i(i,j) = inputs.vw_ref(i)*((outputs.h_inCycle(i,j) + outputs.Rp(i,j)*cos(outputs.beta(i)))/inputs.h_ref)^inputs.windShearExp;
        else
          % Extrapolated from dataset
          outputs.vw_i(i,j) = inputs.vw_ref(i)*interp1(inputs.windProfile_h, inputs.windProfile_vw, (outputs.h_inCycle(i,j) + outputs.Rp(i,j)*cos(outputs.beta(i))), 'linear', 'extrap');
        end
     
        % Wind velocity vector
        outputs.vw_r_i(i,j)     = outputs.vw_i(i,j)*sin(outputs.theta_i(i,j))*cos(outputs.phi_i(i,j));
        outputs.vw_theta_i(i,j) = outputs.vw_i(i,j)*cos(outputs.theta_i(i,j))*cos(outputs.phi_i(i,j));
        outputs.vw_phi_i(i,j)   = -outputs.vw_i(i,j)*sin(outputs.phi_i(i,j)); % 

        % Reel-in factor
        outputs.f_i(i,j) = outputs.vk_r_i(i,j)/outputs.vw_i(i,j);
        
        % Apparent speed vector
        outputs.va_r_i(i,j)     = outputs.vw_r_i(i,j) - (-outputs.vk_r_i(i,j)); % Reel-in speed direction is in the negative radial direction)
        outputs.va_theta_i(i,j) = outputs.vw_theta_i(i,j) - 0;
        outputs.va_phi_i(i,j)   = outputs.vw_phi_i(i,j) - 0;

        % Apparent wind velocity magnitude
        outputs.va_i(i,j) = sqrt(outputs.va_r_i(i,j)^2 + outputs.va_theta_i(i,j)^2); % Wind has components only in r and theta directions

        % Aerodynamic force magnitude
        outputs.CD_k_i(i,j)     = inputs.Cd0+(outputs.CL_i(i,j)- inputs.Cl0_airfoil)^2/(pi()*inputs.AR*inputs.e);
        outputs.CD_i(i,j)       = outputs.CD_k_i(i,j) + outputs.CD_t(i,j);
        outputs.E_i(i,j)        = outputs.CL_i(i,j)/outputs.CD_i(i,j);
        outputs.Fa_i(i,j)       = outputs.halfRhoS*sqrt(outputs.CL_i(i,j)^2+outputs.CD_i(i,j)^2)*outputs.va_i(i,j)^2;

        % Gravitational force vector (kite + tether)
        outputs.Fg_r_i(i,j)       = -outputs.W(i)*cos(outputs.theta_i(i,j));
        outputs.Fg_theta_i(i,j)   = outputs.W(i)*sin(outputs.theta_i(i,j));
        outputs.Fg_phi_i(i,j)     = 0;

        % Aerodynamic force vector
        outputs.Fa_theta_i(i,j) = -outputs.Fg_theta_i(i,j);
        outputs.Fa_phi_i(i,j)   = -outputs.Fg_phi_i(i,j);
        outputs.Fa_r_i(i,j)     = sqrt(outputs.Fa_i(i,j)^2-outputs.Fa_theta_i(i,j)^2-outputs.Fa_phi_i(i,j)^2);

        % Dot product F_a*v_a;
        outputs.F_dot_v_i(i,j) = outputs.Fa_r_i(i,j)*outputs.va_r_i(i,j) + outputs.Fa_theta_i(i,j)*outputs.va_theta_i(i,j) + outputs.Fa_phi_i(i,j)*outputs.va_phi_i(i,j);
        
        % Drag vector from dot product
        outputs.D_r_i(i,j)     = (outputs.F_dot_v_i(i,j)/outputs.va_i(i,j)^2)*outputs.va_r_i(i,j);
        outputs.D_theta_i(i,j) = (outputs.F_dot_v_i(i,j)/outputs.va_i(i,j)^2)*outputs.va_theta_i(i,j);
        outputs.D_phi_i(i,j)   = (outputs.F_dot_v_i(i,j)/outputs.va_i(i,j)^2)*outputs.va_phi_i(i,j);
        
        % Lift vector
        outputs.L_r_i(i,j)     = outputs.Fa_r_i(i,j) - outputs.D_r_i(i,j);
        outputs.L_theta_i(i,j) = outputs.Fa_theta_i(i,j) - outputs.D_theta_i(i,j);
        outputs.L_phi_i(i,j)   = outputs.Fa_phi_i(i,j) - outputs.D_phi_i(i,j);
        
        % Drag magnitude
        outputs.D_i(i,j) = sqrt(outputs.D_r_i(i,j)^2 + outputs.D_theta_i(i,j)^2 + outputs.D_phi_i(i,j)^2);
        
        % Lift magnitude
        outputs.L_i(i,j) = sqrt(outputs.L_r_i(i,j)^2 + outputs.L_theta_i(i,j)^2 + outputs.L_phi_i(i,j)^2);

        % Lift-to-drag ratio that follows from the chosen kinematic ratio
        outputs.E_result_i(i,j) = sqrt(((outputs.Fa_i(i,j)*outputs.va_i(i,j))/outputs.F_dot_v_i(i,j))^2-1);

        % Straight-tether force 
        outputs.Ft_i(i,j) = outputs.Fa_r_i(i,j) + outputs.Fg_r_i(i,j);

        % Effective mechanical reel-out power
        outputs.P_m_i_eff(i,j) = outputs.Ft_i(i,j)*outputs.vk_r_i(i,j); %[W]  

        % Generator efficiency during RI: As a function of RPM/RPM_max, where RPM_max is driven by winch i.e Max VRI
        outputs.etaGen_i(i) = (inputs.etaGen.param(1)*(outputs.vk_r_i(i)/inputs.etaGen.v_max)^3 + ...
                                 inputs.etaGen.param(2)*(outputs.vk_r_i(i)/inputs.etaGen.v_max)^2 + ...
                                  inputs.etaGen.param(3)*(outputs.vk_r_i(i)/inputs.etaGen.v_max)+inputs.etaGen.param(4))^sign(1);
        
        % Effective electrical reel-in power
        outputs.P_e_i_eff(i,j) = outputs.P_m_i_eff(i,j)/inputs.etaGearbox/inputs.etaSto/outputs.etaGen_i(i)/inputs.etaPE;

     end
         
     %% Cycle calculation
     
      % Reel-out time
      outputs.t1(i)       = outputs.vk_r(i,1)/inputs.a_d_max;
      outputs.to_eff(i,:) = outputs.elemDeltaL(i)./outputs.vk_r(i,:);
      outputs.to(i)       = outputs.t1(i) + sum(outputs.to_eff(i,:));
     
      % Reel-out power during transition
      outputs.P1_m_o(i) = outputs.P_m_o_eff(i,1)/2;
      outputs.P1_e_o(i) = outputs.P_e_o_eff(i,1)/2;

      % Reel-out power 
        outputs.P_m_o(i) = (sum(outputs.P_m_o_eff(i,:).*outputs.to_eff(i,:)) + outputs.P1_m_o(i)*outputs.t1(i))/outputs.to(i);
        outputs.P_e_o(i) = (sum(outputs.P_e_o_eff(i,:).*outputs.to_eff(i,:)) + outputs.P1_e_o(i)*outputs.t1(i))/outputs.to(i);

      % Reel-in time
        outputs.t2(i)       = outputs.vk_r_i(i)/inputs.a_d_max;
        outputs.ti_eff(i,:) = outputs.elemDeltaL(i)./outputs.vk_r_i(i,:);
        outputs.ti(i)       = outputs.t2(i) + sum(outputs.ti_eff(i,:));
      
        % Reel-in power duing transition
        outputs.P2_m_i(i)     = outputs.P_m_i_eff(i,1)/2;
        outputs.P2_e_i(i)     = outputs.P_e_i_eff(i,1)/2;
  
        % Reel-in power 
        outputs.P_m_i(i) = (sum(outputs.P_m_i_eff(i,:).*outputs.ti_eff(i,:)) + outputs.P2_m_i(i)*outputs.t2(i))/outputs.ti(i);
        outputs.P_e_i(i) = (sum(outputs.P_e_i_eff(i,:).*outputs.ti_eff(i,:)) + outputs.P2_e_i(i)*outputs.t2(i))/outputs.ti(i);
  
        % Cycle time
        outputs.tCycle(i) = outputs.to(i)+outputs.ti(i);
  
        % Time for one pattern revolution and number of patterns in the cycle
        outputs.tPatt(i,:)     = 2*pi()*outputs.Rp(i,:)./outputs.vk_phi(i,:);
        outputs.numOfPatt(i,:) = outputs.to(i)./outputs.tPatt(i,:);
  
        % Electrical cycle power
         outputs.P_e_avg(i) = (sum(outputs.to_eff(i,:).*outputs.P_e_o_eff(i,:)) + outputs.t1(i)*outputs.P1_e_o(i) - ...
                                     sum(outputs.ti_eff(i,:).*outputs.P_e_i_eff(i,:)) -outputs.t2(i)*outputs.P2_e_i(i))/outputs.tCycle(i);   
  
        % Mechanical cycle power - without drivetrain eff
        outputs.P_m_avg(i) = (sum(outputs.to_eff(i,:).*outputs.P_m_o_eff(i,:)) + outputs.t1(i)*outputs.P1_m_o(i) - ...
                                     sum(outputs.ti_eff(i,:).*outputs.P_m_i_eff(i,:)) - outputs.t2(i)*outputs.P2_m_i(i))/outputs.tCycle(i);

end      