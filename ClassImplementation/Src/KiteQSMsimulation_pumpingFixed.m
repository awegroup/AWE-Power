classdef KiteQSMsimulation_pumpingFixed < KiteQSMsimulationBase
    methods
        function obj = KiteQSMsimulation_pumpingFixed(inputs, loggerIdentifier)
            % Class Constructor for KiteQSMsimulation (pumping fixed-wing)
            %
            % Args:
            %     inputs (struct): A structure containing the required and optional simulation parameters.
            %         Required Parameters:
            %             - windSpeedReference (float): Wind speed at reference height (m/s).
            %             - heightWindReference (float): Height at which wind speed is measured (m).
            %             - areaWing (float): Area of the wing (m^2).
            %             - span (float): Wing span (m).
            %             - forceTether_max (float): Maximum tether force (N).
            %             - P_ratedElec (float): Rated electrical power (W).
            %             - Cl_maxAirfoil (float): Maximum coefficient of lift of the airfoil.
            %             - Cl_eff_F (float): Effective lift coefficient factor.
            %             - Cd0 (float): Zero-lift drag coefficient.
            %
            %         Optional Parameters:
            %             - aspectRatio (float, default = span^2/areaWing): Aspect ratio of the wing.
            %             - numDeltaLelems (int, default = 5): Number of delta elements.
            %             - doIncludeGravity (bool, default = true): Flag indicating whether to include gravity effects.
            %             - doMassOverride (bool, default = false): Flag indicating whether to override the kite mass.
            %             - maxTeLen (float, default = 1000): Maximum tether length (m).
            %             - maxHeight (float, default = 1000): Maximum height (m).
            %             - minGroundClear (float, default = 100): Minimum ground clearance (m).
            %             - safetyFactor_forceTetherMax (float, default = 0.8): Safety factor for maximum tether force.
            %             - peakM2E_F (float, default = 2.5): Peak mechanical-to-electric energy factor.
            %             - maxStrengthTether (float, default = 7e8): Maximum tether strength (Pa).
            %             - densityTether (float, default = 980): Density of the tether material (kg/m^3).
            %             - speedReelout_max (float, default = 20): Maximum reel-out speed (m/s).
            %             - accReel_max (float, default = 5): Maximum reel-out acceleration (m/s^2).
            %             - e (float, default = 0.6): Oswald efficiency factor.
            %             - Cd_cylinder (float, default = 1.2): Drag coefficient for the cylinder.
            %             - etaGen_param (1x4 float array, default = [0.671, -1.4141, 0.9747, 0.7233]): Generator efficiency parameters.
            %             - efficiency_Gearbox (float, default = 0.9): Efficiency of the gearbox.
            %             - efficiency_Storage (float, default = 0.9): Efficiency of the energy storage system.
            %             - efficiency_PowerElectronics (float, default = 0.95): Efficiency of power electronics.
            %             - accGravity (float, default = 9.81): Acceleration due to gravity (m/s^2).
            %             - densityAir (float, default = 1.225): Density of air (kg/m^3).
            %             - doUseReferenceWindProfile (bool, default = false): Flag to use a reference wind profile.
            %             - vwCutOut_patternHeight (float, default = 25): Cut-out wind speed at pattern height.
            %
            %     loggerIdentifyer (string, optional): Identifier for the logger. Defaults to 'loggerQSM'.
            %
            % Returns:
            %     obj (KiteQSMsimulation): Instance of the pumping fixed simulation class.
            obj@KiteQSMsimulationBase(inputs, loggerIdentifier);
        end
        
        %% Implement abstract methods
        function set_optimization_defaults_and_parse(obj, inputs, p3)
            % Setting defaults for p2
            nx = ones(1,inputs.numDeltaLelems);
            % Setting defaults for p3
            addParameter(p3, 'x0', [200, ... % deltaL
                deg2rad(30), ... % avgPattEle
                deg2rad(5), ... % coneAngle
                3*inputs.span, ... % Rp_start
                inputs.speedReelout_max*nx, ... % v_i
                inputs.Cl_maxAirfoil*inputs.Cl_eff_F*nx, ... % CL_i
                0.8*nx, ... % v_o
                90*nx, ... % kinematicRatio
                inputs.Cl_maxAirfoil*inputs.Cl_eff_F*nx], ... % CL
                obj.validationFcn.ArrayFloatPositive);
            addParameter(p3, 'lb', [50, ...
                deg2rad(1), ...
                deg2rad(1), ...
                3*inputs.span, ...
                1*nx, ...
                0.1*nx, ...
                0.8*nx, ...
                1*nx, ...
                0.1*nx], obj.validationFcn.ArrayFloatPositive);
            addParameter(p3, 'ub', [500, ...
                deg2rad(90), ...
                deg2rad(60), ...
                100, ...
                inputs.speedReelout_max*nx, ...
                inputs.Cl_maxAirfoil*inputs.Cl_eff_F*nx,...
                inputs.speedReelout_max*nx, ...
                200*nx, ...
                inputs.Cl_maxAirfoil*inputs.Cl_eff_F*nx], ...
                obj.validationFcn.ArrayFloatPositive);
            parse(p3, inputs);
        end
        
        function obj = runsimulation(obj, inputs, optionsFMINCON)
            persistent outputs; % Without specifying outputs as persistent it works slower
            if nargin<2
                optionsFMINCON                      = optimoptions('fmincon');
                optionsFMINCON.Display              = 'none';
                optionsFMINCON.Algorithm            = 'sqp';
                optionsFMINCON.FiniteDifferenceType = 'forward';
                optionsFMINCON.ScaleProblem         = true;
            else
                if ~isa(optionsFMINCON, 'optim.options.Fmincon')
                    error('KiteQSMsimulation:runSimulation Wrong set of options provided')
                end
            end

            % Initial guess for optimization
            x0 = inputs.x0;
            xSol = zeros(length(inputs.windSpeedReference), length(x0));

            % Bounds
            lb = inputs.lb;
            ub = inputs.ub;

            % Initialize variables to store results
            exitflag = zeros(1, length(inputs.windSpeedReference));
            optHist = cell(1, length(inputs.windSpeedReference));
            lambda = cell(1, length(inputs.windSpeedReference));

            % Preinitialize the struct with all fields set to empty arrays
            % varNames = {'beta', 'gamma', 'CL', 'd_t', 'l_t_inCycle', ...
            %           'pattGrClr', 'm_t', 'm_eff', 'theta', 'phi', 'chi', ...
            %           'CD_k', 'CD_t', 'CD', 'E', 'h_inCycle', 'Rp', 'vw', ...
            %           'vw_r', 'vw_theta', 'vw_phi', 'va', 'Fa', 'lambda', ...
            %           'vk_tau', 'vk_theta', 'vk_phi', 'W', 'Fg_r', 'Fg_theta', ...
            %           'Fg_phi', 'Fa_theta', 'Fa_phi', 'Fa_r', 'rollAngle', 'va_r', ...
            %           'va_theta', 'va_phi', 'F_dot_v', 'D_r', 'D_theta', 'D_phi', ...
            %           'L_r', 'L_theta', 'L_phi', 'D', 'L', 'Ft', 'E_result', ...
            %           'kByE', 'P_m_o_eff', 'zetaMech', 'P_e_o_eff', 'tPatt', ...
            %           'numOfPatt', 'P_e_avg', 'P_m_avg', 'thrust', 'maxThrust'};
            varNames = {'vw_reference', 'deltaL', 'beta', 'gamma', 'Rp_start', 'vk_r_i', ...
                      'CL_i', 'vk_r', 'kRatio', 'CL', 'l_t_min', ...
                      'pattStartGrClr', 'h_cycleStart', 'l_t_max', ...
                      'pattEndGrClr', 'l_t_avg', 'h_cycleAvg', 'h_cycleEnd', ...
                      'd_t', 'elemDeltaL', 'l_t_inCycle', ...
                      'pattGrClr', 'm_t', 'm_eff', 'theta', 'phi', 'chi', ...
                      'CD_k', 'CD_t', 'CD', 'E', 'h_inCycle', 'Rp', 'vw', ...
                      'vw_r', 'vw_theta', 'vw_phi', 'f', 'va', 'Fa', 'lambda', ...
                      'vk_tau', 'vk_theta', 'vk_phi', 'W', 'Fg_r', 'Fg_theta', ...
                      'Fg_phi', 'Fa_theta', 'Fa_phi', 'Fa_r', 'rollAngle', 'va_r', ...
                      'va_theta', 'va_phi', 'F_dot_v', 'D_r', 'D_theta', 'D_phi', ...
                      'L_r', 'L_theta', 'L_phi', 'D', 'L', 'Ft', 'E_result', ...
                      'kByE', 'P_m_o_eff', 'zetaMech', 'etaGen_o', 'P_e_o_eff', ...
                      'theta_i', 'phi_i', 'vw_i', 'vw_r_i', 'vw_theta_i', ...
                      'vw_phi_i', 'f_i', 'va_r_i', 'va_theta_i', 'va_phi_i', ...
                      'va_i', 'CD_k_i', 'CD_i', 'E_i', 'Fa_i', 'Fg_r_i', ...
                      'Fg_theta_i', 'Fg_phi_i', 'Fa_theta_i', 'Fa_phi_i', ...
                      'Fa_r_i', 'F_dot_v_i', 'D_r_i', 'D_theta_i', 'D_phi_i', ...
                      'L_r_i', 'L_theta_i', 'L_phi_i', 'D_i', 'L_i', 'E_result_i', ...
                      'Ft_i', 'P_m_i_eff', 'etaGen_i', 'P_e_i_eff', 't1', 'to_eff', ...
                      'to', 'P1_m_o', 'P1_e_o', 'P_m_o', 'P_e_o', 't2', 'ti_eff', ...
                      'ti', 'P2_m_i', 'P2_e_i', 'P_m_i', 'P_e_i', 'tCycle', ...
                      'tPatt', 'numOfPatt', 'P_e_avg', 'P_m_avg'};

            % Initialize an empty structure
            outputSave = struct();
            
            % Add fields to the empty structure using the variable names
            % for i = 1:length(varNames)
            %     outputSave.(varNames{i})(length(inputs.windSpeedReference),1) = 0;
            % end
            outputs = [];

            halfRhoS = 0.5*inputs.densityAir*inputs.areaWing;

            % Optimize operation for every wind speed
            for i = 1:length(inputs.windSpeedReference)
                windSpeed = inputs.windSpeedReference(i); 
                % Define objective and constraint functions
                conFunc = @(x) constraints(x, inputs);
                objFunc = @(x) objective(x, inputs, obj.kiteMass, halfRhoS, windSpeed);

                % Run optimization
                try
                    if ~isreal(objFunc(x0))
                        x0=inputs.x0;
                    end
                    [x, ~, exitflag(i), optHist{i}, lambda{i}] = fmincon(objFunc, x0, [], [], [], [], lb, ub, conFunc, optionsFMINCON);
    
                    % Store final results
                    [~] = objFunc(x);
                    for var_i = 1:length(varNames)
                        outputSave.(varNames{var_i})(i,:) = outputs.(varNames{var_i});
                    end
    
                    % Update initial guess for the next iteration
                    if abs(mean((outputs.E_result - outputs.E), 2)) > 1 %0.01
                        x0 = inputs.x0;
                        xSol(i,:) = x0;
                    else
                        x0 = x;
                        xSol(i,:) = x;
                    end
                    
                catch err %#ok<*NASGU>
                    x0 = inputs.x0;
                    xSol(i,:) = x0;
                end
            end
            if inputs.windSpeedReference(1)>inputs.windSpeedReference(end)
                for i = 1 %repeat highest wind speed with solution of the wind speed one step lower
                    x0 = xSol(i+1,:);
                    windSpeed = inputs.windSpeedReference(i); 
                    % Define objective and constraint functions
                    conFunc = @(x) constraints(x, inputs);
                    objFunc = @(x) objective(x, inputs, obj.kiteMass, halfRhoS, windSpeed);
    
                    % Run optimization
                    try
                        [x, ~, exitflag(i), optHist{i}, lambda{i}] = fmincon(objFunc, x0, [], [], [], [], lb, ub, conFunc, optionsFMINCON);
        
                        % Store final results
                        [~] = objFunc(x);
                        for var_i = 1:length(varNames)
                            outputSave.(varNames{var_i})(i,:) = outputs.(varNames{var_i});
                        end
        
                        % Update initial guess for the next iteration
                        if abs(mean((outputs.E_result - outputs.E), 2)) > 1 %0.01
                            xSol(i,:) = x0;
                        else
                            xSol(i,:) = x;
                        end
                    catch err %#ok<*NASGU>
                        xSol(i,:) = x0;
                    end
                end
            end

            % Store optimization details
            obj.outputs = outputSave;
            obj.optimDetails.optHist = optHist;
            obj.optimDetails.exitflag = exitflag;
            obj.optimDetails.lambda = lambda;
            
            function fval = objective(x, inputs, kiteMass, halfRhoS, windSpeedReference)
                numDeltaLelems = inputs.numDeltaLelems;
            
                outputs.vw_reference = windSpeedReference;

                outputs.deltaL = x(1);
                outputs.beta = x(2);
                outputs.gamma = x(3);
                outputs.Rp_start = x(4);
                outputs.vk_r_i(1,1:numDeltaLelems) = x(0*numDeltaLelems+5:1*numDeltaLelems+4);
                outputs.CL_i(1,1:numDeltaLelems) = x(1*numDeltaLelems+5:2*numDeltaLelems+4);
                outputs.vk_r(1,1:numDeltaLelems) = x(2*numDeltaLelems+5:3*numDeltaLelems+4);
                outputs.kRatio(1,1:numDeltaLelems) = x(3*numDeltaLelems+5:4*numDeltaLelems+4);
                outputs.CL(1,1:numDeltaLelems) = x(4*numDeltaLelems+5:5*numDeltaLelems+4);
    
                % Constant for the cycle
                outputs.l_t_min = outputs.Rp_start/sin(outputs.gamma);
                outputs.pattStartGrClr = outputs.l_t_min*sin(outputs.beta-outputs.gamma);
                outputs.h_cycleStart = outputs.l_t_min*cos(outputs.gamma)*sin(outputs.beta);
                outputs.l_t_max = outputs.l_t_min+outputs.deltaL/cos(outputs.gamma);
                outputs.pattEndGrClr = outputs.l_t_max*sin(outputs.beta-outputs.gamma);
                outputs.l_t_avg = (outputs.l_t_max+outputs.l_t_min)/2; %[m]
                outputs.h_cycleAvg = outputs.l_t_avg*cos(outputs.gamma)*sin(outputs.beta);
                outputs.h_cycleEnd = outputs.l_t_max*cos(outputs.gamma)*sin(outputs.beta);
                outputs.d_t = sqrt(inputs.forceTether_max/inputs.maxStrengthTether*4/pi); %[m] safety factor could be added (say *1.1)
                
                % Discretizing the reel-out length in chosen number of elements
                % Found to be not highly sensitive to the number of elements
                outputs.elemDeltaL = outputs.deltaL/numDeltaLelems;
                
                % Assigning and evaluating a single flight state equilibrium for each length element for reel-out and reel-in phase
                for j = 1:numDeltaLelems
                    % Reel-out phase:
                    % Tether length at jth element
                    if j == 1
                        outputs.l_t_inCycle(1,j) = outputs.l_t_min + outputs.elemDeltaL/2/cos(outputs.gamma);
                    else
                        outputs.l_t_inCycle(1,j) = outputs.l_t_inCycle(1,j-1) + outputs.elemDeltaL/cos(outputs.gamma);
                    end
                    
                    % Pattern ground clearance
                    outputs.pattGrClr(1,j) = outputs.l_t_inCycle(1,j)*sin(outputs.beta-outputs.gamma);
                
                    % Effective mass lumped at kite point (Kite + tether)
                    outputs.m_t(1,j) = inputs.densityTether*pi/4*outputs.d_t^2*outputs.l_t_inCycle(1,j);
                    outputs.m_eff(1,j) = kiteMass + 0.5*outputs.m_t(1,j);
                
                    % Coordinates of the Kite's position and orientation in the Spherical ref. frame
                    % Center point (representative point)
                    outputs.theta(1,j) = pi/2 - (outputs.beta);
                    outputs.phi(1,j) = asin(4*sin(outputs.gamma)/3/pi); % Centroid of semicircle
                    outputs.chi(1,j) = pi/2;
                
                    % Effective CD
                    % outputs.CD_k(1,j) = inputs.Cd0 + (outputs.CL(1,j))^2/(pi*inputs.aspectRatio*inputs.e);
                    outputs.CD_k(1,j) = inputs.Cd0 + (outputs.CL(1,j)-inputs.Cl0_airfoil)^2/(pi()*inputs.aspectRatio*inputs.e) + inputs.Cd_additional;
                    outputs.CD_t(1,j) = (1/4)*inputs.Cd_cylinder*outputs.d_t*outputs.l_t_inCycle(1,j)/inputs.areaWing;
                    outputs.CD(1,j) = outputs.CD_k(1,j) + outputs.CD_t(1,j);
                
                    % Effective Glide ratio
                    outputs.E(1,j) = outputs.CL(1,j)/outputs.CD(1,j);
                
                    % Average pattern height at point of interest on deltaL
                    if j == 1
                        outputs.h_inCycle(1,j) = outputs.h_cycleStart + outputs.elemDeltaL/2*sin(outputs.beta);
                    else
                        outputs.h_inCycle(1,j) = outputs.h_inCycle(1,j-1) + outputs.elemDeltaL*sin(outputs.beta);
                    end
                
                    % Pattern radius at point of interest on deltaL
                    if j == 1
                        outputs.Rp(1,j) = outputs.Rp_start + outputs.elemDeltaL/2*tan(outputs.gamma);
                    else
                        outputs.Rp(1,j) = outputs.Rp(1,j-1) + outputs.elemDeltaL*tan(outputs.gamma);
                    end
                    
                    if ~inputs.doUseReferenceWindProfile
                        % Modelled
                        outputs.vw(1,j) = windSpeedReference*((outputs.h_inCycle(1,j))/inputs.heightWindReference)^inputs.windShearExp;
                    else
                        % Extrapolated from dataset
                        outputs.vw(1,j) = windSpeedReference*interp1(inputs.windProfile_h, inputs.windProfile_vw, (outputs.h_inCycle(1,j)), 'linear', 'extrap');
                    end
                
                    % Wind velocity vector
                    outputs.vw_r(1,j) = outputs.vw(1,j)*sin(outputs.theta(1,j))*cos(outputs.phi(1,j));
                    outputs.vw_theta(1,j) = outputs.vw(1,j)*cos(outputs.theta(1,j))*cos(outputs.phi(1,j));
                    outputs.vw_phi(1,j) = -outputs.vw(1,j)*sin(outputs.phi(1,j));
                
                    % Reel-out factor
                    outputs.f(1,j) = outputs.vk_r(1,j)/outputs.vw(1,j);
                
                    % Apparent wind velocity magnitude
                    outputs.va(1,j) = (outputs.vw_r(1,j)-outputs.vk_r(1,j))*sqrt(1+outputs.kRatio(1,j)^2);
                
                    % Aerodynamic force magnitude
                    outputs.Fa(1,j) = halfRhoS*sqrt(outputs.CL(1,j)^2+outputs.CD(1,j)^2)*outputs.va(1,j)^2;
                
                    % Tangential kite velocity factor
                    a = cos(outputs.theta(1,j))*cos(outputs.phi(1,j))*cos(outputs.chi(1,j))-sin(outputs.phi(1,j))*sin(outputs.chi(1,j));
                    b = sin(outputs.theta(1,j))*cos(outputs.phi(1,j));
                    outputs.lambda(1,j) = a + sqrt(a^2 + b^2 - 1 + outputs.kRatio(1,j)^2*(b - outputs.f(1,j))^2);
                    % Tangential kite velocity
                    outputs.vk_tau(1,j) = outputs.lambda(1,j)*outputs.vw(1,j);
                    % Tangential kite velocity components in theta and phi directions
                    outputs.vk_theta(1,j) = outputs.vk_tau(1,j)*cos(outputs.chi(1,j));
                    outputs.vk_phi(1,j) = outputs.vk_tau(1,j)*sin(outputs.chi(1,j));
                
                    % Gravity toggle
                    if inputs.doIncludeGravity
                        outputs.W(1,j) = outputs.m_eff(1,j)*inputs.accGravity;
                    else
                        outputs.W(1,j) = 0;
                    end
                
                    % Gravitational force vector (kite + tether)
                    outputs.Fg_r(1,j) = -outputs.W(1,j)*cos(outputs.theta(1,j));
                    outputs.Fg_theta(1,j) = outputs.W(1,j)*sin(outputs.theta(1,j));
                    outputs.Fg_phi(1,j) = 0;
                
                    % Aerodynamic force vector
                    outputs.Fa_theta(1,j) = -outputs.Fg_theta(1,j);
                    outputs.Fa_phi(1,j) = -outputs.Fg_phi(1,j);
                    outputs.Fa_r(1,j) = sqrt(outputs.Fa(1,j)^2-outputs.Fa_theta(1,j)^2-outputs.Fa_phi(1,j)^2);
                
                    % Roll angle in case of Gravity (N.A when CF is included)
                    outputs.rollAngle(1,j) = atan(outputs.Fa_theta(1,j)/outputs.Fa_r(1,j));
                
                    % Apparent wind velocity vector
                    outputs.va_r(1,j) = outputs.vw(1,j)*(sin(outputs.theta(1,j))*cos(outputs.phi(1,j))-outputs.f(1,j));
                    outputs.va_theta(1,j) = outputs.vw(1,j)*(cos(outputs.theta(1,j))*cos(outputs.phi(1,j))-outputs.lambda(1,j)*cos(outputs.chi(1,j)));
                    outputs.va_phi(1,j) = outputs.vw(1,j)*(-sin(outputs.phi(1,j))-outputs.lambda(1,j)*sin(outputs.chi(1,j)));
                
                    % Dot product F_a*v_a;
                    outputs.F_dot_v(1,j) = outputs.Fa_r(1,j)*outputs.va_r(1,j) + outputs.Fa_theta(1,j)*outputs.va_theta(1,j) + outputs.Fa_phi(1,j)*outputs.va_phi(1,j);
                
                    % Drag vector from dot product
                    outputs.D_r(1,j) = (outputs.F_dot_v(1,j)/outputs.va(1,j)^2)*outputs.va_r(1,j);
                    outputs.D_theta(1,j) = (outputs.F_dot_v(1,j)/outputs.va(1,j)^2)*outputs.va_theta(1,j);
                    outputs.D_phi(1,j) = (outputs.F_dot_v(1,j)/outputs.va(1,j)^2)*outputs.va_phi(1,j);
                
                    % Lift vector
                    outputs.L_r(1,j) = outputs.Fa_r(1,j) - outputs.D_r(1,j);
                    outputs.L_theta(1,j) = outputs.Fa_theta(1,j) - outputs.D_theta(1,j);
                    outputs.L_phi(1,j) = outputs.Fa_phi(1,j) - outputs.D_phi(1,j);
                
                    % Drag magnitude
                    outputs.D(1,j) = sqrt(outputs.D_r(1,j)^2 + outputs.D_theta(1,j)^2 + outputs.D_phi(1,j)^2);
                
                    % Lift magnitude
                    outputs.L(1,j) = sqrt(outputs.L_r(1,j)^2 + outputs.L_theta(1,j)^2 + outputs.L_phi(1,j)^2);
                
                    % Straight-tether force
                    outputs.Ft(1,j) = outputs.Fa_r(1,j) + outputs.Fg_r(1,j);
                
                    % Lift-to-drag ratio that follows from the chosen kinematic ratio
                    outputs.E_result(1,j) = sqrt(((outputs.Fa(1,j)*outputs.va(1,j))/outputs.F_dot_v(1,j))^2-1);
                
                    % Ratio of Kinemactic Ratio and Glide Ratio
                    outputs.kByE(1,j) = outputs.kRatio(1,j)/outputs.E_result(1,j);
                
                    % Effective mechanical reel-out power
                    outputs.P_m_o_eff(1,j) = outputs.Ft(1,j)*outputs.vk_r(1,j); %[W]
                
                    outputs.zetaMech(1,j) = outputs.P_m_o_eff(1,j)/(halfRhoS*outputs.vw(1,j)^3);
                
                    % Effective electrical reel-out power
                    % Generator efficiency. As a function of RPM/RPM_max, where RPM_max is driven by winch
                    vkr_ratio = outputs.vk_r(1,j)/inputs.speedReelout_max;
                    outputs.etaGen_o(1,j) = inputs.etaGen_param(1)*vkr_ratio^3 + ...
                        inputs.etaGen_param(2)*vkr_ratio^2 + ...
                        inputs.etaGen_param(3)*vkr_ratio+inputs.etaGen_param(4);
                    outputs.P_e_o_eff(1,j) = outputs.P_m_o_eff(1,j)*inputs.efficiency_Gearbox*outputs.etaGen_o(1,j)*inputs.efficiency_PowerElectronics*inputs.efficiency_Storage*inputs.efficiency_PowerElectronics;
                
                    %% Retraction Phase
                    % Kite position in spherical coordinates
                    % Reel-in is assumed to start from the top of the pattern
                    outputs.theta_i(1,j) = pi/2 - (outputs.beta+outputs.gamma);
                    outputs.phi_i(1,j) = 0;
                
                    if ~inputs.doUseReferenceWindProfile
                        % Modelled
                        outputs.vw_i(1,j) = windSpeedReference*((outputs.h_inCycle(1,j) + outputs.Rp(1,j)*cos(outputs.beta))/inputs.heightWindReference)^inputs.windShearExp;
                    else
                        % Extrapolated from dataset
                        outputs.vw_i(1,j) = windSpeedReference*interp1(inputs.windProfile_h, inputs.windProfile_vw, (outputs.h_inCycle(1,j) + outputs.Rp(1,j)*cos(outputs.beta)), 'linear', 'extrap');
                    end
                
                    % Wind velocity vector
                    outputs.vw_r_i(1,j) = outputs.vw_i(1,j)*sin(outputs.theta_i(1,j))*cos(outputs.phi_i(1,j));
                    outputs.vw_theta_i(1,j) = outputs.vw_i(1,j)*cos(outputs.theta_i(1,j))*cos(outputs.phi_i(1,j));
                    outputs.vw_phi_i(1,j) = -outputs.vw_i(1,j)*sin(outputs.phi_i(1,j)); %
                
                    % Reel-in factor
                    outputs.f_i(1,j) = outputs.vk_r_i(1,j)/outputs.vw_i(1,j);
                
                    % Apparent speed vector
                    outputs.va_r_i(1,j) = outputs.vw_r_i(1,j) - (-outputs.vk_r_i(1,j)); % Reel-in speed direction is in the negative radial direction)
                    outputs.va_theta_i(1,j) = outputs.vw_theta_i(1,j) - 0;
                    outputs.va_phi_i(1,j) = outputs.vw_phi_i(1,j) - 0;
                
                    % Apparent wind velocity magnitude
                    outputs.va_i(1,j) = sqrt(outputs.va_r_i(1,j)^2 + outputs.va_theta_i(1,j)^2); % Wind has components only in r and theta directions
                
                    % Aerodynamic force magnitude
                    % outputs.CD_k_i(1,j) = inputs.Cd0+(outputs.CL_i(1,j))^2/(pi*inputs.aspectRatio*inputs.e);
                    outputs.CD_k_i(1,j) = inputs.Cd0 + (outputs.CL_i(1,j)-inputs.Cl0_airfoil)^2/(pi()*inputs.aspectRatio*inputs.e) + inputs.Cd_additional;
                    outputs.CD_i(1,j) = outputs.CD_k_i(1,j) + outputs.CD_t(1,j);
                    outputs.E_i(1,j) = outputs.CL_i(1,j)/outputs.CD_i(1,j);
                    outputs.Fa_i(1,j) = halfRhoS*sqrt(outputs.CL_i(1,j)^2+outputs.CD_i(1,j)^2)*outputs.va_i(1,j)^2;
                
                    % Gravitational force vector (kite + tether)
                    outputs.Fg_r_i(1,j) = -outputs.W(1,j)*cos(outputs.theta_i(1,j));
                    outputs.Fg_theta_i(1,j) = outputs.W(1,j)*sin(outputs.theta_i(1,j));
                    outputs.Fg_phi_i(1,j) = 0;
                
                    % Aerodynamic force vector
                    outputs.Fa_theta_i(1,j) = -outputs.Fg_theta_i(1,j);
                    outputs.Fa_phi_i(1,j) = -outputs.Fg_phi_i(1,j);
                    outputs.Fa_r_i(1,j) = sqrt(outputs.Fa_i(1,j)^2-outputs.Fa_theta_i(1,j)^2-outputs.Fa_phi_i(1,j)^2);
                
                    % Dot product F_a*v_a;
                    outputs.F_dot_v_i(1,j) = outputs.Fa_r_i(1,j)*outputs.va_r_i(1,j) + outputs.Fa_theta_i(1,j)*outputs.va_theta_i(1,j) + outputs.Fa_phi_i(1,j)*outputs.va_phi_i(1,j);
                
                    % Drag vector from dot product
                    outputs.D_r_i(1,j) = (outputs.F_dot_v_i(1,j)/outputs.va_i(1,j)^2)*outputs.va_r_i(1,j);
                    outputs.D_theta_i(1,j) = (outputs.F_dot_v_i(1,j)/outputs.va_i(1,j)^2)*outputs.va_theta_i(1,j);
                    outputs.D_phi_i(1,j) = (outputs.F_dot_v_i(1,j)/outputs.va_i(1,j)^2)*outputs.va_phi_i(1,j);
                
                    % Lift vector
                    outputs.L_r_i(1,j) = outputs.Fa_r_i(1,j) - outputs.D_r_i(1,j);
                    outputs.L_theta_i(1,j) = outputs.Fa_theta_i(1,j) - outputs.D_theta_i(1,j);
                    outputs.L_phi_i(1,j) = outputs.Fa_phi_i(1,j) - outputs.D_phi_i(1,j);
                
                    % Drag magnitude
                    outputs.D_i(1,j) = sqrt(outputs.D_r_i(1,j)^2 + outputs.D_theta_i(1,j)^2 + outputs.D_phi_i(1,j)^2);
                
                    % Lift magnitude
                    outputs.L_i(1,j) = sqrt(outputs.L_r_i(1,j)^2 + outputs.L_theta_i(1,j)^2 + outputs.L_phi_i(1,j)^2);
                
                    % Lift-to-drag ratio that follows from the chosen kinematic ratio
                    outputs.E_result_i(1,j) = sqrt(((outputs.Fa_i(1,j)*outputs.va_i(1,j))/outputs.F_dot_v_i(1,j))^2-1);
                
                    % Straight-tether force
                    outputs.Ft_i(1,j) = outputs.Fa_r_i(1,j) + outputs.Fg_r_i(1,j);
                
                    % Effective mechanical reel-out power
                    outputs.P_m_i_eff(1,j) = outputs.Ft_i(1,j)*outputs.vk_r_i(1,j); %[W]
                
                    % Generator efficiency during RI: As a function of RPM/RPM_max, where RPM_max is driven by winch i.e Max VRI
                    vkri_ratio = outputs.vk_r_i(1,j)/inputs.speedReelout_max;
                    outputs.etaGen_i = inputs.etaGen_param(1)*vkri_ratio^3 + ...
                        inputs.etaGen_param(2)*vkri_ratio^2 + ...
                        inputs.etaGen_param(3)*vkri_ratio+inputs.etaGen_param(4);
                
                    % Effective electrical reel-in power
                    outputs.P_e_i_eff(1,j) = outputs.P_m_i_eff(1,j)/inputs.efficiency_Gearbox/inputs.efficiency_Storage/outputs.etaGen_i/inputs.efficiency_PowerElectronics;
                
                end
    
                
                %% Cycle calculation
                % Reel-out time
                outputs.t1       = outputs.vk_r(1)/inputs.accReel_max;
                outputs.to_eff = outputs.elemDeltaL./outputs.vk_r;
                outputs.to       = outputs.t1 + sum(outputs.to_eff);
                
                % Reel-out power during transition
                outputs.P1_m_o = outputs.P_m_o_eff(1)/2;
                outputs.P1_e_o = outputs.P_e_o_eff(1)/2;
                
                % Reel-out power 
                outputs.P_m_o = (sum(outputs.P_m_o_eff.*outputs.to_eff) + outputs.P1_m_o*outputs.t1)/outputs.to;
                outputs.P_e_o = (sum(outputs.P_e_o_eff.*outputs.to_eff) + outputs.P1_e_o*outputs.t1)/outputs.to;
                
                % Reel-in time
                outputs.t2       = outputs.vk_r_i(1)/inputs.accReel_max;
                outputs.ti_eff = outputs.elemDeltaL./outputs.vk_r_i;
                outputs.ti       = outputs.t2 + sum(outputs.ti_eff);
                
                % Reel-in power during transition
                outputs.P2_m_i     = outputs.P_m_i_eff(1)/2;
                outputs.P2_e_i     = outputs.P_e_i_eff(1)/2;
                
                % Reel-in power 
                outputs.P_m_i = (sum(outputs.P_m_i_eff.*outputs.ti_eff) + outputs.P2_m_i*outputs.t2)/outputs.ti;
                outputs.P_e_i = (sum(outputs.P_e_i_eff.*outputs.ti_eff) + outputs.P2_e_i*outputs.t2)/outputs.ti;
                
                % Cycle time
                outputs.tCycle = outputs.to+outputs.ti;
                
                % Time for one pattern revolution and number of patterns in the cycle
                outputs.tPatt = 2*pi*outputs.Rp./outputs.vk_phi;
                outputs.numOfPatt = outputs.to./outputs.tPatt;
                
                % Electrical cycle power
                outputs.P_e_avg = (sum(outputs.to_eff.*outputs.P_e_o_eff) + outputs.t1*outputs.P1_e_o - ...
                                           sum(outputs.ti_eff.*outputs.P_e_i_eff) - outputs.t2*outputs.P2_e_i)/outputs.tCycle;
                
                % Mechanical cycle power - without drivetrain eff
                outputs.P_m_avg = (sum(outputs.to_eff.*outputs.P_m_o_eff) + outputs.t1*outputs.P1_m_o - ...
                                           sum(outputs.ti_eff.*outputs.P_m_i_eff) - outputs.t2*outputs.P2_m_i)/outputs.tCycle;
    
                fval = -outputs.P_e_avg;

            end
            function [c, ceq] = constraints(x, inputs) %#ok<INUSD>
                numDeltaLelems = inputs.numDeltaLelems;
                % Inequality constraints
    
                % Min clearance between ground and bottom point of pattern
                c(1)   = (inputs.minGroundClear - outputs.pattStartGrClr);
                c(2)   = (inputs.minGroundClear - outputs.pattEndGrClr);
    
                % Capping for requested electrical rated power
                c(3)   = (outputs.P_e_avg - inputs.P_ratedElec); 
    
                % Tether length limit
                c(4) = (outputs.l_t_max - inputs.maxTeLen); 
    
                % Min number of patterns to get into transition 
                c(5) = (1 - mean(outputs.numOfPatt));
    
                % Max. cycle avg height
                c(6) = (outputs.h_cycleEnd - inputs.maxHeight);
    
                % Peak mechanical power limit
                c(7:numDeltaLelems+6) = (outputs.P_m_o_eff - inputs.peakM2E_F*inputs.P_ratedElec);
    
                % Maximum tether force
                c(7+numDeltaLelems:2*numDeltaLelems+6) = (outputs.Ft - inputs.forceTether_max*inputs.safetyFactor_forceTetherMax);
    
                % Apparent speed cannot be negative
                c(7+2*numDeltaLelems:3*numDeltaLelems+6) = (0 - outputs.va);
    
                % Tangential velocity factor cannot be negative
                c(7+3*numDeltaLelems:4*numDeltaLelems+6) = (0 - outputs.lambda);
    
                % Tether force during reel-in cannot be negative
                c(7+4*numDeltaLelems:5*numDeltaLelems+6) = (0 - outputs.Ft_i);  
    
                %% Equality constraints
                % Glide ratio during reel-out
                ceq(0*numDeltaLelems+1:1*numDeltaLelems) = (outputs.E_result - outputs.E);
    
                % Glide ratio during reel-in
                ceq(1*numDeltaLelems+1:2*numDeltaLelems) = (outputs.E_result_i - outputs.E_i);
    
            end
        end
        
        function obj = processoutputs(obj)

            inputs = obj.inputs;
            outputs = obj.outputs;
            
            processedOutputs = struct();
            %% Post processing
            % vw = inputs.windSpeedReference; % Wind speed at ref. height
            vw = outputs.vw_reference; % Wind speed at ref. height
            processedOutputs.vw_reference = vw;
            % On exit flag convergence is too harsch, glide ratio is a good
            % indicator for convergence
            converged_idx = find(abs(mean((outputs.E_result - outputs.E),2)) < 0.01);
            vw_converged = vw(converged_idx);
            %% Cut-in wind speed
            % Glide ratio constraint violation acceptance for feasible solution
            % temp1                  = vw(abs(mean((outputs.E_result - outputs.E),2)) < 0.01); %0.01
            processedOutputs.cutIn = min(vw_converged(outputs.P_e_avg(converged_idx)>0)); %#ok<*PROPLC>
            
            %% Rated wind and power
            vw_above_rated              = vw(round(outputs.P_e_avg./max(outputs.P_e_avg),2)>=0.99);
            processedOutputs.ratedWind  = vw_above_rated(1); % At Reference height 
            processedOutputs.ratedPower = outputs.P_e_avg(vw==processedOutputs.ratedWind);
            
            %% Cut-out wind speed
            processedOutputs.cutOut   = inputs.vwCutOut_patternHeight; % maximum at pattern height
            % Operating range traslated to wind speed at ref. height
            processedOutputs.vw_100m_operRange = vw((vw>=processedOutputs.cutIn) & mean(outputs.vw,2)<=processedOutputs.cutOut);
            processedOutputs.cutOut_referenceHeight = max(processedOutputs.vw_100m_operRange);
            %% Extract feasible results in the operational range
            processedOutputs.Dia_te = outputs.d_t(1); % always constant

            variableNamesToCopy = {'vw','vw_i','P1_m_o','P2_m_i','P_m_o_eff',...
                        'P_m_i_eff','P1_e_o', 'P2_e_i','P_e_o_eff','P_e_i_eff',...
                        'P_m_o','P_m_i','P_e_o','P_e_i','deltaL','t1','t2',...
                        'to_eff','ti_eff','to','ti','tCycle','Ft','Ft_i','Fa',...
                        'Fa_i','W','vk_r','vk_tau','vk_r_i','P_e_avg','P_m_avg','Rp',...
                        'h_cycleStart','h_cycleAvg','l_t_max','l_t_min','l_t_inCycle',...
                        'h_inCycle','va','beta','rollAngle','gamma','CL','CD',...
                        'CD_k','CD_k_i','CD_t','CL_i','CD_i','E','E_i','numOfPatt',...
                        'lambda','f','f_i','tPatt','zetaMech','d_t'};

            variableNamesToInit_scalar = {'dutyCycle', 'cycleEff_elec', 'cycleEff_mech', ...
                'storageExchange'};
            variableNamesToInit_array = {'sweptArea', 'Cp_m_o', 'Cp_e_avg'};

            for i=1:length(vw)
                %vw is reference height so cut_out at reference height should be taken
                if vw(i)>=processedOutputs.cutIn && vw(i)<=processedOutputs.cutOut
                    % Define the list of variable names
                    processedOutputs = KiteQSMsimulationBase.copy_outputs_to_processedoutputs(...
                        processedOutputs, outputs, variableNamesToCopy, i);

                    %convert angles to degrees
                    processedOutputs.beta(i)           = rad2deg(processedOutputs.beta(i));
                    processedOutputs.rollAngle(i,:)    = rad2deg(-processedOutputs.rollAngle(i,:));
                    processedOutputs.gamma(i)          = rad2deg(processedOutputs.gamma(i));

                    % calculate dutycycle
                    processedOutputs.dutyCycle(i)      = processedOutputs.to(i)/processedOutputs.tCycle(i);

                    % Cycle efficiency
                    processedOutputs.cycleEff_elec(i)  = (processedOutputs.P_e_avg(i)*processedOutputs.tCycle(i))/...
                                                            (processedOutputs.P_e_o(i)*processedOutputs.to(i));
                    processedOutputs.cycleEff_mech(i)  = (processedOutputs.P_m_avg(i)*processedOutputs.tCycle(i))/...
                                                            (processedOutputs.P_m_o(i)*processedOutputs.to(i));
                    
                    processedOutputs.storageExchange(i) = (processedOutputs.P_m_avg(i)+processedOutputs.P_m_i(i))*processedOutputs.ti(i)/3.6e3;  %[Wh]

                    % Coefficient of power (Cp) as defined for HAWTs
                    processedOutputs.sweptArea(i,:)     = pi.*((outputs.Rp(i,:)+inputs.span/2).^2 - (outputs.Rp(i,:)-inputs.span/2).^2);
                    processedOutputs.Cp_m_o(i,:)        = outputs.P_m_o(i)./(0.5.*inputs.densityAir.*processedOutputs.sweptArea(i,:).*outputs.vw(i,:).^3);
                    processedOutputs.Cp_e_avg(i,:)      = outputs.P_e_avg(i)./(0.5.*inputs.densityAir.*processedOutputs.sweptArea(i,:).*outputs.vw(i,:).^3);
                else
                    processedOutputs = obj.copy_zeros_to_processedoutputs_based_on_Size(...
                        processedOutputs, [1,1], variableNamesToInit_scalar, i);
                    processedOutputs = obj.copy_zeros_to_processedoutputs_based_on_Size(...
                        processedOutputs, size(outputs.Rp(i,:)), variableNamesToInit_array, i);
                    processedOutputs = obj.copy_zeros_to_processedoutputs(...
                        processedOutputs, outputs, variableNamesToCopy, i);
                end
            end

            
            
            if processedOutputs.vw_reference(1)>processedOutputs.vw_reference(end)
                variableNamesProcessed = fields(processedOutputs);
                for var_i = 1:numel(variableNamesProcessed)
                    variableName = variableNamesProcessed{var_i};
                    if size(processedOutputs.(variableName),1)==1
                        processedOutputs.(variableName) = fliplr(processedOutputs.(variableName));
                    elseif size(processedOutputs.(variableName),2)==1
                        processedOutputs.(variableName) = flipud(processedOutputs.(variableName));
                    elseif size(processedOutputs.(variableName),1)==inputs.numDeltaLelems
                        processedOutputs.(variableName) = fliplr(processedOutputs.(variableName));
                    elseif size(processedOutputs.(variableName),2)==inputs.numDeltaLelems
                        processedOutputs.(variableName) = flipud(processedOutputs.(variableName));
                    end
                end
            end

            obj.processedOutputs = processedOutputs;

            
        end
    end
    
    methods (Static)
        %% Implement abstract static methods
        function requiredParams = get_required_params()
            requiredParams = {'windSpeedReference', 'heightWindReference', ...
                'areaWing', 'span', 'forceTether_max', 'P_ratedElec', ...
                'Cl_maxAirfoil', 'Cl_eff_F', ...
                'Cd0', 'Cl0_airfoil'};
        end
        
        function parametersWithDefaults = get_parameters_withdefaults(inputs,validationFcn)
            parametersWithDefaults = {
                'aspectRatio', (inputs.span^2/inputs.areaWing), validationFcn.FloatPositive;
                'numDeltaLelems', 5, validationFcn.Int;
                'doIncludeGravity', true, validationFcn.Bool;
                'doMassOverride', false, validationFcn.Bool;
                'maxTeLen', 1000, validationFcn.FloatPositive;
                'maxHeight', 1000, validationFcn.FloatPositive;
                'minGroundClear', 100, validationFcn.FloatPositive;
                'safetyFactor_forceTetherMax', 0.8, validationFcn.FloatPositive;
                'peakM2E_F', 2, validationFcn.FloatPositive;
                'maxStrengthTether', 7e8, validationFcn.FloatPositive;
                'densityTether', 980, validationFcn.FloatPositive;
                'speedReelout_max', 20, validationFcn.FloatPositive;
                'accReel_max', 5, validationFcn.FloatPositive;
                'e', 0.6, validationFcn.FloatPositive;
                'Cd_cylinder', 1.2, validationFcn.FloatPositive;
                'Cd_additional',0.0, validationFcn.FloatPositiveOrZero;
                'etaGen_param', [0.671, -1.4141, 0.9747, 0.7233], validationFcn.etaGen_param;
                'efficiency_Gearbox', 0.9, validationFcn.FloatPositive;
                'efficiency_Storage', 0.9, validationFcn.FloatPositive;
                'efficiency_PowerElectronics', 0.95, validationFcn.FloatPositive;
                'accGravity', 9.81, validationFcn.FloatPositive;
                'densityAir', 1.225, validationFcn.FloatPositive;
                'doUseReferenceWindProfile', false, validationFcn.Bool;
                'vwCutOut_patternHeight',25,validationFcn.FloatPositive;
                'name','simulationQSM_PumpingFixed',validationFcn.String;
                'kiteMass',0,validationFcn.FloatPositiveOrZero;
            };
        end
    end

end
