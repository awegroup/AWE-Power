classdef KiteQSMsimulation_flygen < KiteQSMsimulationBase
    methods
        function obj = KiteQSMsimulation_flygen(inputs, loggerIdentifier)
            % Class Constructor for KiteQSMsimulation_flygen
            %
            % Args:
            %     inputs (struct): A structure containing the required and optional simulation parameters.
            %         Required Parameters:
            %             - windSpeedReference (float): Wind speed at reference height (m/s).
            %             - heightWindReference (float): Height at which wind speed is measured (m).
            %             - areaWing (float): Area of the wing (m²).
            %             - span (float): Wing span (m).
            %             - forceTether_max (float): Maximum tether force (N).
            %             - P_ratedElec (float): Rated electrical power (W).
            %             - Cl_maxAirfoil (float): Maximum coefficient of lift of the airfoil.
            %             - Cl_eff_F (float): Effective lift coefficient factor.
            %             - Cd0 (float): Zero-lift drag coefficient.
            %             - numTurbines (int): Number of turbines.
            %             - areaTurbine (float): Area of each turbine (m²).
            %
            %         Optional Parameters:
            %             - aspectRatio (float, default = span^2/areaWing): Aspect ratio of the wing.
            %             - doIncludeGravity (bool, default = true): Flag to include gravity effect.
            %             - doMassOverride (bool, default = false): Flag to override the kite mass.
            %             - maxTeLen (float, default = 1000): Maximum tether length (m).
            %             - minRp (float, default = 3 * span): Minimum radius of pattern (m).
            %             - maxHeight (float, default = 1000): Maximum height (m).
            %             - minGroundClear (float, default = 100): Minimum ground clearance (m).
            %             - safetyFactor_forceTetherMax (float, default = 0.8): Safety factor for max tether force.
            %             - peakM2E_F (float, default = 2.5): Peak mechanical-to-electric energy factor.
            %             - maxStrengthTether (float, default = 7e8): Maximum strength of the tether (Pa).
            %             - densityTether (float, default = 980): Density of the tether material (kg/m³).
            %             - e (float, default = 0.6): Oswald efficiency factor.
            %             - Cd_cylinder (float, default = 1.1): Drag coefficient of tether.
            %             - efficiency_PowerElectronics (float, default = 0.95): Efficiency of power electronics.
            %             - accGravity (float, default = 9.81): Acceleration due to gravity (m/s²).
            %             - densityAir (float, default = 1.225): Density of air (kg/m³).
            %             - doUseReferenceWindProfile (bool, default = false): Flag to use reference wind profile.
            %             - vwCutOut_patternHeight (float, default = 25): Cut-out wind speed at pattern height.
            %             - Cd_additional (float, default = 0): Additional drag coefficient.
            %
            %     loggerIdentifyer (string, optional): Identifier for the logger. Defaults to 'loggerQSM'.
            %
            % Returns:
            %     obj (KiteQSMsimulation_flygen): Instance of the class.
            obj@KiteQSMsimulationBase(inputs, loggerIdentifier);
        end
        
        %% Implement abstract methods
        function set_optimization_defaults_and_parse(obj, inputs, p3)
            % Use the flygen-specific optimization defaults.
            addParameter(p3, 'x0', [deg2rad(30), deg2rad(5), 100, 90, ...
                inputs.Cl_maxAirfoil*inputs.Cl_eff_F, 0.5], obj.validationFcn.ArrayFloatPositive);
            % Add lower and upper bounds, then parse:
            addParameter(p3, 'lb', [deg2rad(1), deg2rad(1), 50, 1, 0.1, 0.1], obj.validationFcn.ArrayFloatPositive);
            addParameter(p3, 'ub', [deg2rad(90), deg2rad(60), 1000, 200, inputs.Cl_maxAirfoil*inputs.Cl_eff_F, 3], obj.validationFcn.ArrayFloatPositive);
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
            x0 = obj.inputs.x0;
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
            varNames = {'vw_reference', 'beta', 'gamma', 'Rp', 'lambda', 'CL', 'l_t_inCycle', ...
                'pattGrClr', 'h_inCycle', 'd_t', 'm_t', 'm_eff', 'theta', ...
                'phi', 'chi', 'CD_k', 'CD_t', 'CD', 'CT', 'E', 'vw', 'vw_r', ...
                'vw_theta', 'vw_phi', 'vk_tau', 'vk_theta', 'vk_phi', 'va', ...
                'Fthrust', 'maxThrust', 'maxPower', 'axialInduction', 'Fa', 'W', 'Fg_r', 'Fg_theta', 'Fg_phi', ...
                'Fa_theta', 'Fa_phi', 'Fa_r', 'rollAngle', 'Ft', 'va_r', 'va_theta', ...
                'va_phi', 'F_dot_v', 'D_r', 'D_theta', 'D_phi', 'L_r', 'L_theta', ...
                'L_phi', 'D', 'L', 'E_result', 'P_m_eff', 'zetaMech', 'P_e_eff', ...
                'tPatt','turbineGain'};

            % Initialize an empty structure
            outputSave = struct();
            
            % Add fields to the empty structure using the variable names
            for i = 1:length(varNames)
                outputSave.(varNames{i})(length(inputs.windSpeedReference),1) = 0;
            end
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
                    [x, ~, exitflag(i), optHist{i}, lambda{i}] = fmincon(objFunc, x0, [], [], [], [], lb, ub, conFunc, optionsFMINCON);
    
                    % Store final results
                    [~] = objFunc(x);
                    for var_i = 1:length(varNames)
                        outputSave.(varNames{var_i})(i) = outputs.(varNames{var_i});
                    end
    
                    % Update initial guess for the next iteration
                    if abs(mean((outputs.E_result - outputs.E), 2)) > 0.2 %0.01
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
                            outputSave.(varNames{var_i})(i) = outputs.(varNames{var_i});
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

                outputs.vw_reference = windSpeedReference;

                outputs.beta = x(1);
                outputs.gamma = x(2);
                outputs.l_t_inCycle = x(3);
                outputs.lambda = x(4);
                outputs.CL = x(5);
                outputs.turbineGain = x(6);
                
                % Constant for the cycle
                outputs.Rp = outputs.l_t_inCycle*sin(outputs.gamma);
                outputs.pattGrClr = outputs.l_t_inCycle*sin(outputs.beta-outputs.gamma);
                outputs.h_inCycle = outputs.l_t_inCycle*cos(outputs.gamma)*sin(outputs.beta);
                outputs.d_t = sqrt(inputs.forceTether_max/inputs.maxStrengthTether*4/pi); %[m] safety factor could be added (say *1.1)
                    
                % Tether core mass
                outputs.m_t = inputs.densityTether*pi/4*outputs.d_t^2*outputs.l_t_inCycle;

                % increase diameter by a factor to account for electrical
                % wires.
                outputs.d_t = outputs.d_t*inputs.tetherCoreIncreaseFactor;

                % simplified model given Bauer, F. PhD dissertation
                % core is only 60% of total mass.
                outputs.m_t = outputs.m_t/0.6; 
                
                % Effective mass lumped at kite point (Kite + tether)
                outputs.m_eff = kiteMass + 0.5*outputs.m_t;
            
                % Coordinates of the Kite's position and orientation in the Spherical ref. frame
                % Center point (representative point)
                outputs.theta = pi/2 - (outputs.beta);
                outputs.phi = asin(4*sin(outputs.gamma)/3/pi); % Centroid of semicircle
                outputs.chi = pi/2;
            
                % Effective CD
                outputs.CD_k = inputs.Cd0 + (outputs.CL-inputs.Cl0_airfoil)^2/(pi*inputs.aspectRatio*inputs.e) + inputs.Cd_additional;
                outputs.CD_t = (1/4)*inputs.Cd_cylinder*outputs.d_t*outputs.l_t_inCycle/inputs.areaWing;
            
                outputs.CD = outputs.CD_k + outputs.CD_t;
                
                % turbine thrust
                outputs.CT = outputs.turbineGain*outputs.CD_k;
                
                % Effective Glide ratio with turbine thrust
                outputs.E = outputs.CL/(outputs.CD+outputs.CT);
            
                if ~inputs.doUseReferenceWindProfile
                    % Modelled
                    outputs.vw = windSpeedReference*(outputs.h_inCycle/inputs.heightWindReference)^inputs.windShearExp;
                else
                    % Extrapolated from dataset
                    outputs.vw = windSpeedReference*interp1(inputs.windProfile_h, inputs.windProfile_vw, outputs.h_inCycle, 'linear', 'extrap');
                end
            
                % Wind velocity vector
                outputs.vw_r = outputs.vw*sin(outputs.theta)*cos(outputs.phi);
                outputs.vw_theta = outputs.vw*cos(outputs.theta)*cos(outputs.phi);
                outputs.vw_phi = -outputs.vw*sin(outputs.phi);
                
                % Tangential kite velocity
                outputs.vk_tau = outputs.lambda*outputs.vw;
                % Tangential kite velocity components in theta and phi directions
                outputs.vk_theta = outputs.vk_tau*cos(outputs.chi);
                outputs.vk_phi = outputs.vk_tau*sin(outputs.chi);
            
                % Apparent wind velocity magnitude
                outputs.va = outputs.vw*sqrt(1+outputs.lambda^2);
            
                outputs.Fthrust = outputs.CT*halfRhoS*outputs.va^2;
            
                % Aerodynamic force magnitude + turbine thrust
                outputs.Fa = halfRhoS*sqrt(outputs.CL^2+(outputs.CD+outputs.CT)^2)*outputs.va^2;
                
                % Gravity toggle
                if inputs.doIncludeGravity
                    outputs.W = outputs.m_eff*inputs.accGravity;
                else
                    outputs.W = 0;
                end
            
                % Gravitational force vector (kite + tether)
                outputs.Fg_r = -outputs.W*cos(outputs.theta);
                outputs.Fg_theta = outputs.W*sin(outputs.theta);
                outputs.Fg_phi = 0;
            
                % Aerodynamic force vector
                outputs.Fa_theta = -outputs.Fg_theta;
                outputs.Fa_phi = -outputs.Fg_phi;
                outputs.Fa_r = sqrt(outputs.Fa^2-outputs.Fa_theta^2-outputs.Fa_phi^2);
            
                % Roll angle in case of Gravity (N.A when CF is included)
                outputs.rollAngle = atan(outputs.Fa_theta/outputs.Fa_r);
            
                % Straight-tether force
                outputs.Ft = outputs.Fa_r + outputs.Fg_r;
            
                % Apparent wind velocity vector
                outputs.va_r = outputs.vw_r;
                outputs.va_theta = outputs.vw*(cos(outputs.theta)*cos(outputs.phi)-outputs.lambda*cos(outputs.chi));
                outputs.va_phi = outputs.vw*(-sin(outputs.phi)-outputs.lambda*sin(outputs.chi));
            
                % Dot product F_a*v_a;
                outputs.F_dot_v = outputs.Fa_r*outputs.va_r + outputs.Fa_theta*outputs.va_theta + outputs.Fa_phi*outputs.va_phi;
            
                % Drag vector from dot product
                outputs.D_r = (outputs.F_dot_v/outputs.va^2)*outputs.va_r;
                outputs.D_theta = (outputs.F_dot_v/outputs.va^2)*outputs.va_theta;
                outputs.D_phi = (outputs.F_dot_v/outputs.va^2)*outputs.va_phi;
            
                % Lift vector
                outputs.L_r = outputs.Fa_r - outputs.D_r;
                outputs.L_theta = outputs.Fa_theta - outputs.D_theta;
                outputs.L_phi = outputs.Fa_phi - outputs.D_phi;
            
                % Drag magnitude
                outputs.D = sqrt(outputs.D_r^2 + outputs.D_theta^2 + outputs.D_phi^2);
            
                % Lift magnitude
                outputs.L = sqrt(outputs.L_r^2 + outputs.L_theta^2 + outputs.L_phi^2);
            
                % Lift-to-drag ratio that follows from the chosen kinematic ratio
                outputs.E_result = sqrt(((outputs.Fa*outputs.va)/outputs.F_dot_v)^2-1);
            
                % LIMITS
                % Assume CT = 1 (or a = 0.5), different CT than above
                outputs.maxThrust = inputs.numTurbines*inputs.areaTurbine*0.5*inputs.densityAir*outputs.va^2;
                % Assume Cp = 0.593 (or a = 0.33)
                Cp_betz_limit = 0.593;
                outputs.maxPower = inputs.numTurbines*inputs.areaTurbine*Cp_betz_limit*0.5*inputs.densityAir*outputs.va^3;

                % Effective mechanical turbine power
                % outputs.P_m_eff = sqrt(outputs.Fthrust^3/2/inputs.densityAir/inputs.turbineArea);
                CTturbine = outputs.Fthrust/outputs.maxThrust;
                outputs.axialInduction = fzero(@(x)4*x*(1-x)-CTturbine,0.5);
                outputs.P_m_eff = 4*outputs.axialInduction*(1-outputs.axialInduction)^2*outputs.maxThrust*outputs.va;
            
                outputs.zetaMech = outputs.P_m_eff/(halfRhoS*outputs.vw^3);
            
                outputs.P_e_eff = outputs.P_m_eff*inputs.efficiency_PowerElectronics;
            
                % Time for one pattern revolution and number of patterns in the cycle
                outputs.tPatt = 2*pi*outputs.Rp./outputs.vk_phi;

                

                fval = -outputs.P_e_eff;
            
            end
            function [c, ceq] = constraints(x, inputs) %#ok<INUSD>
                % Inequality constraints
    
                % Min clearance between ground and bottom point of pattern
                c(1)   = (inputs.minGroundClear - outputs.pattGrClr);
    
                % Capping for requested electrical rated power
                c(2)   = (outputs.P_e_eff - inputs.P_ratedElec); 
    
                % Tether length limit
                c(3) = (outputs.l_t_inCycle - inputs.maxTeLen);

                % Maximum thrust
                c(4) = (outputs.Fthrust - outputs.maxThrust);
    
                % Max. cycle avg height
                c(5) = (outputs.h_inCycle - inputs.maxHeight);
    
                % Peak mechanical power limit
                c(6) = (outputs.P_m_eff - inputs.peakM2E_F*inputs.P_ratedElec);
    
                % Maximum tether force
                c(7) = (outputs.Ft - inputs.forceTether_max*inputs.safetyFactor_forceTetherMax);
    
                % Apparent speed cannot be negative
                c(8) = (0 - outputs.va);
                
                c(9) = (0 - outputs.Fthrust);

                c(10) = (inputs.minRp - outputs.Rp);

                c(11) = (outputs.P_m_eff - outputs.maxPower);
                
                %% Equality constraints
                % Glide ratio
                ceq(1) = (outputs.E_result - outputs.E);

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
            % processedOutputs.vw_reference = sort(vw,'ascend');
            % On exit flag convergence is too harsch, glide ratio is a good
            % indicator for convergence
            converged_idx = find(abs(mean((outputs.E_result - outputs.E),2)) < 0.01);
            vw_converged = vw(converged_idx);
            %% Cut-in wind speed
            % Glide ratio constraint violation acceptance for feasible solution
            % temp1                  = vw(abs(mean((outputs.E_result - outputs.E),2)) < 0.01); %0.01
            processedOutputs.cutIn = min(vw_converged(outputs.P_e_eff(converged_idx)>0)); %#ok<*PROPLC>
            
            %% Rated wind and power
            vw_above_rated              = vw(round(outputs.P_e_eff./max(outputs.P_e_eff),2)>=0.99);
            processedOutputs.ratedWind  = min(vw_above_rated); % At Reference height 
            processedOutputs.ratedPower = outputs.P_e_eff(vw==processedOutputs.ratedWind);
            
            %% Cut-out wind speed
            processedOutputs.cutOut   = inputs.vwCutOut_patternHeight; % maximum at pattern height
            % Operating range traslated to wind speed at ref. height
            processedOutputs.vw_100m_operRange = vw((vw>=processedOutputs.cutIn) & mean(outputs.vw,2)<=processedOutputs.cutOut);
            % processedOutputs.vw_100m_operRange = sort(processedOutputs.vw_100m_operRange,'ascend');
            processedOutputs.cutOut_referenceHeight = max(processedOutputs.vw_100m_operRange);
            %% Extract feasible results in the operational range
            processedOutputs.Dia_te = outputs.d_t(1); % always constant

            variableNamesToCopy = {'vw','P_m_eff','P_e_eff','Ft','Fa','W',...
                'vk_tau','Rp','l_t_inCycle','h_inCycle','va','beta',...
                'rollAngle','gamma','CL','CD','CD_k','CD_t','CT','E',...
                'lambda','tPatt','zetaMech','d_t','turbineGain',...
                'Fthrust','axialInduction'};

            variableNamesToInit_scalar = {'sweptArea', 'Cp_m', 'Cp_e', 'storageExchange'};

            for i=1:length(vw)
                %vw is reference height so cut_out at reference height should be taken
                if vw(i)>=processedOutputs.cutIn && vw(i)<=processedOutputs.cutOut_referenceHeight
                    % Define the list of variable names
                    processedOutputs = KiteQSMsimulationBase.copy_outputs_to_processedoutputs(...
                        processedOutputs, outputs, variableNamesToCopy, i);

                    %convert angles to degrees
                    processedOutputs.beta(i)           = rad2deg(processedOutputs.beta(i));
                    processedOutputs.rollAngle(i)    = rad2deg(-processedOutputs.rollAngle(i));
                    processedOutputs.gamma(i)          = rad2deg(processedOutputs.gamma(i));
                    
                    processedOutputs.storageExchange(i) = obj.kiteMass*inputs.accGravity*outputs.Rp(i)/3.6e3;  %[Wh]

                    % Coefficient of power (Cp) as defined for HAWTs
                    processedOutputs.sweptArea(i)     = pi.*((outputs.Rp(i)+inputs.span/2).^2 - (outputs.Rp(i)-inputs.span/2).^2);
                    processedOutputs.Cp_m(i)        = outputs.P_m_eff(i)./(0.5.*inputs.densityAir.*processedOutputs.sweptArea(i).*outputs.vw(i).^3);
                    processedOutputs.Cp_e(i)      = outputs.P_e_eff(i)./(0.5.*inputs.densityAir.*processedOutputs.sweptArea(i).*outputs.vw(i).^3);
                else
                    processedOutputs = obj.copy_zeros_to_processedoutputs_based_on_Size(...
                        processedOutputs, [1,1], variableNamesToInit_scalar, i);
                    processedOutputs = KiteQSMsimulationBase.copy_zeros_to_processedoutputs(...
                        processedOutputs, outputs, variableNamesToCopy, i);
                end
            end

            ind_vw_notzero = find(processedOutputs.vw>0);
            [~, sort_ind] = sort(processedOutputs.vw(ind_vw_notzero),'ascend');
            sorted_ind_vw_notzero = ind_vw_notzero(sort_ind);
            
            if ind_vw_notzero(1)~=sorted_ind_vw_notzero(1)
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
                'Cd0', 'Cl0_airfoil','numTurbines','areaTurbine'};
        end
        % 
        function parametersWithDefaults = get_parameters_withdefaults(inputs,validationFcn)
            parametersWithDefaults = {
                'aspectRatio', (inputs.span^2/inputs.areaWing), validationFcn.FloatPositive;
                'doIncludeGravity', true, validationFcn.Bool;
                'doMassOverride', false, validationFcn.Bool;
                'maxTeLen', 1000, validationFcn.FloatPositive;
                'minRp', 3*inputs.span, validationFcn.FloatPositive;
                'maxHeight', 1000, validationFcn.FloatPositive;
                'minGroundClear', 100, validationFcn.FloatPositive;
                'safetyFactor_forceTetherMax', 0.8, validationFcn.FloatPositive;
                'peakM2E_F', 2, validationFcn.FloatPositive;
                'maxStrengthTether', 7e8, validationFcn.FloatPositive;
                'densityTether', 980, validationFcn.FloatPositive;
                'e', 0.6, validationFcn.FloatPositive;
                'Cd_cylinder', 1.1, validationFcn.FloatPositive;
                'efficiency_PowerElectronics', 0.95, validationFcn.FloatPositive;
                'accGravity', 9.81, validationFcn.FloatPositive;
                'densityAir', 1.225, validationFcn.FloatPositive;
                'doUseReferenceWindProfile', false, validationFcn.Bool;
                'vwCutOut_patternHeight',25,validationFcn.FloatPositive;
                'Cd_additional',0,validationFcn.FloatPositiveOrZero;
                'name','simulationQSM_Flygen',validationFcn.String;
                'kiteMass',0,validationFcn.FloatPositiveOrZero;
                'tetherCoreIncreaseFactor',1.3,validationFcn.FloatPositive;
            };
        end
    end

end
