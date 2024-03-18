classdef KiteQSMsimulation < handle
    properties
        inputs
        outputs
        optimDetails
        processedOutputs
        kiteMass
    end
    
    methods
        function obj = KiteQSMsimulation(inputs)
            % Constructor
            p = inputParser;
            p.KeepUnmatched = true;
            validationFcnFloatPositive = @(x) assert(isnumeric(x) && isscalar(x) ...
                && (x > 0), 'Value must be positive, scalar, and numeric.');
            validationFcnInt = @(x) assert(isnumeric(x) && isscalar(x) ...
                && (x > 0) && mod(x,1)==0, 'Value must be positive, integer.');
            validationFcnBool = @(x) assert(islogical(x) || x==0 ...
                || x==1, 'Value must be boolean or 0/1.');
            validationFcnArrayFloatPositive = @(x) assert(isnumeric(x) ...
                && all(x > 0), 'Value must be positive, scalar, and numeric.');

            % List of required parameters
            requiredParams = {'windSpeedReference', 'heightWindReference', 'areaWing', 'span', 'forceTether_max', ...
                'P_ratedElec', 'Cl_maxAirfoil', 'Cl_eff_F', 'Cl0_airfoil', ...
                'Cd0'}; % Add all required parameters
            
            % Loop through required parameters and add them
            for i = 1:numel(requiredParams)
                param = requiredParams{i};
                if isfield(inputs, param)
                    addParameter(p, param, inputs.(param)); % Add parameter with value from input structure
                else
                    error(['KiteQSMsimulation:constructor Input parameter ', param, ' is required.']);
                end
            end

            if isfield(inputs,'doMassOverride')
                if inputs.doMassOverride
                    if isfield(inputs, 'kiteMass')
                        addParameter(p, 'kiteMass', 0, validationFcnFloatPositive);
                    else
                        error('KiteQSMsimulation:constructor Input parameter kiteMass is required when doMassOverride is true.');
                    end
                end
            end

            % Optional parameters with defaults if not provided
            addParameter(p, 'aspectRatio', (inputs.span^2/inputs.areaWing), validationFcnFloatPositive);
            addParameter(p, 'numDeltaLelems', 5, validationFcnInt); 
            addParameter(p, 'doIncludeGravity', true, validationFcnBool); % false = No, true = Yes
            addParameter(p, 'doMassOverride', false, validationFcnBool);
            addParameter(p, 'maxTeLen', 1000, validationFcnFloatPositive); %[m]
            addParameter(p, 'maxHeight', 1000, validationFcnFloatPositive); %[m]
            addParameter(p, 'minGroundClear', 100, validationFcnFloatPositive); %[m]
            addParameter(p, 'safetyFactor_forceTetherMax', 0.8, validationFcnFloatPositive); % [-] 0.8 for gust margin
            addParameter(p, 'peakM2E_F', 2.5, validationFcnFloatPositive); % [-]
            addParameter(p, 'maxStrengthTether', 7e8, validationFcnFloatPositive); % [Pa]
            addParameter(p, 'densityTether', 980, validationFcnFloatPositive); %[kg/m^3] 
            addParameter(p, 'speedReelout_max', 25, validationFcnFloatPositive); %[m/s]
            addParameter(p, 'accReel_max', 5, validationFcnFloatPositive); %[m/s^2]
            addParameter(p, 'e', 0.6, validationFcnFloatPositive); %[-] % oswald efficiency or lifting line span efficiency?
            addParameter(p, 'Cd_cylinder', 1.2, validationFcnFloatPositive); %[-]
            addParameter(p, 'etaGen_param', [0.671, -1.4141, 0.9747, 0.7233], ...
                @(x) assert(isnumeric(x) && numel(x) == 4, ...
                'Parameter must be a numeric array with 4 elements.')); %[-]
            addParameter(p, 'speedGeneratorReelout_max', 25, validationFcnFloatPositive); % Assuming speedReelout_max is the default
            addParameter(p, 'efficiency_Gearbox', 0.9, validationFcnFloatPositive); %[-]
            addParameter(p, 'efficiency_Storage', 0.9, validationFcnFloatPositive); %[-]
            addParameter(p, 'efficiency_PowerElectronics', 0.95, validationFcnFloatPositive); %[-]
            addParameter(p, 'accGravity', 9.81, validationFcnFloatPositive); %[m/s^2]
            addParameter(p, 'densityAir', 1.225, validationFcnFloatPositive); %[kg/m^3]  
            
            addParameter(p, 'doUseReferenceWindProfile', false, validationFcnBool); % false = Modelled, true = From dataset 
            if isfield(inputs,'doUseReferenceWindProfile')
                if inputs.doUseReferenceWindProfile
                    % Sample vertical wind profile datasets (i.e. for inputs.doUseReferenceWindProfile = true)
                    % Profiles from: https://doi.org/10.1016/j.renene.2022.06.094. Wind speeds normalized with value at 100m. 
                    h_windDataset          = [10,20,40,60,80,100,120,140,150,160,180,200,220,250,300,500,600]; % [m]
                    v_windDataset_Cabauw   = [0.541219682843206,0.607355091566827,0.768630154201962,0.868484406441142,0.941395360902529,1,1.04810058627160,1.08638854381156,1.10277338731106,1.11715868927737,1.14412258234309,1.16573551308321,1.18394938534465,1.20653423381438,1.23266397972046,1.26662287360302,1.26414483994687];
                    v_windDataset_Ijmuiden = [0.847612611633547,0.870603040595613,0.927240267828556,0.959346286990695,0.982291573490674,1,1.01377720773809,1.02356771954493,1.02766760602000,1.03079423355205,1.03659625208888,1.04025827758100,1.04284618416620,1.04496440015282,1.04461712713371,1.02473617783789,1.01076976884552]; %#ok<NASGU>
    
                    addParameter(p, 'windProfile_h', h_windDataset, validationFcnArrayFloatPositive);
                    addParameter(p, 'windProfile_vw', v_windDataset_Cabauw, validationFcnArrayFloatPositive);
                else
                    addParameter(p, 'windShearExp', 0.143, validationFcnFloatPositive); %[-] % 0.143 over land, 0.11 over sea
                end
            else
                addParameter(p, 'windShearExp', 0.143, validationFcnFloatPositive); %[-] % 0.143 over land, 0.11 over sea
            end

            parse(p, inputs);
            
            
            p2 = inputParser;
            p2.KeepUnmatched = true;
            % Optimisation problem data
            addParameter(p2, 'nx', ones(1,p.Results.numDeltaLelems), validationFcnArrayFloatPositive);
            parse(p2, inputs);

            p3 = inputParser;
            p3.KeepUnmatched = true;
            addParameter(p3, 'x0', [200, ... % deltaL
                deg2rad(30), ... % avgPattEle
                deg2rad(5), ... % coneAngle
                50, ... % Rp_start
                p.Results.speedReelout_max*p2.Results.nx, ... % v_i
                p.Results.Cl_maxAirfoil*p.Results.Cl_eff_F*p2.Results.nx, ... % CL_i
                0.8*p2.Results.nx, ... % v_o
                90*p2.Results.nx, ... % kinematicRatio
                p.Results.Cl_maxAirfoil*p.Results.Cl_eff_F*p2.Results.nx], ... % CL
                validationFcnArrayFloatPositive);
            addParameter(p3, 'lb', [50, ...
                deg2rad(1), ...
                deg2rad(1), ...
                50, ...
                1*p2.Results.nx, ...
                0.1*p2.Results.nx, ...
                0.8*p2.Results.nx, ...
                1*p2.Results.nx, ...
                0.1*p2.Results.nx], validationFcnArrayFloatPositive);
            addParameter(p3, 'ub', [500, ...
                deg2rad(90), ...
                deg2rad(60), ...
                100, ...
                p.Results.speedReelout_max*p2.Results.nx, ...
                p.Results.Cl_maxAirfoil*p.Results.Cl_eff_F*p2.Results.nx,...
                p.Results.speedReelout_max*p2.Results.nx, ...
                200*p2.Results.nx, ...
                p.Results.Cl_maxAirfoil*p.Results.Cl_eff_F*p2.Results.nx], ...
                validationFcnArrayFloatPositive);
            parse(p3, inputs);

            % write inputs to class property
            f = fieldnames(p.Results);
            for i = 1:length(f)
                obj.inputs = setfield(obj.inputs, f{i}, p.Results.(f{i}));
            end
            f = fieldnames(p2.Results);
            for i = 1:length(f)
                obj.inputs = setfield(obj.inputs, f{i}, p2.Results.(f{i}));
            end
            f = fieldnames(p3.Results);
            for i = 1:length(f)
                obj.inputs = setfield(obj.inputs, f{i}, p3.Results.(f{i}));
            end

        end
        
        function obj = runsimulation(obj, inputs, optionsFMINCON)
            persistent outputs;
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

            % Bounds
            lb = inputs.lb;
            ub = inputs.ub;

            % Initialize variables to store results
            exitflag = zeros(1, length(inputs.windSpeedReference));
            optHist = cell(1, length(inputs.windSpeedReference));
            lambda = cell(1, length(inputs.windSpeedReference));

            % Preinitialize the struct with all fields set to empty arrays
            varNames = {'deltaL', 'beta', 'gamma', 'Rp_start', 'vk_r_i', ...
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
            for i = 1:length(varNames)
                outputSave.(varNames{i})(length(inputs.windSpeedReference),1) = 0;
            end
            outputs = [];

            % obj.initialise_outputs_to_zero();
            halfRhoS = 0.5*inputs.densityAir*inputs.areaWing;

            % Kite mass 
            if inputs.doMassOverride == 1
                obj.kiteMass = inputs.kiteMass;
            else
                % Vincent Bonnin's simple mass model developed at Ampyx Power. Based on AP3 data and projected data for larger systems (AP4-AP5)      
                a1     = 0.002415;       a2     = 0.0090239;       b1     = 0.17025;       b2     = 3.2493;
                k1     = 5;              c1     = 0.46608;         d1     = 0.65962;       k2     = 1.1935;
                AR_ref = 12;
                a_massCalc = a1*(inputs.forceTether_max/1000/inputs.areaWing) + a2;
                b_massCalc = b1*(inputs.forceTether_max/1000/inputs.areaWing) + b2;
                obj.kiteMass = 10*(a_massCalc*inputs.areaWing^2 +b_massCalc*inputs.areaWing-k1)...
                    *(c1*(inputs.aspectRatio/AR_ref)^2-d1*(inputs.aspectRatio/AR_ref)+k2); 
            end

            % Optimize operation for every wind speed
            for i = 1:length(inputs.windSpeedReference)
                windSpeed = inputs.windSpeedReference(i); 
                % Define objective and constraint functions
                conFunc = @(x) constraints(x, inputs);
                objFunc = @(x) objective(x, inputs, obj.kiteMass, halfRhoS, windSpeed);

                % Run optimization

                [x, ~, exitflag(i), optHist{i}, lambda{i}] = fmincon(objFunc, x0, [], [], [], [], lb, ub, conFunc, optionsFMINCON);

                % Store final results
                [~] = objective(x, inputs, obj.kiteMass, halfRhoS, windSpeed);
                for var_i = 1:length(varNames)
                    tempOut = outputs.(varNames{var_i});
                    outputSave.(varNames{var_i})(i,1:numel(tempOut)) = tempOut;
                end

                % Update initial guess for the next iteration
                if abs(mean((outputs.E_result - outputs.E), 2)) > 1 %0.01
                    x0 = inputs.x0;
                else
                    x0 = x;
                end
            end

            % Store optimization details
            obj.outputs = outputSave;
            obj.optimDetails.optHist = optHist;
            obj.optimDetails.exitflag = exitflag;
            obj.optimDetails.lambda = lambda;
            
            function fval = objective(x, inputs, kiteMass, halfRhoS, windSpeedReference)
                numDeltaLelems = inputs.numDeltaLelems;
            
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
                    outputs.phi(1,j) = 0;
                    outputs.chi(1,j) = deg2rad(90);
                
                    % Effective CD
                    outputs.CD_k(1,j) = inputs.Cd0 + (outputs.CL(1,j)-inputs.Cl0_airfoil)^2/(pi*inputs.aspectRatio*inputs.e);
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
                    vkr_ratio = outputs.vk_r(1,j)/inputs.speedGeneratorReelout_max;
                    outputs.etaGen_o(1,j) = inputs.etaGen_param(1)*vkr_ratio^3 + ...
                        inputs.etaGen_param(2)*vkr_ratio^2 + ...
                        inputs.etaGen_param(3)*vkr_ratio+inputs.etaGen_param(4);
                    outputs.P_e_o_eff(1,j) = outputs.P_m_o_eff(1,j)*inputs.efficiency_Gearbox*outputs.etaGen_o(1,j)*inputs.efficiency_PowerElectronics;
                
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
                    outputs.CD_k_i(1,j) = inputs.Cd0+(outputs.CL_i(1,j)- inputs.Cl0_airfoil)^2/(pi*inputs.aspectRatio*inputs.e);
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
                    vkri_ratio = outputs.vk_r_i(1,j)/inputs.speedGeneratorReelout_max;
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
                outputs.t2       = outputs.vk_r_i/inputs.accReel_max;
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
        
        function [inputs, outputs, optimDetails, processedOutputs] = getresult(obj)
            % Post processing and saving outputs
            % Define processedOutputs and other processing tasks here
            % ...

            % Return results
            inputs = obj.inputs;
            outputs = obj.outputs;
            optimDetails = obj.optimDetails;
            processedOutputs = obj.processedOutputs;
        end
        
        function obj = plotresults(obj, inputs, processedOutputs)
            if nargin<2
                inputs = obj.inputs;
                processedOutputs = obj.processedOutputs;
            elseif nargin<3
                processedOutputs = obj.processedOutputs;
            end
            %% Plot settings
            vw = inputs.windSpeedReference;
            x_axis_limit = [0 processedOutputs.vw_100m_operRange(end)];
            
            newcolors = [ % 0.25, 0.25, 0.25
                0 0.4470 0.7410
                0.8500 0.3250 0.0980
                0.4660, 0.6740, 0.1880
                0.9290, 0.6940, 0.1250
                0.4940, 0.1840, 0.5560
                0.6350 0.0780 0.1840
                0.3010 0.7450 0.9330];

            % % Plots

            % Plots showing the per element results for rated wind speed
            ws = processedOutputs.ratedWind;
            fig = figure();
            colororder(newcolors)
            % Lengths
            subplot(2,2,1)
            hold on; grid on; box on
            plot(processedOutputs.Rp(ws,:),':o','linewidth',1,'markersize',3);
            plot(processedOutputs.l_t_inCycle(ws,:),':^','linewidth',1,'markersize',3);
            plot(processedOutputs.h_inCycle(ws,:),':s','linewidth',1,'markersize',3);
            ylabel('(m)');
            legend('R_{p}','l_{t}','h_{p}','location','northwest');
            hold off
            % Speeds
            subplot(2,2,2); hold on; grid on; box on; 
            yyaxis left
            plot(processedOutputs.lambda(ws,:),':o','linewidth',1,'markersize',3);
            plot(processedOutputs.f(ws,:),':s','linewidth',1,'markersize',3);
            plot(processedOutputs.f_i(ws,:),':v','linewidth',1,'markersize',3);
            ylabel('(-)');
            yyaxis right
            plot(processedOutputs.vw(ws,:),':s','linewidth',1,'markersize',3);
            plot(processedOutputs.vw_i(ws,:),':v','linewidth',1,'markersize',3);
            ylabel('(m/s)');
            legend('λ','f_o','f_i','v_{w,o}','v_{w,i}','location','northwest');
            hold off
            % Glide ratios
            subplot(2,2,3); hold on; grid on; box on
            plot(processedOutputs.E(ws,:),':s','linewidth',1,'markersize',3);
            plot(processedOutputs.E_i(ws,:),':v','linewidth',1,'markersize',3);
            % ylim([0 2.5]);
            ylabel('(-)');
            legend('E_o','E_i','location','northwest');
            hold off
            % Forces
            subplot(2,2,4); hold on; grid on; box on
            plot(processedOutputs.Fa(ws,:)/1e3,':^','linewidth',1,'markersize',3);
            plot(processedOutputs.Fa_i(ws,:)/1e3,':o','linewidth',1,'markersize',3);
            plot(processedOutputs.Ft(ws,:)/1e3,':s','linewidth',1,'markersize',3);
            plot(processedOutputs.Ft_i(ws,:)/1e3,':v','linewidth',1,'markersize',3);
            plot(processedOutputs.W(ws,:)/1e3,':x','linewidth',1,'markersize',3);
            legend('F_{a,o}','F_{a,i}','F_{t,o}','F_{t,i}','W','location','northwest');
            ylabel('(kN)');
            hold off
            han=axes(fig,'visible','off'); 
            han.Title.Visible='on'; han.XLabel.Visible='on'; han.YLabel.Visible='off';
            ylabel(han,'yourYLabel');
            xlabel(han,'Discretized reeling distance over cycle');
            title(han,['Wind speed at 100m = ' num2str(ws) 'm/s']);
            
            % Speeds
            figure('units','inch','Position', [0.2 0.5 3.5 2.2])
            colororder(newcolors)
            hold on
            grid on
            box on
            plot(vw, mean(processedOutputs.lambda,2),':o','linewidth',1,'markersize',3);
            plot(vw, mean(processedOutputs.f,2),':s','linewidth',1,'markersize',3);
            plot(vw, mean(processedOutputs.f_i,2),':v','linewidth',1,'markersize',3);
            ylabel('(-)');
            legend('λ','f_{o}','f_{i}','location','northwest');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            %ylim([0 160]);
            hold off
        
            % CL
            figure('units','inch','Position', [4.2 0.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            plot(vw, mean(processedOutputs.CL,2),'o:','linewidth',1,'markersize',3);
            plot(vw, mean(processedOutputs.CL_i,2),'d:','linewidth',1,'markersize',3);
            ylabel('(-)');
            legend('C_{L,o}','C_{L,i}','location','northwest','Orientation','vertical');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            %ylim([0 160]);
            hold off
        
            % CD
            figure('units','inch','Position', [8.2 0.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            plot(vw, mean(processedOutputs.CD,2),'o:','linewidth',1,'markersize',3);
            plot(vw, mean(processedOutputs.CD_i,2),'d:','linewidth',1,'markersize',3);
            plot(vw, mean(processedOutputs.CD_k,2),'^:','linewidth',1,'markersize',3);
            plot(vw, mean(processedOutputs.CD_k_i,2),'v:','linewidth',1,'markersize',3);
            plot(vw, mean(processedOutputs.CD_t,2),'s:','linewidth',1,'markersize',3);
            ylabel('(-)');
            legend('C_{D,o}','C_{D,i}','C_{D,k,o}','C_{D,k,i}','C_{D,t}','location','northwest','Orientation','vertical');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            %ylim([0 160]);
            hold off
            
            % Lengths
            figure('units','inch','Position', [12.2 0.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            plot(vw, processedOutputs.h_cycleAvg,'d:','linewidth',1,'markersize',3);
            plot(vw, mean(processedOutputs.Rp,2),'o:','linewidth',1,'markersize',3);
            plot(vw, processedOutputs.deltaL,'^:','linewidth',1,'markersize',3);
            plot(vw, processedOutputs.l_t_max,'s:','linewidth',1,'markersize',3);
            plot(vw, processedOutputs.l_t_min,'x:','linewidth',1,'markersize',3);
            ylabel('(m)');
            legend('h_{p,avg}','R_{p,avg}','Δl','l_{t,max}','l_{t,min}','location','northwest','Orientation','vertical');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            %ylim([0 160]);
            hold off
        
            % Roll angle, avg patt elevation
            figure('units','inch','Position', [16.2 0.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            plot(vw, mean(processedOutputs.rollAngle,2),'s:','linewidth',1,'markersize',3);
            plot(vw, processedOutputs.beta,'o:','linewidth',1,'markersize',3);
            plot(vw, processedOutputs.gamma,'d:','linewidth',1,'markersize',3);
            ylabel('(deg)');
            legend('$$\Psi$$','$$\overline{\beta}$$','$$\gamma$$','location','northwest','Orientation','vertical','Interpreter','latex');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            %ylim([0 160]);
            hold off
        
            % Reel-out forces
            figure('units','inch','Position', [0.2 3.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            plot(vw, mean(processedOutputs.W,2)./10^3,':o','markersize',3)
            plot(vw, mean(processedOutputs.Fa,2)./10^3,':s','linewidth',1,'markersize',3);
            plot(vw, mean(processedOutputs.Ft,2)./10^3,':^','linewidth',1,'markersize',3);
            ylabel('(kN)');
            legend('W','F_{a,o}','F_{t,o}','location','northwest');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            %ylim([0 160]);
            hold off
        
            % Reel-in forces
            figure('units','inch','Position', [4.2 3.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            plot(vw, mean(processedOutputs.W,2)./10^3,':o','markersize',3)
            plot(vw, mean(processedOutputs.Fa_i,2)./10^3,':s','linewidth',1,'markersize',3);
            plot(vw, mean(processedOutputs.Ft_i,2)./10^3,':^','linewidth',1,'markersize',3);
            ylabel('(kN)');
            legend('W','F_{a,i}','F_{t,i}','location','northwest');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            hold off
        
            % Time, num of patterns
            figure('units','inch','Position', [8.2 3.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            yyaxis left
            plot(vw, processedOutputs.to,':s','linewidth',1,'markersize',3);
            plot(vw, processedOutputs.ti,':d','linewidth',1,'markersize',3);
            plot(vw, mean(processedOutputs.tPatt,2),':o','linewidth',1,'markersize',3);
            ylabel('(s)');
            yyaxis right
            plot(vw, mean(processedOutputs.numOfPatt,2),':^','linewidth',1,'markersize',3);
            ylabel('(-)');
            xlabel('Wind speed at 100m height (m/s)');
            legend('t_{o}','t_{i}','t_{patt,avg}','N_{p}','location','northwest');
            xlim(x_axis_limit);
            hold off
        
            % Reel-out power
            figure('units','inch','Position', [12.2 3.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            plot(vw, processedOutputs.P_m_o./10^3,':o','linewidth',1,'markersize',3);
            plot(vw, processedOutputs.P_e_o./10^3,':s','linewidth',1,'markersize',3);
            plot(vw, processedOutputs.P_e_avg./10^3,':x','linewidth',1,'markersize',5);
            ylabel('(kW)');
            %title('Cycle averages');
            legend('P_{m,o}','P_{e,o}','P_{e,avg}','location','northwest');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            hold off
        
            % Reel-in power
            figure('units','inch','Position', [16.2 3.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            plot(vw, processedOutputs.P_m_i./10^3,':o','linewidth',1,'markersize',3);
            plot(vw, processedOutputs.P_e_i./10^3,':s','linewidth',1,'markersize',3);
            ylabel('(kW)');
            %title('Cycle averages');
            legend('P_{m,i}','P_{e,i}','location','northwest');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            hold off

            % Power curve comparison plot
            % Loyd
            P_Loyd = zeros(1,length(vw));
            CL_loyd = inputs.Cl_maxAirfoil * inputs.Cl_eff_F;
            CD_loyd = inputs.Cd0 + (CL_loyd - inputs.Cl0_airfoil)^2 / (pi * inputs.aspectRatio * inputs.e);
            for i = 1:length(vw)
                P_Loyd(i) = (4/27) * (CL_loyd^3 / CD_loyd^2) * (1/2) * inputs.densityAir * inputs.areaWing * (vw(i)^3);
                if P_Loyd(i) > inputs.P_ratedElec
                    P_Loyd(i) = inputs.P_ratedElec;
                end 
            end
            
            % AP3 6DoF simulation results
            AP3_PC_ws = [7.32E+00, 8.02E+00, 9.05E+00, 1.00E+01, 1.10E+01, 1.20E+01, ...
              1.30E+01, 1.41E+01, 1.50E+01, 1.60E+01, 1.70E+01, 1.80E+01, 1.90E+01]; %[m/s]
            AP3_PC_power = [0.00E+00, 7.01E+03, 2.37E+04, 4.31E+04, 6.47E+04, 8.46E+04, ...
              1.02E+05, 1.20E+05, 1.34E+05, 1.49E+05, 1.50E+05, 1.50E+05, 1.50E+05]./10^3; %[kW]
            
            figure('units','inch','Position', [0.2 6.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            plot(vw, P_Loyd./10^3,'-','linewidth',1);
            plot(vw, processedOutputs.P_e_avg./10^3,':s','linewidth',1,'MarkerSize',3);
            plot(AP3_PC_ws, AP3_PC_power,':^k','linewidth',1,'MarkerSize',3);
            ylabel('Power (kW)');
            legend('Ideal Reel-out','QSM','6-DOF','location','southeast');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            hold off

            % Cycle timeseries plots: Pattern averages
            windSpeeds = [processedOutputs.cutIn, processedOutputs.ratedWind, processedOutputs.cutOut];
            for i = windSpeeds
                idx = find(vw == i);
                figure('units','inch','Position', [4.2 6.5 3.5 2.2])
                hold on
                grid on
                box on
                % Plot your data
                plot(processedOutputs.cyclePowerRep(i).t_inst, processedOutputs.cyclePowerRep(i).P_e_inst,'linewidth',1.5, 'DisplayName', 'P_{e}');
                plot(processedOutputs.cyclePowerRep(i).t_inst, processedOutputs.cyclePowerRep(i).P_m_inst,'--','linewidth',1.5, 'DisplayName', 'P_{m}');
                yline(processedOutputs.P_e_avg(idx)/10^3, ':', 'linewidth',1.2, 'DisplayName', 'P_{e,avg}');
                % Set yline positions within the range of your data
                yline(0, '-', 'linewidth', 1, 'HandleVisibility', 'off');      
                ylabel('(kW)');
                xlabel('Time (s)');      
                % Create a custom legend with 'Data1' and 'Data2' only
                legend('Location', 'Best');      
                % Customize legend colors
                legend('TextColor', 'k', 'Color', 'w'); % 'TextColor' sets the text color, 'Color' sets the background color      
                xlim([0 processedOutputs.cyclePowerRep(i).t_inst(end)]);
                ylim([1.5*min(processedOutputs.cyclePowerRep(i).P_e_inst) 1.05*max(processedOutputs.cyclePowerRep(i).P_m_inst)]);
                title(strcat('Wind speed at 100m:', num2str(processedOutputs.cyclePowerRep(idx).ws), ' m/s'));      
                hold off
            end
            
            % Wind profile
            if ~inputs.doUseReferenceWindProfile
                Vref = 10; % m/s
                z = 10:10:600; % m
                V = Vref * (z/inputs.heightWindReference).^inputs.windShearExp/Vref;
                %   MegAWES Onshore location Cabauw. Wind speeds normalized with value at 100m. Profiles from: https://doi.org/10.1016/j.renene.2022.06.094
                z_MegAWES = [10,20,40,60,80,100,120,140,150,160,180,200,220,250,300,500,600];
                V_MegAWES_Cabauw = [0.541219682843206,0.607355091566827,0.768630154201962,0.868484406441142,0.941395360902529,1,1.04810058627160,1.08638854381156,1.10277338731106,1.11715868927737,1.14412258234309,1.16573551308321,1.18394938534465,1.20653423381438,1.23266397972046,1.26662287360302,1.26414483994687];
                V_MegAWES_Ijmuiden = [0.847612611633547,0.870603040595613,0.927240267828556,0.959346286990695,0.982291573490674,1,1.01377720773809,1.02356771954493,1.02766760602000,1.03079423355205,1.03659625208888,1.04025827758100,1.04284618416620,1.04496440015282,1.04461712713371,1.02473617783789,1.01076976884552];
                figure('units','inch','Position', [8.2 6.5 2 2.2])
                hold on
                box on
                grid on
                plot(V,z,'linewidth',1.2)
                plot(V_MegAWES_Cabauw,z_MegAWES,'--','linewidth',1.2)
                plot(V_MegAWES_Ijmuiden,z_MegAWES,'-.','linewidth',1.2)
                legend('α = 0.143','Cabauw,NL','Ijmuiden,NL');
                xlim([0.5 1.5])
                xlabel('Wind Speed (-)')
                ylabel('Height (m)')
                % saveas(gcf, fullfile(directory, ['vert_wind_prof', '.fig']));
                % print(gcf, fullfile(directory, ['vert_wind_prof', '.eps']), '-depsc2');
                hold off
            else
                % Continuous function from vertical wind profile
                % Define the heights at which you want to interpolate
                heights = (0:10:1000); % Adjust the desired height values as needed
                % Interpolate at the specified heights
                V_interpolated = interp1(inputs.windProfile_h, inputs.windProfile_vw, heights, 'linear', 'extrap');
                % Plot the interpolated function
                figure('units','inch','Position', [15 6 2 2.2])
                hold on 
                box on
                grid on
                plot(inputs.windProfile_vw, inputs.windProfile_h, 'o', V_interpolated, heights, '-');
                xlabel('Wind speed (m/s)');
                ylabel('Height (m)');
                xlim([0.5 1.5])
                title('Wind profile');
                legend('Data', 'Interpolated', 'Location', 'Best');
                hold off
            end


        end
                
        function obj = processoutputs(obj, inputs, outputs)
            if nargin<2
                inputs = obj.inputs;
                outputs = obj.outputs;
            elseif nargin<3
                outputs = obj.outputs;
            end
            %% Post processing
            vw = inputs.windSpeedReference; % Wind speed at ref. height
            
            %% Cut-in wind speed
            % Glide ratio constraint violation acceptance for feasible solution
            temp1                  = vw(abs(mean((outputs.E_result - outputs.E),2)) < 0.01); %0.01
            processedOutputs.cutIn = max(temp1(1)); %#ok<*PROPLC>
            
            %% Rated wind and power
            temp3                       = round(outputs.P_e_avg./max(outputs.P_e_avg),2);
            temp4                       = vw(temp3>=0.99);
            processedOutputs.ratedWind  = temp4(1); % At Reference height 
            processedOutputs.ratedPower = outputs.P_e_avg(vw==temp4(1));
            
            %% Cut-out wind speed
            processedOutputs.cutOut   = 25; % At operational height (Assumption)
            % Operating range traslated to wind speed at ref. height
            processedOutputs.vw_100m_operRange = vw(mean(outputs.vw,2)<=processedOutputs.cutOut);
            
            %% Extract feasible results in the operational range
            processedOutputs.Dia_te = outputs.d_t;
            for i=1:length(vw)
                if vw(i)>=processedOutputs.cutIn && vw(i)<= processedOutputs.cutOut
                    processedOutputs.vw(i,:)           = outputs.vw(i,:);
                    processedOutputs.vw_i(i,:)         = outputs.vw_i(i,:);
                    processedOutputs.P1_m_o(i)         = outputs.P1_m_o(i);
                    processedOutputs.P2_m_i(i)         = outputs.P2_m_i(i); 
                    processedOutputs.P_m_o_eff(i,:)    = outputs.P_m_o_eff(i,:);
                    processedOutputs.P_m_i_eff(i,:)    = outputs.P_m_i_eff(i,:);  
                    processedOutputs.P1_e_o(i)         = outputs.P1_e_o(i);
                    processedOutputs.P2_e_i(i)         = outputs.P2_e_i(i); 
                    processedOutputs.P_e_o_eff(i,:)    = outputs.P_e_o_eff(i,:);
                    processedOutputs.P_e_i_eff(i,:)    = outputs.P_e_i_eff(i,:);
                    processedOutputs.P_m_o(i)          = outputs.P_m_o(i);
                    processedOutputs.P_m_i(i)          = outputs.P_m_i(i);
                    processedOutputs.P_e_o(i)          = outputs.P_e_o(i);
                    processedOutputs.P_e_i(i)          = outputs.P_e_i(i);
                    processedOutputs.deltaL(i)         = outputs.deltaL(i);
                    processedOutputs.t1(i)             = outputs.t1(i);
                    processedOutputs.t2(i)             = outputs.t2(i);
                    processedOutputs.to_eff(i,:)       = outputs.to_eff(i,:);
                    processedOutputs.ti_eff(i,:)       = outputs.ti_eff(i,:);
                    processedOutputs.to(i)             = outputs.to(i);
                    processedOutputs.ti(i)             = outputs.ti(i);
                    processedOutputs.tCycle(i)         = outputs.tCycle(i);
                    processedOutputs.Ft(i,:)           = outputs.Ft(i,:);
                    processedOutputs.Ft_i(i,:)         = outputs.Ft_i(i,:);
                    processedOutputs.Fa(i,:)           = outputs.Fa(i,:);
                    processedOutputs.Fa_i(i,:)         = outputs.Fa_i(i,:);
                    processedOutputs.W(i,:)            = outputs.W(i,:);
                    processedOutputs.vo(i,:)           = outputs.vk_r(i,:);
                    processedOutputs.vk_tau(i,:)       = outputs.vk_tau(i,:);
                    processedOutputs.vi(i,:)           = outputs.vk_r_i(i,:);
                    processedOutputs.P_e_avg(i)        = outputs.P_e_avg(i);
                    processedOutputs.P_m_avg(i)        = outputs.P_m_avg(i);
                    processedOutputs.Rp(i,:)           = outputs.Rp(i,:);
                    processedOutputs.h_cycleStart(i)   = outputs.h_cycleStart(i);
                    processedOutputs.h_cycleAvg(i)     = outputs.h_cycleAvg(i);
                    processedOutputs.l_t_max(i)        = outputs.l_t_max(i);
                    processedOutputs.l_t_min(i)        = outputs.l_t_min(i);
                    processedOutputs.l_t_inCycle(i,:)  = outputs.l_t_inCycle(i,:);
                    processedOutputs.h_inCycle(i,:)    = outputs.h_inCycle(i,:); 
                    processedOutputs.va(i,:)           = outputs.va(i,:);
                    processedOutputs.beta(i)           = rad2deg(outputs.beta(i));
                    processedOutputs.rollAngle(i,:)    = rad2deg(-outputs.rollAngle(i,:));
                    processedOutputs.gamma(i)          = rad2deg(outputs.gamma(i));
                    processedOutputs.CL(i,:)           = outputs.CL(i,:);
                    processedOutputs.CD(i,:)           = outputs.CD(i,:);
                    processedOutputs.CD_k(i,:)         = outputs.CD_k(i,:);
                    processedOutputs.CD_k_i(i,:)       = outputs.CD_k_i(i,:);
                    processedOutputs.CD_t(i,:)         = outputs.CD_t(i,:);
                    processedOutputs.CL_i(i,:)         = outputs.CL_i(i,:);
                    processedOutputs.CD_i(i,:)         = outputs.CD_i(i,:);
                    processedOutputs.E(i,:)            = outputs.E(i,:);
                    processedOutputs.E_i(i,:)          = outputs.E_i(i,:);
                    processedOutputs.numOfPatt(i,:)    = outputs.numOfPatt(i,:);
                    processedOutputs.lambda(i,:)       = outputs.lambda(i,:); 
                    processedOutputs.f(i,:)            = outputs.f(i,:); 
                    processedOutputs.f_i(i,:)          = outputs.f_i(i,:); 
                    processedOutputs.tPatt(i,:)        = outputs.tPatt(i,:);
                    processedOutputs.dutyCycle(i)      = processedOutputs.to(i)/processedOutputs.tCycle(i);
                    processedOutputs.zetaMech(i,:)     = outputs.zetaMech(i,:);  
                    processedOutputs.d_te              = outputs.d_t;  
                    
                    % Cycle efficiency
                    processedOutputs.cycleEff_elec(i)  = (processedOutputs.P_e_avg(i)*processedOutputs.tCycle(i))/...
                                                            (processedOutputs.P_e_o(i)*processedOutputs.to(i));
                    processedOutputs.cycleEff_mech(i)  = (processedOutputs.P_m_avg(i)*processedOutputs.tCycle(i))/...
                                                            (processedOutputs.P_m_o(i)*processedOutputs.to(i));
                    
                    % Coefficient of power (Cp) as defined for HAWTs
                    processedOutputs.sweptArea(i,:)     = pi.*((outputs.Rp(i,:)+inputs.span/2).^2 - (outputs.Rp(i,:)-inputs.span/2).^2);
                    processedOutputs.Cp_m_o(i,:)        = outputs.P_m_o(i)./(0.5.*inputs.densityAir.*processedOutputs.sweptArea(i,:).*outputs.vw(i,:).^3);
                    processedOutputs.Cp_e_avg(i,:)      = outputs.P_e_avg(i)./(0.5.*inputs.densityAir.*processedOutputs.sweptArea(i,:).*outputs.vw(i,:).^3);
                end
            end
            
            %% Cycle power representation for wind speeds in the operational range
            for vw_i = vw(vw>=processedOutputs.cutIn&vw<=processedOutputs.cutOut)
                processedOutputs.cyclePowerRep(vw_i) = createCyclePowerRep(vw_i,...
                    inputs.windSpeedReference, processedOutputs.t1, ...
                    processedOutputs.to_eff, processedOutputs.t2,...
                    processedOutputs.ti_eff, processedOutputs.P_e_o_eff,...
                    processedOutputs.P_e_i_eff, processedOutputs.P_m_o_eff, ...
                    processedOutputs.P_m_i_eff); 
            end

            obj.processedOutputs = processedOutputs;

            function cyclePowerRep = createCyclePowerRep(vw_i, ...
                    windSpeedReference, t1, to_eff, t2, ti_eff, P_e_o_eff, ...
                    P_e_i_eff, P_m_o_eff, P_m_i_eff)
                cyclePowerRep.ws     = vw_i;
                cyclePowerRep.idx    = find(windSpeedReference==vw_i);
                cyclePowerRep.t1     = round(t1(cyclePowerRep.idx),2);
                cyclePowerRep.to_eff = round(sum(to_eff(cyclePowerRep.idx,:),2),2);
                cyclePowerRep.t2     = round(t2(cyclePowerRep.idx),2);
                cyclePowerRep.ti_eff = round(sum(ti_eff(cyclePowerRep.idx,:),2),2);
                cyclePowerRep.to     = cyclePowerRep.t1 + cyclePowerRep.to_eff; %[s]
                cyclePowerRep.ti     = cyclePowerRep.t2 + cyclePowerRep.ti_eff; %[s]
                cyclePowerRep.tCycle = cyclePowerRep.to + cyclePowerRep.to; %[s]
                
                cyclePowerRep.t_inst   = cumsum([0 cyclePowerRep.t1 (cyclePowerRep.to_eff-cyclePowerRep.t1) cyclePowerRep.t1 ...
                                      cyclePowerRep.t2 (cyclePowerRep.ti_eff-cyclePowerRep.t2) cyclePowerRep.t2]);
                
                cyclePowerRep.P_e_inst = [0, P_e_o_eff(cyclePowerRep.idx,1), P_e_o_eff(cyclePowerRep.idx,end), 0, ...
                                      -P_e_i_eff(cyclePowerRep.idx,end), -P_e_i_eff(cyclePowerRep.idx,1), 0]./10^3;
                
                cyclePowerRep.P_m_inst = [0 P_m_o_eff(cyclePowerRep.idx,1), P_m_o_eff(cyclePowerRep.idx,end), 0, ...
                                      -P_m_i_eff(cyclePowerRep.idx,end), -P_m_i_eff(cyclePowerRep.idx,1), 0]./10^3;
            end
        end

        

        function obj = initialise_outputs_to_zero(obj)
            % Initialize all obj.outputs fields with zeros
            num_ref = length(obj.inputs.windSpeedReference);
            % num_elems = obj.inputs.numDeltaLelems;

            obj.kiteMass = 0;
            
            % Single-column variables
            outputs(num_ref).l_t_min = []; %#ok<*PROP>
            outputs(num_ref).pattStartGrClr = [];
            outputs(num_ref).h_cycleStart = [];
            outputs(num_ref).l_t_max = [];
            outputs(num_ref).pattEndGrClr = [];
            outputs(num_ref).l_t_avg = [];
            outputs(num_ref).h_cycleAvg = [];
            outputs(num_ref).h_cycleEnd = [];
            outputs(num_ref).d_t = [];
            outputs(num_ref).deltaLelems = [];
            outputs(num_ref).elemDeltaL = [];
            outputs(num_ref).t1 = [];
            outputs(num_ref).P1_m_o = [];
            outputs(num_ref).P1_e_o = [];
            outputs(num_ref).P_m_o = [];
            outputs(num_ref).P_e_o = [];
            outputs(num_ref).t2 = [];
            outputs(num_ref).P2_m_i = [];
            outputs(num_ref).P2_e_i = [];
            outputs(num_ref).P_m_i = [];
            outputs(num_ref).P_e_i = [];
            outputs(num_ref).tCycle = [];
            outputs(num_ref).tPatt = [];
            outputs(num_ref).numOfPatt = [];
            outputs(num_ref).P_e_avg = [];
            outputs(num_ref).P_m_avg = [];
            outputs(num_ref).deltaL = [];
            outputs(num_ref).beta = [];
            outputs(num_ref).gamma = [];
            outputs(num_ref).Rp_start = [];
            outputs(num_ref).etaGen_i = [];
            outputs(num_ref).to = [];
            outputs(num_ref).ti = [];
            
            % Variables with multiple columns
            outputs(num_ref).l_t_inCycle = [];
            outputs(num_ref).pattGrClr = [];
            outputs(num_ref).m_t = [];
            outputs(num_ref).m_eff = [];
            outputs(num_ref).theta = [];
            outputs(num_ref).phi = [];
            outputs(num_ref).chi = [];
            outputs(num_ref).CD_k = [];
            outputs(num_ref).CD_t = [];
            outputs(num_ref).CD = [];
            outputs(num_ref).E = [];
            outputs(num_ref).h_inCycle = [];
            outputs(num_ref).Rp = [];
            outputs(num_ref).vw = [];
            outputs(num_ref).vw_r = [];
            outputs(num_ref).vw_theta = [];
            outputs(num_ref).vw_phi = [];
            outputs(num_ref).f = [];
            outputs(num_ref).va = [];
            outputs(num_ref).lambda = [];
            outputs(num_ref).vk_tau = [];
            outputs(num_ref).vk_theta = [];
            outputs(num_ref).vk_phi = [];
            outputs(num_ref).W = [];
            outputs(num_ref).Fg_r = [];
            outputs(num_ref).Fg_theta = [];
            outputs(num_ref).Fg_phi = [];
            outputs(num_ref).Fa_theta = [];
            outputs(num_ref).Fa_phi = [];
            outputs(num_ref).Fa_r = [];
            outputs(num_ref).Fa = [];
            outputs(num_ref).va_r = [];
            outputs(num_ref).va_theta = [];
            outputs(num_ref).va_phi = [];
            outputs(num_ref).F_dot_v = [];
            outputs(num_ref).D_r = [];
            outputs(num_ref).D_theta = [];
            outputs(num_ref).D_phi = [];
            outputs(num_ref).L_r = [];
            outputs(num_ref).L_theta = [];
            outputs(num_ref).L_phi = [];
            outputs(num_ref).D = [];
            outputs(num_ref).L = [];
            outputs(num_ref).Ft = [];
            outputs(num_ref).E_result = [];
            outputs(num_ref).kByE = [];
            outputs(num_ref).P_m_o_eff = [];
            outputs(num_ref).zetaMech = [];
            outputs(num_ref).etaGen_o = [];
            outputs(num_ref).P_e_o_eff = [];
            outputs(num_ref).theta_i = [];
            outputs(num_ref).phi_i = [];
            outputs(num_ref).vw_i = [];
            outputs(num_ref).vw_r_i = [];
            outputs(num_ref).vw_theta_i = [];
            outputs(num_ref).vw_phi_i = [];
            outputs(num_ref).f_i = [];
            outputs(num_ref).va_r_i = [];
            outputs(num_ref).va_theta_i = [];
            outputs(num_ref).va_phi_i = [];
            outputs(num_ref).va_i = [];
            outputs(num_ref).CD_k_i = [];
            outputs(num_ref).CD_i = [];
            outputs(num_ref).E_i = [];
            outputs(num_ref).Fa_i = [];
            outputs(num_ref).Fg_r_i = [];
            outputs(num_ref).Fg_theta_i = [];
            outputs(num_ref).Fg_phi_i = [];
            outputs(num_ref).Fa_theta_i = [];
            outputs(num_ref).Fa_phi_i = [];
            outputs(num_ref).Fa_r_i = [];
            outputs(num_ref).F_dot_v_i = [];
            outputs(num_ref).D_r_i = [];
            outputs(num_ref).D_theta_i = [];
            outputs(num_ref).D_phi_i = [];
            outputs(num_ref).L_r_i = [];
            outputs(num_ref).L_theta_i = [];
            outputs(num_ref).L_phi_i = [];
            outputs(num_ref).D_i = [];
            outputs(num_ref).L_i = [];
            outputs(num_ref).Ft_i = [];
            outputs(num_ref).E_result_i = [];
            outputs(num_ref).P_m_i_eff = [];
            outputs(num_ref).P_e_i_eff = [];
            outputs(num_ref).vk_r_i = [];
            outputs(num_ref).CL_i = [];
            outputs(num_ref).vk_r = [];
            outputs(num_ref).kRatio = [];
            outputs(num_ref).CL = [];
            outputs(num_ref).ti_eff = [];
            outputs(num_ref).to_eff = [];
            outputs(num_ref).rollAngle = [];

            obj.outputs = outputs;
        end
    end
end