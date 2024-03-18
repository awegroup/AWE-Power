classdef KiteQSMsimulation < handle
    properties
        inputs
        outputs
        optimDetails
        processedOutputs
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
        
        function obj = runsimulation(obj, optionsFMINCON)

            if nargin<2
                optionsFMINCON                      = optimoptions('fmincon');
                optionsFMINCON.Display              = 'none';
                optionsFMINCON.Algorithm            = 'sqp';
                optionsFMINCON.FiniteDifferenceType = 'central';
                optionsFMINCON.ScaleProblem         = true;
            else
                if ~isa(optionsFMINCON, 'optim.options.Fmincon')
                    error('KiteQSMsimulation:runSimulation Wrong set of options provided')
                end
            end

            % Initial guess for optimization
            x0 = obj.inputs.x0;

            % Bounds
            lb = obj.inputs.lb;
            ub = obj.inputs.ub;

            % Initialize variables to store results
            exitflag = zeros(1, length(obj.inputs.windSpeedReference));
            optHist = cell(1, length(obj.inputs.windSpeedReference));
            lambda = cell(1, length(obj.inputs.windSpeedReference));

            obj.initialise_outputs_to_zero();
            obj.outputs.halfRhoS = 0.5*obj.inputs.densityAir*obj.inputs.areaWing;

            % Kite mass 
            if obj.inputs.doMassOverride == 1
                obj.outputs.m_k = obj.inputs.kiteMass;
            else
                % Vincent Bonnin's simple mass model developed at Ampyx Power. Based on AP3 data and projected data for larger systems (AP4-AP5)      
                a1     = 0.002415;       a2     = 0.0090239;       b1     = 0.17025;       b2     = 3.2493;
                k1     = 5;              c1     = 0.46608;         d1     = 0.65962;       k2     = 1.1935;
                AR_ref = 12;
                a = a1*(obj.inputs.forceTether_max/1000/obj.inputs.areaWing) + a2;
                b = b1*(obj.inputs.forceTether_max/1000/obj.inputs.areaWing) + b2;
                obj.outputs.m_k = 10*(a*obj.inputs.areaWing^2 +b*obj.inputs.areaWing-k1)*(c1*(obj.inputs.aspectRatio/AR_ref)^2-d1*(obj.inputs.aspectRatio/AR_ref)+k2); 
            end

            % Optimize operation for every wind speed
            for i = 1:length(obj.inputs.windSpeedReference)

                % Define objective and constraint functions
                conFunc = @(x) constraints(obj, i);
                objFunc = @(x) objective(obj, x, i);

                % Run optimization
                [x, ~, exitflag(i), optHist{i}, lambda{i}] = fmincon(objFunc, x0, [], [], [], [], lb, ub, conFunc, optionsFMINCON);

                % Store final results
                obj.objective(x, i);

                % Update initial guess for the next iteration
                if abs(mean((obj.outputs.E_result - obj.outputs.E), 2)) > 1 %0.01
                    x0 = obj.inputs.x0;
                else
                    x0 = x;
                end
            end

            % Store optimization details
            obj.optimDetails.optHist = optHist;
            obj.optimDetails.exitflag = exitflag;
            obj.optimDetails.lambda = lambda;
        end

        function [c, ceq] = constraints(obj, i)
            %% Inequality constraints
            % Access outputs through obj.outputs

            % Min clearance between ground and bottom point of pattern
            c(1)   = (obj.inputs.minGroundClear - obj.outputs.pattStartGrClr(i));
            c(2)   = (obj.inputs.minGroundClear - obj.outputs.pattEndGrClr(i));

            % Capping for requested electrical rated power
            c(3)   = (obj.outputs.P_e_avg(i) - obj.inputs.P_ratedElec); 

            % Tether length limit
            c(4) = (obj.outputs.l_t_max(i) - obj.inputs.maxTeLen); 

            % Min number of patterns to get into transition 
            c(5) = (1 - mean(obj.outputs.numOfPatt(i,:)));

            % Max. cycle avg height
            c(6) = (obj.outputs.h_cycleEnd(i) - obj.inputs.maxHeight);

            % Peak mechanical power limit
            c(7:obj.inputs.numDeltaLelems+6) = (obj.outputs.P_m_o_eff(i,:) - obj.inputs.peakM2E_F*obj.inputs.P_ratedElec);

            % Maximum tether force
            c(7+obj.inputs.numDeltaLelems:2*obj.inputs.numDeltaLelems+6) = (obj.outputs.Ft(i,:) - obj.inputs.forceTether_max*obj.inputs.safetyFactor_forceTetherMax);

            % Apparent speed cannot be negative
            c(7+2*obj.inputs.numDeltaLelems:3*obj.inputs.numDeltaLelems+6) = (0 - obj.outputs.va(i,:));

            % Tangential velocity factor cannot be negative
            c(7+3*obj.inputs.numDeltaLelems:4*obj.inputs.numDeltaLelems+6) = (0 - obj.outputs.lambda(i,:));

            % Tether force during reel-in cannot be negative
            c(7+4*obj.inputs.numDeltaLelems:5*obj.inputs.numDeltaLelems+6) = (0 - obj.outputs.Ft_i(i,:));  

            %% Equality constraints
            % Glide ratio during reel-out
            ceq(0*obj.inputs.numDeltaLelems+1:1*obj.inputs.numDeltaLelems) = (obj.outputs.E_result(i,:) - obj.outputs.E(i,:));

            % Glide ratio during reel-in
            ceq(1*obj.inputs.numDeltaLelems+1:2*obj.inputs.numDeltaLelems) = (obj.outputs.E_result_i(i,:) - obj.outputs.E_i(i,:));
        end

        function [fval, obj] = objective(obj, x, i)
            % Assign variable values
            obj.outputs.deltaL(i) = x(1);
            obj.outputs.beta(i) = x(2);
            obj.outputs.gamma(i) = x(3);
            obj.outputs.Rp_start(i) = x(4);
            obj.outputs.vk_r_i(i, 1:obj.inputs.numDeltaLelems) = x(0*obj.inputs.numDeltaLelems+5:1*obj.inputs.numDeltaLelems+4);
            obj.outputs.CL_i(i, 1:obj.inputs.numDeltaLelems) = x(1*obj.inputs.numDeltaLelems+5:2*obj.inputs.numDeltaLelems+4);
            obj.outputs.vk_r(i, 1:obj.inputs.numDeltaLelems) = x(2*obj.inputs.numDeltaLelems+5:3*obj.inputs.numDeltaLelems+4);
            obj.outputs.kRatio(i, 1:obj.inputs.numDeltaLelems) = x(3*obj.inputs.numDeltaLelems+5:4*obj.inputs.numDeltaLelems+4);
            obj.outputs.CL(i, 1:obj.inputs.numDeltaLelems) = x(4*obj.inputs.numDeltaLelems+5:5*obj.inputs.numDeltaLelems+4);

            % Main computation
            obj.compute(i);

            % Objective: Maximise electrical cycle average power
            fval = -obj.outputs.P_e_avg(i);

            % Objective: Maximise mechanical reel-out power
%             fval = -obj.outputs.P_m_o(i);
        end

        function obj = compute(obj, i)
            
            
            % Constant for the cycle
            obj.outputs.l_t_min(i) = obj.outputs.Rp_start(i)/sin(obj.outputs.gamma(i));
            obj.outputs.pattStartGrClr(i) = obj.outputs.l_t_min(i)*sin(obj.outputs.beta(i)-obj.outputs.gamma(i));
            obj.outputs.h_cycleStart(i) = obj.outputs.l_t_min(i)*cos(obj.outputs.gamma(i))*sin(obj.outputs.beta(i));
            obj.outputs.l_t_max(i) = obj.outputs.l_t_min(i)+obj.outputs.deltaL(i)/cos(obj.outputs.gamma(i));
            obj.outputs.pattEndGrClr(i) = obj.outputs.l_t_max(i)*sin(obj.outputs.beta(i)-obj.outputs.gamma(i));
            obj.outputs.l_t_avg(i) = (obj.outputs.l_t_max(i)+obj.outputs.l_t_min(i))/2; %[m]
            obj.outputs.h_cycleAvg(i) = obj.outputs.l_t_avg(i)*cos(obj.outputs.gamma(i))*sin(obj.outputs.beta(i));
            obj.outputs.h_cycleEnd(i) = obj.outputs.l_t_max(i)*cos(obj.outputs.gamma(i))*sin(obj.outputs.beta(i));
            obj.outputs.d_t = sqrt(obj.inputs.forceTether_max/obj.inputs.maxStrengthTether*4/pi); %[m] safety factor could be added (say *1.1)
            
            % Discretizing the reel-out length in chosen number of elements
            % Found to be not highly sensitive to the number of elements
            obj.outputs.deltaLelems = obj.inputs.numDeltaLelems;
            obj.outputs.elemDeltaL(i) = obj.outputs.deltaL(i)/obj.outputs.deltaLelems;
            
            % Assigning and evaluating a single flight state equilibrium for each length element for reel-out and reel-in phase
            for j = 1:obj.outputs.deltaLelems
                % Reel-out phase:
                % Tether length at jth element
                if j == 1
                    obj.outputs.l_t_inCycle(i,j) = obj.outputs.l_t_min(i) + obj.outputs.elemDeltaL(i)/2/cos(obj.outputs.gamma(i));
                else
                    obj.outputs.l_t_inCycle(i,j) = obj.outputs.l_t_inCycle(i,j-1) + obj.outputs.elemDeltaL(i)/cos(obj.outputs.gamma(i));
                end
                
                % Pattern ground clearance
                obj.outputs.pattGrClr(i,j) = obj.outputs.l_t_inCycle(i,j)*sin(obj.outputs.beta(i)-obj.outputs.gamma(i));
            
                % Effective mass lumped at kite point (Kite + tether)
                obj.outputs.m_t(i,j) = obj.inputs.densityTether*pi/4*obj.outputs.d_t^2*obj.outputs.l_t_inCycle(i,j);
                obj.outputs.m_eff(i,j) = obj.outputs.m_k + 0.5*obj.outputs.m_t(i,j);
            
                % Coordinates of the Kite's position and orientation in the Spherical ref. frame
                % Center point (representative point)
                obj.outputs.theta(i,j) = pi/2 - (obj.outputs.beta(i));
                obj.outputs.phi(i,j) = 0;
                obj.outputs.chi(i,j) = deg2rad(90);
            
                % Effective CD
                obj.outputs.CD_k(i,j) = obj.inputs.Cd0 + (obj.outputs.CL(i,j)-obj.inputs.Cl0_airfoil)^2/(pi*obj.inputs.aspectRatio*obj.inputs.e);
                obj.outputs.CD_t(i,j) = (1/4)*obj.inputs.Cd_cylinder*obj.outputs.d_t*obj.outputs.l_t_inCycle(i,j)/obj.inputs.areaWing;
                obj.outputs.CD(i,j) = obj.outputs.CD_k(i,j) + obj.outputs.CD_t(i,j);
            
                % Effective Glide ratio
                obj.outputs.E(i,j) = obj.outputs.CL(i,j)/obj.outputs.CD(i,j);
            
                % Average pattern height at point of interest on deltaL
                if j == 1
                    obj.outputs.h_inCycle(i,j) = obj.outputs.h_cycleStart(i) + obj.outputs.elemDeltaL(i)/2*sin(obj.outputs.beta(i));
                else
                    obj.outputs.h_inCycle(i,j) = obj.outputs.h_inCycle(i,j-1) + obj.outputs.elemDeltaL(i)*sin(obj.outputs.beta(i));
                end
            
                % Pattern radius at point of interest on deltaL
                if j == 1
                    obj.outputs.Rp(i,j) = obj.outputs.Rp_start(i) + obj.outputs.elemDeltaL(i)/2*tan(obj.outputs.gamma(i));
                else
                    obj.outputs.Rp(i,j) = obj.outputs.Rp(i,j-1) + obj.outputs.elemDeltaL(i)*tan(obj.outputs.gamma(i));
                end
                
                if ~obj.inputs.doUseReferenceWindProfile
                    % Modelled
                    obj.outputs.vw(i,j) = obj.inputs.windSpeedReference(i)*((obj.outputs.h_inCycle(i,j))/obj.inputs.heightWindReference)^obj.inputs.windShearExp;
                else
                    % Extrapolated from dataset
                    obj.outputs.vw(i,j) = obj.inputs.windSpeedReference(i)*interp1(obj.inputs.windProfile_h, obj.inputs.windProfile_vw, (obj.outputs.h_inCycle(i,j)), 'linear', 'extrap');
                end
            
                % Wind velocity vector
                obj.outputs.vw_r(i,j) = obj.outputs.vw(i,j)*sin(obj.outputs.theta(i,j))*cos(obj.outputs.phi(i,j));
                obj.outputs.vw_theta(i,j) = obj.outputs.vw(i,j)*cos(obj.outputs.theta(i,j))*cos(obj.outputs.phi(i,j));
                obj.outputs.vw_phi(i,j) = -obj.outputs.vw(i,j)*sin(obj.outputs.phi(i,j));
            
                % Reel-out factor
                obj.outputs.f(i,j) = obj.outputs.vk_r(i,j)/obj.outputs.vw(i,j);
            
                % Apparent wind velocity magnitude
                obj.outputs.va(i,j) = (obj.outputs.vw_r(i,j)-obj.outputs.vk_r(i,j))*sqrt(1+obj.outputs.kRatio(i,j)^2);
            
                % Aerodynamic force magnitude
                obj.outputs.Fa(i,j) = obj.outputs.halfRhoS*sqrt(obj.outputs.CL(i,j)^2+obj.outputs.CD(i,j)^2)*obj.outputs.va(i,j)^2;
            
                % Tangential kite velocity factor
                a = cos(obj.outputs.theta(i,j))*cos(obj.outputs.phi(i,j))*cos(obj.outputs.chi(i,j))-sin(obj.outputs.phi(i,j))*sin(obj.outputs.chi(i,j));
                b = sin(obj.outputs.theta(i,j))*cos(obj.outputs.phi(i,j));
                obj.outputs.lambda(i,j) = a + sqrt(a^2 + b^2 - 1 + obj.outputs.kRatio(i,j)^2*(b - obj.outputs.f(i,j))^2);
                % Tangential kite velocity
                obj.outputs.vk_tau(i,j) = obj.outputs.lambda(i,j)*obj.outputs.vw(i,j);
                % Tangential kite velocity components in theta and phi directions
                obj.outputs.vk_theta(i,j) = obj.outputs.vk_tau(i,j)*cos(obj.outputs.chi(i,j));
                obj.outputs.vk_phi(i,j) = obj.outputs.vk_tau(i,j)*sin(obj.outputs.chi(i,j));
            
                % Gravity toggle
                if obj.inputs.doIncludeGravity
                    obj.outputs.W(i,j) = obj.outputs.m_eff(i,j)*obj.inputs.accGravity;
                else
                    obj.outputs.W(i,j) = 0;
                end
            
                % Gravitational force vector (kite + tether)
                obj.outputs.Fg_r(i,j) = -obj.outputs.W(i,j)*cos(obj.outputs.theta(i,j));
                obj.outputs.Fg_theta(i,j) = obj.outputs.W(i,j)*sin(obj.outputs.theta(i,j));
                obj.outputs.Fg_phi(i,j) = 0;
            
                % Aerodynamic force vector
                obj.outputs.Fa_theta(i,j) = -obj.outputs.Fg_theta(i,j);
                obj.outputs.Fa_phi(i,j) = -obj.outputs.Fg_phi(i,j);
                obj.outputs.Fa_r(i,j) = sqrt(obj.outputs.Fa(i,j)^2-obj.outputs.Fa_theta(i,j)^2-obj.outputs.Fa_phi(i,j)^2);
            
                % Roll angle in case of Gravity (N.A when CF is included)
                obj.outputs.rollAngle(i,j) = atan(obj.outputs.Fa_theta(i,j)/obj.outputs.Fa_r(i,j));
            
                % Apparent wind velocity vector
                obj.outputs.va_r(i,j) = obj.outputs.vw(i,j)*(sin(obj.outputs.theta(i,j))*cos(obj.outputs.phi(i,j))-obj.outputs.f(i,j));
                obj.outputs.va_theta(i,j) = obj.outputs.vw(i,j)*(cos(obj.outputs.theta(i,j))*cos(obj.outputs.phi(i,j))-obj.outputs.lambda(i,j)*cos(obj.outputs.chi(i,j)));
                obj.outputs.va_phi(i,j) = obj.outputs.vw(i,j)*(-sin(obj.outputs.phi(i,j))-obj.outputs.lambda(i,j)*sin(obj.outputs.chi(i,j)));
            
                % Dot product F_a*v_a;
                obj.outputs.F_dot_v(i,j) = obj.outputs.Fa_r(i,j)*obj.outputs.va_r(i,j) + obj.outputs.Fa_theta(i,j)*obj.outputs.va_theta(i,j) + obj.outputs.Fa_phi(i,j)*obj.outputs.va_phi(i,j);
            
                % Drag vector from dot product
                obj.outputs.D_r(i,j) = (obj.outputs.F_dot_v(i,j)/obj.outputs.va(i,j)^2)*obj.outputs.va_r(i,j);
                obj.outputs.D_theta(i,j) = (obj.outputs.F_dot_v(i,j)/obj.outputs.va(i,j)^2)*obj.outputs.va_theta(i,j);
                obj.outputs.D_phi(i,j) = (obj.outputs.F_dot_v(i,j)/obj.outputs.va(i,j)^2)*obj.outputs.va_phi(i,j);
            
                % Lift vector
                obj.outputs.L_r(i,j) = obj.outputs.Fa_r(i,j) - obj.outputs.D_r(i,j);
                obj.outputs.L_theta(i,j) = obj.outputs.Fa_theta(i,j) - obj.outputs.D_theta(i,j);
                obj.outputs.L_phi(i,j) = obj.outputs.Fa_phi(i,j) - obj.outputs.D_phi(i,j);
            
                % Drag magnitude
                obj.outputs.D(i,j) = sqrt(obj.outputs.D_r(i,j)^2 + obj.outputs.D_theta(i,j)^2 + obj.outputs.D_phi(i,j)^2);
            
                % Lift magnitude
                obj.outputs.L(i,j) = sqrt(obj.outputs.L_r(i,j)^2 + obj.outputs.L_theta(i,j)^2 + obj.outputs.L_phi(i,j)^2);
            
                % Straight-tether force
                obj.outputs.Ft(i,j) = obj.outputs.Fa_r(i,j) + obj.outputs.Fg_r(i,j);
            
                % Lift-to-drag ratio that follows from the chosen kinematic ratio
                obj.outputs.E_result(i,j) = sqrt(((obj.outputs.Fa(i,j)*obj.outputs.va(i,j))/obj.outputs.F_dot_v(i,j))^2-1);
            
                % Ratio of Kinemactic Ratio and Glide Ratio
                obj.outputs.kByE(i,j) = obj.outputs.kRatio(i,j)/obj.outputs.E_result(i,j);
            
                % Effective mechanical reel-out power
                obj.outputs.P_m_o_eff(i,j) = obj.outputs.Ft(i,j)*obj.outputs.vk_r(i,j); %[W]
            
                obj.outputs.zetaMech(i,j) = obj.outputs.P_m_o_eff(i,j)/(obj.outputs.halfRhoS*obj.outputs.vw(i,j)^3);
            
                % Effective electrical reel-out power
                % Generator efficiency. As a function of RPM/RPM_max, where RPM_max is driven by winch
                vkr_ratio = obj.outputs.vk_r(i,j)/obj.inputs.speedGeneratorReelout_max;
                obj.outputs.etaGen_o(i,j) = obj.inputs.etaGen_param(1)*vkr_ratio^3 + ...
                    obj.inputs.etaGen_param(2)*vkr_ratio^2 + ...
                    obj.inputs.etaGen_param(3)*vkr_ratio+obj.inputs.etaGen_param(4);
                obj.outputs.P_e_o_eff(i,j) = obj.outputs.P_m_o_eff(i,j)*obj.inputs.efficiency_Gearbox*obj.outputs.etaGen_o(i,j)*obj.inputs.efficiency_PowerElectronics;
            
                %% Retraction Phase
                % Kite position in spherical coordinates
                % Reel-in is assumed to start from the top of the pattern
                obj.outputs.theta_i(i,j) = pi/2 - (obj.outputs.beta(i)+obj.outputs.gamma(i));
                obj.outputs.phi_i(i,j) = 0;
            
                if ~obj.inputs.doUseReferenceWindProfile
                    % Modelled
                    obj.outputs.vw_i(i,j) = obj.inputs.windSpeedReference(i)*((obj.outputs.h_inCycle(i,j) + obj.outputs.Rp(i,j)*cos(obj.outputs.beta(i)))/obj.inputs.heightWindReference)^obj.inputs.windShearExp;
                else
                    % Extrapolated from dataset
                    obj.outputs.vw_i(i,j) = obj.inputs.windSpeedReference(i)*interp1(obj.inputs.windProfile_h, obj.inputs.windProfile_vw, (obj.outputs.h_inCycle(i,j) + obj.outputs.Rp(i,j)*cos(obj.outputs.beta(i))), 'linear', 'extrap');
                end
            
                % Wind velocity vector
                obj.outputs.vw_r_i(i,j) = obj.outputs.vw_i(i,j)*sin(obj.outputs.theta_i(i,j))*cos(obj.outputs.phi_i(i,j));
                obj.outputs.vw_theta_i(i,j) = obj.outputs.vw_i(i,j)*cos(obj.outputs.theta_i(i,j))*cos(obj.outputs.phi_i(i,j));
                obj.outputs.vw_phi_i(i,j) = -obj.outputs.vw_i(i,j)*sin(obj.outputs.phi_i(i,j)); %
            
                % Reel-in factor
                obj.outputs.f_i(i,j) = obj.outputs.vk_r_i(i,j)/obj.outputs.vw_i(i,j);
            
                % Apparent speed vector
                obj.outputs.va_r_i(i,j) = obj.outputs.vw_r_i(i,j) - (-obj.outputs.vk_r_i(i,j)); % Reel-in speed direction is in the negative radial direction)
                obj.outputs.va_theta_i(i,j) = obj.outputs.vw_theta_i(i,j) - 0;
                obj.outputs.va_phi_i(i,j) = obj.outputs.vw_phi_i(i,j) - 0;
            
                % Apparent wind velocity magnitude
                obj.outputs.va_i(i,j) = sqrt(obj.outputs.va_r_i(i,j)^2 + obj.outputs.va_theta_i(i,j)^2); % Wind has components only in r and theta directions
            
                % Aerodynamic force magnitude
                obj.outputs.CD_k_i(i,j) = obj.inputs.Cd0+(obj.outputs.CL_i(i,j)- obj.inputs.Cl0_airfoil)^2/(pi*obj.inputs.aspectRatio*obj.inputs.e);
                obj.outputs.CD_i(i,j) = obj.outputs.CD_k_i(i,j) + obj.outputs.CD_t(i,j);
                obj.outputs.E_i(i,j) = obj.outputs.CL_i(i,j)/obj.outputs.CD_i(i,j);
                obj.outputs.Fa_i(i,j) = obj.outputs.halfRhoS*sqrt(obj.outputs.CL_i(i,j)^2+obj.outputs.CD_i(i,j)^2)*obj.outputs.va_i(i,j)^2;
            
                % Gravitational force vector (kite + tether)
                obj.outputs.Fg_r_i(i,j) = -obj.outputs.W(i,j)*cos(obj.outputs.theta_i(i,j));
                obj.outputs.Fg_theta_i(i,j) = obj.outputs.W(i,j)*sin(obj.outputs.theta_i(i,j));
                obj.outputs.Fg_phi_i(i,j) = 0;
            
                % Aerodynamic force vector
                obj.outputs.Fa_theta_i(i,j) = -obj.outputs.Fg_theta_i(i,j);
                obj.outputs.Fa_phi_i(i,j) = -obj.outputs.Fg_phi_i(i,j);
                obj.outputs.Fa_r_i(i,j) = sqrt(obj.outputs.Fa_i(i,j)^2-obj.outputs.Fa_theta_i(i,j)^2-obj.outputs.Fa_phi_i(i,j)^2);
            
                % Dot product F_a*v_a;
                obj.outputs.F_dot_v_i(i,j) = obj.outputs.Fa_r_i(i,j)*obj.outputs.va_r_i(i,j) + obj.outputs.Fa_theta_i(i,j)*obj.outputs.va_theta_i(i,j) + obj.outputs.Fa_phi_i(i,j)*obj.outputs.va_phi_i(i,j);
            
                % Drag vector from dot product
                obj.outputs.D_r_i(i,j) = (obj.outputs.F_dot_v_i(i,j)/obj.outputs.va_i(i,j)^2)*obj.outputs.va_r_i(i,j);
                obj.outputs.D_theta_i(i,j) = (obj.outputs.F_dot_v_i(i,j)/obj.outputs.va_i(i,j)^2)*obj.outputs.va_theta_i(i,j);
                obj.outputs.D_phi_i(i,j) = (obj.outputs.F_dot_v_i(i,j)/obj.outputs.va_i(i,j)^2)*obj.outputs.va_phi_i(i,j);
            
                % Lift vector
                obj.outputs.L_r_i(i,j) = obj.outputs.Fa_r_i(i,j) - obj.outputs.D_r_i(i,j);
                obj.outputs.L_theta_i(i,j) = obj.outputs.Fa_theta_i(i,j) - obj.outputs.D_theta_i(i,j);
                obj.outputs.L_phi_i(i,j) = obj.outputs.Fa_phi_i(i,j) - obj.outputs.D_phi_i(i,j);
            
                % Drag magnitude
                obj.outputs.D_i(i,j) = sqrt(obj.outputs.D_r_i(i,j)^2 + obj.outputs.D_theta_i(i,j)^2 + obj.outputs.D_phi_i(i,j)^2);
            
                % Lift magnitude
                obj.outputs.L_i(i,j) = sqrt(obj.outputs.L_r_i(i,j)^2 + obj.outputs.L_theta_i(i,j)^2 + obj.outputs.L_phi_i(i,j)^2);
            
                % Lift-to-drag ratio that follows from the chosen kinematic ratio
                obj.outputs.E_result_i(i,j) = sqrt(((obj.outputs.Fa_i(i,j)*obj.outputs.va_i(i,j))/obj.outputs.F_dot_v_i(i,j))^2-1);
            
                % Straight-tether force
                obj.outputs.Ft_i(i,j) = obj.outputs.Fa_r_i(i,j) + obj.outputs.Fg_r_i(i,j);
            
                % Effective mechanical reel-out power
                obj.outputs.P_m_i_eff(i,j) = obj.outputs.Ft_i(i,j)*obj.outputs.vk_r_i(i,j); %[W]
            
                % Generator efficiency during RI: As a function of RPM/RPM_max, where RPM_max is driven by winch i.e Max VRI
                vkri_ratio = obj.outputs.vk_r_i(i)/obj.inputs.speedGeneratorReelout_max;
                obj.outputs.etaGen_i(i) = obj.inputs.etaGen_param(1)*vkri_ratio^3 + ...
                    obj.inputs.etaGen_param(2)*vkri_ratio^2 + ...
                    obj.inputs.etaGen_param(3)*vkri_ratio+obj.inputs.etaGen_param(4);
            
                % Effective electrical reel-in power
                obj.outputs.P_e_i_eff(i,j) = obj.outputs.P_m_i_eff(i,j)/obj.inputs.efficiency_Gearbox/obj.inputs.efficiency_Storage/obj.outputs.etaGen_i(i)/obj.inputs.efficiency_PowerElectronics;
            
            end

            
            %% Cycle calculation
            % Reel-out time
            obj.outputs.t1(i)       = obj.outputs.vk_r(i,1)/obj.inputs.accReel_max;
            obj.outputs.to_eff(i,:) = obj.outputs.elemDeltaL(i)./obj.outputs.vk_r(i,:);
            obj.outputs.to(i)       = obj.outputs.t1(i) + sum(obj.outputs.to_eff(i,:));
            
            % Reel-out power during transition
            obj.outputs.P1_m_o(i) = obj.outputs.P_m_o_eff(i,1)/2;
            obj.outputs.P1_e_o(i) = obj.outputs.P_e_o_eff(i,1)/2;
            
            % Reel-out power 
            obj.outputs.P_m_o(i) = (sum(obj.outputs.P_m_o_eff(i,:).*obj.outputs.to_eff(i,:)) + obj.outputs.P1_m_o(i)*obj.outputs.t1(i))/obj.outputs.to(i);
            obj.outputs.P_e_o(i) = (sum(obj.outputs.P_e_o_eff(i,:).*obj.outputs.to_eff(i,:)) + obj.outputs.P1_e_o(i)*obj.outputs.t1(i))/obj.outputs.to(i);
            
            % Reel-in time
            obj.outputs.t2(i)       = obj.outputs.vk_r_i(i)/obj.inputs.accReel_max;
            obj.outputs.ti_eff(i,:) = obj.outputs.elemDeltaL(i)./obj.outputs.vk_r_i(i,:);
            obj.outputs.ti(i)       = obj.outputs.t2(i) + sum(obj.outputs.ti_eff(i,:));
            
            % Reel-in power during transition
            obj.outputs.P2_m_i(i)     = obj.outputs.P_m_i_eff(i,1)/2;
            obj.outputs.P2_e_i(i)     = obj.outputs.P_e_i_eff(i,1)/2;
            
            % Reel-in power 
            obj.outputs.P_m_i(i) = (sum(obj.outputs.P_m_i_eff(i,:).*obj.outputs.ti_eff(i,:)) + obj.outputs.P2_m_i(i)*obj.outputs.t2(i))/obj.outputs.ti(i);
            obj.outputs.P_e_i(i) = (sum(obj.outputs.P_e_i_eff(i,:).*obj.outputs.ti_eff(i,:)) + obj.outputs.P2_e_i(i)*obj.outputs.t2(i))/obj.outputs.ti(i);
            
            % Cycle time
            obj.outputs.tCycle(i) = obj.outputs.to(i)+obj.outputs.ti(i);
            
            % Time for one pattern revolution and number of patterns in the cycle
            obj.outputs.tPatt(i,:) = 2*pi*obj.outputs.Rp(i,:)./obj.outputs.vk_phi(i,:);
            obj.outputs.numOfPatt(i,:) = obj.outputs.to(i)./obj.outputs.tPatt(i,:);
            
            % Electrical cycle power
            obj.outputs.P_e_avg(i) = (sum(obj.outputs.to_eff(i,:).*obj.outputs.P_e_o_eff(i,:)) + obj.outputs.t1(i)*obj.outputs.P1_e_o(i) - ...
                                       sum(obj.outputs.ti_eff(i,:).*obj.outputs.P_e_i_eff(i,:)) - obj.outputs.t2(i)*obj.outputs.P2_e_i(i))/obj.outputs.tCycle(i);
            
            % Mechanical cycle power - without drivetrain eff
            obj.outputs.P_m_avg(i) = (sum(obj.outputs.to_eff(i,:).*obj.outputs.P_m_o_eff(i,:)) + obj.outputs.t1(i)*obj.outputs.P1_m_o(i) - ...
                                       sum(obj.outputs.ti_eff(i,:).*obj.outputs.P_m_i_eff(i,:)) - obj.outputs.t2(i)*obj.outputs.P2_m_i(i))/obj.outputs.tCycle(i);

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
        
        function obj = plotresults(obj)
                        
            %% Plot settings
            vw = obj.inputs.windSpeedReference;
            x_axis_limit = [0 obj.processedOutputs.vw_100m_operRange(end)];
            
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
            ws = obj.processedOutputs.ratedWind;
            fig = figure();
            colororder(newcolors)
            % Lengths
            subplot(2,2,1)
            hold on; grid on; box on
            plot(obj.processedOutputs.Rp(ws,:),':o','linewidth',1,'markersize',3);
            plot(obj.processedOutputs.l_t_inCycle(ws,:),':^','linewidth',1,'markersize',3);
            plot(obj.processedOutputs.h_inCycle(ws,:),':s','linewidth',1,'markersize',3);
            ylabel('(m)');
            legend('R_{p}','l_{t}','h_{p}','location','northwest');
            hold off
            % Speeds
            subplot(2,2,2); hold on; grid on; box on; 
            yyaxis left
            plot(obj.processedOutputs.lambda(ws,:),':o','linewidth',1,'markersize',3);
            plot(obj.processedOutputs.f(ws,:),':s','linewidth',1,'markersize',3);
            plot(obj.processedOutputs.f_i(ws,:),':v','linewidth',1,'markersize',3);
            ylabel('(-)');
            yyaxis right
            plot(obj.processedOutputs.vw(ws,:),':s','linewidth',1,'markersize',3);
            plot(obj.processedOutputs.vw_i(ws,:),':v','linewidth',1,'markersize',3);
            ylabel('(m/s)');
            legend('λ','f_o','f_i','v_{w,o}','v_{w,i}','location','northwest');
            hold off
            % Glide ratios
            subplot(2,2,3); hold on; grid on; box on
            plot(obj.processedOutputs.E(ws,:),':s','linewidth',1,'markersize',3);
            plot(obj.processedOutputs.E_i(ws,:),':v','linewidth',1,'markersize',3);
            % ylim([0 2.5]);
            ylabel('(-)');
            legend('E_o','E_i','location','northwest');
            hold off
            % Forces
            subplot(2,2,4); hold on; grid on; box on
            plot(obj.processedOutputs.Fa(ws,:)/1e3,':^','linewidth',1,'markersize',3);
            plot(obj.processedOutputs.Fa_i(ws,:)/1e3,':o','linewidth',1,'markersize',3);
            plot(obj.processedOutputs.Ft(ws,:)/1e3,':s','linewidth',1,'markersize',3);
            plot(obj.processedOutputs.Ft_i(ws,:)/1e3,':v','linewidth',1,'markersize',3);
            plot(obj.processedOutputs.W(ws,:)/1e3,':x','linewidth',1,'markersize',3);
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
            plot(vw, mean(obj.processedOutputs.lambda,2),':o','linewidth',1,'markersize',3);
            plot(vw, mean(obj.processedOutputs.f,2),':s','linewidth',1,'markersize',3);
            plot(vw, mean(obj.processedOutputs.f_i,2),':v','linewidth',1,'markersize',3);
            ylabel('(-)');
            legend('λ','f_{o}','f_{i}','location','northwest');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            %ylim([0 160]);
            hold off
        
            % CL
            figure('units','inch','Position', [4.2 0.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            plot(vw, mean(obj.processedOutputs.CL,2),'o:','linewidth',1,'markersize',3);
            plot(vw, mean(obj.processedOutputs.CL_i,2),'d:','linewidth',1,'markersize',3);
            ylabel('(-)');
            legend('C_{L,o}','C_{L,i}','location','northwest','Orientation','vertical');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            %ylim([0 160]);
            hold off
        
            % CD
            figure('units','inch','Position', [8.2 0.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            plot(vw, mean(obj.processedOutputs.CD,2),'o:','linewidth',1,'markersize',3);
            plot(vw, mean(obj.processedOutputs.CD_i,2),'d:','linewidth',1,'markersize',3);
            plot(vw, mean(obj.processedOutputs.CD_k,2),'^:','linewidth',1,'markersize',3);
            plot(vw, mean(obj.processedOutputs.CD_k_i,2),'v:','linewidth',1,'markersize',3);
            plot(vw, mean(obj.processedOutputs.CD_t,2),'s:','linewidth',1,'markersize',3);
            ylabel('(-)');
            legend('C_{D,o}','C_{D,i}','C_{D,k,o}','C_{D,k,i}','C_{D,t}','location','northwest','Orientation','vertical');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            %ylim([0 160]);
            hold off
            
            % Lengths
            figure('units','inch','Position', [12.2 0.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            plot(vw, obj.processedOutputs.h_cycleAvg,'d:','linewidth',1,'markersize',3);
            plot(vw, mean(obj.processedOutputs.Rp,2),'o:','linewidth',1,'markersize',3);
            plot(vw, obj.processedOutputs.deltaL,'^:','linewidth',1,'markersize',3);
            plot(vw, obj.processedOutputs.l_t_max,'s:','linewidth',1,'markersize',3);
            plot(vw, obj.processedOutputs.l_t_min,'x:','linewidth',1,'markersize',3);
            ylabel('(m)');
            legend('h_{p,avg}','R_{p,avg}','Δl','l_{t,max}','l_{t,min}','location','northwest','Orientation','vertical');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            %ylim([0 160]);
            hold off
        
            % Roll angle, avg patt elevation
            figure('units','inch','Position', [16.2 0.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            plot(vw, mean(obj.processedOutputs.rollAngle,2),'s:','linewidth',1,'markersize',3);
            plot(vw, obj.processedOutputs.beta,'o:','linewidth',1,'markersize',3);
            plot(vw, obj.processedOutputs.gamma,'d:','linewidth',1,'markersize',3);
            ylabel('(deg)');
            legend('$$\Psi$$','$$\overline{\beta}$$','$$\gamma$$','location','northwest','Orientation','vertical','Interpreter','latex');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            %ylim([0 160]);
            hold off
        
            % Reel-out forces
            figure('units','inch','Position', [0.2 3.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            plot(vw, mean(obj.processedOutputs.W,2)./10^3,':o','markersize',3)
            plot(vw, mean(obj.processedOutputs.Fa,2)./10^3,':s','linewidth',1,'markersize',3);
            plot(vw, mean(obj.processedOutputs.Ft,2)./10^3,':^','linewidth',1,'markersize',3);
            ylabel('(kN)');
            legend('W','F_{a,o}','F_{t,o}','location','northwest');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            %ylim([0 160]);
            hold off
        
            % Reel-in forces
            figure('units','inch','Position', [4.2 3.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            plot(vw, mean(obj.processedOutputs.W,2)./10^3,':o','markersize',3)
            plot(vw, mean(obj.processedOutputs.Fa_i,2)./10^3,':s','linewidth',1,'markersize',3);
            plot(vw, mean(obj.processedOutputs.Ft_i,2)./10^3,':^','linewidth',1,'markersize',3);
            ylabel('(kN)');
            legend('W','F_{a,i}','F_{t,i}','location','northwest');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            hold off
        
            % Time, num of patterns
            figure('units','inch','Position', [8.2 3.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            yyaxis left
            plot(vw, obj.processedOutputs.to,':s','linewidth',1,'markersize',3);
            plot(vw, obj.processedOutputs.ti,':d','linewidth',1,'markersize',3);
            plot(vw, mean(obj.processedOutputs.tPatt,2),':o','linewidth',1,'markersize',3);
            ylabel('(s)');
            yyaxis right
            plot(vw, mean(obj.processedOutputs.numOfPatt,2),':^','linewidth',1,'markersize',3);
            ylabel('(-)');
            xlabel('Wind speed at 100m height (m/s)');
            legend('t_{o}','t_{i}','t_{patt,avg}','N_{p}','location','northwest');
            xlim(x_axis_limit);
            hold off
        
            % Reel-out power
            figure('units','inch','Position', [12.2 3.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            plot(vw, obj.processedOutputs.P_m_o./10^3,':o','linewidth',1,'markersize',3);
            plot(vw, obj.processedOutputs.P_e_o./10^3,':s','linewidth',1,'markersize',3);
            plot(vw, obj.processedOutputs.P_e_avg./10^3,':x','linewidth',1,'markersize',5);
            ylabel('(kW)');
            %title('Cycle averages');
            legend('P_{m,o}','P_{e,o}','P_{e,avg}','location','northwest');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            hold off
        
            % Reel-in power
            figure('units','inch','Position', [16.2 3.5 3.5 2.2]); hold on; grid on; box on
            colororder(newcolors)
            plot(vw, obj.processedOutputs.P_m_i./10^3,':o','linewidth',1,'markersize',3);
            plot(vw, obj.processedOutputs.P_e_i./10^3,':s','linewidth',1,'markersize',3);
            ylabel('(kW)');
            %title('Cycle averages');
            legend('P_{m,i}','P_{e,i}','location','northwest');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            hold off

            % Power curve comparison plot
            % Loyd
            P_Loyd = zeros(1,length(vw));
            CL_loyd = obj.inputs.Cl_maxAirfoil * obj.inputs.Cl_eff_F;
            CD_loyd = obj.inputs.Cd0 + (CL_loyd - obj.inputs.Cl0_airfoil)^2 / (pi * obj.inputs.aspectRatio * obj.inputs.e);
            for i = 1:length(vw)
                P_Loyd(i) = (4/27) * (CL_loyd^3 / CD_loyd^2) * (1/2) * obj.inputs.densityAir * obj.inputs.areaWing * (vw(i)^3);
                if P_Loyd(i) > obj.inputs.P_ratedElec
                    P_Loyd(i) = obj.inputs.P_ratedElec;
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
            plot(vw, obj.processedOutputs.P_e_avg./10^3,':s','linewidth',1,'MarkerSize',3);
            plot(AP3_PC_ws, AP3_PC_power,':^k','linewidth',1,'MarkerSize',3);
            ylabel('Power (kW)');
            legend('Ideal Reel-out','QSM','6-DOF','location','southeast');
            xlabel('Wind speed at 100m height (m/s)');
            xlim(x_axis_limit);
            hold off

            % Cycle timeseries plots: Pattern averages
            windSpeeds = [obj.processedOutputs.cutIn, obj.processedOutputs.ratedWind, obj.processedOutputs.cutOut];
            for i = windSpeeds
                idx = find(vw == i);
                figure('units','inch','Position', [4.2 6.5 3.5 2.2])
                hold on
                grid on
                box on
                % Plot your data
                plot(obj.processedOutputs.cyclePowerRep(i).t_inst, obj.processedOutputs.cyclePowerRep(i).P_e_inst,'linewidth',1.5, 'DisplayName', 'P_{e}');
                plot(obj.processedOutputs.cyclePowerRep(i).t_inst, obj.processedOutputs.cyclePowerRep(i).P_m_inst,'--','linewidth',1.5, 'DisplayName', 'P_{m}');
                yline(obj.processedOutputs.P_e_avg(idx)/10^3, ':', 'linewidth',1.2, 'DisplayName', 'P_{e,avg}');
                % Set yline positions within the range of your data
                yline(0, '-', 'linewidth', 1, 'HandleVisibility', 'off');      
                ylabel('(kW)');
                xlabel('Time (s)');      
                % Create a custom legend with 'Data1' and 'Data2' only
                legend('Location', 'Best');      
                % Customize legend colors
                legend('TextColor', 'k', 'Color', 'w'); % 'TextColor' sets the text color, 'Color' sets the background color      
                xlim([0 obj.processedOutputs.cyclePowerRep(i).t_inst(end)]);
                ylim([1.5*min(obj.processedOutputs.cyclePowerRep(i).P_e_inst) 1.05*max(obj.processedOutputs.cyclePowerRep(i).P_m_inst)]);
                title(strcat('Wind speed at 100m:', num2str(obj.processedOutputs.cyclePowerRep(idx).ws), ' m/s'));      
                hold off
            end
            
            % Wind profile
            if ~obj.inputs.doUseReferenceWindProfile
                Vref = 10; % m/s
                z = 10:10:600; % m
                V = Vref * (z/obj.inputs.heightWindReference).^obj.inputs.windShearExp/Vref;
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
                V_interpolated = interp1(obj.inputs.windProfile_h, obj.inputs.windProfile_vw, heights, 'linear', 'extrap');
                % Plot the interpolated function
                figure('units','inch','Position', [15 6 2 2.2])
                hold on 
                box on
                grid on
                plot(obj.inputs.windProfile_vw, obj.inputs.windProfile_h, 'o', V_interpolated, heights, '-');
                xlabel('Wind speed (m/s)');
                ylabel('Height (m)');
                xlim([0.5 1.5])
                title('Wind profile');
                legend('Data', 'Interpolated', 'Location', 'Best');
                hold off
            end


        end
        
        function obj = processoutputs(obj)
            %% Post processing
            vw = obj.inputs.windSpeedReference; % Wind speed at ref. height
            
            %% Cut-in wind speed
            % Glide ratio constraint violation acceptance for feasible solution
            temp1                  = vw(abs(mean((obj.outputs.E_result - obj.outputs.E),2)) < 0.01); %0.01
            obj.processedOutputs.cutIn = max(temp1(1));
            
            %% Rated wind and power
            temp3                       = round(obj.outputs.P_e_avg./max(obj.outputs.P_e_avg),2);
            temp4                       = vw(temp3>=0.99);
            obj.processedOutputs.ratedWind  = temp4(1); % At Reference height 
            obj.processedOutputs.ratedPower = obj.outputs.P_e_avg(vw==temp4(1));
            
            %% Cut-out wind speed
            obj.processedOutputs.cutOut   = 25; % At operational height (Assumption)
            % Operating range traslated to wind speed at ref. height
            obj.processedOutputs.vw_100m_operRange = vw(mean(obj.outputs.vw,2)<=obj.processedOutputs.cutOut);
            
            %% Extract feasible results in the operational range
            obj.processedOutputs.Dia_te = obj.outputs.d_t;
            for i=1:length(vw)
                if vw(i)>=obj.processedOutputs.cutIn && vw(i)<= obj.processedOutputs.cutOut
                    obj.processedOutputs.vw(i,:)           = obj.outputs.vw(i,:);
                    obj.processedOutputs.vw_i(i,:)         = obj.outputs.vw_i(i,:);
                    obj.processedOutputs.m_k               = obj.outputs.m_k;  
                    obj.processedOutputs.P1_m_o(i)         = obj.outputs.P1_m_o(i);
                    obj.processedOutputs.P2_m_i(i)         = obj.outputs.P2_m_i(i); 
                    obj.processedOutputs.P_m_o_eff(i,:)    = obj.outputs.P_m_o_eff(i,:);
                    obj.processedOutputs.P_m_i_eff(i,:)    = obj.outputs.P_m_i_eff(i,:);  
                    obj.processedOutputs.P1_e_o(i)         = obj.outputs.P1_e_o(i);
                    obj.processedOutputs.P2_e_i(i)         = obj.outputs.P2_e_i(i); 
                    obj.processedOutputs.P_e_o_eff(i,:)    = obj.outputs.P_e_o_eff(i,:);
                    obj.processedOutputs.P_e_i_eff(i,:)    = obj.outputs.P_e_i_eff(i,:);
                    obj.processedOutputs.P_m_o(i)          = obj.outputs.P_m_o(i);
                    obj.processedOutputs.P_m_i(i)          = obj.outputs.P_m_i(i);
                    obj.processedOutputs.P_e_o(i)          = obj.outputs.P_e_o(i);
                    obj.processedOutputs.P_e_i(i)          = obj.outputs.P_e_i(i);
                    obj.processedOutputs.deltaL(i)         = obj.outputs.deltaL(i);
                    obj.processedOutputs.t1(i)             = obj.outputs.t1(i);
                    obj.processedOutputs.t2(i)             = obj.outputs.t2(i);
                    obj.processedOutputs.to_eff(i,:)       = obj.outputs.to_eff(i,:);
                    obj.processedOutputs.ti_eff(i,:)       = obj.outputs.ti_eff(i,:);
                    obj.processedOutputs.to(i)             = obj.outputs.to(i);
                    obj.processedOutputs.ti(i)             = obj.outputs.ti(i);
                    obj.processedOutputs.tCycle(i)         = obj.outputs.tCycle(i);
                    obj.processedOutputs.Ft(i,:)           = obj.outputs.Ft(i,:);
                    obj.processedOutputs.Ft_i(i,:)         = obj.outputs.Ft_i(i,:);
                    obj.processedOutputs.Fa(i,:)           = obj.outputs.Fa(i,:);
                    obj.processedOutputs.Fa_i(i,:)         = obj.outputs.Fa_i(i,:);
                    obj.processedOutputs.W(i,:)            = obj.outputs.W(i,:);
                    obj.processedOutputs.vo(i,:)           = obj.outputs.vk_r(i,:);
                    obj.processedOutputs.vk_tau(i,:)       = obj.outputs.vk_tau(i,:);
                    obj.processedOutputs.vi(i,:)           = obj.outputs.vk_r_i(i,:);
                    obj.processedOutputs.P_e_avg(i)        = obj.outputs.P_e_avg(i);
                    obj.processedOutputs.P_m_avg(i)        = obj.outputs.P_m_avg(i);
                    obj.processedOutputs.Rp(i,:)           = obj.outputs.Rp(i,:);
                    obj.processedOutputs.h_cycleStart(i)   = obj.outputs.h_cycleStart(i);
                    obj.processedOutputs.h_cycleAvg(i)     = obj.outputs.h_cycleAvg(i);
                    obj.processedOutputs.l_t_max(i)        = obj.outputs.l_t_max(i);
                    obj.processedOutputs.l_t_min(i)        = obj.outputs.l_t_min(i);
                    obj.processedOutputs.l_t_inCycle(i,:)  = obj.outputs.l_t_inCycle(i,:);
                    obj.processedOutputs.h_inCycle(i,:)    = obj.outputs.h_inCycle(i,:); 
                    obj.processedOutputs.va(i,:)           = obj.outputs.va(i,:);
                    obj.processedOutputs.beta(i)           = rad2deg(obj.outputs.beta(i));
                    obj.processedOutputs.rollAngle(i,:)    = rad2deg(-obj.outputs.rollAngle(i,:));
                    obj.processedOutputs.gamma(i)          = rad2deg(obj.outputs.gamma(i));
                    obj.processedOutputs.CL(i,:)           = obj.outputs.CL(i,:);
                    obj.processedOutputs.CD(i,:)           = obj.outputs.CD(i,:);
                    obj.processedOutputs.CD_k(i,:)         = obj.outputs.CD_k(i,:);
                    obj.processedOutputs.CD_k_i(i,:)       = obj.outputs.CD_k_i(i,:);
                    obj.processedOutputs.CD_t(i,:)         = obj.outputs.CD_t(i,:);
                    obj.processedOutputs.CL_i(i,:)         = obj.outputs.CL_i(i,:);
                    obj.processedOutputs.CD_i(i,:)         = obj.outputs.CD_i(i,:);
                    obj.processedOutputs.E(i,:)            = obj.outputs.E(i,:);
                    obj.processedOutputs.E_i(i,:)          = obj.outputs.E_i(i,:);
                    obj.processedOutputs.numOfPatt(i,:)    = obj.outputs.numOfPatt(i,:);
                    obj.processedOutputs.lambda(i,:)       = obj.outputs.lambda(i,:); 
                    obj.processedOutputs.f(i,:)            = obj.outputs.f(i,:); 
                    obj.processedOutputs.f_i(i,:)          = obj.outputs.f_i(i,:); 
                    obj.processedOutputs.tPatt(i,:)        = obj.outputs.tPatt(i,:);
                    obj.processedOutputs.dutyCycle(i)      = obj.processedOutputs.to(i)/obj.processedOutputs.tCycle(i);
                    obj.processedOutputs.zetaMech(i,:)     = obj.outputs.zetaMech(i,:);  
                    obj.processedOutputs.d_te              = obj.outputs.d_t;  
                    
                    % Cycle efficiency
                    obj.processedOutputs.cycleEff_elec(i)  = (obj.processedOutputs.P_e_avg(i)*obj.processedOutputs.tCycle(i))/...
                                                            (obj.processedOutputs.P_e_o(i)*obj.processedOutputs.to(i));
                    obj.processedOutputs.cycleEff_mech(i)  = (obj.processedOutputs.P_m_avg(i)*obj.processedOutputs.tCycle(i))/...
                                                            (obj.processedOutputs.P_m_o(i)*obj.processedOutputs.to(i));
                    
                    % Coefficient of power (Cp) as defined for HAWTs
                    obj.processedOutputs.sweptArea(i,:)     = pi.*((obj.outputs.Rp(i,:)+obj.inputs.span/2).^2 - (obj.outputs.Rp(i,:)-obj.inputs.span/2).^2);
                    obj.processedOutputs.Cp_m_o(i,:)        = obj.outputs.P_m_o(i)./(0.5.*obj.inputs.densityAir.*obj.processedOutputs.sweptArea(i,:).*obj.outputs.vw(i,:).^3);
                    obj.processedOutputs.Cp_e_avg(i,:)      = obj.outputs.P_e_avg(i)./(0.5.*obj.inputs.densityAir.*obj.processedOutputs.sweptArea(i,:).*obj.outputs.vw(i,:).^3);
                end
            end
            
            %% Cycle power representation for wind speeds in the operational range
            for vw_i = vw(vw>=obj.processedOutputs.cutIn&vw<=obj.processedOutputs.cutOut)
                obj.processedOutputs.cyclePowerRep(vw_i) = obj.createCyclePowerRep(vw_i); 
            end
        end

        function cyclePowerRep = createCyclePowerRep(obj, vw_i)
            cyclePowerRep.ws     = vw_i;
            cyclePowerRep.idx    = find(obj.inputs.windSpeedReference==vw_i);
            cyclePowerRep.t1     = round(obj.processedOutputs.t1(cyclePowerRep.idx),2);
            cyclePowerRep.to_eff = round(sum(obj.processedOutputs.to_eff(cyclePowerRep.idx,:),2),2);
            cyclePowerRep.t2     = round(obj.processedOutputs.t2(cyclePowerRep.idx),2);
            cyclePowerRep.ti_eff = round(sum(obj.processedOutputs.ti_eff(cyclePowerRep.idx,:),2),2);
            cyclePowerRep.to     = cyclePowerRep.t1 + cyclePowerRep.to_eff; %[s]
            cyclePowerRep.ti     = cyclePowerRep.t2 + cyclePowerRep.ti_eff; %[s]
            cyclePowerRep.tCycle = cyclePowerRep.to + cyclePowerRep.to; %[s]
            
            cyclePowerRep.t_inst   = cumsum([0 cyclePowerRep.t1 (cyclePowerRep.to_eff-cyclePowerRep.t1) cyclePowerRep.t1 ...
                                  cyclePowerRep.t2 (cyclePowerRep.ti_eff-cyclePowerRep.t2) cyclePowerRep.t2]);
            
            cyclePowerRep.P_e_inst = [0 obj.processedOutputs.P_e_o_eff(cyclePowerRep.idx,1) obj.processedOutputs.P_e_o_eff(cyclePowerRep.idx,end) 0 ...
                                  -obj.processedOutputs.P_e_i_eff(cyclePowerRep.idx,end) -obj.processedOutputs.P_e_i_eff(cyclePowerRep.idx,1) 0]./10^3;
            
            cyclePowerRep.P_m_inst = [0 obj.processedOutputs.P_m_o_eff(cyclePowerRep.idx,1) obj.processedOutputs.P_m_o_eff(cyclePowerRep.idx,end) 0 ...
                                  -obj.processedOutputs.P_m_i_eff(cyclePowerRep.idx,end) -obj.processedOutputs.P_m_i_eff(cyclePowerRep.idx,1) 0]./10^3;
        end

        function obj = initialise_outputs_to_zero(obj)
            % Initialize all obj.outputs fields with zeros
            num_ref = length(obj.inputs.windSpeedReference);
            num_elems = obj.inputs.numDeltaLelems;

            obj.outputs.halfRhoS = 0;
            obj.outputs.m_k = 0;
            
            % Single-column variables
            obj.outputs.l_t_min = zeros(1, num_ref);
            obj.outputs.pattStartGrClr = zeros(1, num_ref);
            obj.outputs.h_cycleStart = zeros(1, num_ref);
            obj.outputs.l_t_max = zeros(1, num_ref);
            obj.outputs.pattEndGrClr = zeros(1, num_ref);
            obj.outputs.l_t_avg = zeros(1, num_ref);
            obj.outputs.h_cycleAvg = zeros(1, num_ref);
            obj.outputs.h_cycleEnd = zeros(1, num_ref);
            obj.outputs.d_t = zeros(1, num_ref);
            obj.outputs.deltaLelems = zeros(1, num_ref);
            obj.outputs.elemDeltaL = zeros(1, num_ref);
            obj.outputs.t1 = zeros(1, num_ref);
            obj.outputs.P1_m_o = zeros(1, num_ref);
            obj.outputs.P1_e_o = zeros(1, num_ref);
            obj.outputs.P_m_o = zeros(1, num_ref);
            obj.outputs.P_e_o = zeros(1, num_ref);
            obj.outputs.t2 = zeros(1, num_ref);
            obj.outputs.P2_m_i = zeros(1, num_ref);
            obj.outputs.P2_e_i = zeros(1, num_ref);
            obj.outputs.P_m_i = zeros(1, num_ref);
            obj.outputs.P_e_i = zeros(1, num_ref);
            obj.outputs.tCycle = zeros(1, num_ref);
            obj.outputs.tPatt = zeros(num_ref, num_elems);
            obj.outputs.numOfPatt = zeros(num_ref, num_elems);
            obj.outputs.P_e_avg = zeros(1, num_ref);
            obj.outputs.P_m_avg = zeros(1, num_ref);
            obj.outputs.deltaL = zeros(1, num_ref);
            obj.outputs.beta = zeros(1, num_ref);
            obj.outputs.gamma = zeros(1, num_ref);
            obj.outputs.Rp_start = zeros(1, num_ref);
            obj.outputs.etaGen_i = zeros(1, num_ref);
            obj.outputs.to = zeros(1, num_ref);
            obj.outputs.ti = zeros(1, num_ref);
            
            % Variables with multiple columns
            obj.outputs.l_t_inCycle = zeros(num_ref, num_elems);
            obj.outputs.pattGrClr = zeros(num_ref, num_elems);
            obj.outputs.m_t = zeros(num_ref, num_elems);
            obj.outputs.m_eff = zeros(num_ref, num_elems);
            obj.outputs.theta = zeros(num_ref, num_elems);
            obj.outputs.phi = zeros(num_ref, num_elems);
            obj.outputs.chi = zeros(num_ref, num_elems);
            obj.outputs.CD_k = zeros(num_ref, num_elems);
            obj.outputs.CD_t = zeros(num_ref, num_elems);
            obj.outputs.CD = zeros(num_ref, num_elems);
            obj.outputs.E = zeros(num_ref, num_elems);
            obj.outputs.h_inCycle = zeros(num_ref, num_elems);
            obj.outputs.Rp = zeros(num_ref, num_elems);
            obj.outputs.vw = zeros(num_ref, num_elems);
            obj.outputs.vw_r = zeros(num_ref, num_elems);
            obj.outputs.vw_theta = zeros(num_ref, num_elems);
            obj.outputs.vw_phi = zeros(num_ref, num_elems);
            obj.outputs.f = zeros(num_ref, num_elems);
            obj.outputs.va = zeros(num_ref, num_elems);
            obj.outputs.lambda = zeros(num_ref, num_elems);
            obj.outputs.vk_tau = zeros(num_ref, num_elems);
            obj.outputs.vk_theta = zeros(num_ref, num_elems);
            obj.outputs.vk_phi = zeros(num_ref, num_elems);
            obj.outputs.W = zeros(num_ref, num_elems);
            obj.outputs.Fg_r = zeros(num_ref, num_elems);
            obj.outputs.Fg_theta = zeros(num_ref, num_elems);
            obj.outputs.Fg_phi = zeros(num_ref, num_elems);
            obj.outputs.Fa_theta = zeros(num_ref, num_elems);
            obj.outputs.Fa_phi = zeros(num_ref, num_elems);
            obj.outputs.Fa_r = zeros(num_ref, num_elems);
            obj.outputs.Fa = zeros(num_ref, num_elems);
            obj.outputs.va_r = zeros(num_ref, num_elems);
            obj.outputs.va_theta = zeros(num_ref, num_elems);
            obj.outputs.va_phi = zeros(num_ref, num_elems);
            obj.outputs.F_dot_v = zeros(num_ref, num_elems);
            obj.outputs.D_r = zeros(num_ref, num_elems);
            obj.outputs.D_theta = zeros(num_ref, num_elems);
            obj.outputs.D_phi = zeros(num_ref, num_elems);
            obj.outputs.L_r = zeros(num_ref, num_elems);
            obj.outputs.L_theta = zeros(num_ref, num_elems);
            obj.outputs.L_phi = zeros(num_ref, num_elems);
            obj.outputs.D = zeros(num_ref, num_elems);
            obj.outputs.L = zeros(num_ref, num_elems);
            obj.outputs.Ft = zeros(num_ref, num_elems);
            obj.outputs.E_result = zeros(num_ref, num_elems);
            obj.outputs.kByE = zeros(num_ref, num_elems);
            obj.outputs.P_m_o_eff = zeros(num_ref, num_elems);
            obj.outputs.zetaMech = zeros(num_ref, num_elems);
            obj.outputs.etaGen_o = zeros(num_ref, num_elems);
            obj.outputs.P_e_o_eff = zeros(num_ref, num_elems);
            obj.outputs.theta_i = zeros(num_ref, num_elems);
            obj.outputs.phi_i = zeros(num_ref, num_elems);
            obj.outputs.vw_i = zeros(num_ref, num_elems);
            obj.outputs.vw_r_i = zeros(num_ref, num_elems);
            obj.outputs.vw_theta_i = zeros(num_ref, num_elems);
            obj.outputs.vw_phi_i = zeros(num_ref, num_elems);
            obj.outputs.f_i = zeros(num_ref, num_elems);
            obj.outputs.va_r_i = zeros(num_ref, num_elems);
            obj.outputs.va_theta_i = zeros(num_ref, num_elems);
            obj.outputs.va_phi_i = zeros(num_ref, num_elems);
            obj.outputs.va_i = zeros(num_ref, num_elems);
            obj.outputs.CD_k_i = zeros(num_ref, num_elems);
            obj.outputs.CD_i = zeros(num_ref, num_elems);
            obj.outputs.E_i = zeros(num_ref, num_elems);
            obj.outputs.Fa_i = zeros(num_ref, num_elems);
            obj.outputs.Fg_r_i = zeros(num_ref, num_elems);
            obj.outputs.Fg_theta_i = zeros(num_ref, num_elems);
            obj.outputs.Fg_phi_i = zeros(num_ref, num_elems);
            obj.outputs.Fa_theta_i = zeros(num_ref, num_elems);
            obj.outputs.Fa_phi_i = zeros(num_ref, num_elems);
            obj.outputs.Fa_r_i = zeros(num_ref, num_elems);
            obj.outputs.F_dot_v_i = zeros(num_ref, num_elems);
            obj.outputs.D_r_i = zeros(num_ref, num_elems);
            obj.outputs.D_theta_i = zeros(num_ref, num_elems);
            obj.outputs.D_phi_i = zeros(num_ref, num_elems);
            obj.outputs.L_r_i = zeros(num_ref, num_elems);
            obj.outputs.L_theta_i = zeros(num_ref, num_elems);
            obj.outputs.L_phi_i = zeros(num_ref, num_elems);
            obj.outputs.D_i = zeros(num_ref, num_elems);
            obj.outputs.L_i = zeros(num_ref, num_elems);
            obj.outputs.Ft_i = zeros(num_ref, num_elems);
            obj.outputs.E_result_i = zeros(num_ref, num_elems);
            obj.outputs.P_m_i_eff = zeros(num_ref, num_elems);
            obj.outputs.P_e_i_eff = zeros(num_ref, num_elems);
            obj.outputs.vk_r_i = zeros(num_ref, num_elems);
            obj.outputs.CL_i = zeros(num_ref, num_elems);
            obj.outputs.vk_r = zeros(num_ref, num_elems);
            obj.outputs.kRatio = zeros(num_ref, num_elems);
            obj.outputs.CL = zeros(num_ref, num_elems);
            obj.outputs.ti_eff = zeros(num_ref, num_elems);
            obj.outputs.to_eff = zeros(num_ref, num_elems);
            obj.outputs.rollAngle = zeros(num_ref, num_elems);
        end
    end
end