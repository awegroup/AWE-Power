classdef (Abstract) KiteQSMsimulationBase < handle
    % KiteQSMsimulationBase is an abstract class that contains
    % all properties and methods common to all kite simulation types.
    
    properties
        inputs
        outputs
        optimDetails
        processedOutputs
        kiteMass
        logger
        validationFcn
    end
    
    methods
        function obj = KiteQSMsimulationBase(inputs, loggerIdentifier)
            if nargin<2
                obj.logger = mlog.Logger('loggerQSM');
            else
                obj.logger = mlog.Logger(loggerIdentifier);
            end
            
            cls = class(obj);   % Gets the actual (sub)class name.
            
            obj.set_inputvalidation_functions();

            p = inputParser; p.KeepUnmatched = true;
            % Constructor: common initialization
            get_required_params = str2func([cls '.get_required_params']);
            requiredParams = get_required_params();
            obj.add_required_parameters(p, inputs, requiredParams);
            
            get_parameters_withdefaults = str2func([cls '.get_parameters_withdefaults']);
            parametersWithDefaults = get_parameters_withdefaults(inputs, obj.validationFcn);
            obj.add_parameters_withdefaults(p, parametersWithDefaults);

            obj.add_windprofile_withdefaults(p, inputs);

            parse(p, inputs);
            obj.add_inputparserresults_to_property(p);
            
            obj.add_mass_parameter(obj.inputs);

            p3 = inputParser; p3.KeepUnmatched = true;
            obj.set_optimization_defaults_and_parse(inputs, p3);
            obj.add_inputparserresults_to_property(p3);
            
            % Log default values if used
            obj.write_inputparserdefaults_to_logger(obj.logger, p, obj.inputs, inputs, obj.kiteMass)
            obj.write_inputparserdefaults_to_logger(obj.logger, p3, obj.inputs, inputs)
        end
        
        %% Common methods (matching in both subclasses)
        function obj = add_required_parameters(obj, p, inputs, requiredParams)
            for i = 1:numel(requiredParams)
                param = requiredParams{i};
                if isfield(inputs, param)
                    p.addParameter(param, inputs.(param)); % Add parameter with value from input structure
                else
                    try
                        error('KiteQSMsimulation:constructor', ['Input parameter ', param, ' is required.']);
                    catch err
                        obj.logger.write(err);
                    end
                end
            end
        end
        %
        function obj = add_parameters_withdefaults(obj, p, parametersWithDefaults)
            for i = 1:size(parametersWithDefaults, 1)
                param = parametersWithDefaults{i, 1};
                defaultVal = parametersWithDefaults{i, 2};
                validationFcn_i = parametersWithDefaults{i, 3};

                addParameter(p, param, defaultVal, validationFcn_i);
            end
        end
        
        function obj = set_inputvalidation_functions(obj)
            obj.validationFcn.FloatPositive = @(x) assert(isnumeric(x) && isscalar(x) ...
                && (x > 0), 'Value must be positive, scalar, and numeric.');
            obj.validationFcn.FloatPositiveOrZero = @(x) assert(isnumeric(x) && isscalar(x) ...
                && (x >= 0), 'Value must be positive, scalar, and numeric.');
            obj.validationFcn.Float = @(x) assert(isnumeric(x) && isscalar(x), ...
                'Value must be numeric.');
            obj.validationFcn.Int = @(x) assert(isnumeric(x) && isscalar(x) ...
                && (x > 0) && mod(x,1)==0, 'Value must be positive, integer.');
            obj.validationFcn.Bool = @(x) assert(islogical(x) || x==0 ...
                || x==1, 'Value must be boolean or 0/1.');
            obj.validationFcn.String = @(x) assert(ischar(x) || isstring(x), ...
                'Value must be a string or character vector.');
            obj.validationFcn.ArrayFloatPositive = @(x) assert(isnumeric(x) ...
                && all(x > 0), 'Value must be positive, scalar, and numeric.');
            obj.validationFcn.etaGen_param= @(x) assert(isnumeric(x) && numel(x) == 4, ...
                'Parameter must be a numeric array with 4 elements.');
        end
        
        function obj = add_mass_parameter(obj, inputs)
            if isfield(inputs,'doMassOverride')
                if inputs.doMassOverride
                    if isfield(inputs, 'kiteMass')
                        try
                            obj.validationFcn.FloatPositive(inputs.kiteMass);
                        catch err
                            obj.logger.error('error in value of kiteMass');
                            obj.logger.write(err);
                        end
                        obj.kiteMass = inputs.kiteMass;
                    else
                        try
                            error('KiteQSMsimulation:constructor', 'Input parameter kiteMass is required when doMassOverride is true.');
                        catch err
                            obj.logger.write(err);
                        end
                    end
                else
                    obj.logger.warning(['Mass is calculated using Vincent ' ...
                        'Bonnin his simple mass model developed at Ampyx ' ...
                        'Power. Based on AP3 data and projected data for ' ...
                        'larger systems (AP4-AP5)']);
                    a1     = 0.002415;
                    a2     = 0.0090239;
                    b1     = 0.17025;
                    b2     = 3.2493;
                    k1     = 5;
                    c1     = 0.46608;
                    d1     = 0.65962;
                    k2     = 1.1935;
                    AR_ref = 12;
                    a_massCalc = a1*(inputs.forceTether_max/1000/inputs.areaWing) + a2;
                    b_massCalc = b1*(inputs.forceTether_max/1000/inputs.areaWing) + b2;
                    obj.kiteMass = 10*(a_massCalc*inputs.areaWing^2 +b_massCalc*inputs.areaWing-k1)...
                        *(c1*(inputs.aspectRatio/AR_ref)^2-d1*(inputs.aspectRatio/AR_ref)+k2); 
                end
            end
        end
        
        function obj = add_windprofile_withdefaults(obj, p, inputs)
            if isfield(inputs, 'doUseReferenceWindProfile')
                if inputs.doUseReferenceWindProfile
                    % Sample vertical wind profile datasets (i.e. for inputs.doUseReferenceWindProfile = true)
                    % Profiles from: https://doi.org/10.1016/j.renene.2022.06.094. Wind speeds normalized with value at 100m. 
                    windData = yaml.ReadYaml('winddata_Lidar_Avg.yaml',false,true);

                    h_windDataset = windData.h;
                    if isfield(inputs, 'referenceLocation')
                        if strcmpi(inputs.referenceLocation,'ijmuiden')
                            v_windDataset = windData.vw_norm_metmast_ijmuiden;
                        elseif strcmpi(inputs.referenceLocation,'cabauw')
                            v_windDataset = windData.vw_norm_metmast_cabauw;
                        else
                            error('KiteQSMsimulation:runSimulation Wrong reference location provided')
                        end
                    else
                        v_windDataset = windData.vw_norm_metmast_ijmuiden;
                    end
                    
                    vw_to_normalise = interp1(h_windDataset,v_windDataset,inputs.heightWindReference);
                    addParameter(p, 'windProfile_h', h_windDataset, obj.validationFcn.ArrayFloatPositive);
                    % normalise to correct reference height
                    addParameter(p, 'windProfile_vw', v_windDataset./vw_to_normalise, obj.validationFcn.ArrayFloatPositive);
                else
                    addParameter(p, 'windShearExp', 0.143, obj.validationFcn.FloatPositive); %[-] % 0.143 over land, 0.11 over sea
                end
            else
                addParameter(p, 'windShearExp', 0.143, obj.validationFcn.FloatPositive); %[-] % 0.143 over land, 0.11 over sea
            end
        end
        
        function obj = add_inputparserresults_to_property(obj, p)
            f = fieldnames(p.Results);
            for i = 1:length(f)
                obj.inputs = setfield(obj.inputs, f{i}, p.Results.(f{i}));
            end
        end

        function obj = plot_all_results(obj)
            % inputs = obj.inputs;
            processedOutputs = obj.processedOutputs; %#ok<*PROP>

            %% Plot settings
            x_axis_limit = [0 processedOutputs.vw_100m_operRange(end)];
            
            colorOrderList = [ % 0.25, 0.25, 0.25
                0 0.4470 0.7410
                0.8500 0.3250 0.0980
                0.4660, 0.6740, 0.1880
                0.9290, 0.6940, 0.1250
                0.4940, 0.1840, 0.5560
                0.6350 0.0780 0.1840
                0.3010 0.7450 0.9330];

            % Speeds
            obj.plot_speed_results(colorOrderList, x_axis_limit);
        
            % CL
            obj.plot_cl_results(colorOrderList, x_axis_limit);
        
            % CD
            obj.plot_cd_results(colorOrderList, x_axis_limit);
            
            % Lengths
            obj.plot_tether_lengths(colorOrderList, x_axis_limit);
        
            % Roll angle, avg patt elevation
            obj.plot_roll_avg_pattern(colorOrderList, x_axis_limit);

            obj.plot_forces(colorOrderList, x_axis_limit);
        
            % Reel-out power
            obj.plot_power(colorOrderList, x_axis_limit);

            % Power curve comparison plot
            % Loyd
            obj.plot_power_curve_comparisonLoyd(colorOrderList, x_axis_limit);

            % Wind profile
            obj.plot_wind_profile(colorOrderList);

        end
        %  
        function obj=plot_speed_results(obj, colorOrderList, x_axis_limit)
            vw = obj.processedOutputs.vw_reference;
            if nargin<3;x_axis_limit  = [min(vw), max(vw)];end
            yAxesOptions = {struct('xData', vw, ...
                                    'yData', {{'lambda'}}, ...
                                    'yFcn', @(x) mean(x,2) ,...
                                    'yLabel', 'Tangential kite velocity ratio (-)', ...
                                    'xLabel', 'Wind speed at 100m height (m/s)', ...
                                    'xLim', x_axis_limit);};
            if nargin>1;plotOptions.colors = colorOrderList;end
            if nargin>1
                obj.plotResults(obj.processedOutputs, yAxesOptions,[],plotOptions);
            else
                obj.plotResults(obj.processedOutputs, yAxesOptions,[],struct());
            end
        end
        %
        function obj=plot_cl_results(obj, colorOrderList, x_axis_limit)
            vw = obj.processedOutputs.vw_reference;
            if nargin<3;x_axis_limit  = [min(vw), max(vw)];end
            if isfield(obj.inputs,'numTurbines')
                yAxesOptions = {struct('xData', vw, ...
                                        'yData', {{'CL'}}, ...
                                        'yFcn', @(x) mean(x,2) ,...
                                        'yLabel', 'Lift coefficient (-)', ...
                                        'xLabel', 'Wind speed at 100m height (m/s)', ...
                                        'xLim', x_axis_limit);};
            else
                yAxesOptions = {struct('xData', vw, ...
                                        'yData', {{'CL','CL_i'}}, ...
                                        'yFcn', @(x) mean(x,2) ,...
                                        'yLabel', 'Lift coefficient (-)', ...
                                        'xLabel', 'Wind speed at 100m height (m/s)', ...
                                        'Legend', {{'C_{L,o}','C_{L,i}'}}, ...
                                        'xLim', x_axis_limit);};
            end
            if nargin>1;plotOptions.colors = colorOrderList;end
            if nargin>1
                obj.plotResults(obj.processedOutputs, yAxesOptions,[],plotOptions);
            else
                obj.plotResults(obj.processedOutputs, yAxesOptions,[],struct());
            end
        end
        %
        function obj = plot_cd_results(obj, colorOrderList, x_axis_limit)
            vw = obj.processedOutputs.vw_reference;
            if nargin<3;x_axis_limit  = [min(vw), max(vw)];end
            if isfield(obj.inputs,'numTurbines')
                yAxesOptions = {struct('xData', vw, ...
                                    'yData', {{'CD', 'CD_k', 'CD_t', 'CT'}}, ...
                                    'yFcn', @(x) mean(x,2), ...
                                    'yLabel', 'Drag coefficient (-)', ...
                                    'xLabel', 'Wind speed at 100m height (m/s)', ...
                                    'Legend', {{'C_{D}','C_{D,k}','C_{D,t}','C_{T}'}}, ...
                                    'xLim', x_axis_limit);};
            else
                yAxesOptions = {struct('xData', vw, ...
                                    'yData', {{'CD', 'CD_i','CD_k', 'CD_k_i', 'CD_t'}}, ...
                                    'yFcn', @(x) mean(x,2), ...
                                    'yLabel', 'Drag coefficient (-)', ...
                                    'xLabel', 'Wind speed at 100m height (m/s)', ...
                                    'Legend', {{'C_{D}','C{D_i}','C_{D,k}','C_{D,k,i}','C_{D,t}'}}, ...
                                    'xLim', x_axis_limit);};
            end
            if nargin>1;plotOptions.colors = colorOrderList;end
            if nargin>1
                obj.plotResults(obj.processedOutputs, yAxesOptions,[],plotOptions);
            else
                obj.plotResults(obj.processedOutputs, yAxesOptions,[],struct());
            end
        end
        %
        function obj = plot_tether_lengths(obj, colorOrderList, x_axis_limit)
            vw = obj.processedOutputs.vw_reference;
            if nargin<3;x_axis_limit  = [min(vw), max(vw)];end
            if isfield(obj.inputs,'numTurbines')
                yAxesOptions = {struct('xData', vw, ...
                                        'yData', {{'h_inCycle', 'Rp', 'l_t_inCycle'}}, ...
                                        'yFcn', @(x) mean(x,2), ...
                                        'yLabel', '(m)', ...
                                        'xLabel', 'Wind speed at 100m height (m/s)', ...
                                        'Legend', {{'Pattern height','Pattern radius (circle)','Tether length'}}, ...
                                        'xLim', x_axis_limit);};
            else
                yAxesOptions = {struct('xData', vw, ...
                                    'yData', {{'h_cycleAvg', 'Rp', 'deltaL', 'l_t_max', 'l_t_min'}}, ...
                                    'yFcn', @(x) mean(x,2), ...
                                    'yLabel', '(m)', ...
                                    'xLabel', 'Wind speed at 100m height (m/s)', ...
                                    'Legend', {{'Pattern height','Pattern radius (circle)', 'Reel-out length','Min tether length','Max tether length'}}, ...
                                    'xLim', x_axis_limit);};
            end
            if nargin>1;plotOptions.colors = colorOrderList;end
            if nargin>1
                obj.plotResults(obj.processedOutputs, yAxesOptions,[],plotOptions);
            else
                obj.plotResults(obj.processedOutputs, yAxesOptions,[],struct());
            end
        end
        %
        function obj = plot_roll_avg_pattern(obj, colorOrderList, x_axis_limit)
            vw = obj.processedOutputs.vw_reference;
            if nargin<3;x_axis_limit  = [min(vw), max(vw)];end
            yAxesOptions = {struct('xData', vw, ...
                                    'yData', {{'rollAngle', 'beta', 'gamma'}}, ...
                                    'yFcn', @(x) mean(x,2), ...
                                    'yLabel', '(deg)', ...
                                    'xLabel', 'Wind speed at 100m height (m/s)', ...
                                    'Legend', {{'Roll angle','Elevation angle','Pattern opening angle'}}, ...
                                    'xLim', x_axis_limit);};
            if nargin>1;plotOptions.colors = colorOrderList;end
            if nargin>1
                obj.plotResults(obj.processedOutputs, yAxesOptions,[],plotOptions);
            else
                obj.plotResults(obj.processedOutputs, yAxesOptions,[],struct());
            end
        end
        %
        function obj = plot_forces(obj, colorOrderList, x_axis_limit)
            vw = obj.processedOutputs.vw_reference;
            if nargin<3;x_axis_limit  = [min(vw), max(vw)];end
            if isfield(obj.inputs,'numTurbines')
                yAxesOptions = {struct('xData', vw, ...
                                        'yData', {{'W', 'Fa', 'Ft', 'Fthrust'}}, ...
                                        'yFcn', @(x) mean(x,2)./1e3, ...
                                        'yLabel', 'Force (kN)', ...
                                        'xLabel', 'Wind speed at 100m height (m/s)', ...
                                        'Legend', {{'Weight','Aero resultant','Tether', 'Thrust'}}, ...
                                        'xLim', x_axis_limit)};
            else
                yAxesOptions = {struct('xData', vw, ...
                                        'yData', {{'W', 'Fa', 'Ft','Fa_i','Ft_i'}}, ...
                                        'yFcn', @(x) mean(x,2)./1e3, ...
                                        'yLabel', 'Force (kN)', ...
                                        'xLabel', 'Wind speed at 100m height (m/s)', ...
                                        'Legend', {{'Weight','Aero reel-out','Tether reel-out','Aero reel-in','Tether reel-in'}}, ...
                                        'xLim', x_axis_limit)};
            end
            if nargin>1;plotOptions.colors = colorOrderList;end
            if nargin>1
                obj.plotResults(obj.processedOutputs, yAxesOptions,[],plotOptions);
            else
                obj.plotResults(obj.processedOutputs, yAxesOptions,[],struct());
            end
        end
        %
        function obj = plot_power(obj, colorOrderList, x_axis_limit)
            vw = obj.processedOutputs.vw_reference;
            if nargin<3;x_axis_limit  = [min(vw), max(vw)];end
            if isfield(obj.inputs,'numTurbines')
                ydata = {'P_m_eff', 'P_e_eff'};
            else
                ydata = {'P_m_avg', 'P_e_avg'};
            end
            yAxesOptions = {struct('xData', vw, ...
                                    'yData', {ydata}, ...
                                    'yFcn', @(x) x./10^3, ...
                                    'yLabel', 'Power (kW)', ...
                                    'XLabel', 'Wind speed at 100m height (m/s)', ...
                                    'Legend', {{'Mechanical','Electrical'}}, ...
                                    'xLim', x_axis_limit)};
            if nargin>1;plotOptions.colors = colorOrderList;end
            if nargin>1
                obj.plotResults(obj.processedOutputs, yAxesOptions,[],plotOptions);
            else
                obj.plotResults(obj.processedOutputs, yAxesOptions,[],struct());
            end
        end
        %
        function obj = plot_power_curve_comparisonLoyd(obj, colorOrderList, x_axis_limit)
            vw = obj.processedOutputs.vw_reference;
            if nargin<3;x_axis_limit  = [min(vw), max(vw)];end
            % Calculate Loyd power
            P_Loyd = zeros(1,length(vw));
            CL_loyd = obj.inputs.Cl_maxAirfoil * obj.inputs.Cl_eff_F;
            CD_loyd_noTether = obj.inputs.Cd0 + (CL_loyd-obj.inputs.Cl0_airfoil)^2 / (pi * obj.inputs.aspectRatio * obj.inputs.e) + obj.inputs.Cd_additional;
            for i = 1:length(vw)
                CD_loyd = CD_loyd_noTether+obj.outputs.CD_t(i);
                P_Loyd(i) = (4/27) * (CL_loyd^3 / CD_loyd^2) * ...
                    (1/2) * obj.inputs.densityAir * obj.inputs.areaWing * (vw(i)^3);
                if P_Loyd(i) > obj.inputs.P_ratedElec
                    P_Loyd(i) = obj.inputs.P_ratedElec;
                end 
            end
            if isfield(obj.inputs,'numTurbines')
                ydata = {P_Loyd./1e3, 'P_e_eff'};
            else
                ydata = {P_Loyd./1e3, 'P_e_avg'};
            end
            yAxesOptions = {struct('xData', {{vw, vw}}, ...
                                    'yData', {ydata}, ...
                                    'yFcn', @(x) x./1e3, ...
                                    'yLabel', 'Power (kW)', ...
                                    'XLabel', 'Wind speed at 100m height (m/s)', ...
                                    'Legend', {{'Ideal Loyd','QSM'}}, ...
                                    'xLim', x_axis_limit)};
            if nargin>1;plotOptions.colors = colorOrderList;end
            if nargin>1
                obj.plotResults(obj.processedOutputs, yAxesOptions,[],plotOptions);
            else
                obj.plotResults(obj.processedOutputs, yAxesOptions,[],struct());
            end
        end
        %
        function obj = plot_wind_profile(obj, colorOrderList)
            if ~obj.inputs.doUseReferenceWindProfile
                Vref = 10; % m/s
                z = 10:10:600; % m
                V = Vref * (z/obj.inputs.heightWindReference).^obj.inputs.windShearExp/Vref;
                % MegAWES Onshore location Cabauw. Wind speeds normalized with value at 100m. Profiles from: https://doi.org/10.1016/j.renene.2022.06.094
                z_MegAWES = [10,20,40,60,80,100,120,140,150,160,180,200,220,250,300,500,600];
                V_MegAWES_Cabauw = [0.541219682843206,0.607355091566827,0.768630154201962,0.868484406441142,0.941395360902529,1,1.04810058627160,1.08638854381156,1.10277338731106,1.11715868927737,1.14412258234309,1.16573551308321,1.18394938534465,1.20653423381438,1.23266397972046,1.26662287360302,1.26414483994687];
                V_MegAWES_Ijmuiden = [0.847612611633547,0.870603040595613,0.927240267828556,0.959346286990695,0.982291573490674,1,1.01377720773809,1.02356771954493,1.02766760602000,1.03079423355205,1.03659625208888,1.04025827758100,1.04284618416620,1.04496440015282,1.04461712713371,1.02473617783789,1.01076976884552];
        
                % Prepare yAxesOptions for plotResults
                yAxesOptions = {
                    struct('xData', {{V,V_MegAWES_Cabauw,V_MegAWES_Ijmuiden}}, ...
                    'yData', {{z,z_MegAWES, z_MegAWES}}, ...
                    'yLabel', 'Height (m)', ...
                    'Legend', {{'α = 0.143','Cabauw,NL','Ijmuiden,NL'}}, ...
                    'XLabel', 'Wind Speed (-)', ...
                    'xLim', [0.5, 1.5], ...
                    'title', 'Wind profile'),...
                };
            else
                % Continuous function from vertical wind profile
                % Define the heights at which you want to interpolate
                heights = (0:10:1000); % Adjust the desired height values as needed
                % Interpolate at the specified heights
                V_interpolated = interp1(obj.inputs.windProfile_h, obj.inputs.windProfile_vw, heights, 'linear', 'extrap');
                
                % Prepare yAxesOptions for plotResults
                yAxesOptions = {
                    struct('xData', {{obj.inputs.windProfile_vw, V_interpolated}}, ...
                    'yData', {{obj.inputs.windProfile_h, heights}}, ...
                    'yLabel', 'Height (m)', ...
                    'Legend', {{'Data','Interpolated'}}, ...
                    'XLabel', 'Wind Speed (-)', ...
                    'xLim', [0.5, 1.5], ...
                    'title', 'Wind profile'),...
                };
            end
        
            if nargin>1;plotOptions.colors = colorOrderList;end
            plotOptions.LineStyle = '-';
            plotOptions.markerList = repmat({'none'},1,3);
            obj.plotResults([], yAxesOptions, [], plotOptions);
          end
    end
    
    methods (Static)
        function write_inputparserdefaults_to_logger(logObj, p, inputsProcessed, inputs, kiteMass)
            if ~isempty(p.UsingDefaults)
                logObj.warning('Some QSM inputs are set to defaults; see details below.');
                for i = 1:numel(p.UsingDefaults)
                    if ~ismember(p.UsingDefaults{i},{'aspectRatio','windProfile_h','windProfile_vw','kiteMass'})
                        logObj.info('Parameter %s set to default: %s', p.UsingDefaults{i}, mat2str(inputsProcessed.(p.UsingDefaults{i})));
                    elseif strcmp(p.UsingDefaults{i},'aspectRatio')
                        logObj.info('Parameter %s is calculated to be: %s', p.UsingDefaults{i}, mat2str(inputsProcessed.(p.UsingDefaults{i})));
                    elseif strcmp(p.UsingDefaults{i},'kiteMass')
                        if ~inputsProcessed.doMassOverride
                            logObj.info('Parameter %s is calculated to be: %s', p.UsingDefaults{i}, mat2str(kiteMass));
                        end
                    elseif ismember(p.UsingDefaults{i},{'windProfile_h','windProfile_vw'})
                        if inputsProcessed.doUseReferenceWindProfile
                            if ~isfield(inputs, 'referenceLocation')
                                logObj.info('Parameter %s set to default: %s', p.UsingDefaults{i}, mat2str(inputsProcessed.(p.UsingDefaults{i})));
                            end
                        end
                    end
                end
            end
        end
        function processedOutputs = copy_outputs_to_processedoutputs(processedOutputs, outputs, variableNames, idx)
            for i = 1:numel(variableNames)
                name = variableNames{i};
                if any(size(outputs.(name))==1)
                    processedOutputs.(name)(idx) = outputs.(name)(idx);
                else
                    processedOutputs.(name)(idx,:) = outputs.(name)(idx,:);
                end
            end
        end
        
        function processedOutputs = copy_zeros_to_processedoutputs(processedOutputs, outputs, variableNames, idx)
            for i = 1:numel(variableNames)
                name = variableNames{i};
                if any(size(outputs.(name))==1)
                    processedOutputs.(name)(idx) = 0;
                else
                    processedOutputs.(name)(idx,:) = zeros(size(outputs.(name)(idx,:)));
                end
            end
        end

        function fig = plotResults(processedOutputs, yAxesOptions, subplotLayout, options)
    
            p = inputParser; p.KeepUnmatched = true;
            addParameter(p, 'markerList', {'o','^','s','v','*','+','x','pentagram','diamond'});
            addParameter(p, 'LineWidth', 1);
            addParameter(p, 'MarkerSize', 6);
            addParameter(p, 'LineStyle', ':');
        
            if nargin<4
                parse(p,struct())
            else
                parse(p,options)
            end
            lineOptions = p.Results;
        
            if nargin < 3
                subplotLayout = [];
            end
        
            if isempty(subplotLayout)
                numSubplots = 1;
            else
                numSubplots = prod(subplotLayout);
            end
        
            fig=figure; hold on; grid on; box on
        
            if isfield(options,'colors')
                colororder(options.colors);
            end
            
            for i = 1:numSubplots
                if numSubplots==1
                    subplot(1,1,1);
                else
                    subplot(subplotLayout(1), subplotLayout(2), i);
                end
                plotSinglePlot(processedOutputs, yAxesOptions{i}, lineOptions)
            end
        
            if numSubplots>1
                han=axes(fig,'visible','off'); 
                han.Title.Visible='on'; han.XLabel.Visible='on'; han.YLabel.Visible='off';
                if isfield(options,'yLabel')
                    ylabel(options.yLabel);
                end
                if isfield(options,'xLabel')
                    xlabel(options.xLabel);
                end
                if isfield(options,'title')
                    title(options.title);
                end
            end
        
            function plotSinglePlot(processedOutputs, subplotOptions, lineOptions)
                % Input:
                % processedOutputs: Structure containing processed data
                % subplotOptions: Structure specifying options for a single subplot
                
                hold on; grid on; box on
                % Plot each data series
                for i = 1:size(subplotOptions,1)
                    % Determine yyaxis position
                    if i==1 && size(subplotOptions,1)==2
                        yyaxis left
                    elseif i==2
                        yyaxis right
                    end
                    ylabel(subplotOptions(i).yLabel);
                    % Plot data
                    markerInd = 1;
                    if isfield(subplotOptions(i),'xData')
                        for j = 1:numel(subplotOptions(i).yData)
                            if iscell(subplotOptions(i).xData)
                                if isnumeric(subplotOptions(i).xData{j})
                                    x = subplotOptions(i).xData{j};
                                else
                                    x = processedOutputs.(subplotOptions(i).xData{j});
                                end
                            else
                                if isnumeric(subplotOptions(i).xData)
                                    x = subplotOptions(i).xData;
                                else
                                    x = processedOutputs.(subplotOptions(i).xData);
                                end
                            end
                            if isnumeric(subplotOptions(i).yData{j})
                                y = subplotOptions(i).yData{j};
                            else
                                if isfield(subplotOptions(i),'yFcn')
                                    if size(processedOutputs.(subplotOptions(i).yData{j}),1) == 1
                                        y = subplotOptions(i).yFcn(processedOutputs.(subplotOptions(i).yData{j})');
                                    else
                                        y = subplotOptions(i).yFcn(processedOutputs.(subplotOptions(i).yData{j}));
                                    end
                                else
                                    y = processedOutputs.(subplotOptions(i).yData{j});
                                end
                            end
                            plot(x, y, 'LineStyle', lineOptions.LineStyle, 'Marker', lineOptions.markerList{markerInd}, ...
                                'LineWidth', lineOptions.LineWidth, 'MarkerSize', lineOptions.MarkerSize);
                            markerInd = markerInd+1;
                        end
                    else
                        for j = 1:numel(subplotOptions(i).yData)
                            if isfield(subplotOptions(i),'windIndex')
                                ws = subplotOptions(i).windIndex;
                                if isnumeric(subplotOptions(i).yData{j})
                                    y = subplotOptions(i).yData{j};
                                else
                                    if isfield(subplotOptions(i),'yFcn')
                                        y = subplotOptions(i).yFcn(processedOutputs.(subplotOptions(i).yData{j})(ws,:));
                                    else
                                        y = processedOutputs.(subplotOptions(i).yData{j})(ws,:);
                                    end
                                end
                            else
                                if isnumeric(subplotOptions(i).yData{j})
                                    y = subplotOptions(i).yData{j};
                                else
                                    if isfield(subplotOptions(i),'yFcn')
                                        y = subplotOptions(i).yFcn(processedOutputs.(subplotOptions(i).yData{j}));
                                    else
                                        y = processedOutputs.(subplotOptions(i).yData{j});
                                    end
                                end
                            end
                            plot(y, 'LineStyle', lineOptions.LineStyle, 'Marker', lineOptions.markerList{markerInd}, ...
                                'LineWidth', lineOptions.LineWidth, 'MarkerSize', lineOptions.MarkerSize);
                            markerInd = markerInd+1;
                        end
                    end  
                end
                
                % Set x-axis label and legend
                if isfield(subplotOptions(i),'xLabel')
                    xlabel(subplotOptions(i).xLabel);
                end
                if isfield(subplotOptions(i),'Legend')
                    legendLabels = [subplotOptions.Legend];
                    legend(legendLabels, 'location', 'northwest');
                end
            
                if isfield(subplotOptions(i),'xLim')
                    xlim(subplotOptions(i).xLim);
                end
                if isfield(subplotOptions(i),'yLim')
                    ylim(subplotOptions(i).yLim);
                end
                if isfield(subplotOptions(i),'title')
                    title(subplotOptions(i).title);
                end
                hold off;
            end
        end
    end
    
    %% Abstract methods (non‐matching functions)
    methods (Abstract)
        % Each subclass must set the optimization defaults (using a second inputParser, p3)
        set_optimization_defaults_and_parse(obj, inputs, p3)
        
        % Each subclass must implement its own simulation run routine.
        runsimulation(obj, inputs, optionsFMINCON)
        
        % Each subclass must process the raw simulation outputs.
        processoutputs(obj)
    end
    
    %% Static Abstract methods (non‐matching functions)
    methods (Abstract, Static)
         % Each subclass must tell the base class which required parameters
        % (a cell array of strings) it needs.
        req = get_required_params()
        
        % Each subclass must return its optional parameters with defaults,
        % as a cell array where each row is {parameterName, defaultValue, validationFcn}
        params = get_parameters_withdefaults(inputs,validationFcn)
    end
end
