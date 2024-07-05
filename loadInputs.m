function inputs = loadInputs(yamlFilePath)
    % Check if the file exists
    if exist(yamlFilePath, 'file') ~= 2
        error('File not found: %s', yamlFilePath);
    end
    
    % Read the YAML file
    yamlData = yaml.ReadYaml(yamlFilePath);

    % Get all field names
    fields = fieldnames(yamlData);
    
    % Loop through each field and convert cell arrays to double arrays
    for i = 1:numel(fields)
        field = fields{i};
        if iscell(yamlData.(field))
            yamlData.(field) = cell2mat(yamlData.(field));
        end
    end
    
    % Populate the inputs structure
    inputs = struct();
    
    % Define Input Structure
    inputs.name = yamlData.name;
    
    % Model Settings
    inputs.numDeltaLelems = yamlData.numDeltaLelems;
    inputs.FgToggle = yamlData.FgToggle;
    inputs.vertWindProfile = yamlData.vertWindProfile;
    
    % Wind Parameters
    inputs.vw_ref = yamlData.vw_ref;
    inputs.h_ref = yamlData.h_ref;
    inputs.windShearExp = yamlData.windShearExp;
    
    % Sample vertical wind profile datasets (for inputs.vertWindProfile = 1)
    inputs.h_windDataset = yamlData.h_windDataset;
    inputs.v_windDataset_Cabauw = yamlData.v_windDataset_Cabauw;
    inputs.v_windDataset_Ijmuiden = yamlData.v_windDataset_Ijmuiden;
    inputs.windProfile_h = yamlData.h_windDataset;
    inputs.windProfile_vw = yamlData.v_windDataset_Cabauw;
    
    % System Parameters
    inputs.S = yamlData.S;
    inputs.AR = yamlData.AR;
    inputs.P_ratedElec = yamlData.P_ratedElec;
    inputs.massOverride = yamlData.massOverride;
    inputs.kiteMass = yamlData.kiteMass;
    inputs.peakM2E_F = yamlData.peakM2E_F;
    
    % Tether Parameters
    inputs.Ft_max = yamlData.Ft_max;
    inputs.Ft_max_SF = yamlData.Ft_max_SF;
    inputs.maxTeLen = yamlData.maxTeLen;
    inputs.maxHeight = yamlData.maxHeight;
    inputs.minGroundClear = yamlData.minGroundClear;
    inputs.Te_matStrength = yamlData.Te_matStrength;
    inputs.Te_matDensity = yamlData.Te_matDensity;
    
    % Aerodynamic Parameters
    inputs.Cl_maxAirfoil = yamlData.Cl_maxAirfoil;
    inputs.Cl_eff_F = yamlData.Cl_eff_F;
    inputs.Cl0_airfoil = yamlData.Cl0_airfoil;
    inputs.e = yamlData.e;
    inputs.Cd0 = yamlData.Cd0;
    inputs.Cd_c = yamlData.Cd_c;
    
    % Drum operational Limits
    inputs.v_d_max = yamlData.v_d_max;
    inputs.a_d_max = yamlData.a_d_max;
    
    % Drivetrain Efficiency Parameters
    inputs.etaGen.param = yamlData.etaGen_param;
    inputs.etaGearbox = yamlData.etaGearbox;
    inputs.etaSto = yamlData.etaSto;
    inputs.etaPE = yamlData.etaPE;
    
    % Environmental Constants
    inputs.gravity = yamlData.gravity;
    inputs.airDensity = yamlData.airDensity;
end
