function schema = inputValidators()
% Definition of the data schema 
%
% Each parameter will be validated through a defined function handle.
% Two types of function handles are supported:
%   1) Validate attributes
%   (https://nl.mathworks.com/help/matlab/ref/validateattributes.html)
%   2) Custom validation function that throws an error if false
%
% Returns:
%   struct : data schema
%
schema = struct();

% File details
schema.name = @(x) validateattributes(x, {'char'}, {'nonempty'});

% Model Settings
schema.numDeltaLelems = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.FgToggle = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.vertWindProfile = @mustBeValidNumericBoolean;

% Wind Parameters
schema.vw_ref = @(x) validateattributes(x, {'numeric'}, {'increasing', 'positive'});
schema.h_ref = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.windShearExp = @(x) validateattributes(x, {'numeric'}, {'nonempty'});

% Sample vertical wind profile datasets
schema.h_windDataset = @(x) validateattributes(x, {'numeric'}, {'increasing', 'positive'});
schema.v_windDataset_Cabauw = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.v_windDataset_Ijmuiden = @(x) validateattributes(x, {'numeric'}, {'nonempty'});

% Wind Profile selection
schema.windProfile_h = @(x) validateattributes(x, {'char', 'numeric'}, {'nonempty'});
schema.windProfile_vw = @(x) validateattributes(x, {'char', 'numeric'}, {'nonempty'});

% System Parameters
schema.S = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.AR = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.P_ratedElec = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.massOverride = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.kiteMass = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.peakM2E_F =@(x) validateattributes(x, {'numeric'}, {'nonempty'});

% Tether Parameters
schema.Ft_max = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.Ft_max_SF = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.maxTeLen = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.maxHeight = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.minGroundClear = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.Te_matStrength = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.Te_matDensity = @(x) validateattributes(x, {'numeric'}, {'nonempty'});

% Aerodynamic Parameters
schema.Cl_maxAirfoil = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.Cl_eff_F = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.Cl0_airfoil = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.e = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.Cd0 = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.Cd_c = @(x) validateattributes(x, {'numeric'}, {'nonempty'});

% Drum Operational Limits
schema.v_d_max = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.a_d_max = @(x) validateattributes(x, {'numeric'}, {'nonempty'});

% Drivetrain Efficiency Parameters
schema.etaGen.param = @(x) validateattributes(x, {'numeric'}, {'row', 'numel', 4});
schema.etaGearbox = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.etaSto = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.etaPE = @(x) validateattributes(x, {'numeric'}, {'nonempty'});

% Environmental Constants
schema.gravity = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
schema.airDensity = @(x) validateattributes(x, {'numeric'}, {'nonempty'});

end


%% Validation functions

function mustBeValidNumericBoolean(x)
    if ~(x == 0 || x == 1)
        error('Field vertWindProfile must be 0 or 1.');
    end
end
