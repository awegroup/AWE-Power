function schema = inputValidators()
% Definition of the data schema 
%
% Each parameter will be validated through a defined function handle.
% Three types of function handles are recommended:
%   1) Data type identification
%   (https://nl.mathworks.com/help/matlab/data-type-identification.html)
%   2) Validate attributes
%   (https://nl.mathworks.com/help/matlab/ref/validateattributes.html)
%   3) Custom validation function
%
% Returns:
%   struct : data schema
%
schema = struct();

% File details
schema.name = @ischar;

% Model Settings
schema.numDeltaLelems = @isnumeric;
schema.FgToggle = @islogical;
schema.vertWindProfile = @islogical;

% Wind Parameters
schema.vw_ref = @(x) validateattributes(x, {'numeric'}, {'increasing', 'positive'});
schema.h_ref = @isnumeric;
schema.windShearExp = @mustBeValidWindShearExp;

% Sample vertical wind profile datasets
schema.h_windDataset = @(x) validateattributes(x, {'numeric'}, {'increasing', 'positive'});
schema.v_windDataset_Cabauw = @isnumeric;
schema.v_windDataset_Ijmuiden = @isnumeric;

% System Parameters
schema.S = @isnumeric;
schema.AR = @isnumeric;
schema.P_ratedElec = @isnumeric;
schema.massOverride = @isnumeric;
schema.kiteMass = @isnumeric;
schema.peakM2E_F = @isnumeric;

% Tether Parameters
schema.Ft_max = @isnumeric;
schema.Ft_max_SF = @isnumeric;
schema.maxTeLen = @isnumeric;
schema.maxHeight = @isnumeric;
schema.minGroundClear = @isnumeric;
schema.Te_matStrength = @isnumeric;
schema.Te_matDensity = @isnumeric;

% Aerodynamic Parameters
schema.Cl_maxAirfoil = @isnumeric;
schema.Cl_eff_F = @isnumeric;
schema.Cl0_airfoil = @isnumeric;
schema.e = @isnumeric;
schema.Cd0 = @isnumeric;
schema.Cd_c = @isnumeric;

% Drum Operational Limits
schema.v_d_max = @isnumeric;
schema.a_d_max = @isnumeric;

% Drivetrain Efficiency Parameters
schema.etaGen.param = @(x) validateattributes(x, {'numeric'}, {'row', 'numel', 4});
schema.etaGearbox = @isnumeric;
schema.etaSto = @isnumeric;
schema.etaPE = @isnumeric;

% Environmental Constants
schema.gravity = @isnumeric;
schema.airDensity = @isnumeric;

end


%% Validation functions

function mustBeValidWindShearExp(x)
    if ~(x == 0.143 || x == 0.11)
        error('Field windShearExp must be 0.143 (over land) or 0.11 (over sea).');
    end
end
