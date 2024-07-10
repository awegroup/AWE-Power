function inputs = loadInputs(filepath)
% Load inputs from .yml and .m files 
%
% Args:
%   filepath (char): path to the parameters file
%
% Returns:
%   struct

% Check if the file exists
if exist(filepath, 'file') ~= 2
    error('File not found: %s', filepath);
end

[~, name, ext] = fileparts(filepath);

switch ext
    case {'.yml', '.yaml'}        
        inputs = yaml.ReadYaml(filepath);        
    
        % Get all field names
        fields = fieldnames(inputs);
        
        % Loop through each field and convert cell arrays to double arrays
        % TODO: replace with recursive function
        for i = 1:numel(fields)
            field = fields{i};
            if isstruct(inputs.(field))
                subfields = fieldnames(inputs.(field));
                for j = 1:numel(subfields)
                    subfield = subfields{j};
                    if iscell(inputs.(field).(subfield))
                        inputs.(field).(subfield) = cell2mat(inputs.(field).(subfield));
                    end
                end
            elseif iscell(inputs.(field))
                inputs.(field) = cell2mat(inputs.(field));
            end
        end

    case '.m'
        inputs = feval(name);

    otherwise
        error('Expected input file of type .yml or .m, instead received %s', ext);
end

