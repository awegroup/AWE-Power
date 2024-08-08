function validateInput(data, schema)
% Validate the input structure for the data schema 
%
% Args:
%   data (struct): input data structure to validate 
%   schema (struct): data schema 
%
% Raises:
%   error when validation fails

% Initialize variables
fields = fieldnames(schema);
errorMessages = {};

for i = 1:numel(fields)
    field = fields{i};

    % Check if field is present
    if ~isfield(data, field)
        error('Missing required field: %s', field);
    end

    % Perform data validations
    % TODO: add recursive validation for nested schema
    try
        if isa(schema.(field), 'function_handle')
            validationFunc = schema.(field);
            validationFunc(data.(field));          
        elseif isstruct(schema.(field))
            subfields = fieldnames(schema.(field));
            for j = 1:numel(subfields)
                subfield = subfields{j};
                if isa(schema.(field).(subfield), 'function_handle')
                    validationFunc = schema.(field).(subfield);
                    validationFunc(data.(field).(subfield));   
                end
            end
        else
            warning('Field %s has no recognized validation function handle', field);        
        end        
    catch ME
        errorMessages{end+1} = [sprintf('Validation error for %s: ', field), ME.message];
    end

end

% Report error messages
if ~isempty(errorMessages)
    errorStr = strjoin(errorMessages, [newline '-----' newline]);
    error(['Input validation failed with the following errors:', [newline '-----' newline], errorStr]);
else
   % disp('All validations passed.');
end

end
