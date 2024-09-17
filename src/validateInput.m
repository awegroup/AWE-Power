function validateInput(data, schema)
  % VALIDATEINPUT Validates an input structure against a given schema
  %
  % This function checks if the input data structure meets the criteria defined
  % in the schema. It verifies the presence of required fields and applies
  % validation functions to ensure the data conforms to expected formats or
  % constraints. If validation fails, an error is raised with detailed messages.
  %
  % Args:
  %   data (struct): The input data structure to validate. This structure must
  %       contain all fields specified in the schema.
  %
  %   schema (struct): The schema defining the validation rules. Each field in
  %       the schema corresponds to a field in the data. The schema can include
  %       function handles for validation or nested structures with further
  %       validation rules.
  %
  % Raises:
  %   error: If any required field is missing or if any validation function
  %          fails, an error is thrown with a message detailing the validation
  %          errors.
  %
  % Example:
  %   % Define a schema with validation functions
  %   schema = struct('field1', @validateField1, 'field2', @validateField2);
  %
  %   % Validate input data
  %   validateInput(data, schema);
  %
  % Notes:
  %   - The schema can include nested structures for validating subfields.
  %   - Validation functions should throw an error if validation fails.
  %   - The function currently does not handle recursive validation for deeply
  %     nested schemas but can be extended for such use cases.


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
