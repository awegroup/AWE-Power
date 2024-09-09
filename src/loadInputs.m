function inputs = loadInputs(filepath)
  % loadInputs: Load inputs from .yml or .m files
  %
  % This function reads input parameters from either a YAML (.yml or .yaml) file or a MATLAB script (.m) file.
  % It processes the content and ensures any cell arrays are converted to double arrays, so the data is usable
  % in subsequent operations.
  %
  % Args:
  %   filepath (char): Full path to the parameters file, which can be either a YAML file or a MATLAB .m file.
  %
  % Returns:
  %   struct: A structure containing the loaded inputs.
  %
  % Raises:
  %   error: If the file does not exist or the file type is not supported.


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

