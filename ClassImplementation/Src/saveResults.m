function saveResults(inputs, optimDetails, outputs, processedOutputs)
  % SAVERESULTS Saves the results of the airborne wind energy system to MAT-files
  %
  % This function saves the details of the optimization, raw outputs, and
  % processed outputs into separate MAT-files in the 'outputFiles' directory.
  % The files are named based on the provided input file name, ensuring that
  % the saved files are associated with specific input configurations.
  %
  % Inputs:
  %   inputs - A structure containing system parameters and reference values.
  %     .name: The name of the input file, used to generate the output file names.
  %
  %   optimDetails - A structure containing details of the optimization process.
  %
  %   outputs - A structure containing raw output data from the system.
  %
  %   processedOutputs - A structure containing processed data from the system.
  %
  % Files Saved:
  %   - 'outputFiles/[inputs.name]_optimDetails.mat': Contains the optimization details.
  %   - 'outputFiles/[inputs.name]_outputs.mat': Contains the raw output data.
  %   - 'outputFiles/[inputs.name]_processedOutputs.mat': Contains the processed data.
  %
  % Example:
  %   saveResults(inputs, optimDetails, outputs, processedOutputs);
  %
  % This will save the results in the 'outputFiles' directory with filenames
  % based on the 'name' field in the 'inputs' structure.


  if not(isfolder('outputFiles'))
    mkdir 'outputFiles';
  end
  % Change names to associate with specific input file
  save(['outputFiles/' inputs.name '_' 'optimDetails' '.mat'], 'optimDetails');
  save(['outputFiles/' inputs.name '_' 'outputs' '.mat'], 'outputs');
  save(['outputFiles/' inputs.name '_' 'processedOutputs' '.mat'], 'processedOutputs');

end