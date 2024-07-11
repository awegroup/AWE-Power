function saveResults(inputs, optimDetails, outputs, processedOutputs)

  if not(isfolder('outputFiles'))
    mkdir 'outputFiles';
  end
  % Change names to associate with specific input file
  save(['outputFiles/' inputs.name '_' 'optimDetails' '.mat'], 'optimDetails');
  save(['outputFiles/' inputs.name '_' 'outputs' '.mat'], 'outputs');
  save(['outputFiles/' inputs.name '_' 'processedOutputs' '.mat'], 'processedOutputs');
  
end