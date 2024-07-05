function [testResult] = test_new_input_method()

  inputFile_example_awePower;

  yamlInputs = loadInputs('inputFile_example_awePower.yml');
  yamlInputs = appendInputs(yamlInputs);

  testResult = isequal(inputs, yamlInputs);

end