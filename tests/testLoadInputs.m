classdef testLoadInputs < matlab.unittest.TestCase

    properties
        InputFilePathYAML
        InputFilePathM
    end

    methods (TestMethodSetup)
        function setUpPath(testCase)
            % Get the full path to the current file
            currentFilePath = mfilename('fullpath');
            
            % Get the directory of the current file
            [currentDir, ~, ~] = fileparts(currentFilePath);
            
            % Add data folder to path
            addpath([currentDir, '/data']);

            % Construct the full path to the input file
            testCase.InputFilePathM = fullfile(currentDir,  'data', 'inputFile.m');
            testCase.InputFilePathYAML = fullfile(currentDir, 'data', 'inputFile.yml');
        end
    end

    methods(Test)
        
        function loadYaml(testCase)
            inputFile = testCase.InputFilePathYAML;
            inputs = loadInputs(inputFile);
            try
                validateInput(inputs, inputValidators());
                pass = true;
            catch ME
                pass = false;
            end
            testCase.verifyTrue(pass, 'Data validation for .yml input failed')
        end

        function UpdateYamlData(testCase)
            inputFile = testCase.InputFilePathYAML;
            inputs = loadInputs(inputFile);
            validateInput(inputs, inputValidators());
            actual_inputs = appendInputs(inputs);
            expected_inputs = loadInputs(testCase.InputFilePathM);                        
            testCase.verifyEqual(actual_inputs, expected_inputs, 'Updated YAML data is not equal to MATLAB file')

        end

        function loadmFile(testCase)
            inputFile = testCase.InputFilePathM;
            inputs = loadInputs(inputFile);
            try
                validateInput(inputs, inputValidators());
                pass = true;
            catch ME
                pass = false;
            end
            testCase.verifyTrue(pass, 'Data validation for .m input failed')
        end
    end    
end