% Run this script after successful execution of 'run_awePower.m'

% Add folders and subfolders to path
addpath(genpath([pwd '/inputFiles']));
addpath(genpath([pwd '/outputFiles'])); 
addpath(genpath([pwd '/src']));

% Filepath to simulation parameters
inputFile = './inputFiles/inputFile_example_awePower.yml';

% Load and validate parameters
inputs = loadInputs(inputFile);
validateInput(inputs, inputValidators());

% Calculate additional parameters
inputs = appendInputs(inputs);

% Run 
plotResults_awePower(inputs);

