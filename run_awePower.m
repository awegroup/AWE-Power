% Clean workspace
clc; clearvars;

% Add folders and subfolders to path
addpath(genpath([pwd '/inputFiles']));
addpath(genpath([pwd '/outputFiles'])); 
addpath(genpath([pwd '/src']));
addpath(genpath([pwd '/lib']));

% Filepath to simulation parameters
inputFile = './inputFiles/inputFile_scalingEffects_awePower.yml';

% Load and validate parameters
inputs = loadInputs(inputFile);
validateInput(inputs, inputValidators());

% Calculate additional parameters
inputs = appendInputs(inputs);

% Run
[outputs, optimDetails, processedOutputs] = main_awePower(inputs);
