% Run this script after successful execution of 'run_awePower.m'

% Add folders and subfolders to path
addpath(genpath([pwd '/inputFiles']));
addpath(genpath([pwd '/outputFiles'])); 
addpath(genpath([pwd '/src']));

% Load input file
inputs = loadInputs('inputFile_example_awePower.yml');
% Append dependent inputs
inputs = appendInputs(inputs);

% Run 
plotResults_awePower(inputs);

