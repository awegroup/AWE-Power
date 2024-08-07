% Clean workspace
clc; clearvars;

% Add folders and subfolders to path
addpath(genpath([pwd '/inputFiles']));
addpath(genpath([pwd '/outputFiles'])); 
addpath(genpath([pwd '/src']));
addpath(genpath([pwd '/lib']));

% Load input file
inputs = loadInputs('inputFile_example_awePower.yml');

% Run
[inputs, outputs, optimDetails, processedOutputs] = main_awePower(inputs);
