% Clean workspace
clc; clearvars;

% Add folders and subfolders to path
addpath(genpath([pwd '/inputFiles']));
addpath(genpath([pwd '/outputFiles'])); 
addpath(genpath([pwd '/src']));
addpath(genpath([pwd '/lib']));
addpath(genpath([pwd '/tests']));

% Load defined input file
inputs = loadInputs('inputFile_example_awePower.yml');


% Test to check inputs
test_new_input_method()


% Run
[outputs, optimDetails, processedOutputs] = main_awePower(inputs);

