% Clean workspace
clc; clearvars;

% Add folders and subfolders to path
addpath(genpath([pwd '/inputFiles']));
addpath(genpath([pwd '/outputFiles'])); 
addpath(genpath([pwd '/src']));

% Load defined input file
inputFile_example_awePower;

% Run
[inputs, outputs, optimDetails, processedOutputs] = main_awePower(inputs);

