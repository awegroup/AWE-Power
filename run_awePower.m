% Clean workspace
clc; clearvars;

% Add folders and subfolders to path
addpath(genpath([pwd '/inputFiles']));
addpath(genpath([pwd '/outputFiles'])); 
addpath(genpath([pwd '/src']));
addpath(genpath([pwd '/lib']));
addpath(genpath([pwd '/tests']));

% Load defined input file
% inputFile_example_awePower;

newInputs = loadInputs('inputs.yml');
newInputs = appendInputs(newInputs);

test_new_input_method(newInputs)


% Run
[outputs, optimDetails, processedOutputs] = main_awePower(newInputs);

