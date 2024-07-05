% Run this script after successful execution of 'run_awePower.m'

% Add folders and subfolders to path
addpath(genpath([pwd '/inputFiles']));
addpath(genpath([pwd '/outputFiles'])); 
addpath(genpath([pwd '/src']));

% Load input sheet
% inputFile_example_awePower;

newInputs = loadInputs('inputs.yml');
newInputs = appendInputs(newInputs);

% Run 
plotResults_awePower(newInputs);

