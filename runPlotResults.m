% Run the below function after successful execution of 'runSimulation.m' to plot results

clear global

% Add folders and subfolders to path
addpath(genpath([pwd '\inputSheets']));
addpath(genpath([pwd '\outputFiles'])); 
addpath(genpath([pwd '\src']));

% Load input sheet
% inputSheet_Example;
inputSheet_MW_scale;
% inputSheet_AWEC2024_ParamStudy;

% Load processed outputs
% load("outputFiles\inputSheet_MW_scale_processedOutputs.mat");
load("outputFiles\inputSheet_MW_scale_processedOutputs.mat");
% load("outputFiles\inputSheet_AWEC2024_ParamStudy_processedOutputs.mat");

% Run 
plotResults(inputs, processedOutputs);

