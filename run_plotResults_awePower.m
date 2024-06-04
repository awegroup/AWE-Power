% Run the below function after successful execution of 'runSimulation.m' to plot results

% Add folders and subfolders to path
addpath(genpath([pwd '\inputSheets']));
addpath(genpath([pwd '\outputFiles'])); 
addpath(genpath([pwd '\src']));

% Load input sheet
inputSheet_Example;

% Load processed outputs
load("outputFiles\inputSheet_Example_processedOutputs.mat");

% Run 
plotResults_awePower(inputs, processedOutputs);

