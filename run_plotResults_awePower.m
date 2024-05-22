% Run the below function after successful execution of 'runSimulation.m' to plot results

clear global

% Add folders and subfolders to path
addpath(genpath([pwd '\inputSheets']));
addpath(genpath([pwd '\outputFiles'])); 
addpath(genpath([pwd '\src']));

% Load input sheet
inputSheet_Example;
% inputSheet_MW_scale;
% inputSheet_MW_scale_EcoModel;

% Load processed outputs
load("outputFiles\inputSheet_Example_processedOutputs.mat");
% load("outputFiles\inputSheet_MW_scale_processedOutputs.mat");
% load("outputFiles\inputSheet_MW_scale_EcoModel_processedOutputs.mat");

% Run 
plotResults_awePower(inputs, processedOutputs);

