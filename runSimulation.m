% Make sure to clean the variables from the workspace
clc
clearvars
clear global

% Add folders and subfolders to path
addpath(genpath([pwd '\inputSheets']));
addpath(genpath([pwd '\outputFiles'])); 
addpath(genpath([pwd '\src']));

% Load defined input sheet
inputSheet_Example;
% inputSheet_MW_scale;
% inputSheet_MW_scale_EcoModel;


% Enter your defined input sheet name as an argument 
[inputs, outputs, optimDetails, processedOutputs] = main(inputs);

