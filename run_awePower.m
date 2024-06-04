% Clean workspace
clc; clearvars;

% Add folders and subfolders to path
addpath(genpath([pwd '\inputSheets']));
addpath(genpath([pwd '\outputFiles'])); 
addpath(genpath([pwd '\src']));

% Load defined input sheet
inputSheet_Example;

% Enter your defined input sheet name as an argument 
[inputs, outputs, optimDetails, processedOutputs] = main_awePower(inputs);

