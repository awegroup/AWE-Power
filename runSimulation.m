% Make sure to clean the variables from the workspace
clc
clearvars
clear global

tic % log runtim

% Add folders and subfolders to path
addpath(genpath([pwd '\InputSheets']));
addpath(genpath([pwd '\OutputFiles'])); 
addpath(genpath([pwd '\Src']));

% Load defined input sheet
inputSheet_Dylan;

% Enter your defined input sheet name as an argument 
[inputs, outputs, optimDetails, processedOutputs] = main(inputs, 'inputSheet_Dylan');

toc

