% Run the below function after successful execution of 'runSimulation.m' to plot results

clear global

% Load input sheet
inputSheet_MW_scale;

% Load processed outputs
load("outputFiles\inputSheet_MW_scale_processedOutputs.mat");

% Run 
plotResults(inputs, processedOutputs);
