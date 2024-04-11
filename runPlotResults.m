% Run the below function after successful execution of 'runSimulation.m' to plot results

clear global

% Load input sheet
inputSheet_Example;

% Load processed outputs
load("outputFiles\inputSheet_Example_processedOutputs.mat");

% Run 
plotResults(inputs, processedOutputs);
