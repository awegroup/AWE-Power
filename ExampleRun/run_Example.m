%% Run this script to save outputs of the power computation model
clc
clearvars
clear global

%% Defined input sheet
inputSheet_Example;

%% Run Model
[optData,outputs,processedOutputs] = main(inputs);

%% Save outputs
save(fullfile('Outputs_Example','inputs.mat'))
save(fullfile('Outputs_Example','optData.mat'))
save(fullfile('Outputs_Example','outputs.mat'))
save(fullfile('Outputs_Example','processedOutputs.mat'))

%% run the script 'plotResultsExample.m' to visulaize some of the relevant results


