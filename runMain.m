%% Run this script to save outputs of the power computation model.
clc
clearvars
clear global

%% Call your defined input sheet here. Comment out or delete any other input file name.
inputSheet_AP3;

%% Run Model
[optData,outputs,processedOutputs] = main(inputs);

%% Save outputs
save(fullfile('Outputs','inputs.mat'))
save(fullfile('Outputs','optData.mat'))
save(fullfile('Outputs','outputs.mat'))
save(fullfile('Outputs','processedOutputs.mat'))

%% run the script 'plotResults.m' to visulaize some of the relevant results.


