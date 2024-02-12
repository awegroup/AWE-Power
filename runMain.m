%% Run this script to save outputs of the power computation model
clc
clearvars
clear global

%% Defined input sheet
inputSheet_AP3;
%inputSheet;

%% Run Model
[optData,outputs,processedOutputs] = main(inputs);

%% Save outputs
save(fullfile('Outputs','inputs.mat'))
save(fullfile('Outputs','optData.mat'))
save(fullfile('Outputs','outputs.mat'))
save(fullfile('Outputs','processedOutputs.mat'))


