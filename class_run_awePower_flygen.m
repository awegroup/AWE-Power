% Clean workspace
clc; clearvars;

% Add folders and subfolders to path
addpath(genpath([pwd '/inputFiles']));
addpath(genpath([pwd '/outputFiles'])); 
addpath(genpath([pwd '/src']));
addpath(genpath([pwd '/lib']));

% Load input file
inputs = loadInputs('inputFile_example_Class_flygen.yml');

% Run
[inputs, outputs, optimDetails, processedOutputs, simulationQSM] = main_Class_flygen(inputs);

% Plotting
% simulationQSM.plot_all_results();

% simulationQSM.plot_speed_results();
% simulationQSM.plot_cl_results();
% simulationQSM.plot_cd_results();
% simulationQSM.plot_tether_lengths();
% simulationQSM.plot_roll_avg_pattern();
% simulationQSM.plot_forces();
% simulationQSM.plot_power();
% simulationQSM.plot_power_curve_comparisonLoyd();
% simulationQSM.plot_wind_profile();