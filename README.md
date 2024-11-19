# AWE-Power
![GitHub License](https://img.shields.io/github/license/awegroup/AWE-Power)
![Static Badge](https://img.shields.io/badge/MATLAB-R2021b-blue)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13842297.svg)](https://doi.org/10.5281/zenodo.13842297)

A quasi-steady model simulating the reel-out and reel-in phases of fixed-wing ground-generation airborne wind energy systems. The model is formulated as an optimisation problem where the net electrical cycle power of a system is maximised for given wind conditions.

1. The model accounts for the effects of pattern elevation, mass (gravity), vertical wind shear, and drivetrain losses.
2. It is suitable for sensitivity and scalability studies, making it a valuable tool for evaluating design and innovation trade-offs.
3. It is suitable for integration with cost models and systems engineering tools. This enhances the applicability of the proposed model in exploring the potential of airborne wind energy in the energy system.


## How to cite

If you use the model then please cite the corresponding paper:

Joshi, R., Schmehl, R., and Kruijff, M.: Power curve modelling and scaling of fixed-wing ground-generation airborne wind energy systems, Wind Energ. Sci., 9, 2195–2215, https://doi.org/10.5194/wes-9-2195-2024, 2024.

## Dependencies

The model is built and tested in MATLAB R2021b (without additional add-ons). Try installing this version if your version of MATLAB does not execute the code successfully.


## Installation and execution 

Please Clone or Download the repository to start using the model.


## Overview of the Repository

The Repository consists following folders:

1. `src`: Contains the source code of the model.
2. `inputFiles`: Contains some pre-defined input files.
3. `outputFiles`: Contains some generated output files using the pre-defined input files.
4. `lib`: Contains external functions required to run the model.
5. `scripts`: Contains some scripts used to generate plots used in the associated journal publication.
6. `tests`: Contains some tests used during code developement.


## Pre-defined Example Simulation

A pre-defined input file named `inputFile_example_awePower.yml` is stored in the folder `inputFiles`. This input file simulates a system with rated electrical power of 150 kW.
1. Execute `run_awePower.m` to generate and save model outputs. The outputs will be stored in the folder `outputFiles`. Name of every output file is prefixed with the respective input file name.
2. Execute `run_plotResults_awePower.m` to visualize relevant outputs.


## To Run with User-defined Inputs

1. A data structure named 'inputs' needs to be created with all necessary inputs as defined in `inputFile_example_awePower.yml` from the `inputFiles` folder.
2. Results can be generated by executing the script `run_awePower.m`. Within the script, ensure to pass as an argument, the name of your newly defined input sheet.
3. Results can be visualized by executing the script `run_plotResults_awePower.m`. Within the script, ensure to load your specific input file.

### Generated Output Files

The generated output files have a prefix as the name of the used input file. Following .mat files get stored in the `outputFiles` folder.

1. `optimDetails` has the details regarding the optimization.
2. `outputs` has all the raw outputs.
3. `processedOutputs` has post-processed relevant outputs for better visualization.

## Licence
This project is licensed under the MIT License. Please see the below WAIVER in association with the license.

## Acknowledgement
The project was supported by the Digital Competence Centre, Delft University of Technology.

## WAIVER
Technische Universiteit Delft hereby disclaims all copyright interest in the program “AWE-Power” (a fast cycle-power computation model) written by the Author(s).

Prof.dr. H.G.C. (Henri) Werij, Dean of Aerospace Engineering

Copyright (c) 2024 Rishikesh Joshi






