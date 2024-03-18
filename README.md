# FixedWingGG-Power

A quasi-steady model simulating the reel-out and reel-in phases of fixed-wing ground-generation airborne wind energy systems. The model is formulated as an optimisation problem where the net electrical cycle power of a system is maximised for given wind conditions.

1. The model accounts for the effects of pattern elevation, mass (gravity), vertical wind shear, and drivetrain losses.
2. It is suitable for sensitivity and scalability studies, making it a valuable tool for evaluating design and innovation trade-offs.
3. It is suitable for integration with cost models and systems engineering tools. This enhances the applicability of the proposed model in exploring the potential of airborne wind energy in the energy system.

[Link to the associated scientific article: ...]

---

## Dependencies

The model is built and tested in MATLAB R2023a. Try installing this version if your version of MATLAB does not execute the code successfully.

---

## Overview of the Repository

The Repository does not need any of 3 present Folders (old code)

---

## Pre-defined Example Simulation

1. Execute `Class_runSimulation.m` to generate model outputs. The outputs will be stored in the class property `processedOutputs`.
2. Plotting method is implemented and shown at the end of `Class_runSimulation.m`

---

## To Run with User-defined Inputs

1. A data structure named 'inputs' needs to be created with all necessary inputs as defined in Class_runSimulation.m`, optional inputs can also be given which are defined in class constructor of `KiteQSMsimulation.m`.
2. Results can be generated similar to executing the script `Class_runSimulation.m`.

### Generated Output Files

1. `optimDetails` has the details regarding the internal optimization.
2. `outputs` has all the raw outputs.
3. `processedOutputs` has post-processed relevant outputs for better visualization.

---

## LICENSE

MIT License : Copyright (c) 2024 TU Delft

---

## Author Details

Name: Rishikesh Joshi  
Email: [r.joshi@tudelft.nl](mailto:r.joshi@tudelft.nl), [rishikeshsjoshi@gmail.com](mailto:rishikeshsjoshi@gmail.com)  
Affiliation: Wind Energy Section, Delft University of Technology

---

## Citation Details

---


