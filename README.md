# FixedWingGG-Power

This code is an implementation of a quasi-steady model estimating the power curve of fixed-wing ground-generation airborne wind energy systems.

1. It is a fast cycle-power computation model based on quasi-steady equilibrium, which accounts for the effects of pattern elevation, mass (gravity), vertical wind shear, and drivetrain losses.
2. It is suitable for sensitivity and scalability studies, making it a valuable tool for evaluating design and innovation trade-offs.
3. It is suitable for integration with cost models and systems engineering tools. This enhances the applicability of the proposed model in exploring the potential of airborne wind energy in the energy system.

[Link to the associated scientific article: ...]

---

## Overview of Scripts

The code is divided into 5 scripts:

1. `runMain.m`
2. `main.m`
3. `objective.m`
4. `constraints.m`
5. `compute.m`

The model requires definition of an input file. Some pre-defined input files can be found in the folder named `Input-sheets`. 

---

## Example Run

As an example, the model can be run using pre-defined inputs simulating a 150kW system by simply executing `runMain.m` in the working directory. This will also plot all the relevant outputs.

---

## To Run with User-defined Inputs

To run the model with user-defined inputs, a data structure named 'inputs' needs to be created with all necessary inputs as defined in `inputSheet_AP3.m` from the `Input-sheets` folder. 
The model can then be run by calling the newly defined inputs script instead of the pre-defined `InputSheet_AP3.m` at the start of `runMain.m`.

### Generated Output Files

1. `optData` has the details regarding the optimization.
2. `outputs` has all the raw outputs.
3. `processedOutputs` has post-processed outputs with relevant results for better visualization.

---

## Author Details

Name: Rishikesh Joshi  
Email: [r.joshi@tudelft.nl](mailto:r.joshi@tudelft.nl), [rishikeshsjoshi@gmail.com](mailto:rishikeshsjoshi@gmail.com)  
Affiliation: Wind Energy Section, Delft University of Technology


