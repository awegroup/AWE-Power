Author details:
Rishikesh Joshi
r.joshi@tudelft.nl
Wind Energy Section, Delft University of Technology.

This code is an implementation of a Quasi-steady model estimating the power curve of fixed-wing ground-generation airborne wind energy systems.
Following is the abstract of the paper attached to this code:

Paper abstract:
We present a fast cycle-power computation model by capturing the interactions between the critical high-level design and operational parameters of the fixed-wing ground-generation airborne wind energy systems. The kite is modelled as a point mass and uses a quasi-steady approach to model the reel-out and reel-in phases of the system. The reel-out phase is discretized in a number of length elements, and a single flight state is approximated per element. The model accounts for the effects of pattern elevation, mass (gravity), vertical wind shear and drivetrain losses. The system design parameters are the kite wing area, aspect ratio, aerodynamic properties, tether properties and dimensions, drivetrain properties, etc. The operational design parameters are the reel-out length, pattern elevation angle, opening cone angle, starting pattern radius, and kite speed. For a defined system, cycle power is maximised using gradient-based optimisation, resulting in optimal operational set points for the system. The model provides fast power curve computations, mainly suitable for sensitivity and scalability studies to support design and innovation trade-offs. It is also suitable for integration with cost models and systems engineering tools. The model results are compared against 6-DoF simulation results of a $\SI{150}{kW}$ system. The interdependency between key system design parameters and their impact on upscaling is presented.

The model is divided into 4 scripts:
2) main.m
3) objective.m
4) constraints.m
5) compute.m

(Please ignore the rest of the files when trying to understand the model and the implementation for the first time)

As an example, the model can be run using a pre-initialised input file for a 150kW system by directly running 'runMain.m'. This will also plot all the relevant outputs.

The model requires a data structure named 'inputs' and can be run using the command '[optData,outputs,processedOutputs] = main(inputs);'.

'optData' stores the optimisation results
'outputs' stores the detailed inputs, and
'processedOutputs' post processess and stores some outputs

'inputSheet_AP3' can be run as an example to generate the required 'inputs' structure. 
The user can then create a similar file or edit the same file to generate a new 'inputs' structure.




