This code is an implementation of a quasi-steady model estimating the power curve of fixed-wing ground-generation airborne wind energy systems. 

Link to the associated scientific article: ...

Following are the key points describing the model and its purpose:
1) It is a fast cycle-power computation model which accounts for the effects of pattern elevation, mass (gravity), vertical wind shear and drivetrain losses.
2) It is suitable for sensitivity and scalability studies making it a valuable tool for evaluating design and innovation trade-offs.
3) It is suitable for integration with cost models and systems engineering tools. This enhances the applicability of the proposed model in exploring potential of airborne wind energy to the energy system.

The model is divided into 4 scripts:
2) main.m
3) objective.m
4) constraints.m
5) compute.m

(Please ignore the rest of the files when trying to understand the model and the implementation for the first time)

As an example, the model can be run using pre-defined inputs simulating a 150kW system by simply executing 'runMain.m'. This will also plot all the relevant outputs.

To run the model with user-defined inputs, a data structure named 'inputs' needs to be created with all necessary inputs as defined in 'inpusSheet_AP3.m'. 
The model can then be run using the command '[optData,outputs,processedOutputs] = main(inputs);' by parsing the newly defined inputs.

Following are the generated output files:
1) 'optData' stores the details regarding the optimisation itself.
2) 'outputs' stores all the raw outputs from the model.
3) 'processedOutputs' post processess and stores relevant outputs for better visualisation.

Note: The pre-defined 'runMain.m' could be used as an example for visualizing results.

Author details:
Rishikesh Joshi
r.joshi@tudelft.nl, rishikeshsjoshi@gmail.com
Wind Energy Section, Delft University of Technology.


