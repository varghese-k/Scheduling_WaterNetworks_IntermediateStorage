
## Script to generate optimal schedules for operation of water distribution netwoks with intermediate storage facilities

+ Implementation of the method described in the paper https://doi.org/10.1016/j.compchemeng.2018.04.017

+ Software required: MATLAB 

+ The folder `epanet` has the example network file, the epanet network simulator file, and the functions to connect epanet with MATLAB

+ The folder `out` has the outputs generated after each run. The folders are arranged based on the date and time of invoking the solver.

+ The folder `src` has the script to run the optimizer: `Scheduling.m`
    + The first section invokes the hydraulic solver to estimate the flow rates under different valve settings. Please read the note before the definition of `pipeIndices` if using a different network. EPANET works a little counterintuitively on pipe indexing. 
    + The second section formulates and solves a linear program to identify a set of efficeint states 
    + The third section formulates and solves an integer linear program to identify the optimal sequence of states which ensures that the tank levels are within the allowable minimum and maximum values