# interaction_model

 This repository contains the interaction model and fitting process. 
 
 To generate simulations run: "generate_simulations.R"
 To fit a simulation run: "run_simulation_fit.R"
 
 There is an example file "Data/Simulation1_1" as an example, so that the fit can be demonstrated without running simulations
 
 The data file contains a spreadsheet with all the details of each simulation

 Other files sourced are:
 - Setup.R: Loads the required packages
 - Flex_age_model_trickle.cpp: The model itself, written in Rcpp
 - Mixing_UK_ages: generates the contact matrix
 - SIPR_1: the MCMC fitting process for the model (first chain, for second chain change seed to 302). This is the version used for running on a cluster.
 - mcmcMH_S: the mcmc algorithm
 - Params_and_summary_for_MCMC: Loads the fixed parameters and some functions to run the MCMC
 - beta_priors_r0_calc: used in the calculation of the priors for the transmission parameters

 
