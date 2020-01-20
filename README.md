# interaction_model

 This repository contains the interaction model and fitting process. 

 Files are:
 - Setup: Loads the required packages
 - Flex_age_model_trickle.cpp: The model itself, written in Rcpp
 - Generate_simulations_trickle: used to generate simulations for the model
 - Mixing_UK_ages: generates the contact matrix
 - SIPR_1: the MCMC fitting process for the model (first chain, for second chain change seed to 302)
 - mcmcMH_S: the mcmc algorithm
 - Params_and_summary_for_MCMC: Loads the fixed parameters and some functions to run the MCMC
 - beta_priors_r0_calc: used in the calculation of the priors for the transmission parameters

 The data file contains a spreadsheet with all the details of each simulation
 
