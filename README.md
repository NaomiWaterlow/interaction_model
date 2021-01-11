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

 
Code was run on R version 3.6.0, and checks where subsequently made in R version 4.0.0

Session information is: 
R version 4.0.0 (2020-04-24)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Catalina 10.15.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] tmvtnorm_1.4-10   gmm_1.6-5         sandwich_2.5-1    Matrix_1.2-18     mvtnorm_1.1-1    
 [6] coda_0.19-4       bbmle_1.0.23.1    MASS_7.3-53       socialmixr_0.1.8  reshape2_1.4.4   
[11] ggplot2_3.3.3     data.table_1.13.6 plyr_1.8.6        gridExtra_2.3     deSolve_1.28     
[16] Rcpp_1.0.5.4     

loaded via a namespace (and not attached):
 [1] pillar_1.4.7        compiler_4.0.0      tools_4.0.0         lattice_0.20-41    
 [5] jsonlite_1.7.2      lubridate_1.7.9.2   lifecycle_0.2.0     tibble_3.0.4       
 [9] gtable_0.3.0        pkgconfig_2.0.3     rlang_0.4.10        rstudioapi_0.13    
[13] curl_4.3            withr_2.3.0         dplyr_1.0.2         stringr_1.4.0      
[17] httr_1.4.2          xml2_1.3.2          generics_0.1.0      vctrs_0.3.6        
[21] grid_4.0.0          tidyselect_1.1.0    glue_1.4.2          R6_2.5.0           
[25] bdsmatrix_1.3-4     oai_0.3.0           XML_3.99-0.5        purrr_0.3.4        
[29] magrittr_2.0.1      scales_1.1.1        ellipsis_0.3.1      countrycode_1.2.0  
[33] colorspace_2.0-0    wpp2015_1.1-2       numDeriv_2016.8-1.1 stringi_1.5.3      
[37] munsell_0.5.0       crayon_1.3.4        zoo_1.8-8          

