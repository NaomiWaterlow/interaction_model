# Run the MCMC for a given simulation

###### CONDITIONS  #######
#Define number of agegroups
Age_groups <- c(0,2,5,16,65)
#If using other age groups need to change summary_stats function in "Params_and_summary"
setwd("~/Documents/GitHub/interaction_model")
source("setup.R")
###############Load required functions #############

setwd("~/Documents/GitHub/interaction_model")
#Load in mcmcMH_function
source("mcmcMH_S.R")
#Load paramaters and summary statistics function
source("Params_and_summary_for_MCMC.R")
#Load in C++ model
setwd("~/Documents/GitHub/VietnamModel/Model")
sourceCpp("Flex_age_model_trickle.cpp")


########## LOAD "REAL" DATA ##################

setwd("~/Documents/GitHub/Data")
observed.data <-read.csv("Simulation_correct_long_1")

############# SIPR MODEL DETAILS ############

SIPR <- list()

SIPR$name <-"SIPR influenza and RSV"

#compartment names
Names <- c("time")
for(i in 1:length(Age_groups)){
  x<-c(paste0("SS", i), paste0("IS", i),paste0("PS", i),paste0("RS", i),
       paste0("SI", i),paste0("II", i),paste0("PI", i),paste0("SP", i),
       paste0("IP", i), paste0("SR", i),paste0("RR", i),paste0("Rcases", i),
       paste0("Icases", i))
  Names <- c(Names, x)
}
SIPR$state.names <- Names

# parameters to estimate
SIPR$theta.names <- c("bR", "bI", "l_sig",'Rdetect2', 'Rdetect5',
                      'Idetect', 'l_rho', 'seedR', 'seedI')


############# MODEL SIMULATION ##############

# simulate the model
SIPR$simulate <- function (theta, init.state, times)
{
  outall <- as.data.table(ode(y=init.state,
                              t=times,
                              func = derivatives,
                              parms = theta,
                              method="ode45"))
  # change to incidence etc. 
  trajectory <- summary_stats(outall)
  #define detection rates
  p_Rdetect2 <-  exp(as.numeric(theta['Rdetect2']))/(1+exp(as.numeric(theta['Rdetect2'])))
  p_Rdetect5 <- exp(as.numeric(theta['Rdetect5']))/(1+exp(as.numeric(theta['Rdetect5'])))
  p_Idetect <- exp(as.numeric(theta['Idetect']))/(1+exp(as.numeric(theta['Idetect'])))
  # work out num detected in youngest two age groups
  trajectory[, Hos_RSV_0 := RSV_incidence_1*p_Rdetect2]
  trajectory[, Hos_INF_0 := INF_incidence_1*p_Idetect]
  trajectory[, Hos_RSV_1 := RSV_incidence_2*p_Rdetect5]
  trajectory[, Hos_INF_1 := INF_incidence_2*p_Idetect]
  
  # return the simulation
  return(trajectory)
}


############### PRIORS ################

SIPR$dprior<-function (theta, log = FALSE)
{# priors based on R0 values
  log.prior.bR <- dnorm((exp(as.numeric(theta["bR"]))*as.numeric(theta["R0ratR"])), mean=3, sd=0.6, log=T)
  log.prior.bI <- dlnorm(((exp(as.numeric(theta["bI"]))*as.numeric(theta["R0ratI"])-1)), meanlog=0.3, sdlog=0.5,log=T)
  # uniform priors
  log.prior.Rdetect2 <- dunif(as.numeric(theta["Rdetect2"]), min = -100000, max = 100000,
                              log = TRUE)
  log.prior.Rdetect5 <- dunif(as.numeric(theta["Rdetect5"]), min = -100000, max = 100000,
                              log = TRUE)
  log.prior.Idetect <- dunif(as.numeric(theta["Idetect"]), min = -100000, max = 100000,
                             log = TRUE)
  log.prior.seedI <- dunif(as.numeric(theta["seedI"]), min = 0, max = 60, log = TRUE)
  log.prior.seedR <- dunif(as.numeric(theta["seedR"]), min = 0, max = 60, log = TRUE)
  log.prior.sigma <- dunif(as.numeric(theta["l_sig"]), min = 0, max = 1, log = TRUE)
  log.prior.rho <- dunif(as.numeric(theta["l_rho"]), min = 0.01, max = 1, log = TRUE)
  
  # sum each prior (as on log scale)
  log.sum <- log.prior.bR + log.prior.bI +
    log.prior.Rdetect2 + log.prior.Rdetect5 + log.prior.Idetect +
    log.prior.seedR + log.prior.seedI +log.prior.sigma + log.prior.rho
  return(ifelse(log, log.sum, exp(log.sum)))
}


########### LIKLIHOOD OF EACH DATA POINT ##########


SIPR$dPointObs_RSV<- function (data.point, model.point, theta, log = FALSE)
{
  # calculate likelihood for RSV cases in youngest age group
  l<- dpois(x = data.point[["Hos_RSV_0"]], lambda = model.point[["Hos_RSV_0"]],
            log = log)
  # calculate likelihood for RSV cases in second youngest age group
  m<- dpois(x = data.point[["Hos_RSV_1"]], lambda = model.point[["Hos_RSV_1"]],
            log = log)
  # return the sum of the likelihoods for the two age groups
  return(l+m)
}

SIPR$dPointObs_INF<- function (data.point, model.point, theta, log = FALSE)
{
  # calculate likelihood for INF cases in youngest age group
  l<- dpois(x = data.point[["Hos_INF_0"]], lambda = model.point[["Hos_INF_0"]],
            log = log)
  # calculate likelihood for RSV cases in second youngest age group
  m<- dpois(x = data.point[["Hos_INF_1"]], lambda = model.point[["Hos_INF_1"]],
            log = log)
  # return the sum of the likelihoods for the two age groups
  return(l+m)
}


################# LIKLIHOOD OF THE TRAJECTORY ################

dTrajObs<- function (fitmodel, theta, init.state, data, log = FALSE)
{
  # time for simulation to run
  times <- seq(from= 1, to=365, by = 1)
  # simulate the model
  traj <- fitmodel$simulate(theta, init.state, times)
  # set up densities
  dens_INF <- 0
  dens_RSV <- 0
  # for each time step
  for (i in 1:(nrow(data))) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[i, ])
    # calculate the influenza likelihood and add to previous
    dens_INF <- dens_INF + fitmodel$dPointObs_INF(data.point = data.point,
                                                  model.point = model.point, theta = theta, log = TRUE)
    #calculate the RSV likelihood and add to previous
    dens_RSV <- dens_RSV + fitmodel$dPointObs_RSV(data.point = data.point,
                                                  model.point = model.point, theta = theta, log = TRUE)
    #sum the rsv and inf densities
    dens <- dens_RSV + dens_INF
  }
  
  return(ifelse(log, dens, exp(dens)))
}

############## CALCULATES LOG POSTERIOR ################

my_dLogPosterior <- function(fitmodel, theta, init.state, data) {
  # convert to real value (instead of sampled)
  theta["l_rho"] <- exp(as.numeric(theta['l_rho']))
  
  # calculate the prior for parameter vector theta
  log.prior <- fitmodel$dprior(theta, log = TRUE)
  # calculate the log-likelihood of `theta`
  log.likelihood <- dTrajObs(fitmodel, theta, init.state, data, log = TRUE)
  # calulate the log-posterior using the log-prior and log-likelihood
  log.posterior <- log.prior + log.likelihood
  
  return(log.posterior)
}

########### SET UP AND RUN MCMC ##################

Log.posterior <- function(input_params) {
  
  # create parameteter list using fixed parameters and estimated parameters (input_params)  
  theta <- list(Rdetect2=unname(input_params["Rdetect2"]), Rdetect5= unname(input_params["Rdetect5"]),
                Idetect = unname(input_params["Idetect"]),
                l_sig=unname(input_params["sigma"]), gammaR = gammaR, gammaI = gammaI,  RSV_Sus = RSV_Sus,
                Contact_Structure = Contact_structure,
                num_grps = length(Age_groups),
                trickleI=0.1,
                l_rho=log(unname(input_params["rho"])),
                bR=unname(input_params["betaR"]), bI=unname(input_params["betaI"]),
                seedI = unname(input_params["seedI"]), seedR =unname(input_params["seedR"]),
                R0ratI = R0ratI, R0ratR= R0ratR # ratio between R0 and beta
  )
  
  return(my_dLogPosterior(fitmodel = SIPR,
                          theta = theta,
                          init.state = initials,
                          data = observed.data))
}

# initial conditions
init.theta<-c(betaR=-3.148707, betaI=-2.751825,Rdetect2= -5.48409, Rdetect5=-6.878253, Idetect=-6.261598
              ,rho=-2, seedR=200, seedI= 13.00155, sigma=0.8
)
# proposal distribution
proposal.sd <-c(betaR=0.001, betaI=0.001 , Rdetect2= 0.01, Rdetect5=0.001, Idetect=0.01
                , seedR=0.001, seedI=0.1, sigma=0.05, rho=0.1
)
# lower limits
limits.lower<- c(betaR=-Inf, betaI=-Inf, Rdetect2= -Inf, Rdetect5=-Inf, Idetect=-Inf
                 , seedR =0, seedI=0, sigma=0, rho=-Inf
)
# upper limits
limits.upper <-c(betaR=Inf, betaI=Inf, Rdetect2= Inf, Rdetect5=Inf, Idetect=Inf
                 ,  seedR= 90, seedI=90, sigma=1, rho=0
)


#run MCM 
# mcmcMH from FitR package
trace <- mcmcMH(target = Log.posterior, # target distribution
                init.theta = init.theta , 
                proposal.sd = proposal.sd, 
                limits = list(lower=limits.lower ,
                              upper=limits.upper),
                adapt.shape.start = 5000,
                adapt.shape.stop = 50000,
                n.iterations = 100000
)

# save as my_trace
my_trace <- mcmc(trace[[1]])

# save the conditions
trace[["init.theta"]] <- init.theta
trace[["proposal.sd"]] <- proposal.sd
trace[["limits.lower"]] <- limits.lower
trace[["limits.upper"]] <- limits.upper


#Store the trace for later analysis
setwd("~/Documents/GitHub/interaction_model")
save(trace, file = "test2")

