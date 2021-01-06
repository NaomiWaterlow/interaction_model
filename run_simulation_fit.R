# Run fitting mechanism on a laptop


### Manual Setup #####

Sim_no <- 1 # Choose which simulation to fit to (from table)
# Choose the fitting parameters
Adapt_start <- 1000 # start adaption
Adapt_stop <- 1000 # stop adaption
Length_run <- 100 # length of run



###### CONDITIONS  #######

simulation_overview <- read.csv("~/Documents/GitHub/interaction_model/Data/Multiple_sims_nov.csv")
Sig_value <- simulation_overview[simulation_ID, "Sigma"] # true value 
Rho_value <- simulation_overview[simulation_ID, "Rho"] # true value 

seed <- 301
#Define number of agegroups
Age_groups <- c(0,2,5,16,65)
#If using other age groups need to change summary_stats function in "Params_and_summary" file

###############Load required functions #############

#print(getwd())
#setwd("home/lsh1402815")#/Core_INF_RSV_MODEL")
#load packages
source("Setup.R")
#Load in mcmcMH_function
source("mcmcMH_S.R")
#Load paramaters and summary statistics function
source("Params_and_summary_for_MCMC.R")

set.seed(seed)
########## LOAD "REAL" DATA ##################
#setwd("home/lsh1402815/fitting_si/mulations")

observed.data <-read.csv(paste0("~/Documents/Github/interaction_model/Data/Simulation",Sim_no,"_",seed_no))

############# SIPR MODEL DETAILS ############


SIPR<-list()

SIPR$name <-"SIPR influenza and RSV"

Names <- c("time")
for(i in 1:length(Age_groups)){
  x<-c(paste0("SS", i), paste0("IS", i),paste0("PS", i),paste0("RS", i),
       paste0("SI", i),paste0("II", i),paste0("PI", i),paste0("SP", i),
       paste0("IP", i), paste0("SR", i),paste0("RR", i),paste0("Rcases", i),paste0("Icases", i))
  Names <- c(Names, x)
}

SIPR$state.names <- Names

SIPR$theta.names <- c("bR", "bI", "l_sig",'Rdetect2', 'Rdetect5', 'Idetect', 'l_rho', 'seedR', 'seedI')


############# MODEL SIMULATION ##############

SIPR$simulate <- function (theta, init.state, times)
{
  outall <- as.data.table(ode(y=init.state,
                              t=times,
                              func = derivatives,
                              parms = theta,
                              method="ode23"))
  
  trajectory <- summary_stats(outall)
  
  trajectory[, Hos_RSV_0 := RSV_incidence_1*as.numeric(theta["Rdetect2"])]
  trajectory[, Hos_INF_0 := INF_incidence_1*as.numeric(theta["Idetect"])]
  trajectory[, Hos_RSV_1 := RSV_incidence_2*as.numeric(theta["Rdetect5"])]
  trajectory[, Hos_INF_1 := INF_incidence_2*as.numeric(theta["Idetect"])]
  #Sum of children
  # trajectory[,Hos_c_RSV := (Hos_RSV_0 + Hos_RSV_1)]
  # trajectory[,Hos_c_INF := (Hos_INF_0 + Hos_INF_1)]
  # trajectory[,RSV_incidence_child := (RSV_incidence_0 + RSV_incidence_1)]
  # trajectory[,INF_incidence_child := (INF_incidence_0 + INF_incidence_1)]
  return(trajectory)
}

################ GENERATES RANDOM SAMPLE ##################

# SIPR$rPointObs<- function (model.point, theta)
# {
#   #
#   obs.point$RSV <- rpois(n = 1, lambda = model.point[["RSV_incidence_child"]])
#   obs.point$INF <- rpois(n = 1, lambda = model.point[["INF_incidence_child"]])
#   return(c(obs = obs.point))
# }



############### PRIORS ################

SIPR$dprior<-function (theta, log = FALSE)
{
  # log.prior.bR <- dunif(as.numeric(theta["bR"]), min = 0, max = 1,
  #                      log = TRUE)
  # log.prior.bI <- dunif(as.numeric(theta["bI"]), min = 0, max = 1,
  #            log = TRUE)
  log.prior.bR <- dnorm((exp(as.numeric(theta["bR"]))*as.numeric(theta["R0ratR"])), mean=3, sd=0.6, log=T)
  log.prior.bI <- dlnorm(((exp(as.numeric(theta["bI"]))*as.numeric(theta["R0ratI"])-1)), meanlog=0.3, sdlog=0.5,log=T)
  log.prior.Rdetect2 <- dunif(as.numeric(theta["Rdetect2"]), min=0, max=1, log=TRUE)
  log.prior.Rdetect5 <- dunif(as.numeric(theta["Rdetect5"]), min=0, max=1, log=TRUE)
  log.prior.Idetect <- dunif(as.numeric(theta["Idetect"]), min=0, max=1, log=TRUE)
  log.prior.seedI <- dunif(as.numeric(theta["seedI"]), min=0, max=60,log=T)
  log.prior.seedR <- dunif(as.numeric(theta["seedR"]), min=0, max=60,log=T)
  log.prior.sigma <- dunif(as.numeric(theta["l_sig"]), min=0, max=1,log=T)
  log.prior.rho <- dunif(as.numeric(theta["l_rho"]), min=0.01, max = 1,log=T)
  
  
  log.sum <- log.prior.bR + log.prior.bI +
    log.prior.Rdetect2 + log.prior.Rdetect5 + log.prior.Idetect +
    log.prior.seedR + log.prior.seedI +log.prior.sigma + log.prior.rho
  return(ifelse(log, log.sum, exp(log.sum)))
}



########### LIKLIHOOD OF EACH DATA POINT ##########


SIPR$dPointObs_RSV<- function (data.point, model.point, theta, log = FALSE)
{
  l<- dpois(x = as.numeric(data.point[["Hos_RSV_0"]]), lambda = model.point[["Hos_RSV_0"]],
            log = log)
  m<- dpois(x = as.numeric(data.point[["Hos_RSV_1"]]), lambda = model.point[["Hos_RSV_1"]],
            log = log)
  
  return(l+m)
}

SIPR$dPointObs_INF<- function (data.point, model.point, theta, log = FALSE)
{

  l<- dpois(x = as.numeric(data.point[["Hos_INF_0"]]), lambda = model.point[["Hos_INF_0"]],
            log = log)
  
  
  m<- dpois(x = as.numeric(data.point[["Hos_INF_1"]]), lambda = model.point[["Hos_INF_1"]],
            log = log)
  
  return(l+m)
}




################# LIKLIHOOD OF THE TRAJECTORY ################

#Combine the  liklehood of each of the data points

dTrajObs<- function (fitmodel, theta, init.state, data, log = FALSE)
{
  times <- seq(from= 1, to=365, by = 1)
  
  traj <- fitmodel$simulate(theta, init.state, times)
  dens_INF <- 0
  dens_RSV <- 0
  ######
  #####
  #####
  ##### / 10 added below
  for (i in 1:(nrow(data))) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj[i, ])
    dens_INF <- dens_INF + fitmodel$dPointObs_INF(data.point = data.point,
                                                  model.point = model.point, theta = theta, log = TRUE)
    dens_RSV <- dens_RSV + fitmodel$dPointObs_RSV(data.point = data.point,
                                                  model.point = model.point, theta = theta, log = TRUE)
    dens <- dens_RSV + dens_INF
    
  }
  return(ifelse(log, dens, exp(dens)))
}

############## CALCULATES LOG POSTERIOR ################

my_dLogPosterior <- function(fitmodel, theta, init.state, data) {
  
  # theta["l_sig"] <- as.numeric(theta["l_sig"])/2+as.numeric(theta["Idetect"])
  # theta["seedI"]<- as.numeric(theta["seedI"])/2 + as.numeric(theta["l_rho"])
  
  theta["l_rho"] <- exp(as.numeric(theta['l_rho']))
  theta["Rdetect2"] <-  exp(as.numeric(theta['Rdetect2']))/(1+exp(as.numeric(theta['Rdetect2'])))
  theta["Rdetect5"] <- exp(as.numeric(theta['Rdetect5']))/(1+exp(as.numeric(theta['Rdetect5'])))
  theta["Idetect"] <- exp(as.numeric(theta['Idetect']))/(1+exp(as.numeric(theta['Idetect'])))
  
  
  # calculate the  prior for parameter vector theta and assign to variable log.prior
  log.prior <- fitmodel$dprior(theta, log = TRUE)
  
  # calculate the log-likelihood of `theta`
  # and `init.state` with respect to the data using `dTrajObs`
  # and assign to a variable `log.likelihood`
  log.likelihood <- dTrajObs(fitmodel, theta, init.state, data, log = TRUE)
  # calulate the log-posterior using the log-prior and log-likelihood
  log.posterior <- log.prior + log.likelihood
  
  return(log.posterior)
  
}

########### SET UP AND RUN MCMC ##################

#Set up initial conditions
Log.posterior.test <- function(input_params) {
  
  
  theta <- list(Rdetect2=unname(input_params["Rdetect2"]), Rdetect5= unname(input_params["Rdetect5"]),
                Idetect = unname(input_params["Idetect"]),
                l_sig=unname(input_params["sigma"]), gammaR = gammaR, gammaI = gammaI,  RSV_Sus = RSV_Sus,
                trickleI = 0.1,
                Contact_Structure = Contact_structure,
                num_grps = length(Age_groups),
                l_rho=unname(input_params["rho"]),
                bR=unname(input_params["betaR"]), bI=unname(input_params["betaI"]),
                seedI = unname(input_params["seedI"]), seedR =unname(input_params["seedR"]),
                R0ratI = R0ratI, R0ratR= R0ratR # ratio between R0 and beta
  )
  
  return(my_dLogPosterior(fitmodel = SIPR,
                          theta = theta,
                          init.state = initials,
                          data = observed.data))
}

# ##Try one sample
# test<-c(betaR=-3.148541, betaI=-2.738075,Rdetect2= -5.506541, Rdetect5=-6.919733, Idetect=-6.246802
#         , seedR=1.038593, seedI= 13.00431, sigma=0.8318222, rho=0.03144764
# )
# Log.posterior.test(input_params = test)


init.theta<-c(betaR=-3.15, betaI=-2.7,Rdetect2= -5.5, Rdetect5=-6.9, Idetect=-6.2
              ,rho=-2, seedR=1, seedI= 10, sigma=0.7
)
proposal.sd <-c(betaR=0.001, betaI=0.001 , Rdetect2= 0.01, Rdetect5=0.001, Idetect=0.01
                , seedR=0.001, seedI=0.1, sigma=0.05, rho=0.1
)
limits.lower<- c(betaR=-Inf, betaI=-Inf, Rdetect2= -Inf, Rdetect5=-Inf, Idetect=-Inf
                 , seedR =0, seedI=0, sigma=0, rho=log(1/365)
)
limits.upper <-c(betaR=Inf, betaI=Inf, Rdetect2= Inf, Rdetect5=Inf, Idetect=Inf
                 ,  seedR= 90, seedI=90, sigma=1, rho=0
)


#run MCM
trace <- mcmcMH(target = Log.posterior.test, # target distribution
                init.theta = init.theta , # intial parameter guess
                proposal.sd = proposal.sd, # standard deviation of
                limits = list(lower=limits.lower ,
                              upper=limits.upper),
                #    adapt.size.start = 100,
                #    adapt.size.cooling = 0.999,
                adapt.shape.start = as.numeric(Adapt_start),
                adapt.shape.stop=as.numeric(Adapt_stop),
                n.iterations = as.numeric(Length_run)
                # covmat = covmat_1
)

print(warnings())

my_trace <- mcmc(trace[[1]])
#plot traces and densities
# require(Cairo)
# Cairo(2000,2000,file="test1.png", type="png", bg="white")
# plot(my_trace, type = "l")
# dev.off()

#covmat_1 <- trace$covmat.empirical


trace[["init.theta"]] <- init.theta
trace[["proposal.sd"]] <- proposal.sd
trace[["limits.lower"]] <- limits.lower
trace[["limits.upper"]] <- limits.upper


#Store the data for fitting to later
#setwd("home/lsh1402815/fitting_simulations")
save(trace, file = paste0("rerun_nov_",Sim_no,"_",seed_no, "_Seed_", seed, ".Rdata"))

