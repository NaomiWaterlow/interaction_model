
# ########### SETUP & PARAMETERSS ########

# Choose which simulation to create
simulation_ID <- 1 # edit manually

simulation_overview <- read.csv("~/Documents/GitHub/interaction_model/Data/Multiple_sims_nov.csv")
set.seed(simulation_overview[simulation_ID, "Seed"])
#Define age groups wanted, Need to include at least 0,2,5
Age_groups <- c(0, 2, 5, 16, 65)

# Load the packages
setwd("~/Documents/GitHub/interaction_model")
source("Setup.R")

#Make mixing matrix using UK Polymod contact matrix
setwd("~/Documents/GitHub/interaction_model")
source("Mixing_UK_Ages.R")
# Source the model
setwd("~/Documents/GitHub/interaction_model")
sourceCpp("Flex_age_model_trickle.cpp")

#Age group specific variables
RSV_Sus <- c(1, 0.75, rep(0.65, (length(Age_groups) - 2)))
Pop_sizes <- unname(unlist(UK_structure$demography[, "population"]))
Start_Sus <- c(1, 0.688, rep(0.525, (length(Age_groups) - 2)))

#Fixed parameters
times <- seq(from = 1, to = 365, by = 1)

Parameters <- list(
  bR = log(0.043),
  bI = log(0.063),
  l_sig = simulation_overview[simulation_ID, "Sigma"],
  gammaR =  0.1111,
  gammaI = 0.2632,
  l_rho = simulation_overview[simulation_ID, "Rho"],
  Rdetect2 = 0.004,
  Rdetect5 = 0.001,
  Idetect = 0.002,
  RSV_Sus = RSV_Sus,
  Contact_Structure = Contact_structure,
  num_grps = length(Age_groups),
  seedI = 10,
  seedR = 1, 
  trickleI = 0.1
)

# set up the initial donditions
initials <- c()
for (i in 1:length(Age_groups)) {
  temp = c((Pop_sizes[i] * Start_Sus[i] - 1),
           #SS
           1,0, 0, 0,0,0, 0,0,(Pop_sizes[i] * (1 - Start_Sus[i])),
           #IS,PS,RS,SI,II,PI,SP,IP,SR
           0, 0, 0
           #    ,0,0  #for tracking intros
  )
  #,RR,R_cases, I_cases
  initials <- c(initials, temp)
}

######### RUN THE MODEL #######

outall <- as.data.table(ode(
  y = initials,
  t = times,
  func = derivatives,
  parms = Parameters,
  method = "ode23"
))


###### MANIPULATE OUTPUT ######

#Name columns
Names <- c("time")
for (i in 1:length(Age_groups)) {
  x <- c(
    paste0("SS", i),
    paste0("IS", i),
    paste0("PS", i),
    paste0("RS", i),
    paste0("SI", i),
    paste0("II", i),
    paste0("PI", i),
    paste0("SP", i),
    paste0("IP", i),
    paste0("SR", i),
    paste0("RR", i),
    paste0("Rcases", i),
    paste0("Icases", i) 
  )
  Names <- c(Names, x)
}
colnames(outall) <- Names

#Calculate incidence per time step, from cumulative incidence.
for (i in 1:length(Age_groups)) {
  outall[, paste0("RSV_incidence_", i) :=
           (get(paste0("Rcases", i)) - shift(get(paste0("Rcases", i)), 1L, type = "lag"))]
  outall[, paste0("INF_incidence_", i) := 
           (get(paste0("Icases", i)) - shift(get(paste0("Icases", i)), 1L, type = "lag"))]
  outall[1, paste0("RSV_incidence_", i)] = 0
  outall[1, paste0("INF_incidence_", i)] = 0
}

#Convert to weekly
outall$time <- as.Date(outall$time, origin = "2014-10-01")
outall$week_begin <- as.Date(cut(outall$time, breaks = "week"))
outall_weekly2 <- melt(outall, id = "week_begin")
outall_weekly2$value <- as.double(outall_weekly2$value)
out_weekly <- as.data.table(dcast(outall_weekly2, week_begin ~ variable, sum))
out_weekly <- out_weekly[-(dim(out_weekly))[1], ]

#Stocahstic observation process - only youngest two age categories
out_weekly[, Hos_RSV_0 := rbinom(dim(out_weekly)[1],round(RSV_incidence_1),
                                 as.numeric(Parameters["Rdetect2"]))]
out_weekly[, Hos_INF_0 := rbinom(dim(out_weekly)[1],round(INF_incidence_1),
                                 as.numeric(Parameters["Idetect"]))]

out_weekly[, Hos_RSV_1 := rbinom(dim(out_weekly)[1],round(RSV_incidence_2),
                                 as.numeric(Parameters["Rdetect5"]))]
out_weekly[, Hos_INF_1 := rbinom(dim(out_weekly)[1],round(INF_incidence_2),
                                 as.numeric(Parameters["Idetect"]))]

#Save the data
out_weekly[,(2:(15 * length(Age_groups) + 2)) := NULL]
simulation.data <- out_weekly
write.table((simulation.data), file = paste0("~/Documents/GitHub/interaction_model/Data/Simulation",
                                             simulation_overview[simulation_ID,"Run"],
                                             "_",
                                             simulation_overview[simulation_ID,"Seed"]),
            sep = ",")
