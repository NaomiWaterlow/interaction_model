# This script contains the fixed parameters and some function for MCMC fitting 

#Load contact matrix
setwd("~/Documents/GitHub/interaction_model")
#UK Polymod contact matrix
source("Mixing_UK_Ages.R") 

# Contact structue
RSV_Sus <- c(1, 0.75, rep(0.65, (length(Age_groups) - 2)))
Pop_sizes <- unname(unlist(UK_structure$demography[, "population"]))
Start_Sus <- c(1, 0.688, rep(0.525, (length(Age_groups) - 2)))

#Contract matrix - from socialmixr UK polymod
Contactmatrix <- Contact_structure
# recovery rates
gammaR <-0.1111
gammaI <- 0.2632

# to calculate the R0s
setwd("~/Documents/GitHub/interaction_model")
source("beta_priors_r0_calc.R")

initials <- c()
for (i in 1:length(Age_groups)) {
  temp = c((Pop_sizes[i] * Start_Sus[i]-1),
           #SS
           1,0, 0,0,0,0, 0,0,(Pop_sizes[i] * (1 - Start_Sus[i])),
           #IS,PS,RS,SI,II,PI,SP,IP,SR
           0, 0, 0)
  #,RR,R_cases, I_cases
  initials <- c(initials, temp)
}

# column names
Names <- c("time")
for (i in 1:length(Age_groups)) {
  x <-c(
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


#  calculate summary statistics with 5 ages groups
summary_stats <- function(outall) {
  
  #name the columns
  colnames(outall) <- Names
  # calculate the incidence from cumulative cases
  outall[, paste0("RSV_incidence_", 1) := (get(paste0("Rcases", 1)) - shift(get(paste0("Rcases", 1)), 1L, type = "lag"))]
  outall[, paste0("INF_incidence_", 1) := (get(paste0("Icases", 1)) - shift(get(paste0("Icases", 1)), 1L, type = "lag"))]
  outall[1, paste0("RSV_incidence_", 1)] = 0
  outall[1, paste0("INF_incidence_", 1)] = 0
  outall[, paste0("RSV_incidence_", 2) := (get(paste0("Rcases", 2)) - shift(get(paste0("Rcases", 2)), 1L, type = "lag"))]
  outall[, paste0("INF_incidence_", 2) := (get(paste0("Icases", 2)) - shift(get(paste0("Icases", 2)), 1L, type = "lag"))]
  outall[1, paste0("RSV_incidence_", 2)] = 0
  outall[1, paste0("INF_incidence_", 2)] = 0
  outall[, paste0("RSV_incidence_", 3) := (get(paste0("Rcases", 3)) - shift(get(paste0("Rcases", 3)), 1L, type = "lag"))]
  outall[, paste0("INF_incidence_", 3) := (get(paste0("Icases", 3)) - shift(get(paste0("Icases", 3)), 1L, type = "lag"))]
  outall[1, paste0("RSV_incidence_", 3)] = 0
  outall[1, paste0("INF_incidence_", 3)] = 0
  
  #Convert to weekly
  outall$time <- as.Date(outall$time, origin = "2014-09-01")
  outall$week_begin <- as.Date(cut(outall$time, breaks = "week"))
  outall_weekly2 <- melt(outall, id = "week_begin")
  outall_weekly2$value <- as.double(outall_weekly2$value)
  out_weekly <-as.data.table(dcast(outall_weekly2, week_begin ~ variable, sum))
  out_weekly <- out_weekly[-(dim(out_weekly))[1], ]
  
  #return as data table
  return(as.data.table(out_weekly))
  
}