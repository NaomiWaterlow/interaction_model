# R0 Calculations - assuming no interaction

# set up betas for calculating R0 ratio for priors
bI_calc_R0 <-1
bR_calc_R0 <- 1

# Influenza transmission matrix
ITransmission <- matrix(nrow = length(Pop_sizes), ncol = length(Pop_sizes))

for (i in 1:length(Pop_sizes)){
  for (j in 1:length(Pop_sizes)) {
    ITransmission[i, j] <- length(Pop_sizes) * Contactmatrix[i, j]*Pop_sizes[j]
  }
}

#Transition matrix
ITransition <- matrix(0, nrow = length(Pop_sizes), ncol = length(Pop_sizes) )
diag(ITransition) <- gammaI

#Inverse 
ITransition_inverse <- ginv(ITransition)
ITransition_inverse

#NGM
INGM <- ITransmission %*% ITransition_inverse
IEigen <- unlist((eigen(INGM)[1]))
R0_INF <- max(IEigen)

# ratio
R0ratI <-  max(IEigen)


# RSV 
# Transmission matrix
RTransmission <- matrix(nrow = length(Pop_sizes), ncol = length(Pop_sizes))

for (i in 1:length(Pop_sizes)){
  for (j in 1:length(Pop_sizes)) {
    RTransmission[i, j] <- RSV_Sus[j] * bR_calc_R0 * Contactmatrix[i, j] * Pop_sizes[j]
  }
}

#Transition matrix
RTransition <- matrix(0,nrow = length(Pop_sizes), ncol = length(Pop_sizes) )
diag(RTransition) <- gammaR

#Inverse 
RTransition_inverse <- ginv(RTransition)

#NGM
RNGM <- RTransmission %*% RTransition_inverse
REigen<- unlist((eigen(RNGM)[1]))
# ratio
R0ratR <- max(REigen)
