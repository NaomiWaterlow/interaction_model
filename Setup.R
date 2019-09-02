# set location to install packages
storage <- "~/Documents/R"

#list of packages to load
packages <- c(
  "Rcpp"
  ,"deSolve"
  , "gridExtra"
  ,"plyr"
  , "data.table"
  , "ggplot2"
  , "reshape2"
  , "socialmixr"
  ,"MASS"
  ,"stats4"
  ,"bbmle"
  , "coda"
  ,"tmvtnorm"
  ,"fitR"
  #,"GGally"
)

#install an individual package
# install.packages("fitR", lib = storage)

# #Install packages (only needed once!)
# for (i in packages){
#   install.packages("Rcpp", lib = storage)
# }

#load the required packages
for (i in packages){
  library(i, character.only = T, lib.loc = storage)
}

