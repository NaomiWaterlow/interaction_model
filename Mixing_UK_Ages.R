data(polymod)
#Get desired contact structure
UK_structure <- contact_matrix(polymod,
                               countries = "United Kingdom",
                               age.limits = c(Age_groups),
                               symmetric = T)
# Symmetrical contact matrix
Contact_structure_temp <- as.matrix(UK_structure[[1]])
for(i in 1:length(Age_groups)){
  Contact_structure_temp[, i] <- Contact_structure_temp[, i] /
    as.numeric(UK_structure[[2]][i, "population"])
}
# rename for later use
Contact_structure <- Contact_structure_temp
