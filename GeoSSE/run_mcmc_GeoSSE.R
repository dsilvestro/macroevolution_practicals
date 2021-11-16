# Joint Bayesian estimation of speciation, dispersal and extinction rates: the GeoSSE model
# The GeoSSE model jointly infers area-specific speciation rates, dispersal rates and local extinction rates (Goldberg et al. 2011). 
# Download mcmc-diversitree here: https://github.com/dsilvestro/mcmc-diversitree
# The mcmc-diversitree program implements a Bayesian algorithm to infer these rates (Silvestro et al. 2014).

library("diversitree")
library("picante")
library("phytools")
library("vioplot")

# Load the mcmc-SEE-library
source("path_to_code/mcmc-diversitree/mcmc-SSE-lib.R")

setwd("practicals")

# Load phylogenetic tree
# Can contain multiple trees, i.e. class “multiPhylo” to account for phylogenetic uncertaineties
tree <- read.tree("Xenarthra_phylogeny.tre")

# Define geographic ranges. Geographic ranges are defined as 0 (present in area A), 1 (present in area B), 2 (present in area AB)
# here areas are (+/-) 0: Andes, 1: lowland South and Central America
range.data <- read.table(file = "Xenarthra_georange.txt", h = T, row.names=1)

# State-specific taxon sampling (A, B, AB) to tune depending on fraction of species included in the analyses in all areas
sampling_fractions = c(1,1,1)

# Define output file name
mcmc.log = "Xenarthra_geosse.log"

# Run the analysis (you can increase the number of iterations to run a longer analysis and improve convergence)
run_mcmc_SSE(tree, range.data, model = "geosse", outfile = mcmc.log, iterations = 50000, 
             rho = sampling_fractions)

# # Run on multiple trees
# run_mcmc_SSE(trees, range.data, model = "geosse", outfile = "geosse.log", iterations = 10000, 
#              nTREES = 5, burnin = 100)

### 3/ Extract outputs ####

post = read.table(mcmc.log, header = T)

# Plot trace of sampled likelihoods
plot(post$likelihood, type = "l")

# Plot speciation, dispersal, and extinction rates per area:
boxplot(post[, c("sA", "sB", "sAB")], col = "blue", main = "speciation rates")
boxplot(post[, c("dA", "dB")], col = "gray", main = "dispersal rates")
boxplot(post[, c("xA", "xB")], col = "red", main = "extinction rates")


# Are dispersals rates significantly different?
delta_dispersal = post$dA - post$dB
hist(delta_dispersal)

# Calculate the posterior probability of dAB > dBA = nb of steps in the posterior distribution in which
# the dispersal rate from A to B (dA) was higher than the dispersal rate from B to A (dB) 
dA_higher_dB <- length(delta_dispersal[delta_dispersal > 0])/length(delta_dispersal) ; dA_higher_dB


# Are extinction rates significantly different?
delta_extinction = post$xA - post$xB
hist(delta_extinction)

# Calculate the posterior probability of xA > xB = nb of steps in the posterior distribution in which
# the extinction rate in A (xA) was higher than the extinction rate in B (dB)
xA_higher_xB <- length(delta_extinction[delta_extinction > 0])/length(delta_extinction) ; xA_higher_xB


# Is allopatric speciation more common than overall sympatric speciation?
delta_speciation_type = post$sAB - (post$sA + post$sB)
hist(delta_speciation_type)

# Calculate the posterior probability of sAB > (sB + sA) = nb of steps in the posterior distribution in which
# the allopatric speciation rate (sAB) was higher than the sympatric speciation rate (sA + sB) 
sAB_higher_sA_sB <- length(delta_speciation_type[delta_speciation_type > 0])/length(delta_speciation_type) ; sAB_higher_sA_sB

# Is allopatric speciation more common than sympatric speciation in either area?

# Calculate the posterior probability of sAB > (sB + sA) = nb of steps in the posterior distribution in which
# the allopatric speciation rate (sAB) was higher than the sympatric speciation rate (sA + sB) 
sAB_higher_sA <- sum(post$sAB  >post$sA) / length(delta_speciation_type) 
sAB_higher_sB <- sum(post$sAB  >post$sB) / length(delta_speciation_type) 



# Are sympatric speciation rates significantly different between the regions?
delta_speciation = post$sA - post$sB
hist(delta_speciation)

# Calculate the posterior probability of sA > sB = nb of steps in the posterior distribution in which
# the speciation rate in A (sA) was higher than the speciation rate in B (sB)
sA_higher_sB <- length(delta_speciation[delta_speciation > 0])/length(delta_speciation) ; sA_higher_sB


# Plot speciation, dispersal, and extinction rates per area:
pdf(file = "GeoSSE_summary_plot.pdf", width = 15, height = 8)
par(mfrow = c(1, 3))
plot(area, col = area$bioregio)
legend(legend =c("Area A", "Area B"), fill = c("red", "limegreen"), x = "bottomleft")
vioplot(post[, c("sA", "sB", "sAB")], col = "blue", main = "speciation rates")
legend(legend = c(paste0("p(sA > sB) = ", sA_higher_sB),  paste0("p(sAB > sA) = ", sAB_higher_sA), paste0("p(sAB > sB) = ", sAB_higher_sB)),
       x = "top", cex = 1.3, bty = "o")
vioplot(post[, c("dA", "dB")], col = "gray", main = "dispersal rates")
legend(legend = paste0("p(dA > dB) = ", dA_higher_dB),
       x = "top", cex = 1.3, bty = "o")
vioplot(post[, c("xA", "xB")], col = "red", main = "extinction rates")
legend(legend = paste0("p(xA > xB) = ", xA_higher_xB),
       x = "top", cex = 1.3, bty = "o")
par(mfrow = c(1, 1))
dev.off()

