# download fossilBM here: https://github.com/dsilvestro/fossilBM
# install dependencies (if necessary)
# install.packages(c("phylobase","adephylo","geiger","phytools","HDInterval"))

# Data from Zhang et al. 2021 Syst Biol (doi: 10.1093/sysbio/syab030)
# define input files
setwd("path") # set working directory
treefile <- "Juglandaceae_trees.tre"
tindex <- 1 # specify which tree should be analyzed if the tree file includes multiple trees
datafile <- "climate_data/Tcold.txt" # estimated temperature of the coldest month

# load FBM library
source("path_to_fossilbm/fossilBM/fossilBM_lib.R")


# set up the fbm input object
fbm_obj <- read_and_transform_data(treefile, 
                                   datafile, 
                                   rescale_trait_data=0.1,  # rescale trait for better convergence
                                   partition_file="Juglandoideae_partition.txt" # define subfamily partition
                                   )
                                   # Optional arguments:
                                   # log_trait_data = 10 # log10 transforme data
                                   # drop_na = FALSE # if set to TRUE NAs will be dropped otherwise they are inputed
                                   # rescale_trait_data = 1 # multiplier to rescale the trait data
                                   # root_calibration = c(0,100) # normal prior (mean,sd) on root state




# run MCMC
output_file <- "Tcold_output.log"
run_mcmc(fbm_obj, 
         logfile=output_file,     
         ngen = 2500, # you will need to increase this number to reach convergence!
         sample = 50, 
         linTrend=TRUE
         )

# plot results
plot_results(fbm_obj, output_file, resfile = "Tcold_summary.pdf")
res <- plot_time_varying_trend(fbm_obj, output_file, resfile="Tcold_trends.pdf")


