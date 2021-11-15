# The Birth-Death-Chronospecies (BDC) model

### For info on how to prepare input files check these tutorials:
* [Generate PyRate input file (option 1)](https://github.com/dsilvestro/PyRate/blob/master/tutorials/pyrate_tutorial_1.md#generate-pyrate-input-file-option-1)  
* [Generate PyRate input file (option 2)](https://github.com/dsilvestro/PyRate/blob/master/tutorials/pyrate_tutorial_1.md#generate-pyrate-input-file-option-2)  



### Set up the PyRate analysis
Fossil and phylogenetic data can be jointly analyzed under the BDC model as described by [Silvestro, Warnock, Gavryushkina & Stadler 2018](https://www.nature.com/articles/s41467-018-07622-y). This analysis requires two input files: a standard PyRate input dataset (that can be generated as explained [here](https://github.com/dsilvestro/PyRate/blob/master/tutorials/pyrate_tutorial_1.md#generate-pyrate-input-file-option-2); see also [examples](https://github.com/dsilvestro/PyRate/tree/master/example_files/BDC_model)) and a tree file in NEXUS format.

To run a joint analysis of fossil and phylogenetic data assuming **independent rates** you can use:

`python PyRate.py example_files/BDC_model/Feliformia.py -tree example_files/BDC_model/Feliformia.tre`

Note that this function requires the [Dendropy library](https://dendropy.org). 
By default the analysis assumes a constant rate birth-death model with independent rate parameters between fossil and phylogenetic data. 
The flag `-bdc` enforces **compatible speciation and extinction rates under the BDC model**

`python PyRate.py example_files/BDC_model/Feliformia.py -tree example_files/BDC_model/Feliformia.tre -bdc`

You can also run an analysis where the rates are constrained to be equal using the flag `-eqr`:

`python PyRate.py example_files/BDC_model/Feliformia.py -tree example_files/BDC_model/Feliformia.tre -bdc`


To run under the **BDC skyline model** you can use the `-fixShift` command as explained [here](https://github.com/dsilvestro/PyRate/blob/master/tutorials/pyrate_tutorial_1.md#speciation-and-extinction-rates-within-fixed-time-bins), for example:

`python PyRate.py example_files/BDC_model/Feliformia.py -tree example_files/BDC_model/Feliformia.tre -fixShift example_files/epochs.txt -bdc`

This command sets up an analysis with rate shifts at the epochs boundaries under the BDC compatible model.


### Skyline model
To run under the **BDC skyline model** you can use the `-fixShift` command as explained [here](https://github.com/dsilvestro/PyRate/blob/master/tutorials/pyrate_tutorial_1.md#speciation-and-extinction-rates-within-fixed-time-bins), for example:

`python PyRate.py example_files/BDC_model/Feliformia.py -tree example_files/BDC_model/Feliformia.tre -fixShift example_files/epochs.txt -bdc`

This command sets up an analysis with rate shifts at the epochs boundaries under the BDC compatible model.

### Plot the results (in R)
Install dependencies and source the [plot\_functions\_BDC_model.R](https://github.com/dsilvestro/PyRate/blob/master/plot_functions_BDC_model.R) script:

```
install.packages("HDInterval")
source("path_to_pyrate/PyRate/plot_functions_BDC_model.R")
```

Provide the time bins used for the skyline model and the log file for plotting:

```
setwd("your_path")
time_bins <- read.table("epochs.txt", h=F)$V1
f <- "feliformia_1_BDC_mcmc.log"
plot_speciation_mode(f,time_bins)
```


