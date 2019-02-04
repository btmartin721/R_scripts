############################################################################

# This R script is meant to be run in conjunction with 
# hybriddetective (https://github.com/bwringe/hybriddetective).

# It assumes that the first two functions (steps) in the hybriddetective 
# pipeline, getTopLoci() and freqbasedsim_AlleleSample(), have already been run.

# The output files from getTopLoci() should be in a directory with nothing
# else in it. Also try not to get too crazy with paths (trust me).

############################################################################

#devtools::install_github("bwringe/parallelnewhybrid")


library("parallel")
library("plyr")
library("stringr")
library("tidyr")
library("parallelnewhybrid")

## Get the file path to the working directory
path.hold <- getwd()

#############################################################################

###*** SETTINGS - SET THESE BEFORE RUNNING ANALYSIS ***###

###*** SET THE DIRECTORY FOR THE ANALYSIS HERE ***###
### Must have forward slashes at beginning and end ###
popDir <- "/simRunsX10/"

###*** SET THE BURNIN AND SWEEPS (AFTER BURN-IN) HERE***###
bIn = 200000
swps = 1000000

## Create an object that is the file path to the folder in which 
## NewHybrids is installed. Note: this folder must be named "newhybrids"
my.NH <- "/home/btm002/hybrid_detective/newhybrids/"

###*** END OF SETTINGS ***###

##############################################################################

## Execute parallelnh. 
parallelnh_LINUX(folder.data = paste0(path.hold, popDir), 
                 where.NH = my.NH, burnin = bIn, sweeps = swps)
