library("hybriddetective")
library("genepopedit")


###*** SET WORKING DIRECTORY HERE ***###

path.hold <- getwd()


##############################################################################

###*** SETTINGS ***###

# Set combined outfile from STEP 6
combined_outfile <- "EAGU_combined.txt"

# Set input file directory
input <- "/all_input/"

ind_file <- paste0("/", combined_outfile, "_individuals.txt")
data <- paste0("/", combined_outfile)

## Copy the individual file to the new folder
file.copy(from = paste0(path.hold, ind_file), to = paste0(path.hold, input))

## Copy the genotype data file to the new folder
file.copy(from = paste0(path.hold, data), to = paste0(path.hold, input))

## Create two copies of the genotype data file to act as technical replicates of the NewHybrids simulation based analysis. This will also serve demonstrate the parallel capabilities of parallelnewhybrid.
file.copy(from = paste0(path.hold, input, data), to = paste0(path.hold, input, data2))
file.copy(from = paste0(path.hold, input, data), to = paste0(path.hold, input, data3))
file.copy(from = paste0(path.hold, input, data), to = paste0(path.hold, input, data4))
file.copy(from = paste0(path.hold, input, data), to = paste0(path.hold, input, data5))

# SET WORKING DIRECTORY HERE #

new_path <- getwd()

nh_plotR(paste0(new_path, "/final_run/postSim_combined/NH.Results/EATT_combined_Zed.txt_Results/EATT_combined_Zed.txt_PofZ.txt"))
