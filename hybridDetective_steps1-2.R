#devtools::install_github("rystanley/genepopedit")
#devtools::install_github("bwringe/hybriddetective")
#devtools::install_github("bwringe/parallelnewhybrid")

library("parallelnewhybrid")
library("hybriddetective")
library("genepopedit")


##*** SET WORKING DIRECTORY HERE ***###

# Store working directory as object.
my_path <- getwd()

###############################################################################

#*** SETTINGS ***###
#*** SET THE FOLLOWING OPTIONS *** #

# Name of GENEPOP input file
GenePopData <- "BOX_genepop_EATT.txt"

# Name of GENEPOP input file containing only two pure popualations.
# Must have .txt extension.
GenePopData2 <- "/EATT_purePops.txt"

# Store GenePopData2 file without as charactervector witout 
# ".txt" extension
without_txt <- substr(GenePopData2, 1, nchar(GenePopData2)-4)

## Select PURE populations to keep. Must be only 2 pops.
PopKeep <- c("EA", "TT")

# Set number of top loci (based on FST and LD) to keep.
panelSize <- 300

# Specify panel file output by getTopLoc()
panelFile <- paste0(without_txt, "_",
                    as.character(panelSize), "_Loci_Panel.txt")

# Get panelFile string without ".txt" extension.
panel_noTXT <- substr(panelFile, 1, nchar(panelFile)-4)

# Set path for plink and PGDspider directories.
# Both must be subdirectories within working directory
plinkDIR <- "/plink_win64_20181202/"
pgdDIR <- "/PGDSpider_2.1.1.5/PGDSpider_2.1.1.5/"

# Set the number of independent simulations
setNumSims <- 10

# Set the number of replicates for each independent simulation
setNumReps <- 3

# Set directory to store simulated files.
# Must be subdirectory within current working directory.
simDIR <- "simRunsX10"

#*** END OF SETTINGS *** #
###############################################################################

# Get population IDs from GENEPOP file using genepopedit. 
# Then print pop list.
pops <- genepop_detective(genepop=GenePopData)
print(pops)

## Subsets GENEPOP file to keep only pur populations.
subset_genepop(genepop = GenePopData, 
                               keep = TRUE, 
                               subs = NULL, 
                               spop = PopKeep, 
                               path = paste0(my_path, GenePopData2))

# Select top loci based on FST. Number of loci set by panel.size
getTopLoc(GPD = paste0(my_path, GenePopData2), panel.size=panelSize, 
          where.PLINK=paste0(my_path, plinkDIR), 
          where.PGDspider=paste0(my_path, pgdDIR))

# Mamma said Visualize the attack. Or some of the results. Either way. ;)
genepop_flatten(genepop = paste0(my_path, panelFile))[1:3, 1:5]

# Make a directory to store the newhybrids input files.
dir.create(file.path(my_path, simDIR))

# Copy panelFile to new directory and append "_1.txt to end of filename.
file.copy(from = paste0(my_path, panelFile), 
          to = paste0(my_path, "/", simDIR, panel_noTXT, "_1.txt"))

# Change to the new directory.
setwd(paste0(my_path, "/", simDIR))

# Remove forward slash from panel_noTXT.
panel_noTXT <- sub('.', '', panel_noTXT)

# Generate the simulation data.
freqbasedsim_AlleleSample(GPD = paste0(panel_noTXT, "_1.txt"), 
                          NumSims = setNumSims, 
                          NumReps = setNumReps)

# Remove panelFile from simDIR. Little bugger!
file.remove(paste0(my_path, "/", simDIR, "/", panel_noTXT, "_1.txt"))

print("You are now ready to run parallelnewhybrid")

###############################################################

###*** 
# YOU ARE NOW READY TO RUN PARALLELNEWHYBRID ON THE
# Simulated datasets.
###***


