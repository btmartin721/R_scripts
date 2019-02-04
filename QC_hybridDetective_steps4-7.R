# Need to install packages first. 
# If you need to install devtools, I hope you have sudo privledges. 
# Good luck.

# library("devtools")
# devtools::install_github("rystanley/genepopedit")
# devtools::install_github("bwringe/parallelnewhybrid")
# devtools::install_github("bwringe/hybriddetective")

library("parallelnewhybrid")
library("genepopedit")
library("hybriddetective")

###*** SET THE WORKING DIRECTORY HERE ***###.

my_path <- getwd()

##############################################################################

###*** SETTINGS - SET THESE BEFORE RUNNING SCRIPT ***###

# Path to parallelNH results file.
results <- paste0(my_path, "/simRunsX10/NH.Results/")

# Path to directory above NH.Results.
simulations <- paste0(my_path, "/simRunsX10/")

# Set file with panel of top N SNPS (based on getTopLoc() function from step 1).
panelFile <- "EATT_purePops_300_Loci_Panel.txt"

# Get panelFile without .txt extension.
panel_noTXT <- substr(panelFile, 1, nchar(panelFile)-4)

# Set initial GENEPOP file (from before running step 1, getTopLoc())
GenePopData <- "BOX_genepop_EATT.txt"

# Get all population IDs.
PopNames <- genepop_detective(GenePopData, "Pops")
print(PopNames)

# Set admixed/unknown populations. This will get
# everything but the two pure populations listed below
unknownPops <- setdiff(PopNames, c("EA", "TT"))

# Set unknown output filename
# for subsetting loci and combining sim/experimental datasets
unknownFile <- "/EATT_unknown_dataset.txt"

# Pick sim and replicate to run through NH with full dataset.
simRep <- "_1_S1R1_NH.txt"

# To be input into nh_analysis_generateR
repToUse <- paste0(panel_noTXT, simRep)

# Filename of combined sim/experimental results.
combined_outfile <- "EATT_combined.txt"

# Get filename without .txt extension.
combined_outfile_noTXT <- substr(combined_outfile, 1, nchar(combined_outfile)-4)

# Get filename of individuals file. Has to end with _individuals.txt for newhybrids to work.
ind_file <- paste0("/", combined_outfile_noTXT, "_individuals.txt")


# Variable name is lower case z. Set file for nh_zcore() to insert known pop 
# column into NH file
# File MUST be CSV that has header line with "Individual" and "Zscore" as 
# column names.
zfile <- "EATT_combined_Zvec.csv"

newhybsDIR <- "/EATT_FINALrun/"

# Check results for convergence
nh_preCheckR(PreDir = results)

# If previous step looks good: Plot the PofZ files output from newhybrids
nh_multiplotR(NHResults = results)

# Generate data and plots to visualize accuracy and efficacy
# of simulations to identify hybrid generations.

hybridPowerComp(dir = results)

# Get locus names for the top N loci (selected in getTopLoci(), step 1)
loci_subset <- genepop_detective(panelFile, variable = "Loci")

# Subset full dataset to get only to N loci from getTopLoci
subset_genepop(genepop = GenePopData, keep = TRUE, subs = loci_subset, 
               path = paste0(my_path, unknownFile))

nh_analysis_generateR(ReferencePopsData = paste0(simulations, repToUse),
                      UnknownIndivs = paste0(my_path, unknownFile), 
                      output.name = combined_outfile)

# Make directory to put combined_outfiles for all analyses.
dir.create(file.path(my_path, newhybsDIR))

# Copy combined_outfile newhybs input directory
file.copy(from = paste0(my_path, "/", combined_outfile), 
          to = paste0(my_path, newhybsDIR, combined_outfile))

file.copy(from = paste0(my_path, ind_file), 
          to = paste0(my_path, newhybsDIR, ind_file))


# Lower case z for zfile. Upper case for Zcore.
nh_Zcore(GetstheZdir = paste0(my_path, newhybsDIR), multiapplyZvec = paste0(my_path, newhybsDIR, zfile))
