library("hybriddetective")
library("genepopedit")


###*** SET WORKING DIRECTORY HERE ***###

path.hold <- getwd()

##############################################################################

###*** SETTINGS ***###

pofz <- "/final_run/postSim_combined/NH.Results/EAON_combined_Zed.txt_Results/EAON_combined_Zed.txt_PofZ.txt"


###*** END OF SETTINGS ***###

##############################################################################


nh_plotR(paste0(path.hold, pofz))
