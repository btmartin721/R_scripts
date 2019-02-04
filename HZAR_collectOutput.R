library("hzar")
library("extrafont")


###*** SET WORKING DIRECTORY HERE ***###

# Get list of .rds filenames in working directory
filenames <- as.list(dir(path = getwd(), pattern="*.rds"))

data <- lapply(filenames, readRDS)

windows()

# Plot cline for first locus
hzar.plot.cline(data[[1]], pch = "")

# Add the other loci to the plot with add = TRUE.
for(i in 2:length(data)){
  hzar.plot.cline(data[[i]], add = TRUE, pch = "")
}


dev.off()
