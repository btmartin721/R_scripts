
get.range <- function(r) {
  
  colrange <- range(r)
  
}

# Sets number of digits for P-values.
# Can change to any integer.
p.digits = 10


# Filenames must end with .txt and .z (for popfiles)
filenames <- list.files(pattern="*.txt")
zfilenames <- list.files(pattern="*.z")


popPvalue=list()
popz=list()

# *.z files produced from CompD
for (i in 1:length(zfilenames)) {
  
  zfiles <- read.delim(zfilenames[i], header=TRUE, sep="\t")
  
  Pvalue <- zfiles[[3]]
  popPvalue[i] <- round(Pvalue, digits = p.digits)
  avgz <- round(zfiles[[2]], digits = 5)
  popz[i] <- avgz
  
}

df_total <- data.frame()

# *.txt file produced from CompD
for (j in 1:length(filenames)) {
  
  test <- tools::file_path_sans_ext(filenames[j])
  
  dataset <- read.delim(filenames[j], header=TRUE, sep="\t")
  rws <- nrow(dataset)
  
  newP <- dataset[[13]]

  zrange <- get.range(dataset[[12]])
  minz <- round(zrange[1], digits = 5)
  maxz <- round(zrange[2], digits = 5)
  Z <- paste0("[", minz, ",", maxz,"]")
  
  stdev_d <- sd(dataset[[8]])
  mean_d <- mean(dataset[[8]])
  d <- round(mean_d, digits = 5)
  STDEV_D <- round(stdev_d, digits = 5)
  D <- paste0(d, " (", STDEV_D, ")")
  
  # Chi-Square No. Significant/Not Significant tests
  sigChi <- sum(dataset[[11]] <= 0.05)
  Chi_sig_notsig <- paste(sigChi, rws, sep="/")
  
  # If P-value <= 0.05
  sigtests <- sum(newP <= 0.05)
  # Adds Z.Significant/NotSignificant tests to DataFrame
  sig_notsig <- paste(sigtests, rws, sep="/")
  
  # Bonferroni Correction from number of tests
  bonferroni <- 0.05 / rws
  
  # P-value <= Bonferroni alpha
  bon_sig <- sum(newP <= bonferroni)
  
  #Adds Sig/NotSig to DataFrame
  bonsig_notsig <- paste(bon_sig, rws, sep="/")
  
  df <- data.frame(test, D, Z, Chi_sig_notsig, sig_notsig, 
                   round(bonferroni, digits = p.digits), 
                   bonsig_notsig, popz[j], popPvalue[j], 
                   stringsAsFactors = FALSE)
  
  colnames(df) <- c("Test", "Mean.D (STDEV.D)", "Z-Range", 
                    "Chi-Square.Sig", "Z.Significant", 
                    "Bonferroni.Alpha", "Bonferroni.Significant", 
                    "Population-Z", "Pop_P-value")
  
  row.names(df_total) <- NULL
  
  df_total <- rbind(df_total, df)
  
}

# Outfile to write to. Will overwrite if script is run again.
ofile <- "summary.tsv"

write.table(format(df_total, digits = p.digits, scientific=F), 
            file=ofile, 
            quote=FALSE, 
            sep="\t", 
            col.names=TRUE, 
            row.names=FALSE)
