
get.p.values <- function(z) {
  
  pvalue2sided <- 2*pnorm(-abs(z))
  
}

get.range <- function(r) {
  
  colrange <- range(r)
  
}


filenames <- list.files(pattern="*.txt")
zfilenames <- list.files(pattern="*.z")


popPvalue=list()
popz=list()

for (i in 1:length(zfilenames)) {
  
  zfiles <- read.delim(zfilenames[i], header=TRUE, sep="\t")
  
  Pvalue <- get.p.values(zfiles[[2]])
  popPvalue[i] <- Pvalue
  avgz <- zfiles[[2]]
  popz[i] <- avgz
  
}

df_total <- data.frame()

for (j in 1:length(filenames)) {
  
  test <- tools::file_path_sans_ext(filenames[j])
  
  dataset <- read.delim(filenames[j], header=TRUE, sep="\t")
  rws <- nrow(dataset)
  
  newP <- get.p.values(dataset[[12]])

  zrange <- get.range(dataset[[12]])
  Z <- paste0("[", zrange[1], ",", zrange[2],"]")
  
  drange <- get.range(dataset[[8]])
  D <- paste0("[", drange[1], ",", drange[2], "]")
  
  
  sigtests <- sum(newP <= 0.05)
  sig_notsig <- paste(sigtests, rws, sep="/")
  
  bonferroni <- 0.05 / rws
  
  bon_sig <- sum(newP <= bonferroni)
  bonsig_notsig <- paste(bon_sig, rws, sep="/")
  
  df <- data.frame(test, D, Z, sig_notsig, bonsig_notsig, popz[j], popPvalue[j], stringsAsFactors = FALSE)
  
  colnames(df) <- c("Test", "D-Range", "Z-Range", "No.Significant", "Bonferroni", "Population-Z", "Pop_P-value")
  row.names(df_total) <- NULL
  
  df_total <- rbind(df_total, df)
  
}

ofile <- "summary.tsv"

write.table(df_total, file=ofile, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
