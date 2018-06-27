
get.p.values <- function(z) {
  
  pvalue2sided <- 2*pnorm(-abs(z))
  
}

get.range <- function(r) {
  
  colrange <- range(r)
  
}

setwd("C:/Users/evobi/Documents/work/BOX_Turtle_Genomics/Analyses/Dtest/lane4/all_output/")

filenames <- list.files(pattern="*.txt")

testdata <- read.delim(filenames[1], header=TRUE, sep="\t")

testrange <- get.range(testdata[[12]])

#print(paste(testrange[2], testrange[1], sep="\t"))

testP <- get.p.values(testdata[[12]])
testsig <- sum(testP < 0.05)

print(testsig)
rs <- nrow(testdata)

paste(testsig, rws, sep="/")

bon <- 0.05 / rs

print(bon)

for (files in filenames) {
  
  data <- read.delim(files, header=TRUE, sep="\t")
  rws <- nrow(data)
  
 #print(data[,12])
  newP <- get.p.values(data[[12]])
  #write.table(newP, file=paste(files, ".testP", sep=""), sep="\t")
  
  zrange <- get.range(data[[12]])
  
  sigtests <- sum(newP <= 0.05)
  sig_notsig <- paste(sigtests, rws, sep="/")
  
  bonferronni <- 0.05 / rws
  
  bonsig <- sum(newP <= bonferronni)
  bonsig_notsig <- paste(bonsig, rws, sep="/")
  
  
  
  
  
#}

#for (i in 1:length(temp)){
  
#  assign(temp[i], read.delim(temp[i], header=TRUE, sep="\t"))
  
#}
#lof <- lapply(filenames, read.delim, header=TRUE, sep="\t")


#for (i in 1:length(lof)){
  
  #newP <- lapply(lof,get.p.values(lof[[i]][,12])
  
  #newdf <- cbind(lof[[i]],newP)
  #print(lof[[i]][,12])
  
}

#lof[[1]]$v12

#lapply(lof, "[",12,drop=FALSE)