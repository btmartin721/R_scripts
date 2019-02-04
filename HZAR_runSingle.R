#!/usr/bin/env Rscript

# Get command-line arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Set j, which corresponds to the current locus
j <- as.integer(args[1])
#print(j)
##### Script for fitting clines using Hzar
##### This script uses output from my vcf2hzar script 
##### from my GitHub: https://github.com/btmartin721/vcf2hzar.git
##### Four files are produced from vcf2hzar: *.locinames.txt, *.nsamples.txt, *.refAllele.txt, *.dist.txt

library(hzar)
library(doMC)
library(foreach)
library(iterators)
library(parallel)
library(extrafont)

#j <- 1
# Set the working directory here.
setwd("./")

## Set chain length. This value is the default setting in the package.
chainLength=1e5;

## Make each model run off a separate seed
mainSeed=
  list(A=c(as.integer(runif(6, 1, 999))), 
       B=c(as.integer(runif(6, 1, 999))), 
       C=c(as.integer(runif(6, 1, 999))))

if(require(doMC)){
  ## If you have doMC, use foreach in parallel mode
  ## to speed up computation.
  registerDoMC()
} else {
  ## Use foreach in sequential mode
  registerDoSEQ();
}

## Import vector of nsamples
nsamples <- read.csv("BOX_HZAR_nsamples.csv", header=TRUE, stringsAsFactors = FALSE)
nsamples <- nsamples[3:length(nsamples)]

# Read names of loci from file
loci <- read.csv("BOX_HZAR_locinames.csv", header=FALSE, stringsAsFactors=FALSE)

# Read reference allele frequencies from file
refAllele <- read.csv("FINAL_HZAR_DATA_refAlleles.csv", header=TRUE, stringsAsFactors = FALSE)
distances <- refAllele[2]
refAllele <- refAllele[3:length(refAllele)]

distances <- read.table("FIANL_HZAR_DATA_distances.csv", header=FALSE, stringsAsFactors = FALSE)


### Submit a separate job file for each locus being fit (setting i for locus of interest)
  #j=1 # THIS IS SET IN JOB FILE
locnames <- loci[j] #column 1, row j - this is the loci ID being used in the loop
locnsamples <- nsamples[j] #column 1, row j - this is the nsample ID being used in the loop

## Save all plots in a series of png files
#png(width=900, height=900, res=200, family="Arial", filename=paste(locnames,"_rawdata",".png",sep=""),pointsize=8)

## Blank out space in memory to hold molecular analysis
#gd stands for genetic data
if(length(apropos("^gd$",ignore.case=FALSE)) == 0 ||
   !is.list(gd) ) gd <- list()

## doing just the one allele at one locus, but it is good to stay organized.
gd$cur <- list();

## Space to hold the observed data
gd$cur$obs <- list();

## Space to hold the models to fit
gd$cur$models <- list();

## Space to hold the compiled fit requests
gd$cur$fitRs <- list();

## Space to hold the output data chains
gd$cur$runs <- list();

## Space to hold the analysed data
gd$cur$analysis <- list();

gd$cur$obs <-
  hzar.doMolecularData1DPops(distances$V1,
                             refAllele[,j],
                             nsamples[,j]);

## Look at a graph of the observed data
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_rawdata",".png"),pointsize=8)
hzar.plot.obsData(gd$cur$obs);
dev.off()

## Make a helper function
gd.load_cur_model <- function(scaling,tails,
                              id=paste(scaling,tails,sep="."))
  gd$cur$models[[id]] <<- hzar.makeCline1DFreq(gd$cur$obs, scaling, tails)

gd.load_cur_model("none","none","model1");
gd.load_cur_model("none" ,"left","model2");
gd.load_cur_model("none" ,"right","model3");
gd.load_cur_model("none" ,"mirror","model4");
gd.load_cur_model("none" ,"both","model5");
gd.load_cur_model("fixed" ,"none","model6");
gd.load_cur_model("fixed" ,"left","model7");
gd.load_cur_model("fixed" ,"right","model8");
gd.load_cur_model("fixed" ,"mirror","model9");
gd.load_cur_model("fixed" ,"both","model10");
gd.load_cur_model("free" ,"none","model11");
gd.load_cur_model("free" ,"left","model12");
gd.load_cur_model("free" ,"right","model13");
gd.load_cur_model("free" ,"mirror","model14");
gd.load_cur_model("free" ,"both","model15");

## Check the default settings
#print(gd$cur$models)

#save output
checkdefaultsettings<-capture.output(print(gd$cur$models))

cat(checkdefaultsettings,file=paste0(locnames,"_checkdefaultsettings.txt"),sep="",append=TRUE)

## Modify all models to focus on the region where the observed data were collected.
## Observations were between 0 and 681.3946 km
gd$cur$models <- sapply(gd$cur$models,
                        hzar.model.addBoxReq,
                        -30 , 600,
                        simplify=FALSE)

## Check the updated settings
#print(gd$cur$models)

#save output
checkupdatedsettings<-capture.output(print(gd$cur$models))

cat(checkupdatedsettings,file=paste0(locnames,"_checkupdatedsettings.txt"),sep="",append=TRUE)

## Compile each of the models to prepare for fitting
gd$cur$fitRs$init <- sapply(gd$cur$models,
                            hzar.first.fitRequest.old.ML,
                            obsData=gd$cur$obs,
                            verbose=FALSE,
                            simplify=FALSE)

gd$cur$fitRs$init$model1$mcmcParam$chainLength <-
  chainLength; #1e5

gd$cur$fitRs$init$model1$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

gd$cur$fitRs$init$model1$mcmcParam$seed[[1]] <-
  mainSeed$A

gd$cur$fitRs$init$model2$mcmcParam$chainLength <-
  chainLength; #1e5

gd$cur$fitRs$init$model2$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

gd$cur$fitRs$init$model2$mcmcParam$seed[[1]] <-
  mainSeed$B

gd$cur$fitRs$init$model3$mcmcParam$chainLength <-
  chainLength; #1e5

gd$cur$fitRs$init$model3$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

gd$cur$fitRs$init$model3$mcmcParam$seed[[1]] <-
  mainSeed$C

gd$cur$fitRs$init$model4$mcmcParam$chainLength <-
  chainLength; #1e5

gd$cur$fitRs$init$model4$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

gd$cur$fitRs$init$model4$mcmcParam$seed[[1]] <-
  mainSeed$C  

gd$cur$fitRs$init$model5$mcmcParam$chainLength <-
  chainLength; #1e5

gd$cur$fitRs$init$model5$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

gd$cur$fitRs$init$model5$mcmcParam$seed[[1]] <-
  mainSeed$C

gd$cur$fitRs$init$model6$mcmcParam$chainLength <-
  chainLength; #1e5

gd$cur$fitRs$init$model6$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

gd$cur$fitRs$init$model6$mcmcParam$seed[[1]] <-
  mainSeed$C

gd$cur$fitRs$init$model7$mcmcParam$chainLength <-
  chainLength; #1e5

gd$cur$fitRs$init$model7$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

gd$cur$fitRs$init$model7$mcmcParam$seed[[1]] <-
  mainSeed$C

gd$cur$fitRs$init$model8$mcmcParam$chainLength <-
  chainLength; #1e5

gd$cur$fitRs$init$model8$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

gd$cur$fitRs$init$model8$mcmcParam$seed[[1]] <-
  mainSeed$C

gd$cur$fitRs$init$model9$mcmcParam$chainLength <-
  chainLength; #1e5

gd$cur$fitRs$init$model9$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

gd$cur$fitRs$init$model9$mcmcParam$seed[[1]] <-
  mainSeed$C

gd$cur$fitRs$init$model10$mcmcParam$chainLength <-
  chainLength; #1e5

gd$cur$fitRs$init$model10$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

gd$cur$fitRs$init$model10$mcmcParam$seed[[1]] <-
  mainSeed$C

gd$cur$fitRs$init$model11$mcmcParam$chainLength <-
  chainLength; #1e5

gd$cur$fitRs$init$model11$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

gd$cur$fitRs$init$model11$mcmcParam$seed[[1]] <-
  mainSeed$C

gd$cur$fitRs$init$model12$mcmcParam$chainLength <-
  chainLength; #1e5

gd$cur$fitRs$init$model12$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

gd$cur$fitRs$init$model12$mcmcParam$seed[[1]] <-
  mainSeed$C

gd$cur$fitRs$init$model13$mcmcParam$chainLength <-
  chainLength; #1e5

gd$cur$fitRs$init$model13$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

gd$cur$fitRs$init$model13$mcmcParam$seed[[1]] <-
  mainSeed$C

gd$cur$fitRs$init$model4$mcmcParam$chainLength <-
  chainLength; #1e5

gd$cur$fitRs$init$model4$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

gd$cur$fitRs$init$model4$mcmcParam$seed[[1]] <-
  mainSeed$C

gd$cur$fitRs$init$model15$mcmcParam$chainLength <-
  chainLength; #1e5

gd$cur$fitRs$init$model15$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

gd$cur$fitRs$init$model15$mcmcParam$seed[[1]] <-
  mainSeed$C

## Check fit request settings
#print(gd$cur$fitRs$init)

checkfitrequestsettings<-capture.output(print(gd$cur$fitRs$init))

cat(checkfitrequestsettings,file=paste0(locnames,"_checkfitrequestsettings.txt"),sep="",append=TRUE)

## Run just one of the models for an initial chain
gd$cur$runs$init <- 
  list()

gd$cur$runs$init$model1 <-
  hzar.doFit(gd$cur$fitRs$init$model1)

## Plot the trace
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_tracemodel_1",".png"),pointsize=8)
plot(hzar.mcmc.bindLL(gd$cur$runs$init$model1))
dev.off()

## Run another model for an initial chain
gd$cur$runs$init$model2 <-
  hzar.doFit(gd$cur$fitRs$init$model2)

## Plot the trace
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_tracemodel_2",".png"),pointsize=8)
plot(hzar.mcmc.bindLL(gd$cur$runs$init$model2))
dev.off()

## Run another model for an initial chain
gd$cur$runs$init$model3 <-
  hzar.doFit(gd$cur$fitRs$init$model3)

## Plot the trace
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_tracemodel_3",".png"),pointsize=8)
plot(hzar.mcmc.bindLL(gd$cur$runs$init$model3))
dev.off()

## Run another model for an initial chain
gd$cur$runs$init$model4 <-
  hzar.doFit(gd$cur$fitRs$init$model4)

## Plot the trace
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_tracemodel_4",".png"),pointsize=8)
plot(hzar.mcmc.bindLL(gd$cur$runs$init$model4))
dev.off()

## Run another model for an initial chain
gd$cur$runs$init$model5 <-
  hzar.doFit(gd$cur$fitRs$init$model5)

## Plot the trace
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_tracemodel_5",".png"),pointsize=8)
plot(hzar.mcmc.bindLL(gd$cur$runs$init$model5))
dev.off()

## Run another model for an initial chain
gd$cur$runs$init$model6 <-
  hzar.doFit(gd$cur$fitRs$init$model6)

## Plot the trace
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_tracemodel_6",".png"),pointsize=8)
plot(hzar.mcmc.bindLL(gd$cur$runs$init$model6))
dev.off()

## Run another model for an initial chain
gd$cur$runs$init$model7 <-
  hzar.doFit(gd$cur$fitRs$init$model7)

## Plot the trace
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_tracemodel_7",".png"),pointsize=8)
plot(hzar.mcmc.bindLL(gd$cur$runs$init$model7))
dev.off()

## Run another model for an initial chain
gd$cur$runs$init$model8 <-
  hzar.doFit(gd$cur$fitRs$init$model8)

## Plot the trace
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_tracemodel_8",".png"),pointsize=8)
plot(hzar.mcmc.bindLL(gd$cur$runs$init$model8))
dev.off()

## Run another model for an initial chain
gd$cur$runs$init$model9 <-
  hzar.doFit(gd$cur$fitRs$init$model9)

## Plot the trace
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_tracemodel_9",".png"),pointsize=8)
plot(hzar.mcmc.bindLL(gd$cur$runs$init$model9))
dev.off()

## Run another model for an initial chain
gd$cur$runs$init$model10 <-
  hzar.doFit(gd$cur$fitRs$init$model10)

## Plot the trace
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_tracemodel_10",".png"),pointsize=8)
plot(hzar.mcmc.bindLL(gd$cur$runs$init$model10))
dev.off()

## Run another model for an initial chain
gd$cur$runs$init$model11 <-
  hzar.doFit(gd$cur$fitRs$init$model11)

## Plot the trace
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_tracemodel_11",".png"),pointsize=8)
plot(hzar.mcmc.bindLL(gd$cur$runs$init$model11))
dev.off()

## Run another model for an initial chain
gd$cur$runs$init$model12 <-
  hzar.doFit(gd$cur$fitRs$init$model12)

## Plot the trace
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_tracemodel_12",".png"),pointsize=8)
plot(hzar.mcmc.bindLL(gd$cur$runs$init$model12))
dev.off()

## Run another model for an initial chain
gd$cur$runs$init$model13 <-
  hzar.doFit(gd$cur$fitRs$init$model13)

## Plot the trace
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_tracemodel_13",".png"),pointsize=8)
plot(hzar.mcmc.bindLL(gd$cur$runs$init$model13))
dev.off()

## Run another model for an initial chain
gd$cur$runs$init$model14 <-
  hzar.doFit(gd$cur$fitRs$init$model14)

## Plot the trace
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_tracemodel_14",".png"),pointsize=8)
plot(hzar.mcmc.bindLL(gd$cur$runs$init$model14))
dev.off()

## Run another model for an initial chain
gd$cur$runs$init$model15 <-
  hzar.doFit(gd$cur$fitRs$init$model15)

## Plot the trace
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_tracemodel_15",".png"),pointsize=8)
plot(hzar.mcmc.bindLL(gd$cur$runs$init$model15))
dev.off()

## Compile a new set of fit requests using the initial chains
gd$cur$fitRs$chains <-
  lapply(gd$cur$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
gd$cur$fitRs$chains <-
  hzar.multiFitRequest(gd$cur$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit

for(i in seq(4, 6)){
  print(gd$cur$fitRs$chains[[i]]$modelParam$init["width"])
  
}

# Get parameters in each model.


## center
random_center <- runif(45, -30, 600) # center for all models
for(i in seq(1, 45)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["center"] <-  random_center[i]
}

## width
random_width <- runif(45, 0, 630)
for(i in seq(1, 45)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["width"] <- random_width[i]
}

## pMin
for(i in seq(31, 45)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["pMin"]=  runif(1, 0, 1)
}

## pMax
for(i in seq(31, 45)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["pMax"]=  runif(1, 0, 1)
  
}

# deltaL
for(i in seq(4, 6)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["deltaL"]=  runif(1, 0, 630)
}
for(i in seq(13, 15)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["deltaL"]=  runif(1, 0, 630)
}
for(i in seq(19, 21)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["deltaL"]=  runif(1, 0, 630)
}
for(i in seq(28, 30)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["deltaL"]=  runif(1, 0, 630)
}
for(i in seq(34, 36)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["deltaL"]=  runif(1, 0, 630)
}
for(i in seq(43, 45)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["deltaL"]=  runif(1, 0, 630)
}

# tauL
for(i in seq(4, 6)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["tauL"]=  runif(1, 0, 1)
}
for(i in seq(13, 15)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["tauL"]=  runif(1, 0, 1)
}
for(i in seq(19, 21)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["tauL"]=  runif(1, 0, 1)
}
for(i in seq(28, 30)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["tauL"]=  runif(1, 0, 1)
}
for(i in seq(34, 36)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["tauL"]=  runif(1, 0, 1)
}
for(i in seq(43, 45)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["tauL"]=  runif(1, 0, 1)
}


## deltaR

for(i in seq(7, 9)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["deltaR"]=  runif(1, 0, 630)
}
for(i in seq(13, 15)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["deltaR"]=  runif(1, 0, 630)
}
for(i in seq(22, 24)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["deltaR"]=  runif(1, 0, 630)
}
for(i in seq(28, 30)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["deltaR"]=  runif(1, 0, 630)
}
for(i in seq(37, 39)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["deltaR"]=  runif(1, 0, 630)
}
for(i in seq(43, 45)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["deltaR"]=  runif(1, 0, 630)
}

## tauR

for(i in seq(7, 9)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["tauR"]=  runif(1, 0, 1)
}
for(i in seq(13, 15)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["tauR"]=  runif(1, 0, 1)
}
for(i in seq(22, 24)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["tauR"]=  runif(1, 0, 1)
}
for(i in seq(28, 30)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["tauR"]=  runif(1, 0, 1)
}
for(i in seq(37, 39)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["tauR"]=  runif(1, 0, 1)
}
for(i in seq(43, 45)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["tauR"]=  runif(1, 0, 1)
}

## deltaM

for(i in seq(10, 12)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["deltaM"]=  runif(1, 0, 630)
}
for(i in seq(25, 27)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["deltaM"]=  runif(1, 0, 630)
}
for(i in seq(40, 42)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["deltaM"]=  runif(1, 0, 630)
}

# tauM

for(i in seq(10, 12)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["tauM"]=  runif(1, 0, 1)
}
for(i in seq(25, 27)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["tauM"]=  runif(1, 0, 1)
}
for(i in seq(40, 42)){
  gd$cur$fitRs$chains[[i]]$modelParam$init["tauM"]=  runif(1, 0, 1)
}

## Go ahead and run a chain of 3 runs for every fit request
gd$cur$runs$chains <- hzar.doChain.multi(gd$cur$fitRs$chains,
                                         doPar=TRUE,
                                         inOrder=FALSE,
                                         count=3)

## Did model1 converge?
#summary(do.call(mcmc.list,
#                    lapply(gd$cur$runs$chains[1:3],
#                            function(x) hzar.mcmc.bindLL(x[[3]]) )) )

#save output
Checkmodel_1_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                               lapply(gd$cur$runs$chains[1:3],
                                                                      function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_1_convergence,file=paste0(locnames,"_Check_model_1_convergence.txt"),sep="",append=TRUE)



## Did model2 converge?
#summary(do.call(mcmc.list,
#                  lapply(gd$cur$runs$chains[4:6],
#                           function(x) hzar.mcmc.bindLL(x[[3]]) )) )

#save output
Checkmodel_2_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                               lapply(gd$cur$runs$chains[4:6],
                                                                      function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_2_convergence,file=paste0(locnames,"_Check_model_2_convergence.txt"),sep="",append=TRUE)



## Did model3 converge?
#summary(do.call(mcmc.list,
#                lapply(gd$cur$runs$chains[7:9],
#                         function(x) hzar.mcmc.bindLL(x[[3]]) )) )

#save output
Checkmodel_3_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                               lapply(gd$cur$runs$chains[7:9],
                                                                      function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_3_convergence,file=paste0(locnames,"_Check_model_3_convergence.txt"),sep="",append=TRUE)



## Did model4 converge?
#summary(do.call(mcmc.list,
#                lapply(gd$cur$runs$chains[10:12],
#                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

#save output
Checkmodel_4_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                               lapply(gd$cur$runs$chains[10:12],
                                                                      function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_4_convergence,file=paste0(locnames,"_Check_model_4_convergence.txt"),sep="",append=TRUE)



## Did model5 converge?
#summary(do.call(mcmc.list,
#                lapply(gd$cur$runs$chains[13:15],
#                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

#save output
Checkmodel_5_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                               lapply(gd$cur$runs$chains[13:15],
                                                                      function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_5_convergence,file=paste0(locnames,"_Check_model_5_convergence.txt"),sep="",append=TRUE)



## Did model6 converge?
#summary(do.call(mcmc.list,
#                lapply(gd$cur$runs$chains[16:18],
#                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

#save output
Checkmodel_6_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                               lapply(gd$cur$runs$chains[16:18],
                                                                      function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_6_convergence,file=paste0(locnames,"_Check_model_6_convergence.txt"),sep="",append=TRUE)



## Did model7 converge?
#summary(do.call(mcmc.list,
#                lapply(gd$cur$runs$chains[19:21],
#                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

#save output
Checkmodel_7_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                               lapply(gd$cur$runs$chains[19:21],
                                                                      function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_7_convergence,file=paste0(locnames,"_Check_model_7_convergence.txt"),sep="",append=TRUE)



## Did model8 converge?
#summary(do.call(mcmc.list,
#                lapply(gd$cur$runs$chains[22:24],
#                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

#save output
Checkmodel_8_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                               lapply(gd$cur$runs$chains[22:24],
                                                                      function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_8_convergence,file=paste0(locnames,"_Check_model_8_convergence.txt"),sep="",append=TRUE)



## Did model9 converge?
#summary(do.call(mcmc.list,
#                lapply(gd$cur$runs$chains[25:27],
#                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

#save output
Checkmodel_9_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                               lapply(gd$cur$runs$chains[25:27],
                                                                      function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_9_convergence,file=paste0(locnames,"_Check_model_9_convergence.txt"),sep="",append=TRUE)


## Did model10 converge?
#summary(do.call(mcmc.list,
#                lapply(gd$cur$runs$chains[28:30],
#                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

#save output
Checkmodel_10_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                                lapply(gd$cur$runs$chains[28:30],
                                                                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_10_convergence,file=paste0(locnames,"_Check_model_10_convergence.txt"),sep="",append=TRUE)


## Did model11 converge?
#summary(do.call(mcmc.list,
#                lapply(gd$cur$runs$chains[31:33],
 #                      function(x) hzar.mcmc.bindLL(x[[3]]) )) )

#save output
Checkmodel_11_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                                lapply(gd$cur$runs$chains[31:33],
                                                                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_11_convergence,file=paste0(locnames,"_Check_model_11_convergence.txt"),sep="",append=TRUE)


## Did model12 converge?
#summary(do.call(mcmc.list,
#                lapply(gd$cur$runs$chains[34:36],
#                      function(x) hzar.mcmc.bindLL(x[[3]]) )) )

#save output
Checkmodel_12_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                                lapply(gd$cur$runs$chains[34:36],
                                                                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_12_convergence,file=paste0(locnames,"_Check_model_12_convergence.txt"),sep="",append=TRUE)



## Did model13 converge?
#summary(do.call(mcmc.list,
#                lapply(gd$cur$runs$chains[37:39],
#                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

#save output
Checkmodel_13_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                                lapply(gd$cur$runs$chains[37:39],
                                                                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_13_convergence,file=paste0(locnames,"_Check_model_13_convergence.txt"),sep="",append=TRUE)



## Did model14 converge?
#summary(do.call(mcmc.list,
#                lapply(gd$cur$runs$chains[40:42],
#                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

#save output
Checkmodel_14_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                                lapply(gd$cur$runs$chains[40:42],
                                                                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_14_convergence,file=paste0(locnames,"_Check_model_14_convergence.txt"),sep="",append=TRUE)


## Did model15 converge?
#print(summary(do.call(mcmc.list,
#                lapply(gd$cur$runs$chains[43:45],
#                       function(x) hzar.mcmc.bindLL(x[[3]]) )) ))

#save output
Checkmodel_15_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                                lapply(gd$cur$runs$chains[43:45],
                                                                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_15_convergence,file=paste0(locnames,"_Check_model_15_convergence.txt"),sep="",append=TRUE)

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
gd$cur$analysis$initDGs <- list(
  nullModel = hzar.dataGroup.null(gd$cur$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
gd$cur$analysis$initDGs$model1 <-
  hzar.dataGroup.add(gd$cur$runs$init$model1)

gd$cur$analysis$initDGs$model2 <-
  hzar.dataGroup.add(gd$cur$runs$init$model2)

gd$cur$analysis$initDGs$model3 <-
  hzar.dataGroup.add(gd$cur$runs$init$model3)

gd$cur$analysis$initDGs$model4 <-
  hzar.dataGroup.add(gd$cur$runs$init$model4)

gd$cur$analysis$initDGs$model5 <-
  hzar.dataGroup.add(gd$cur$runs$init$model5)

gd$cur$analysis$initDGs$model6 <-
  hzar.dataGroup.add(gd$cur$runs$init$model6)

gd$cur$analysis$initDGs$model7 <-
  hzar.dataGroup.add(gd$cur$runs$init$model7)

gd$cur$analysis$initDGs$model8 <-
  hzar.dataGroup.add(gd$cur$runs$init$model8)

gd$cur$analysis$initDGs$model9 <-
  hzar.dataGroup.add(gd$cur$runs$init$model9)

gd$cur$analysis$initDGs$model10 <-
  hzar.dataGroup.add(gd$cur$runs$init$model10)

gd$cur$analysis$initDGs$model11 <-
  hzar.dataGroup.add(gd$cur$runs$init$model11)

gd$cur$analysis$initDGs$model12 <-
  hzar.dataGroup.add(gd$cur$runs$init$model12)

gd$cur$analysis$initDGs$model13 <-
  hzar.dataGroup.add(gd$cur$runs$init$model13)

gd$cur$analysis$initDGs$model14 <-
  hzar.dataGroup.add(gd$cur$runs$init$model14)

gd$cur$analysis$initDGs$model15 <-
  hzar.dataGroup.add(gd$cur$runs$init$model15)

##Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, model1......to model15).
gd$cur$analysis$oDG <-
  hzar.make.obsDataGroup(gd$cur$analysis$initDGs)

gd$cur$analysis$oDG <-
  hzar.copyModelLabels(gd$cur$analysis$initDGs,
                       gd$cur$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
gd$cur$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(gd$cur$runs$chains,
                                hzar.dataGroup.add),
                         gd$cur$analysis$oDG);

## Check to make sure that there are only 16 hzar.dataGroup objects in the hzar.obsDataGroup object.
#print(summary(gd$cur$analysis$oDG$data.groups))

#save output
Checkdatagroupobjs<-capture.output(print(summary(gd$cur$analysis$oDG$data.groups)))

cat(Checkdatagroupobjs,file=paste0(locnames,"_Check_dataGroup_objs.txt"),sep="",append=TRUE)

## Compare the 15 cline models to the null model graphically
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_comparing16models",".png"),pointsize=8)
hzar.plot.cline(gd$cur$analysis$oDG);
dev.off()

## Do model selection based on the AICc scores
print(gd$cur$analysis$AICcTable <-
            hzar.AICc.hzar.obsDataGroup(gd$cur$analysis$oDG));


#save output
AICctableforallmodels<-capture.output(print(gd$cur$analysis$AICcTable <-
                                              hzar.AICc.hzar.obsDataGroup(gd$cur$analysis$oDG)))

cat(AICctableforallmodels,file=paste0(locnames,"_AICc_table_for_all_models.txt"),sep="",append=TRUE)

## Print out the model with the minimum AICc score
print(gd$cur$analysis$model.name <-
            rownames(gd$cur$analysis$AICcTable
                       )[[ which.min(gd$cur$analysis$AICcTable$AICc )]])

#save output
SelectedModel<-capture.output(print(gd$cur$analysis$model.name <-
                                      rownames(gd$cur$analysis$AICcTable
                                      )[[ which.min(gd$cur$analysis$AICcTable$AICc )]]))

cat(SelectedModel,file=paste0(locnames,"_selected_model.txt"),sep="",append=TRUE)

## Extract the hzar.dataGroup object for the selected model
gd$cur$analysis$model.selected <-
  gd$cur$analysis$oDG$data.groups[[gd$cur$analysis$model.name]]

## Look at the variation in parameters for the selected model
#print(hzar.getLLCutParam(gd$cur$analysis$model.selected,
 #                            names(gd$cur$analysis$model.selected$data.param)));

#save the oupt
Var_params<-capture.output(print(hzar.getLLCutParam(gd$cur$analysis$model.selected,
                                                    names(gd$cur$analysis$model.selected$data.param))))
cat(Var_params,file=paste0(locnames,"_MaxLL_var_params_for_selected_model.txt"),sep="",append=TRUE)

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(gd$cur$analysis$model.selected))

#save the oupt
MaxLL_params<-capture.output(print(hzar.get.ML.cline(gd$cur$analysis$model.selected)))
cat(MaxLL_params,file=paste0(locnames,"_MaxLL_params_for_selected_model.txt"),sep="",append=TRUE)

## Plot the maximum likelihood cline for the selected model
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_maxLL_selectedmodel",".png"),pointsize=8)
hzar.plot.cline(gd$cur$analysis$model.selected);
dev.off()

## Plot the 95% credible cline region for the selected model
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_fuzzycline_selectedmodel",".png"),pointsize=8)
hzar.plot.fzCline(gd$cur$analysis$model.selected);
dev.off()

#Rename loc <- replace "cur" with loc name
#names(gd)[[j]]=as.character(locnames)

model.selected <- gd$cur$analysis$model.selected

# Save selected.model object to .rds file for later analysis
loc_obj_rds <- paste0(locnames, ".rds")
saveRDS(model.selected, loc_obj_rds)
