library("genepopedit")
library("hybriddetective")

###*** SET WORKING DIRECTORY HERE ***###

my_path <- getwd()

#################################################################################

###*** SETTINGS ***###


genepopFile <- paste0(my_path, "/BOX_hzar_EAGU_genepop.txt")
fixedGenepop <- paste0(my_path, "/BOX_hzar_EAGU_fixedGenepop.txt")

distancesFile <- paste0(my_path, "/BOX_hzar_distPop.txt")

outputFile <- paste0(my_path, "/BOX_hzar_EAGU_finalInput.csv")

hzar_clinalPops <- paste0(my_path, "/BOX_hzar_EAGU_finalInput.csv")
hzar_clinalDist <- paste0(my_path, "/BOX_hzar_distPop_clinal.txt")

plink <- paste0(my_path, "/plink_win64_20181202/")
pgd <- paste0(my_path, "/PGDSpider_2.1.1.5/PGDSpider_2.1.1.5/")

################################################################################

Distance <- read.delim(file = distancesFile, header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
clinalDist <- read.delim(file = hzar_clinalDist, header = TRUE, sep="\t", row.names = NULL, stringsAsFactors = FALSE)
names(clinalDist) <- c("Pop", "Distance")

names(Distance) <- c("Pop", "Distance")

genepop_ID(genepop = genepopFile, path = fixedGenepop)

pops <- genepop_detective(genepop = genepopFile, variable = "Pops")
print(pops)

genepop_hzar(genepop = fixedGenepop, distances = clinalDist, path = hzar_clinalPops)

hzarFile <- read.csv(file = hzar_clinalPops, header = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)

popCol <- hzarFile$Population
distCol <- hzarFile$Distance

lociNames <- colnames(hzarFile)
lociNames <- data.frame(lociNames[3:length(lociNames)], stringsAsFactors = FALSE, row.names = NULL)

colnames(lociNames) <- c("loci")

locVector <- colnames(hzarFile)
locVector <- locVector[3:length(lociNames)]

Maxpops <- nrow(hzarFile)

noPopdistHZAR <- hzarFile[3:length(hzarFile)]

numLoci <- length(noPopdistHZAR) / 2
print(length(lociNames))

loci_df <- data.frame(lociNames[1:numLoci, 1], stringsAsFactors = FALSE)
colnames(loci_df) <- c("loci")
loci_df_tr <- data.frame(t(loci_df), stringsAsFactors = FALSE)

diagnosticloci <- character(numLoci)

print(Maxpops)

diagnosticlocitokeep <- (noPopdistHZAR[1,1:numLoci] <= 0.2 & noPopdistHZAR[Maxpops,1:numLoci] >= 0.8 | noPopdistHZAR[1,1:numLoci] >= 0.8 & noPopdistHZAR[Maxpops,1:numLoci] <= 0.2)

diagnosticloci <- loci_df_tr[1,diagnosticlocitokeep]

diagnosticAlleleFreq <- noPopdistHZAR[1:Maxpops, diagnosticlocitokeep]

final_df <- cbind(Population = popCol, Distance = distCol, diagnosticAlleleFreq, stringsAsFactors = FALSE)

write.table(final_df, file="FINAL_HZAR_DATA.csv", quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

write.table(diagnosticloci, file="FINALHZAR_DATA_locinames.csv", quote=FALSE, sep=",", row.names=FALSE, col.names=FALSE)

