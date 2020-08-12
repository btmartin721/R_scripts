
library("ggplot2") # Needed for plotting
library("dplyr") # Needed for count() and table()
library("stringr") # Needed for str_split()

# Function to caluclate P-values from Z-scores.
pvalue <- function(z){ 2*pnorm(-abs(z)) }

# Calculate percentage of tests that were significant
perc <- function(df, c){
  df$perc <- (as.numeric(as.character(df$Freq)) / max(as.numeric(as.character(c$n)))) * 100
  return(df)
}


# Function to make plots
makeFig <- function(df, my.ylab){
  p <-
    ggplot() + 
    geom_bar(aes(x = vec, y = perc, fill = P3), data = df, stat = "identity") +
    coord_flip(ylim = c(0, 100)) +
    xlab(label = "Individual") + ylab(label = my.ylab) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.background = element_blank(),
      axis.text.y = element_text(colour = "black", size = rel(0.9)),
      axis.text.x = element_text(colour = "black", size = rel(1.5)),
      axis.ticks = element_line(colour = "black"),
      axis.title = element_text(colour = "black", size = rel(1.8)),
      legend.text = element_text(size = rel(1.8)),
      legend.title = element_text(size = rel(1.8))
    )
  return(p)
}

# Tabulate - Requires dplyr
makeTable <- function(vec, P3vec){
  # Split population ID by string before first underscore.
  P3 <- sapply(str_split(P3vec, "_", n = 2), `[`, 1)
  df <- as.data.frame(table(vec, P3), stringsAsFactors = F)
  return(df)
}

# Read the input files
MO <- read.table(file = "alltests_MO_header.txt", header = T, sep = "\t", stringsAsFactors = T, strip.white = T, blank.lines.skip = T)
ED <- read.table(file = "alltests_ED_header.txt", header = T, sep = "\t", stringsAsFactors = T, strip.white = T, blank.lines.skip = T)

# Re-calculate P-values from Z-scores
# I do this to get additional digits out of the P-values
MO$Z.pval.recalc <- pvalue(MO$Z.score)
ED$Z.pval.recalc <- pvalue(ED$Z.score)

# Count total number of tests for each individual in P1 and P2 groups
MO.P1.counts <- count(MO, vars = P1)
ED.P1.counts <- count(ED, vars = P1)
MO.P2.counts <- count(MO, vars = P2)
ED.P2.counts <- count(ED, vars = P2)

# Get total unique individuals
MO.P1.total <- max(MO.P1.counts$n)
ED.P1.total <- max(ED.P1.counts$n)
MO.P2.total <- max(MO.P2.counts$n)
ED.P2.total <- max(ED.P2.counts$n)

# Do p-value adjustments
MO$bonferroni <- p.adjust(p = MO$Z.pval.recalc, method = "bonferroni")
MO$BY <- p.adjust(p = MO$Z.pval.recalc, method = "BY")
ED$bonferroni <- p.adjust(p = ED$Z.pval.recalc, method = "bonferroni")
ED$BY <- p.adjust(p = ED$Z.pval.recalc, method = "BY")

# Subset original dataframes to significant tests
# Negative Z-scores correspond to ABBA site patterns (introgression from P3 <-> P2)
MO.sig.P2.noCorr <- MO[which(MO$Z.pval.recalc <= 0.05 & MO$Z.score < 0), ]
MO.sig.P2.bonf <- MO[which(MO$bonferroni <= 0.05 & MO$Z.score < 0), ]
MO.sig.P2.BY <- MO[which(MO$BY <= 0.05 & MO$Z.score < 0), ]
ED.sig.P2.noCorr <- ED[which(ED$Z.pval.recalc <= 0.05 & ED$Z.score < 0), ]
ED.sig.P2.bonf <- ED[which(ED$bonferroni <= 0.05 & ED$Z.score < 0), ]
ED.sig.P2.BY <- ED[which(ED$BY <= 0.05 & ED$Z.score < 0), ]

# Whereas positive Z-scores correspond to BABA site patterns (introgressino from P3 <-> P1)
MO.sig.P1.noCorr <- MO[which(MO$Z.pval.recalc <= 0.05 & MO$Z.score >= 0), ]
MO.sig.P1.bonf <- MO[which(MO$bonferroni <= 0.05 & MO$Z.score >= 0), ]
MO.sig.P1.BY <- MO[which(MO$BY <= 0.05 & MO$Z.score >= 0), ]
ED.sig.P1.noCorr <- ED[which(ED$Z.pval.recalc <= 0.05 & ED$Z.score >= 0), ]
ED.sig.P1.bonf <- ED[which(ED$bonferroni <= 0.05 & ED$Z.score >= 0), ]
ED.sig.P1.BY <- ED[which(ED$BY <= 0.05 & ED$Z.score >= 0), ]

# Get frequency of each significant individual in P1 group
MO.tab.P2.noCorr <- makeTable(MO.sig.P2.noCorr$P2, MO.sig.P2.noCorr$P3)
MO.tab.P2.bonf <- makeTable(MO.sig.P2.bonf$P2, MO.sig.P2.bonf$P3)
MO.tab.P2.BY <- makeTable(MO.sig.P2.BY$P2, MO.sig.P2.BY$P3)
ED.tab.P2.noCorr <- makeTable(ED.sig.P2.noCorr$P2, ED.sig.P2.noCorr$P3)
ED.tab.P2.bonf <- makeTable(ED.sig.P2.bonf$P2, ED.sig.P2.bonf$P3)
ED.tab.P2.BY <- makeTable(ED.sig.P2.BY$P2, ED.sig.P2.BY$P3)

# Frequency of each significant individual in P1 group.
MO.tab.P1.noCorr <- makeTable(MO.sig.P1.noCorr$P1, MO.sig.P1.noCorr$P3)
MO.tab.P1.bonf <- makeTable(MO.sig.P1.bonf$P1, MO.sig.P1.bonf$P3)
MO.tab.P1.BY <- makeTable(MO.sig.P1.BY$P1, MO.sig.P1.BY$P3)
ED.tab.P1.noCorr <- makeTable(ED.sig.P1.noCorr$P1, ED.sig.P1.noCorr$P3)
ED.tab.P1.bonf <- makeTable(ED.sig.P1.bonf$P1, ED.sig.P1.bonf$P3)
ED.tab.P1.BY <- makeTable(ED.sig.P1.BY$P1, ED.sig.P1.BY$P3)


# Add percentage column for P2
MO.tab.P2.noCorr <- perc(MO.tab.P2.noCorr, MO.P2.counts)
MO.tab.P2.bonf <- perc(MO.tab.P2.bonf, MO.P2.counts)
MO.tab.P2.BY <- perc(MO.tab.P2.BY, MO.P2.counts)
ED.tab.P2.noCorr <- perc(ED.tab.P2.noCorr, ED.P2.counts)
ED.tab.P2.bonf <- perc(ED.tab.P2.bonf, ED.P2.counts)
ED.tab.P2.BY <- perc(ED.tab.P2.BY, ED.P2.counts)

# Add percentage column for P1
MO.tab.P1.noCorr <- perc(MO.tab.P1.noCorr, MO.P1.counts)
MO.tab.P1.bonf <- perc(MO.tab.P1.bonf, MO.P1.counts)
MO.tab.P1.BY <- perc(MO.tab.P1.BY, MO.P1.counts)
ED.tab.P1.noCorr <- perc(ED.tab.P1.noCorr, ED.P1.counts)
ED.tab.P1.bonf <- perc(ED.tab.P1.bonf, ED.P1.counts)
ED.tab.P1.BY <- perc(ED.tab.P1.BY, ED.P1.counts)

# Combine P2 and P1 dataframes
MO.final.noCorr <- rbind(MO.tab.P2.noCorr, MO.tab.P1.noCorr)
MO.final.bonf <- rbind(MO.tab.P2.bonf, MO.tab.P1.bonf)
MO.final.BY <- rbind(MO.tab.P2.BY, MO.tab.P1.BY)
ED.final.noCorr <- rbind(ED.tab.P2.noCorr, ED.tab.P1.noCorr)
ED.final.bonf <- rbind(ED.tab.P2.bonf, ED.tab.P1.bonf)
ED.final.BY <- rbind(ED.tab.P2.BY, ED.tab.P1.BY)


# Plot it
pdf(file = "mas_dtest_figures.pdf", width = 8.5, height = 11, onefile = T)
makeFig(MO.final.noCorr, "Percentage")
makeFig(MO.final.bonf, "Percentage (Bonferroni Correction)")
makeFig(MO.final.BY, "Percentage (B-Y Correction)")
makeFig(ED.final.noCorr, "Percentage")
makeFig(ED.final.bonf, "Percentage (Bonferroni Correction)")
makeFig(ED.final.BY, "Percentage (B-Y Correction)")
dev.off()

# Write the data to a txt file
# Missouri
write.table(x = MO.final.noCorr, file = "MO.final.noCorr.txt", sep = "\t", col.names = F, quote = F, row.names = F)
write.table(x = MO.final.bonf, file = "MO.final.bonf.txt", sep = "\t", col.names = F, quote = F, row.names = F)
write.table(x = MO.final.BY, file = "MO.final.BY.txt", sep = "\t", col.names = F, quote = F, row.names = F)

# edwardsii
write.table(x = ED.final.noCorr, file = "ED.final.noCorr.txt", sep = "\t", col.names = F, quote = F, row.names = F)
write.table(x = ED.final.bonf, file = "ED.final.bonf.txt", sep = "\t", col.names = F, quote = F, row.names = F)
write.table(x = ED.final.BY, file = "ED.final.BY.txt", sep = "\t", col.names = F, quote = F, row.names = F)



