library(tidyverse)
library(plot.matrix)
# data sets needed from previous scripts: similarity_outputs.csv


#########################################
#### ---- Tests for parallelism ---- ####
#########################################

similarity_outputs <- read.csv("output/similarity_outputs.csv", stringsAsFactors = FALSE, na.strings = c("NA"))

# count the number of QTL in inversions for each mapping population
M1HB <- similarity_outputs %>% filter(MON1_sig & HB) %>% nrow()
M2HB <- similarity_outputs %>% filter(MON2_sig & HB) %>% nrow()
G1HB <- similarity_outputs %>% filter(GSD1_sig & HB) %>% nrow()
G2HB <- similarity_outputs %>% filter(GSD2_sig & HB) %>% nrow()

# count the number of QTL in inversions for each mapping population
M1xHB <- sum(!similarity_outputs[similarity_outputs$MON1_sig, ]$HB)
M2xHB <- sum(!similarity_outputs[similarity_outputs$MON2_sig, ]$HB)
G1xHB <- sum(!similarity_outputs[similarity_outputs$GSD1_sig, ]$HB)
G2xHB <- sum(!similarity_outputs[similarity_outputs$GSD2_sig, ]$HB)

# count the number of QTL shared by pairs of mapping populations
M1M2 <- similarity_outputs %>% filter(MON1_sig & MON2_sig) %>% nrow()
M1G1 <- similarity_outputs %>% filter(MON1_sig & GSD1_sig) %>% nrow()
M1G2 <- similarity_outputs %>% filter(MON1_sig & GSD2_sig) %>% nrow()
M2G1 <- similarity_outputs %>% filter(MON2_sig & GSD1_sig) %>% nrow()
M2G2 <- similarity_outputs %>% filter(MON2_sig & GSD2_sig) %>% nrow()
G1G2 <- similarity_outputs %>% filter(GSD1_sig & GSD2_sig) %>% nrow()

# count the number of QTL shared by pairs of mapping populations within inversions
M1M2_HB <- sum(similarity_outputs[similarity_outputs$MON1_sig & similarity_outputs$MON2_sig, ]$HB)
M1G1_HB <- sum(similarity_outputs[similarity_outputs$MON1_sig & similarity_outputs$GSD1_sig, ]$HB)
M1G2_HB <- sum(similarity_outputs[similarity_outputs$MON1_sig & similarity_outputs$GSD2_sig, ]$HB)
M2G1_HB <- sum(similarity_outputs[similarity_outputs$MON2_sig & similarity_outputs$GSD1_sig, ]$HB)
M2G2_HB <- sum(similarity_outputs[similarity_outputs$MON2_sig & similarity_outputs$GSD2_sig, ]$HB)
G1G2_HB <- sum(similarity_outputs[similarity_outputs$GSD1_sig & similarity_outputs$GSD2_sig, ]$HB)

# count the number of QTL shared by pairs of mapping populations outside inversions
M1M2_xHB <- sum(!similarity_outputs[similarity_outputs$MON1_sig & similarity_outputs$MON2_sig, ]$HB)
M1G1_xHB <- sum(!similarity_outputs[similarity_outputs$MON1_sig & similarity_outputs$GSD1_sig, ]$HB)
M1G2_xHB <- sum(!similarity_outputs[similarity_outputs$MON1_sig & similarity_outputs$GSD2_sig, ]$HB)
M2G1_xHB <- sum(!similarity_outputs[similarity_outputs$MON2_sig & similarity_outputs$GSD1_sig, ]$HB)
M2G2_xHB <- sum(!similarity_outputs[similarity_outputs$MON2_sig & similarity_outputs$GSD2_sig, ]$HB)
G1G2_xHB <- sum(!similarity_outputs[similarity_outputs$GSD1_sig & similarity_outputs$GSD2_sig, ]$HB)

# determine which inversions have QTL shared across pairs of mapping populations (Table 1)
similarity_outputs %>% filter(MON1_sig & MON2_sig & HB) %>% select(CHR) %>% table()
similarity_outputs %>% filter(MON1_sig & GSD1_sig & HB) %>% select(CHR) %>% table()
similarity_outputs %>% filter(MON1_sig & GSD2_sig & HB) %>% select(CHR) %>% table()
similarity_outputs %>% filter(MON2_sig & GSD1_sig & HB) %>% select(CHR) %>% table()
similarity_outputs %>% filter(MON2_sig & GSD2_sig & HB) %>% select(CHR) %>% table()
similarity_outputs %>% filter(GSD1_sig & GSD2_sig & HB) %>% select(CHR) %>% table()


##########################################
#### ---- Complete randomization ---- ####
##########################################

# randomize significant windows and count overlaps
SO_randomized <- similarity_outputs
randomized_overlaps <- data.frame(M1M2=NA, M1G1=NA, M1G2=NA, M2G1=NA, M2G2=NA, G1G2=NA, M1HB=NA, M2HB=NA, G1HB=NA, G2HB=NA, 
                                  M1M2_HB=NA, M1G1_HB=NA, M1G2_HB=NA, M2G1_HB=NA, M2G2_HB=NA, G1G2_HB=NA, 
                                  M1M2_xHB=NA, M1G1_xHB=NA, M1G2_xHB=NA, M2G1_xHB=NA, M2G2_xHB=NA, G1G2_xHB=NA)
set.seed(1)
for (i in 1:10000) {
  # randomization
  SO_randomized$MON1_sig <- sample(SO_randomized$MON1_sig)
  SO_randomized$MON2_sig <- sample(SO_randomized$MON2_sig)
  SO_randomized$GSD1_sig <- sample(SO_randomized$GSD1_sig)
  SO_randomized$GSD2_sig <- sample(SO_randomized$GSD2_sig)
  # count shared QTL between pairs
  randomized_overlaps[i,1] <- SO_randomized %>% filter(MON1_sig & MON2_sig) %>% nrow()
  randomized_overlaps[i,2] <- SO_randomized %>% filter(MON1_sig & GSD1_sig) %>% nrow()
  randomized_overlaps[i,3] <- SO_randomized %>% filter(MON1_sig & GSD2_sig) %>% nrow()
  randomized_overlaps[i,4] <- SO_randomized %>% filter(MON2_sig & GSD1_sig) %>% nrow()
  randomized_overlaps[i,5] <- SO_randomized %>% filter(MON2_sig & GSD2_sig) %>% nrow()
  randomized_overlaps[i,6] <- SO_randomized %>% filter(GSD1_sig & GSD2_sig) %>% nrow()
  # count QTL within inversions
  randomized_overlaps[i,7] <- SO_randomized %>% filter(MON1_sig & HB) %>% nrow()
  randomized_overlaps[i,8] <- SO_randomized %>% filter(MON2_sig & HB) %>% nrow()
  randomized_overlaps[i,9] <- SO_randomized %>% filter(GSD1_sig & HB) %>% nrow()
  randomized_overlaps[i,10] <- SO_randomized %>% filter(GSD2_sig & HB) %>% nrow()
  # count shared QTL within inversions
  randomized_overlaps[i,11] <- sum(SO_randomized[SO_randomized$MON1_sig & SO_randomized$MON2_sig, ]$HB)
  randomized_overlaps[i,12] <- sum(SO_randomized[SO_randomized$MON1_sig & SO_randomized$GSD1_sig, ]$HB)
  randomized_overlaps[i,13] <- sum(SO_randomized[SO_randomized$MON1_sig & SO_randomized$GSD2_sig, ]$HB)
  randomized_overlaps[i,14] <- sum(SO_randomized[SO_randomized$MON2_sig & SO_randomized$GSD1_sig, ]$HB)
  randomized_overlaps[i,15] <- sum(SO_randomized[SO_randomized$MON2_sig & SO_randomized$GSD2_sig, ]$HB)
  randomized_overlaps[i,16] <- sum(SO_randomized[SO_randomized$GSD1_sig & SO_randomized$GSD2_sig, ]$HB)
  # count shared QTL outside inversions
  randomized_overlaps[i,17] <- sum(!SO_randomized[SO_randomized$MON1_sig & SO_randomized$MON2_sig, ]$HB)
  randomized_overlaps[i,18] <- sum(!SO_randomized[SO_randomized$MON1_sig & SO_randomized$GSD1_sig, ]$HB)
  randomized_overlaps[i,19] <- sum(!SO_randomized[SO_randomized$MON1_sig & SO_randomized$GSD2_sig, ]$HB)
  randomized_overlaps[i,20] <- sum(!SO_randomized[SO_randomized$MON2_sig & SO_randomized$GSD1_sig, ]$HB)
  randomized_overlaps[i,21] <- sum(!SO_randomized[SO_randomized$MON2_sig & SO_randomized$GSD2_sig, ]$HB)
  randomized_overlaps[i,22] <- sum(!SO_randomized[SO_randomized$GSD1_sig & SO_randomized$GSD2_sig, ]$HB)
}

write.table(randomized_overlaps, file="output/randomized_overlaps.csv", sep=",", row.names = FALSE, col.names = TRUE, quote=FALSE)
randomized_overlaps <- read.csv("output/randomized_overlaps.csv", stringsAsFactors = FALSE, na.strings = c("NA"))

# plot figure S14
layout(matrix(1:18, nrow=6))
hist(randomized_overlaps$M1HB, xlim=c(0,120), breaks=30, xlab="Number of overlaps", main="QTL Within Inversions\nMON1", las=1); abline(v=M1HB, col="red", lwd=2); mtext("A", 3, at=0)
hist(randomized_overlaps$M2HB, xlim=c(0,120), breaks=30, xlab="Number of overlaps", main="QTL Within Inversions\nMON2", las=1); abline(v=M2HB, col="red", lwd=2); mtext("C", 3, at=0)
hist(randomized_overlaps$G1HB, xlim=c(0,120), breaks=30, xlab="Number of overlaps", main="QTL Within Inversions\nGSD1", las=1); abline(v=G1HB, col="red", lwd=2); mtext("B", 3, at=0)
hist(randomized_overlaps$G2HB, xlim=c(0,120), breaks=30, xlab="Number of overlaps", main="QTL Within Inversions\nGSD2", las=1); abline(v=G2HB, col="red", lwd=2); mtext("D", 3, at=0)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
hist(randomized_overlaps$M1M2, xlim=c(0,125), breaks=30, xlab="Number of overlaps", main="Overall Overlap\nMON1 vs MON2", las=1); abline(v=M1M2, col="red", lwd=2); mtext("E", 3, at=0)
hist(randomized_overlaps$G1G2, xlim=c(0,125), breaks=30, xlab="Number of overlaps", main="Overall Overlap\nGSD1 vs GSD2", las=1); abline(v=G1G2, col="red", lwd=2); mtext("F", 3, at=0)
hist(randomized_overlaps$M1G1, xlim=c(0,125), breaks=30, xlab="Number of overlaps", main="Overall Overlap\nMON1 vs GSD1", las=1); abline(v=M1G1, col="red", lwd=2); mtext("G", 3, at=0)
hist(randomized_overlaps$M1G2, xlim=c(0,125), breaks=30, xlab="Number of overlaps", main="Overall Overlap\nMON1 vs GSD2", las=1); abline(v=M1G2, col="red", lwd=2); mtext("H", 3, at=0)
hist(randomized_overlaps$M2G1, xlim=c(0,125), breaks=30, xlab="Number of overlaps", main="Overall Overlap\nMON2 vs GSD1", las=1); abline(v=M2G1, col="red", lwd=2); mtext("I", 3, at=0)
hist(randomized_overlaps$M2G2, xlim=c(0,125), breaks=30, xlab="Number of overlaps", main="Overall Overlap\nMON2 vs GSD2", las=1); abline(v=M2G2, col="blue", lwd=2); mtext("J", 3, at=0)
text(M2G2, 650, pos=4, cex=1, col="blue", round((sum(randomized_overlaps$M2G2 > M2G2)*2+1)/10001, 3))
hist(randomized_overlaps$M1M2_HB/272-randomized_overlaps$M1M2_xHB/1276, xlim=c(-0.05,0.2), breaks=30, xlab="% overlap within - % overlap outside", main="Within vs. Outside\nMON1 vs MON2", las=1); abline(v=M1M2_HB/272-M1M2_xHB/1276, col="red", lwd=2); mtext("K", 3, at=-0.05)
hist(randomized_overlaps$G1G2_HB/272-randomized_overlaps$G1G2_xHB/1276, xlim=c(-0.05,0.2), breaks=30, xlab="% overlap within - % overlap outside", main="Within vs. Outside\nGSD1 vs GSD2", las=1); abline(v=G1G2_HB/272-G1G2_xHB/1276, col="red", lwd=2); mtext("L", 3, at=-0.05)
hist(randomized_overlaps$M1G1_HB/272-randomized_overlaps$M1G1_xHB/1276, xlim=c(-0.05,0.2), breaks=30, xlab="% overlap within - % overlap outside", main="Within vs. Outside\nMON1 vs GSD1", las=1); abline(v=M1G1_HB/272-M1G1_xHB/1276, col="red", lwd=2); mtext("M", 3, at=-0.05)
hist(randomized_overlaps$M1G2_HB/272-randomized_overlaps$M1G2_xHB/1276, xlim=c(-0.05,0.2), breaks=60, xlab="% overlap within - % overlap outside", main="Within vs. Outside\nMON1 vs GSD2", las=1); abline(v=M1G2_HB/272-M1G2_xHB/1276, col="red", lwd=2); mtext("N", 3, at=-0.05)
hist(randomized_overlaps$M2G1_HB/272-randomized_overlaps$M2G1_xHB/1276, xlim=c(-0.05,0.2), breaks=30, xlab="% overlap within - % overlap outside", main="Within vs. Outside\nMON2 vs GSD1", las=1); abline(v=M2G1_HB/272-M2G1_xHB/1276, col="red", lwd=2); mtext("O", 3, at=-0.05)
hist(randomized_overlaps$M2G2_HB/272-randomized_overlaps$M2G2_xHB/1276, xlim=c(-0.05,0.2), breaks=30, xlab="% overlap within - % overlap outside", main="Within vs. Outside\nMON2 vs GSD2", las=1); abline(v=M2G2_HB/272-M2G2_xHB/1276, col="red", lwd=2); mtext("P", 3, at=-0.05)
# other p-values = 1/10001
# save as 8x10 - Figure S14


###################################################
#### ---- Haploblock aware randomization  ---- ####
###################################################

# randomize significant windows within HB and outside HB and count overlaps
SO_randomized <- similarity_outputs
randomized_overlaps2 <- data.frame(M1M2_HB=NA, M1G1_HB=NA, M1G2_HB=NA, M2G1_HB=NA, M2G2_HB=NA, G1G2_HB=NA, M1M2_xHB=NA, M1G1_xHB=NA, M1G2_xHB=NA, M2G1_xHB=NA, M2G2_xHB=NA, G1G2_xHB=NA)
set.seed(1)
for (i in 1:10000) {
  # randomize QTL within inversions
  SO_randomized[SO_randomized$HB,]$MON1_sig <- sample(SO_randomized[SO_randomized$HB,]$MON1_sig)
  SO_randomized[SO_randomized$HB,]$MON2_sig <- sample(SO_randomized[SO_randomized$HB,]$MON2_sig)
  SO_randomized[SO_randomized$HB,]$GSD1_sig <- sample(SO_randomized[SO_randomized$HB,]$GSD1_sig)
  SO_randomized[SO_randomized$HB,]$GSD2_sig <- sample(SO_randomized[SO_randomized$HB,]$GSD2_sig)
  # randomize QTL outside inversions
  SO_randomized[!SO_randomized$HB,]$MON1_sig <- sample(SO_randomized[!SO_randomized$HB,]$MON1_sig)
  SO_randomized[!SO_randomized$HB,]$MON2_sig <- sample(SO_randomized[!SO_randomized$HB,]$MON2_sig)
  SO_randomized[!SO_randomized$HB,]$GSD1_sig <- sample(SO_randomized[!SO_randomized$HB,]$GSD1_sig)
  SO_randomized[!SO_randomized$HB,]$GSD2_sig <- sample(SO_randomized[!SO_randomized$HB,]$GSD2_sig)
  # count overlaps within inversions
  randomized_overlaps2[i,1] <- sum(SO_randomized[SO_randomized$MON1_sig & SO_randomized$MON2_sig, ]$HB)
  randomized_overlaps2[i,2] <- sum(SO_randomized[SO_randomized$MON1_sig & SO_randomized$GSD1_sig, ]$HB)
  randomized_overlaps2[i,3] <- sum(SO_randomized[SO_randomized$MON1_sig & SO_randomized$GSD2_sig, ]$HB)
  randomized_overlaps2[i,4] <- sum(SO_randomized[SO_randomized$MON2_sig & SO_randomized$GSD1_sig, ]$HB)
  randomized_overlaps2[i,5] <- sum(SO_randomized[SO_randomized$MON2_sig & SO_randomized$GSD2_sig, ]$HB)
  randomized_overlaps2[i,6] <- sum(SO_randomized[SO_randomized$GSD1_sig & SO_randomized$GSD2_sig, ]$HB)
  # count overlaps outside inversions
  randomized_overlaps2[i,7] <- sum(!SO_randomized[SO_randomized$MON1_sig & SO_randomized$MON2_sig, ]$HB)
  randomized_overlaps2[i,8] <- sum(!SO_randomized[SO_randomized$MON1_sig & SO_randomized$GSD1_sig, ]$HB)
  randomized_overlaps2[i,9] <- sum(!SO_randomized[SO_randomized$MON1_sig & SO_randomized$GSD2_sig, ]$HB)
  randomized_overlaps2[i,10] <- sum(!SO_randomized[SO_randomized$MON2_sig & SO_randomized$GSD1_sig, ]$HB)
  randomized_overlaps2[i,11] <- sum(!SO_randomized[SO_randomized$MON2_sig & SO_randomized$GSD2_sig, ]$HB)
  randomized_overlaps2[i,12] <- sum(!SO_randomized[SO_randomized$GSD1_sig & SO_randomized$GSD2_sig, ]$HB)
}

write.table(randomized_overlaps2, file="output/randomized_overlaps2.csv", sep=",", row.names = FALSE, col.names = TRUE, quote=FALSE)
randomized_overlaps2 <- read.csv("output/randomized_overlaps2.csv", stringsAsFactors = FALSE, na.strings = c("NA"))

# plot figure S15
layout(matrix(1:24, nrow=6))
# plot QTL count tables
plot(matrix(c(M1HB, M1xHB, M2HB, M2xHB, M1M2_HB, M1M2_xHB), nrow=2), main="QTL Counts\nMON1 vs MON2", digits = TRUE, key=NULL, fmt.cell='%.0f', axis.col=NULL, axis.row=NULL, xlab='', ylab='', col=NULL)
mtext('M1', side=1, line=1, at=1, cex=0.8); mtext('M2', side=1, line=1, at=2, cex=0.8); mtext('SH', side=1, line=1, at=3, cex=0.8); mtext('W', side=2, line=1, at=2, cex=0.8, las=2); mtext('O', side=2, line=1, at=1, cex=0.8, las=2)
plot(matrix(c(G1HB, G1xHB, G2HB, G2xHB, G1G2_HB, G1G2_xHB), nrow=2), main="QTL Counts\nGSD1 vs GSD2", digits = TRUE, key=NULL, fmt.cell='%.0f', axis.col=NULL, axis.row=NULL, xlab='', ylab='', col=NULL)
mtext('G1', side=1, line=1, at=1, cex=0.8); mtext('G2', side=1, line=1, at=2, cex=0.8); mtext('SH', side=1, line=1, at=3, cex=0.8); mtext('W', side=2, line=1, at=2, cex=0.8, las=2); mtext('O', side=2, line=1, at=1, cex=0.8, las=2)
plot(matrix(c(M1HB, M1xHB, G1HB, G1xHB, M1G1_HB, M1G1_xHB), nrow=2), main="QTL Counts\nMON1 vs GSD1", digits = TRUE, key=NULL, fmt.cell='%.0f', axis.col=NULL, axis.row=NULL, xlab='', ylab='', col=NULL)
mtext('M1', side=1, line=1, at=1, cex=0.8); mtext('G1', side=1, line=1, at=2, cex=0.8); mtext('SH', side=1, line=1, at=3, cex=0.8); mtext('W', side=2, line=1, at=2, cex=0.8, las=2); mtext('O', side=2, line=1, at=1, cex=0.8, las=2)
plot(matrix(c(M1HB, M1xHB, G2HB, G2xHB, M1G2_HB, M1G2_xHB), nrow=2), main="QTL Counts\nMON1 vs GSD2", digits = TRUE, key=NULL, fmt.cell='%.0f', axis.col=NULL, axis.row=NULL, xlab='', ylab='', col=NULL)
mtext('M1', side=1, line=1, at=1, cex=0.8); mtext('G2', side=1, line=1, at=2, cex=0.8); mtext('SH', side=1, line=1, at=3, cex=0.8); mtext('W', side=2, line=1, at=2, cex=0.8, las=2); mtext('O', side=2, line=1, at=1, cex=0.8, las=2)
plot(matrix(c(M2HB, M2xHB, G1HB, G1xHB, M2G1_HB, M2G1_xHB), nrow=2), main="QTL Counts\nMON2 vs GSD1", digits = TRUE, key=NULL, fmt.cell='%.0f', axis.col=NULL, axis.row=NULL, xlab='', ylab='', col=NULL)
mtext('M2', side=1, line=1, at=1, cex=0.8); mtext('G1', side=1, line=1, at=2, cex=0.8); mtext('SH', side=1, line=1, at=3, cex=0.8); mtext('W', side=2, line=1, at=2, cex=0.8, las=2); mtext('O', side=2, line=1, at=1, cex=0.8, las=2)
plot(matrix(c(M2HB, M2xHB, G2HB, G2xHB, M2G2_HB, M2G2_xHB), nrow=2), main="QTL Counts\nMON2 vs GSD2", digits = TRUE, key=NULL, fmt.cell='%.0f', axis.col=NULL, axis.row=NULL, xlab='', ylab='', col=NULL)
mtext('M2', side=1, line=1, at=1, cex=0.8); mtext('G2', side=1, line=1, at=2, cex=0.8); mtext('SH', side=1, line=1, at=3, cex=0.8); mtext('W', side=2, line=1, at=2, cex=0.8, las=2); mtext('O', side=2, line=1, at=1, cex=0.8, las=2)
# overlap inside haploblocks
hist(randomized_overlaps2$M1M2_HB, xlim=c(0,65), breaks=30, xlab="Number of overlaps", main="Within Inversions\nMON1 vs MON2", las=1); abline(v=M1M2_HB, col="blue", lwd=2)
hist(randomized_overlaps2$G1G2_HB, xlim=c(0,65), breaks=30, xlab="Number of overlaps", main="Within Inversions\nGSD1 vs GSD2", las=1); abline(v=G1G2_HB, col="red", lwd=2)
hist(randomized_overlaps2$M1G1_HB, xlim=c(0,65), breaks=30, xlab="Number of overlaps", main="Within Inversions\nMON1 vs GSD1", las=1); abline(v=M1G1_HB, col="red", lwd=2)
hist(randomized_overlaps2$M1G2_HB, xlim=c(0,65), breaks=60, xlab="Number of overlaps", main="Within Inversions\nMON1 vs GSD2", las=1); abline(v=M1G2_HB, col="red", lwd=2)
hist(randomized_overlaps2$M2G1_HB, xlim=c(0,65), breaks=30, xlab="Number of overlaps", main="Within Inversions\nMON2 vs GSD1", las=1); abline(v=M2G1_HB, col="red", lwd=2)
hist(randomized_overlaps2$M2G2_HB, xlim=c(0,65), breaks=30, xlab="Number of overlaps", main="Within Inversions\nMON2 vs GSD2", las=1); abline(v=M2G2_HB, col="blue", lwd=2)
# overlap outside haploblocks
hist(randomized_overlaps2$M1M2_xHB, xlim=c(0,65), breaks=30, xlab="Number of overlaps", main="Outside Inversions\nMON1 vs MON2", las=1); abline(v=M1M2_xHB, col="blue", lwd=2)
hist(randomized_overlaps2$G1G2_xHB, xlim=c(0,65), breaks=30, xlab="Number of overlaps", main="Outside Inversions\nGSD1 vs GSD2", las=1); abline(v=G1G2_xHB, col="red", lwd=2)
hist(randomized_overlaps2$M1G1_xHB, xlim=c(0,65), breaks=30, xlab="Number of overlaps", main="Outside Inversions\nMON1 vs GSD1", las=1); abline(v=M1G1_xHB, col="blue", lwd=2)
hist(randomized_overlaps2$M1G2_xHB, xlim=c(0,65), breaks=60, xlab="Number of overlaps", main="Outside Inversions\nMON1 vs GSD2", las=1); abline(v=M1G2_xHB, col="red", lwd=2)
hist(randomized_overlaps2$M2G1_xHB, xlim=c(0,65), breaks=30, xlab="Number of overlaps", main="Outside Inversions\nMON2 vs GSD1", las=1); abline(v=M2G1_xHB, col="blue", lwd=2)
hist(randomized_overlaps2$M2G2_xHB, xlim=c(0,65), breaks=30, xlab="Number of overlaps", main="Outside Inversions\nMON2 vs GSD2", las=1); abline(v=M2G2_xHB, col="blue", lwd=2)
# difference between overlap within versus outside haploblocks
hist(randomized_overlaps2$M1M2_HB/272*100-randomized_overlaps2$M1M2_xHB/1276*100, xlim=c(-5,20), breaks=30, xlab="% overlap within - % overlap outside", main="Within vs. Outside\nMON1 vs MON2", las=1); abline(v=M1M2_HB/272*100-M1M2_xHB/1276*100, col="blue", lwd=2)
hist(randomized_overlaps2$G1G2_HB/272*100-randomized_overlaps2$G1G2_xHB/1276*100, xlim=c(-5,20), breaks=30, xlab="% overlap within - % overlap outside", main="Within vs. Outside\nGSD1 vs GSD2", las=1); abline(v=G1G2_HB/272*100-G1G2_xHB/1276*100, col="blue", lwd=2)
hist(randomized_overlaps2$M1G1_HB/272*100-randomized_overlaps2$M1G1_xHB/1276*100, xlim=c(-5,20), breaks=30, xlab="% overlap within - % overlap outside", main="Within vs. Outside\nMON1 vs GSD1", las=1); abline(v=M1G1_HB/272*100-M1G1_xHB/1276*100, col="red", lwd=2)
hist(randomized_overlaps2$M1G2_HB/272*100-randomized_overlaps2$M1G2_xHB/1276*100, xlim=c(-5,20), breaks=60, xlab="% overlap within - % overlap outside", main="Within vs. Outside\nMON1 vs GSD2", las=1); abline(v=M1G2_HB/272*100-M1G2_xHB/1276*100, col="red", lwd=2)
hist(randomized_overlaps2$M2G1_HB/272*100-randomized_overlaps2$M2G1_xHB/1276*100, xlim=c(-5,20), breaks=30, xlab="% overlap within - % overlap outside", main="Within vs. Outside\nMON2 vs GSD1", las=1); abline(v=M2G1_HB/272*100-M2G1_xHB/1276*100, col="red", lwd=2)
hist(randomized_overlaps2$M2G2_HB/272*100-randomized_overlaps2$M2G2_xHB/1276*100, xlim=c(-5,20), breaks=30, xlab="% overlap within - % overlap outside", main="Within vs. Outside\nMON2 vs GSD2", las=1); abline(v=M2G2_HB/272*100-M2G2_xHB/1276*100, col="blue", lwd=2)
# save as 8x10 - Figure S14
