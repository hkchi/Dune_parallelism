library(tidyverse)
library(janitor)

# data sets needed from Dryad: haploblocks.petfal.txt, 
# data sets needed from previous scripts: 
# MON1_SNPtable.csv, MON2_SNPtable.csv, GSD1_SNPtable.csv, GSD2_SNPtable.csv


################################
####---- Bin genotypes ---- ####
################################

# define a function that finds the mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# load data
temp1 <- read.csv("output/MON1_SNPtable.csv") %>% t() %>% as.data.frame()
temp1 <- read.csv("output/MON2_SNPtable.csv") %>% t() %>% as.data.frame()
temp1 <- read.csv("output/GSD1_SNPtable.csv") %>% t() %>% as.data.frame()
temp1 <- read.csv("output/GSD2_SNPtable.csv") %>% t() %>% as.data.frame()

# bin markers within 1000 bp and find mode for genotypes
temp2 <- temp1 %>% mutate(chr = as.factor(str_sub(row.names(.), 11, 12)), pos = str_sub(row.names(.), 14, -1)) %>%
  mutate(bin = round(as.numeric(pos), -3))
temp3 <- temp2 %>% group_by(chr, bin) %>% summarise_each_(., funs(getmode), list(quote(c(-pos)))) %>% as.data.frame() %>% filter(!is.na(bin))

# write binned data
write.table(temp3, file="output/MON1_bined_SNPtable.csv", sep=",", row.names = FALSE, col.names = TRUE, quote=FALSE)
write.table(temp3, file="output/MON2_bined_SNPtable.csv", sep=",", row.names = FALSE, col.names = TRUE, quote=FALSE)
write.table(temp3, file="output/GSD1_bined_SNPtable.csv", sep=",", row.names = FALSE, col.names = TRUE, quote=FALSE)
write.table(temp3, file="output/GSD2_bined_SNPtable.csv", sep=",", row.names = FALSE, col.names = TRUE, quote=FALSE)


####################################
####---- Merge data tables ---- ####
####################################

# load all binned data
SS_dat1 <- read.csv("output/MON1_bined_SNPtable.csv", stringsAsFactors = FALSE, na.strings = c("NA")) %>% mutate(bin_pos = paste0(chr, "_", bin)) %>% select(-chr, -bin)
SS_dat2 <- read.csv("output/MON2_bined_SNPtable.csv", stringsAsFactors = FALSE, na.strings = c("NA")) %>% mutate(bin_pos = paste0(chr, "_", bin)) %>% select(-chr, -bin)
SS_dat3 <- read.csv("output/GSD1_bined_SNPtable.csv", stringsAsFactors = FALSE, na.strings = c("NA")) %>% mutate(bin_pos = paste0(chr, "_", bin)) %>% select(-chr, -bin)
SS_dat4 <- read.csv("output/GSD2_bined_SNPtable.csv", stringsAsFactors = FALSE, na.strings = c("NA")) %>% mutate(bin_pos = paste0(chr, "_", bin)) %>% select(-chr, -bin)

# merge all 4 data sets based on bins
temp <- inner_join(SS_dat1, SS_dat2, by="bin_pos") %>% 
  inner_join(., SS_dat3, by="bin_pos") %>%
  inner_join(., SS_dat4, by="bin_pos") %>%
  relocate(bin_pos) %>% t() %>% as.data.frame() %>% row_to_names(1)
temp$pop <- rep(c("MON1", "MON2", "GSD1", "GSD2"), each=48)
temp$size <- as.factor(c(str_sub(row.names(temp)[1:96], -1, -1), rep(rep(c("B", "S"), each=12), 4)))
merged_table <- temp
merged_table[,c(1, 1549:1551)]


#############################################
#### ---- Run Fisher's extact tests ---- ####
#############################################

# Run individual fisher's exact tests on binned data for all mapping populations individually 
similarity_outputs <- data.frame(marker=NA, MON1_p=NA, MON2_p=NA, GSD1_p=NA, GSD2_p=NA, MON1_estimate=NA, MON2_estimate=NA, GSD1_estimate=NA, GSD2_estimate=NA)
for (i in seq(1, ncol(merged_table)-3, by=1)) {
  print(i)
  temp_dat <- merged_table %>% select(i, "pop", "size")
  similarity_outputs[i,1] <- names(temp_dat)[1]
  names(temp_dat) <- c("marker", "pop", "size")
  for (j in 1:4){
    p <- c("MON1", "MON2", "GSD1", "GSD2")[j]
    Tdat_B <- temp_dat %>% filter(pop==p) %>% filter(size=="B") %>% filter(complete.cases(.))
    Tdat_S <- temp_dat %>% filter(pop==p) %>% filter(size=="S") %>% filter(complete.cases(.))
    D_alleles_B <- sum(as.numeric(Tdat_B$marker), na.rm=T)
    D_alleles_S <- sum(as.numeric(Tdat_S$marker), na.rm=T)
    N_alleles_B <- nrow(Tdat_B)*2-D_alleles_B
    N_alleles_S <- nrow(Tdat_S)*2-D_alleles_S
    similarity_outputs[i,j+1] <- fisher.test(matrix(c(D_alleles_B, N_alleles_B, D_alleles_S, N_alleles_S),2))$p.value
    similarity_outputs[i,j+5] <- D_alleles_B/(D_alleles_B+N_alleles_B)-D_alleles_S/(D_alleles_S+N_alleles_S)
  }
}

# add CHR and bp info
similarity_outputs <- similarity_outputs %>% mutate(CHR = str_split_fixed(similarity_outputs$marker, "_", 2)[,1], 
                                                    BP = str_split_fixed(similarity_outputs$marker, "_", 2)[,2],
                                                    MON1_sig = MON1_p < 0.05, MON2_sig = MON2_p < 0.05, 
                                                    GSD1_sig = GSD1_p < 0.05, GSD2_sig = GSD2_p < 0.05) 
similarity_outputs$CHR <- as.numeric(similarity_outputs$CHR)
similarity_outputs$BP <- as.numeric(similarity_outputs$BP)
similarity_outputs$sig_count <- similarity_outputs %>% select(MON1_sig, MON2_sig, GSD1_sig, GSD2_sig) %>% rowSums()

# determine which regions are within haploblocks
haploblocks <- read.table(file="data/haploblocks.petfal.txt", header=T)
haploblock_list <- rep(F, 1548)
for (i in 1:11) {
  temp <- similarity_outputs$CHR==haploblocks$chr[i] & similarity_outputs$BP>haploblocks$start[i] & similarity_outputs$BP<haploblocks$end[i]
  haploblock_list <- mapply("|", temp, haploblock_list)
}
similarity_outputs$HB <- haploblock_list

write.table(similarity_outputs, file="output/similarity_outputs.csv", sep=",", row.names = FALSE, col.names = TRUE, quote=FALSE)

