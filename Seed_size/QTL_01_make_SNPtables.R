library(tidyverse)
library(vcfR)
library(janitor)
# data sets needed from Dryad: MON_SS_GBS_NGM_GATK_total.vcf.gz, GSD_SS_GBS_NGM_GATK_total.vcf.gz

#########################################
####---- Subset and filter vcfs ---- ####
#########################################

# load one vcf and index
# MON data and indexes:
vcf <- read.vcfR("data/MON_SS_GBS_NGM_GATK_total.vcf.gz", verbose = TRUE)
index <- c(3:6,8,10,12,15,17:19,21,22,24,26,27,30,32,33,39,42,45,47:50,52,53,57,59,61:64,66,67,69:71,75,77,80,82,85,87,89,91:93,100) # MON1
index <- setdiff(2:101, index) # MON2
# GSD data and indexes:
vcf <- read.vcfR("data/GSD_SS_GBS_NGM_GATK_total.vcf.gz", verbose = TRUE)
index <- c(4,5, 54:101) # GSD1
index <- c(2,3, 6:53) # GSD2

# subset individuals in the vcf based on the index file
vcf2 <- vcf[,c(1, index)]
head(vcf2)

# filter vcf for biallelic site with <95% heterozygosity, 15% MAF, and 75% of sites represented 
MAF_dat <- maf(vcf2, 2)
GT <- extract.gt(vcf2, return.alleles = FALSE, IDtoRowNames = FALSE, convertNA = FALSE)
heterozygosity <- rowSums(GT == "1/0" | GT == "0/1", na.rm=T)/rowSums(GT != "./.", na.rm=T)
subset_filter <- (heterozygosity<0.95 & MAF_dat[,4]>0.15 & MAF_dat[,1]>75 & is.biallelic(vcf2)) %>% ifelse(is.na(.), FALSE, .)
vcf3 <- vcf2[subset_filter ,]

write.vcf(vcf3, "output/MON1_filt.vcf.gz")
write.vcf(vcf3, "output/MON2_filt.vcf.gz")
write.vcf(vcf3, "output/GSD1_filt.vcf.gz")
write.vcf(vcf3, "output/GSD2_filt.vcf.gz")


######################################
#### ---  Reformat SNPtables ---- ####
######################################

# load data
vcf <- read.vcfR("output/MON1_filt.vcf.gz", verbose = TRUE)
vcf <- read.vcfR("output/MON2_filt.vcf.gz", verbose = TRUE)
vcf <- read.vcfR("output/GSD1_filt.vcf.gz", verbose = TRUE)
vcf <- read.vcfR("output/GSD2_filt.vcf.gz", verbose = TRUE)

# extract information from vcf
info <- getFIX(vcf)
GT <- extract.gt(vcf, return.alleles = FALSE, IDtoRowNames = FALSE, convertNA = FALSE)

# generate a list of genotype codes to flip based on the allele that comes from the dune parent for each mapping population
# MON1 flip list
flip_list <- GT %>% as.data.frame() %>% 
  mutate(flip = case_when(NegSS_3A_DA_B == "1/1" | NegSS_8B_NA_S == "0/0" ~ "no",
                          NegSS_3A_DA_B == "0/0" | NegSS_8B_NA_S == "1/1" ~ "yes",
                          TRUE ~ "NA")) %>% select(flip)
# MON2 flip list
flip_list <- GT %>% as.data.frame() %>% 
  mutate(flip = case_when(NegSS_2A_DB_B == "1/1" | NegSS_9B_NB_S == "0/0" ~ "no",
                          NegSS_2A_DB_B == "0/0" | NegSS_9B_NB_S == "1/1" ~ "yes",
                          TRUE ~ "NA")) %>% select(flip)
# GSD1 flip list
flip_list <- GT %>% as.data.frame() %>% 
  mutate(flip = case_when(PetSS_1701 == "1/1" | PetSS_1791 == "0/0" ~ "no",
                          PetSS_1701 == "0/0" | PetSS_1791 == "1/1" ~ "yes",
                          TRUE ~ "NA")) %>% select(flip)
# GSD2 flip list
flip_list <- GT %>% as.data.frame() %>% 
  mutate(flip = case_when(PetSS_1554 == "1/1" | PetSS_1508 == "0/0" ~ "no",
                          PetSS_1554 == "0/0" | PetSS_1508 == "1/1" ~ "yes",
                          TRUE ~ "NA")) %>% select(flip)

# flip genotype codes based on flip lists 
GT[flip_list=="yes",] <- gsub("0.0", 2, GT[flip_list=="yes",])
GT[flip_list=="no",] <- gsub("0.0", 0, GT[flip_list=="no",])
GT[flip_list=="yes",] <- gsub("1.1", 0, GT[flip_list=="yes",])
GT[flip_list=="no",] <- gsub("1.1", 2, GT[flip_list=="no",])
GT <- gsub("0.1", 1, GT)
GT <- gsub("...", NA, GT)
temp <- cbind(paste0(info[,1], "_", info[,2]), GT) %>% as.data.frame() %>% filter(flip_list!="NA")

# remove parents from data, transpose table, and make columns for seed size and cytoplasm type 
# MON1 (16302 sites)
SNPtable <- temp %>% select(-NegSS_3A_DA_B, -NegSS_8B_NA_S) %>% t() %>% as.data.frame() %>% row_to_names(1) %>% 
  mutate(size = as.factor(str_sub(row.names(.), -1, -1)), cyto = as.factor(str_sub(row.names(.), -4, -4)))
# MON2 (21258 sites)
SNPtable <- temp %>% select(-NegSS_2A_DB_B, -NegSS_9B_NB_S) %>% t() %>% as.data.frame() %>% row_to_names(1) %>% 
  mutate(size = as.factor(str_sub(row.names(.), -1, -1)), cyto = as.factor(str_sub(row.names(.), -4, -4)))
# GSD1 (11969 sites)
SNPtable <- temp %>% select(-PetSS_1701, -PetSS_1791) %>% t() %>% as.data.frame() %>% row_to_names(1) %>% 
  mutate(size = rep(c("B", "S", "B", "S"), each=12), cyto = rep(c("N", "D"), each=24))
# GSD2 (13854 sites)
SNPtable <- temp %>% select(-PetSS_1554, -PetSS_1508) %>% t() %>% as.data.frame() %>% row_to_names(1) %>% 
  mutate(size = rep(c("B", "S", "B", "S"), each=12), cyto = rep(c("N", "D"), each=24))

write.table(SNPtable, file="output/MON1_SNPtable.csv", sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
write.table(SNPtable, file="output/MON2_SNPtable.csv", sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
write.table(SNPtable, file="output/GSD1_SNPtable.csv", sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
write.table(SNPtable, file="output/GSD2_SNPtable.csv", sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)


