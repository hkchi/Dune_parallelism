## Enrichment analysis of the selected regions
## Mojtaba Jahani 2021 Jan, Kaichi Huang 2024 Jul

library(tidyverse)
library(foreach)
library(doParallel)

# Info
my.haploblock <- read_tsv("haploblocks.petfal.txt") %>% select(-name)
my.bed <- my.haploblock[1:11,]
my.bed.HA412 <- "HA412.chr_len.txt"
my.region.file <- "parallel.selected_region.bed.pruning_result"

n.permutation <- 10000
my.out <- matrix(nrow=1, ncol=10)

# Overlap
my.region <- read_tsv(my.region.file, col_names=T) %>%
  select(chr, start, end) %>%
  mutate(size=(end-start)+1)

registerDoParallel(cores=48)
my.overlap <- foreach(i=1:nrow(my.bed), .combine='rbind', .errorhandling='stop') %dopar% {
  mutate(select(filter(my.region,
                       chr == as.character(my.bed[i,1]) &
                         ((start >= as.numeric(my.bed[i,2]) & start <= as.numeric(my.bed[i,3])) | 
                            (end >= as.numeric(my.bed[i,2]) & end <= as.numeric(my.bed[i,3])) | 
                            (start < as.numeric(my.bed[i,2]) & end > as.numeric(my.bed[i,3])))),
                chr,start,end,size),
         inversion_start=as.numeric(my.bed[i,2]),
         inversion_end=as.numeric(my.bed[i,3]))}
my.overlap <- my.overlap %>% 
      mutate(overlap_start=ifelse(start>=inversion_start,start,inversion_start)) %>%
      mutate(overlap_end=ifelse(end>=inversion_end,inversion_end,end)) %>%
      mutate(overlap_size=(overlap_end-overlap_start)+1)

my.out[1,1] <- as.numeric(nrow(my.region))
my.out[1,2] <- as.numeric(nrow(my.overlap))
my.out[1,3] <- my.out[1,2]/my.out[1,1]
my.out[1,6] <- as.numeric(summarise(my.region,overall.size.region=sum(size)))
my.out[1,7] <- as.numeric(summarise(my.overlap,overall.size.overlap=sum(overlap_size)))
my.out[1,8] <- my.out[1,7]/my.out[1,6]

# Permutate and overlap
system("mkdir -p ./tmp.enrichment")
write_tsv(my.bed, "./tmp.enrichment/inv.bed",col_names=F)
permutation_table <- matrix(nrow=n.permutation, ncol=2)
for (n in 1:n.permutation) {
  system(paste0("bedtools shuffle -i ./tmp.enrichment/inv.bed -g ",my.bed.HA412," -noOverlapping > ./tmp.enrichment/perm",n,"_inv.bed")) 
  my.bed.n <- read_tsv(paste0("./tmp.enrichment/perm",n,"_inv.bed"), col_names=F) %>%
    rename(chr=X1, start=X2, end=X3)
  registerDoParallel(cores=48)
  my.overlap.n <- foreach(i=1:nrow(my.bed.n), .combine='rbind', .errorhandling='stop') %dopar% {
    mutate(select(filter(my.region,
                         chr == as.character(my.bed.n[i,1]) &
                           ((start >= as.numeric(my.bed.n[i,2]) & start <= as.numeric(my.bed.n[i,3])) | 
                              (end >= as.numeric(my.bed.n[i,2]) & end <= as.numeric(my.bed.n[i,3])) | 
                              (start < as.numeric(my.bed.n[i,2]) & end > as.numeric(my.bed.n[i,3])))),
                  chr,start,end,size),
           inversion_start=as.numeric(my.bed.n[i,2]),
           inversion_end=as.numeric(my.bed.n[i,3]))}
  my.overlap.n <- my.overlap.n %>% 
    mutate(overlap_start=ifelse(start>=inversion_start,start,inversion_start)) %>%
    mutate(overlap_end=ifelse(end>=inversion_end,inversion_end,end)) %>%
    mutate(overlap_size=(overlap_end-overlap_start)+1)
  
  permutation_table[n,1] <- as.numeric(nrow(my.overlap.n))/my.out[1,1]
  permutation_table[n,2] <- as.numeric(summarise(my.overlap.n,overlap_overall_size=sum(overlap_size)))/my.out[1,6] 
  
  rm(my.bed.n, my.overlap.n)
  system(paste0("rm ./tmp.enrichment/perm",n,"_inv.bed"))
}
system("rm -r ./tmp.enrichment")

my.out[1,4] <- mean(permutation_table[,1])
n.ex.num <- permutation_table %>%
  data.frame() %>%
  filter(X1 > my.out[1,3]) %>%
  nrow()
my.out[1,5] <- n.ex.num / n.permutation
my.out[1,9] <- mean(permutation_table[,2])
n.ex.prop <- permutation_table %>%
  data.frame() %>%
  filter(X2 > my.out[1,8]) %>%
  nrow()
my.out[1,10] <- n.ex.prop / n.permutation

my.out <- my.out %>% data.frame() %>%
  rename(N_SW_window = X1,
         N_overlap_SW_window = X2,
         porportion_N_overlap_SW_window = X3,
         mean_Null_porportion_N_overlap_SW_window = X4,
         P_porportion_N_overlap_SW_window = X5,
         Length_SW_window = X6,
         Length_overlap_SW_window = X7,
         porportion_length_overlap_SW_window = X8,
         mean_Null_porportion_length_overlap_SW_window = X9,
         P_porportion_length_overlap_SW_window = X10)
write_tsv(my.out, "enrichment_cM.txt")
