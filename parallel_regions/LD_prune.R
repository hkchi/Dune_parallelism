## LD pruning of the selected regions
## Mojtaba Jahani 2021 Jan, Kaichi Huang 2024 Jul

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
library(fuzzyjoin)

my.region.file <- "parallel.selected_region.bed"
my.snp.file <- "PetFal.snps_bi.ID_LIST"
my.ld_threshold.file <- "petfal_threshold.95"
my.map.file <- "geneticmap.txt"
my.vcf <- "PetFal.snps_bi.vcf.gz"

my.region <- read_tsv(my.region.file, col_names=F) %>%
  select(chr=X1,start=X2,end=X3) %>%
  mutate(Window=paste0(chr,":",start,"-",end)) %>%
  mutate(pos=(((end-start)/2)+start))

my.snp <- read_tsv(my.snp.file, col_names=F) %>%
  select(snp_id=X1) %>%
  separate(snp_id, c("CHR","POS"), sep=":", remove=F)
my.snp$POS <- as.numeric(my.snp$POS)
# Subset SNPs
registerDoParallel(cores=48)
my.snp.region <- foreach(i=1:nrow(my.region), .combine='rbind', .errorhandling='stop') %dopar% {
  mutate(select(filter(my.snp,
                       CHR == as.character(my.region[i,1]),
                       (POS >= as.numeric(my.region[i,2]) &
                          POS <= as.numeric(my.region[i,3]))),
  snp_id),
  Window=as.character(my.region[i,4]))}

my.ld_threshold <- read_tsv(my.ld_threshold.file, col_names=F) %>%
  select(P_distance=X1, LD_threshold=X2)

my.map <- read_tsv(my.map.file)
my.recombination <- my.map %>%
  group_by(chr) %>%
  arrange(pos, .by_group=T) %>%
  mutate(delta.cM=cM-lag(cM, default=first(cM)),
         delta.bp=pos-lag(pos, default=first(pos)),
         start=lag(pos)) %>% 
  ungroup() %>%
  mutate(recomb_rate=delta.cM/delta.bp) %>%
  select(chr,start,end=pos,recomb_rate) %>%
  filter(!is.na(start)) %>%
  left_join(.,rename(my.map,start=pos))
# Calculate genetic position
registerDoParallel(cores=48) 
my.recombination.region <- foreach(row=1:nrow(my.recombination),.combine='rbind',.errorhandling='stop') %dopar% {
  mutate(select(filter(my.region,
                       chr == as.character(my.recombination[row,1]) &
                        pos >= as.numeric(my.recombination[row,2]) &
                        pos <= as.numeric(my.recombination[row,3])),
                Window, pos),
         start=as.numeric(my.recombination[row,2]), recomb_rate=as.numeric(my.recombination[row,4]), cM=as.numeric(my.recombination[row,5]))}
my.recombination.region <- my.recombination.region %>% 
  mutate(G_position=(((as.numeric(pos)-as.numeric(start))*as.numeric(recomb_rate))+as.numeric(cM))) %>%
  select(Window, G_position)

for (CHROMO in unique(my.region$chr)) {
  my.snp.region.chr <- my.snp.region %>% filter(grepl(CHROMO,snp_id)) %>% select(snp_id)
  write_tsv(my.snp.region.chr, paste0("tmp.my.snp.region.",CHROMO), col_names=F)
  system(paste0("plink --r2 --vcf ",my.vcf," --set-missing-var-ids @:# --extract tmp.my.snp.region.",CHROMO,
                " --ld-window-kb 999999999 --ld-window 999999999 --ld-window-r2 0 --allow-extra-chr --out ",
                "tmp.my.snp.region.",CHROMO))
  my.ld.win <- read_table(paste0("tmp.my.snp.region.",CHROMO,".ld")) %>% select(SNP_A, SNP_B, R2) %>%
    inner_join(., rename(my.snp.region, SNP_A=snp_id, WIN_A=Window)) %>%
    inner_join(., rename(my.snp.region, SNP_B=snp_id, WIN_B=Window)) %>%
    filter(!WIN_A==WIN_B) %>%
    group_by(WIN_A, WIN_B) %>%
    summarise(LD=mean(R2), snp_count=n()) %>%
    ungroup() %>%
    select(WIN_A, WIN_B, LD, snp_count)
  write_tsv(my.ld.win, paste0("tmp.my.snp.region.",CHROMO,".LD_WINDOW"), col_names=T)
  # Compare to the thresholds
  my.ld.win.sig <- my.ld.win %>%
    separate(WIN_A, c("chr_A","start_end_A"), sep=":", remove=F) %>%
    separate(start_end_A, c("start_A","end_A"), sep="-") %>%
    separate(WIN_B, c("chr_B","start_end_B"), sep=":", remove=F) %>%
    separate(start_end_B, c("start_B","end_B"), sep="-") %>%
    mutate(mid_A=(((as.numeric(end_A)-as.numeric(start_A))/2)+as.numeric(start_A))) %>%
    mutate(mid_B=(((as.numeric(end_B)-as.numeric(start_B))/2)+as.numeric(start_B))) %>%
    mutate(P_distance=abs(as.numeric(mid_A)-as.numeric(mid_B))) %>% 
    select(WIN_A, WIN_B, mid_A, LD, P_distance) %>%
    difference_left_join(., my.ld_threshold, max_dist=5e3/2) %>%
    mutate(delta.P_distance=abs(as.numeric(P_distance.x)-as.numeric(P_distance.y))) %>%
    mutate(LD_significant=(LD>=LD_threshold)) %>% 
    left_join(., rename(my.recombination.region, WIN_A=Window, G_position_A=G_position)) %>%
    left_join(., rename(my.recombination.region, WIN_B=Window, G_position_B=G_position)) %>%
    mutate(G_distance=abs(as.numeric(G_position_A)-as.numeric(G_position_B))) %>%
    select(-G_position_A, -G_position_B) %>%
    mutate(GD_significant=(G_distance<=5)) %>%
    mutate(LD_GD_sig = ifelse(LD_significant == "TRUE" & GD_significant == "TRUE",
                              "TRUE",
                              "FALSE"))
  write_tsv(my.ld.win.sig, paste0("tmp.my.snp.region.",CHROMO,".LD_WINDOW_SIG"), col_names=T)
  #system(paste0("rm tmp.my.snp.region.",CHROMO," tmp.my.snp.region.",CHROMO,".log"," tmp.my.snp.region.",CHROMO,".nosex"," tmp.my.snp.region.",CHROMO,".ld"))
}

# Clustering
for (CHROMO in unique(my.region$chr)) {
  my.region.chr <- my.snp.region %>% filter(grepl(CHROMO,snp_id)) %>% select(Window) %>%
    unique() %>%
    separate(Window, c("chrom", "start_end"), sep=":", remove=F) %>%
    separate(start_end, c("start","end"), sep="-") %>%
    arrange(as.numeric(start)) %>%
    select(Window)
  my.ld.win.sig <- read_tsv(paste0("tmp.my.snp.region.",CHROMO,".LD_WINDOW_SIG"), col_names=T) %>%
    select(window_A=WIN_A, window_B=WIN_B, LD_GD_sig)
  
  cluster <- NULL
  # LD clustering steps
  if (nrow(my.region.chr) > 1) {
    for (i in 1:(nrow(my.region.chr)-1)) {
      cluster <- rbind(data.frame(window_A = as.character(my.region.chr[i,1]),
                       window_B = as.character(my.region.chr[i+1,1]),
                       cycle = i),
                       cluster)
    }
  } else {
    cluster <- rbind(data.frame(window_A = my.region.chr$Window[1], 
                                window_B = "-",
                                cycle = 1),
                     cluster)
  }
  cluster <- cluster %>%
    mutate(window_A = as.character(window_A)) %>%
    mutate(window_B = as.character(window_B)) %>%
    left_join(., my.ld.win.sig) %>%
    arrange(cycle) %>%
    mutate(clustering=0)
  for (p in 1:nrow(cluster)) {
    cluster[p,5] <- ifelse(cluster[p,3]==1, 
                           1,
                           cluster[p,5] <- ifelse(cluster[p,5]==0 & cluster[p-1,4]=="TRUE" & cluster[p,4]=="TRUE",
                                                  cluster[p-1,5],
                                                  as.numeric(cluster[p-1,5])+1))
  }
  true_cluster <- cluster %>% filter(LD_GD_sig=="TRUE") %>%
    gather(windows, window_name, 1:2) %>%
    select(window_name, clustering) %>%
    distinct(window_name,.keep_all=T)
  my.cluster.chr <- my.region.chr %>% filter(!Window %in% pull(distinct(true_cluster,window_name),window_name)) %>%
    mutate(clustering=Window) %>%
    select(window_name=Window, clustering) %>%
    rbind(.,true_cluster) %>%
    mutate(cluster=group_indices(.,clustering)) %>%
    select(window_name, cluster)
  write_tsv(my.cluster.chr, paste0("tmp.my.snp.region.",CHROMO,".cluster"), col_names=T)
}

# Output the pruned result
my.cluster.merge <- data.frame()
my.cluster.prune <- data.frame()
for (CHROMO in unique(my.region$chr)) {
  my.cluster.chr <- read_tsv(paste0("tmp.my.snp.region.",CHROMO,".cluster"), col_names=T)
  # Merge
  my.cluster.chr.merge <- my.cluster.chr %>%
    separate(window_name, c("chr","start_end"), sep=":", remove=F) %>%
    separate(start_end, c("start","end"), sep="-") %>%
    mutate(start=as.numeric(start), end=as.numeric(end)) %>%
    group_by(chr, cluster) %>%
    summarise(N_windows=n(), cluster_start=min(start), cluster_end=max(end)) %>%
    ungroup() %>%
    arrange(cluster_start) %>%
    mutate(range=paste0(as.character(cluster_start),":",as.character(cluster_end))) %>%
    mutate(size=as.numeric(cluster_end)-as.numeric(cluster_start)+1) %>%
    select(chr, range, size, N_windows)
  my.cluster.merge <- rbind(my.cluster.merge, my.cluster.chr.merge)
  # Prune (keep the largest window in each cluster)
  my.cluster.chr.prune <- my.cluster.chr %>%
    separate(window_name, c("chr","start_end"), sep=":", remove=F) %>%
    separate(start_end, c("start","end"), sep="-") %>%
    mutate(start=as.numeric(start), end=as.numeric(end)) %>%
    mutate(size=end-start+1) %>%
    group_by(chr, cluster) %>%
    mutate(N_windows=n(), cluster_start=min(start), cluster_end=max(end), max_size=max(size)) %>%
    filter(size==max_size) %>% sample_n(1) %>%
    ungroup() %>%
    mutate(range=paste0(as.character(cluster_start),":",as.character(cluster_end))) %>%
    select(chr, start, end, range, N_windows) %>%
    arrange(start)
  my.cluster.prune <- rbind(my.cluster.prune, my.cluster.chr.prune)
}
#write_tsv(my.cluster.merge, paste0(my.region.file,".clustering_result"), col_names=T)
write_tsv(my.cluster.prune, paste0(my.region.file,".pruning_result"), col_names=T)
#system("rm tmp.my.snp.region.*")
