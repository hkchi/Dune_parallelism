library(tidyverse)
library(qqman)
library(SNPRelate)
library(gridExtra)

# data sets needed from Dryad: haploblocks.petfal.txt, inv_gt.WGS_RF.NEG1.txt, inv_gt.WGS_RF.NEG2.txt, inv_gt.WGS_RF.PET1.txt, inv_gt.WGS_RF.PET2.txt
# data sets needed from previous scripts: 
# MON1_filt.vcf.gz, MON2_filt.vcf.gz, GSD1_filt.vcf.gz, GSD2_filt.vcf.gz

# Load haploblock info
my.haploblock <- read_tsv("data/haploblocks.petfal.txt")
my.haploblock$chr <- factor(my.haploblock$chr)

#############################################
#### ---- Run Fisher's extact tests ---- ####
#############################################

# run tests on each SNP

# load data
SNPtable <- read.csv("output/MON1_SNPtable.csv")
SNPtable <- read.csv("output/MON2_SNPtable.csv")
SNPtable <- read.csv("output/GSD1_SNPtable.csv")
SNPtable <- read.csv("output/GSD2_SNPtable.csv")

# run stats
fisher_outputs <- data.frame(marker=NA, full_P=NA, estimate=NA, diff=NA) 
for (i in seq(1, ncol(SNPtable)-2, by=1)) {
  print(i)
  temp_dat <- SNPtable %>% select(i, "size")
  fisher_outputs[i,1] <- names(temp_dat)[1]
  names(temp_dat) <- c("marker", "size")
  Tdat_B <- temp_dat %>% filter(size=="B") %>% filter(complete.cases(.))
  Tdat_S <- temp_dat %>% filter(size=="S") %>% filter(complete.cases(.))
  D_alleles_B <- sum(as.numeric(Tdat_B$marker), na.rm=T)
  D_alleles_S <- sum(as.numeric(Tdat_S$marker), na.rm=T)
  N_alleles_B <- nrow(Tdat_B)*2-D_alleles_B
  N_alleles_S <- nrow(Tdat_S)*2-D_alleles_S
  fisher_outputs[i,2] <- fisher.test(matrix(c(D_alleles_B, N_alleles_B, D_alleles_S, N_alleles_S),2))$p.value
  fisher_outputs[i,3] <- fisher.test(matrix(c(D_alleles_B, N_alleles_B, D_alleles_S, N_alleles_S),2))$estimate
  fisher_outputs[i,4] <- D_alleles_B/(D_alleles_B+N_alleles_B)-D_alleles_S/(D_alleles_S+N_alleles_S)
}

# add CHR and bp info
fisher_outputs <- fisher_outputs %>% mutate(CHR = str_sub(str_split_fixed(fisher_outputs$marker, "_", 2)[,1], start=-2L, end=-1L),
                                            BP = str_split_fixed(fisher_outputs$marker, "_", 2)[,2]) %>% filter(complete.cases(.))
fisher_outputs$CHR <- as.numeric(fisher_outputs$CHR)
fisher_outputs$BP <- as.numeric(fisher_outputs$BP)

# save data
write.table(fisher_outputs, file="output/MON1_FEtest_outputs.csv", sep=",", row.names = FALSE, col.names = TRUE, quote=FALSE)
write.table(fisher_outputs, file="output/MON2_FEtest_outputs.csv", sep=",", row.names = FALSE, col.names = TRUE, quote=FALSE)
write.table(fisher_outputs, file="output/GSD1_FEtest_outputs.csv", sep=",", row.names = FALSE, col.names = TRUE, quote=FALSE)
write.table(fisher_outputs, file="output/GSD2_FEtest_outputs.csv", sep=",", row.names = FALSE, col.names = TRUE, quote=FALSE)

# run tests on each inversion locus 

my.fisher_outputs <- tibble(dataset=character(), haploblock=character(), p.value=numeric())
for (dataset in c("PET1", "PET2", "NEG1", "NEG2")) {
  # Sample list
  if (dataset=="PET1" | dataset=="PET2") {
    sample_list <- read_tsv(paste0("data/inv_gt.WGS_RF.",dataset,".txt")) %>% filter(haploblock=="pet01.01") %>% select(name)
    sample_list.big <- grep('_B_', sample_list$name, value=T)
    sample_list.small <-  grep('_S_', sample_list$name, value=T)
  } else {
    sample_list <- read_tsv(paste0("data/inv_gt.WGS_RF.",dataset,".txt")) %>% filter(haploblock=="pet01.01") %>% select(name)
    sample_list.big <- grep('F2.*_B', sample_list$name, value=T)
    sample_list.small <- grep('F2.*_S', sample_list$name, value=T)
  }
  # WGS_RF genotypes
  WGS_RF_genotypes <- read_tsv(paste0("data/inv_gt.WGS_RF.",dataset,".txt"))
  # Examine each haploblock
  for (i in 1:nrow(my.haploblock)) {
    name.i <- my.haploblock$name[i]
    genotype.big <- WGS_RF_genotypes %>% filter(haploblock==name.i) %>% filter(name %in% sample_list.big) %>% pull(WGS_RF_genotype)
    genotype.small <- WGS_RF_genotypes %>% filter(haploblock==name.i) %>% filter(name %in% sample_list.small) %>% pull(WGS_RF_genotype)
    big.hap1 <- sum(genotype.big, na.rm=T)
    big.hap0 <- sum(!is.na(genotype.big))*2 - big.hap1
    small.hap1 <- sum(genotype.small, na.rm=T)
    small.hap0 <-  sum(!is.na(genotype.small))*2 - small.hap1
    my.p <- fisher.test(matrix(c(big.hap1,big.hap0,small.hap1,small.hap0),nrow=2))$p.value # some are quite monomorphic and through warnings
    my.fisher_outputs <- rbind(my.fisher_outputs, data.frame(dataset=dataset, haploblock=name.i, p.value=my.p))
  }
  print(dataset)
  mydata_split <- split(WGS_RF_genotypes, WGS_RF_genotypes$haploblock)
  mytable <- lapply(mydata_split, function(x) table(factor(x$WGS_RF_genotype, levels = 0:2)))
  inv_table <- do.call(rbind, mytable) %>% data.frame() %>% mutate(inv_freq = (X0+0.5*X1)/(X0+X1+X2)) %>% mutate(poly = inv_freq < 0.95 & inv_freq > 0.05)
  print(inv_table)
}


################################################################
#### ---  Find thresholds for sigfinicant associations ---- ####
################################################################
# find effective number of tests for bonferroni correction using PCA (Gao et al. 2008) 

# load data
snpgdsVCF2GDS("output/MON1_filt.vcf.gz", "output/MON1.gds", method ="biallelic.only")
snpgdsVCF2GDS("output/MON2_filt.vcf.gz", "output/MON2.gds", method ="biallelic.only")
snpgdsVCF2GDS("output/GSD1_filt.vcf.gz", "output/GSD1.gds", method ="biallelic.only")
snpgdsVCF2GDS("output/GSD1_filt.vcf.gz", "output/GSD2.gds", method ="biallelic.only")
genofile <- snpgdsOpen("output/MON1.gds")
genofile <- snpgdsOpen("output/MON2.gds")
genofile <- snpgdsOpen("output/GSD1.gds")
genofile <- snpgdsOpen("output/GSD2.gds")

# run PCA
pca <- snpgdsPCA(genofile, autosome.only = FALSE, eigen.method = "DSPEV")

# determine how many tests 
cumPVE <- cumsum(pca$varprop)
-log10(0.05/(sum(cumPVE < 0.999)+1)) 
# MON1 threshold is 2.973128 
# all other thresholds are 2.991226 


##########################################
#### ---  Plot data (Figure S13) ---- ####
##########################################

plot.list <- list()
for (cross.dataset in c("PET2", "PET1", "NEG1", "NEG2")) {
  # Sample list
  if (cross.dataset=="PET1") {
    cross.type <- "PetSS_5"
    cross.name <- "GSD2"
    my.threshold <- 2.991226
  } else if (cross.dataset=="PET2") {
    cross.type <- "PetSS_7"
    cross.name <- "GSD1"
    my.threshold <- 2.991226
  } else if (cross.dataset=="NEG1") {
    cross.type <- "NegSS_A"
    cross.name <- "MON1"
    my.threshold <- 2.973128
  } else if (cross.dataset=="NEG2") {
    cross.type <- "NegSS_B"
    cross.name <- "MON2"
    my.threshold <- 2.991226
  }
  
  my.p <- read_csv(paste0("output/",cross.name, "_FEtest_outputs.csv"), col_names=T) %>%
    select(marker, full_P) %>%
    separate(marker, c("chr","pos"),"_") %>%
    filter(chr %in% c(paste0("Ha412HOChr0",1:9),paste0("Ha412HOChr1",0:7)))
  my.e <- read_csv(paste0("output/",cross.name, "_FEtest_outputs.csv"), col_names=T) %>%
    select(marker, diff) %>%
    separate(marker, c("chr","pos"),"_") %>%
    filter(chr %in% c(paste0("Ha412HOChr0",1:9),paste0("Ha412HOChr1",0:7)))
  my.p$chr <- sub("Ha412HOChr","", my.p$chr)
  my.p$chr <- factor(as.numeric(my.p$chr))
  my.p$pos <- as.numeric(my.p$pos)/(10^6)
  
  my.region <- my.haploblock
  my.region$chr <- sub("Ha412HOChr","", my.region$chr)
  my.region$chr <- factor(as.numeric(my.region$chr))
  my.region$start <- my.region$start/(10^6)
  my.region$end <- my.region$end/(10^6)
  
  bp <- my.p$pos[my.p$chr == levels(my.p$chr)[1]]
  mid <- 0+max(bp)/2
  my.region$bp1[my.region$chr == levels(my.p$chr)[1]] = my.region$start[my.region$chr == levels(my.p$chr)[1]]
  my.region$bp2[my.region$chr == levels(my.p$chr)[1]] = my.region$end[my.region$chr == levels(my.p$chr)[1]]
  for (i in 2:length(levels(my.p$chr))) {
    max.now <- max(bp)
    mid <- c(mid, max.now+max(my.p$pos[my.p$chr == i])/2)
    bp <- c(bp, my.p$pos[my.p$chr == i] + max.now)
    my.region$bp1[my.region$chr == i] = my.region$start[my.region$chr == i] + max.now
    my.region$bp2[my.region$chr == i] = my.region$end[my.region$chr == i] + max.now
  }
  my.p$bp <- bp
  my.e$bp <- bp
  
  my.inv.poly <- read_csv(paste0("output/",cross.name, "_inv_f_table.csv"), col_names=T)
  
  my.region <- my.fisher_outputs %>% filter(dataset==cross.dataset) %>% select(-dataset, name=haploblock) %>% 
    inner_join(.,my.region) %>% select(chr, bp1, bp2, p.value)
  
  p <- ggplot() + theme_classic() +
    geom_point(data=my.p, mapping=aes(x=bp, y=-log10(full_P), col=chr), pch=20, na.rm=T, show.legend=FALSE) +
    geom_segment(data=my.region, mapping=aes(x=bp1,xend=bp2,y=-log(p.value),yend=-log(p.value)), col=as.numeric(my.inv.poly$poly)+1, linewidth=2) +
    geom_hline(aes(yintercept=my.threshold), col="red", alpha=1, linetype="dotted") +
    xlab("Chromosome") + ylab(expression(-log(italic(p)))) +
    scale_colour_manual(values=rep(c(alpha("deepskyblue",0.2),alpha("grey70",0.2)),9)) +
    scale_x_continuous(breaks=mid, labels=str_pad(levels(my.p$chr),2,pad="0"), expand=c(0.02,0.02)) +
    scale_y_continuous(expand=c(0,0), limits=c(0,14), breaks=seq(0,12,2)) +
    theme(axis.line.x=element_blank()) +
    ggtitle(cross.name)
  plot.list <- c(plot.list, list(p))
  
  p2 <- ggplot() + theme_classic() +
    geom_point(data=my.e, mapping=aes(x=bp, y=diff, col=chr), pch=20, na.rm=T, show.legend=FALSE) + 
    geom_hline(aes(yintercept=0), col="black", alpha=1, linetype="dotted") +
    xlab("Chromosome") + ylab("D allele freq. diff (B-S)") +
    scale_colour_manual(values=rep(c(alpha("deepskyblue",0.2),alpha("grey70",0.2)),9)) +
    scale_x_continuous(breaks=mid, labels=str_pad(levels(my.p$chr),2,pad="0"), expand=c(0.02,0.02)) +
    scale_y_continuous(expand=c(0,0), limits=c(-1,1)) +
    theme(axis.line.x=element_blank()) +
    ggtitle(cross.name)
  plot.list <- c(plot.list, list(p2))
}

pdf("FigS13_FET_output.pdf", width=9, height=10)
print(
  grid.arrange(
    plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],
    plot.list[[5]],plot.list[[6]],plot.list[[7]],plot.list[[8]],
    nrow=4
  )
)
dev.off()









