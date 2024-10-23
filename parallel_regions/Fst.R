## Format Fst output
## Kaichi Huang 2020 May

library(tidyverse)

my.fst <- read_tsv("PetFal.windowed.weir.fst", col_names=T) %>%
	mutate(POS=(BIN_START+BIN_END)/2, ID=paste(CHROM,BIN_START,BIN_END,sep="_")) %>%
	select(ID,CHROM,POS,WEIGHTED_FST)
my.fst$CHROM <- sub("Ha412HOChr","", my.fst$CHROM)
my.fst$CHROM <- factor(as.numeric(my.fst$CHROM))
my.fst <- my.fst %>% select(CHROM,POS,WEIGHTED_FST)
write_tsv(my.fst, "PetFal.fst.txt")
