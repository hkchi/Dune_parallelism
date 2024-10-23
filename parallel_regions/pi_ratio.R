## Calculate pi_ratio
## Kaichi Huang 2020 Mar

library(tidyverse)

my.pi1 <- read_tsv("PetFal.nonDune.windowed.pi", col_names=T) %>%
	mutate(POS=(BIN_START+BIN_END)/2, ID=paste(CHROM,BIN_START,BIN_END,sep="_")) %>%
	select(ID,CHROM,POS,PI)
my.pi2 <- read_tsv("PetFal.Dune.windowed.pi", col_names=T) %>%
	mutate(POS=(BIN_START+BIN_END)/2, ID=paste(CHROM,BIN_START,BIN_END,sep="_")) %>%
	select(ID,CHROM,POS,PI) %>%
	rename(CHROM2=CHROM,POS2=POS,PI2=PI)
my.rod <- inner_join(my.pi1,my.pi2) %>% mutate(ROD=PI/PI2) %>% select(CHROM,POS,ROD)
my.rod$CHROM <- sub("Ha412HOChr","", my.rod$CHROM)
my.rod$CHROM <- factor(as.numeric(my.rod$CHROM))
my.rod <- my.rod %>% mutate(lg=10*log10(ROD)) %>% select(CHROM,POS,lg)
write_tsv(my.rod, "PetFal.pi_ratio.txt")
