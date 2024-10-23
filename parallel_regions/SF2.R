## Format SweepFinder2 output
## Kaichi Huang 2020 Jun

library(tidyverse)

my.sf2 <- tibble()
for (i in 1:17) {
	chr <- paste0("Ha412HOChr",sprintf("%02d", i))
	tmp <- read_tsv("PetFal_Dune",chr,"SF2.txt",sep="."),col_names=T)
	tmp$chr <- i
	my.sf2 <- rbind(my.sf2,tmp)
}
write_tsv(my.sf2, "PetFal.CLR.txt")
