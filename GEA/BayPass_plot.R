## Plot BayPass GEA results (inv_gt+SNPs), shared variables only
## Kaichi Huang 2022 Dec

library(tidyverse)
library(gridExtra)
library(grid)
select=dplyr::select

my.inv_info <- read_tsv("inv.txt")
my.inv_info$chr <- factor(my.inv_info$chr)

variable <- "percent.cover"
my.data <- read.table("all_summary_betai_reg.out", h=T)
my.pod <- read.table("all_pod_summary_betai_reg.out", h=T)
my.snp <- read.table("MON_GEA.baypass.snp", h=F)
my.cov <- read.table("env.baypass.cov", h=F, sep="\t")
my.inversion <- read.table("inv_all_summary_betai_reg.out"), h=T)
names(my.snp) <- c("CHR", "POS")
my.snp$CHR <- factor(as.numeric(sub("Ha412HOChr","", my.snp$CHR)))
my.snp$POS <- my.snp$POS/(10^6)
my.region <- my.inv_info
my.region$chr <- sub("Ha412HOChr","", my.region$chr)
my.region$chr <- factor(as.numeric(my.region$chr))
my.region$start <- my.region$start/(10^6)
my.region$end <- my.region$end/(10^6)

bp <- my.snp$POS[my.snp$CHR == levels(my.snp$CHR)[1]]
mid <- 0+max(bp)/2
my.region$bp1[my.region$chr == levels(my.snp$CHR)[1]] = my.region$start[my.region$chr == levels(my.snp$CHR)[1]]
my.region$bp2[my.region$chr == levels(my.snp$CHR)[1]] = my.region$end[my.region$chr == levels(my.snp$CHR)[1]]
for (i in 2:length(levels(my.snp$CHR))) {
	max.now <- max(bp)
	mid <- c(mid, max.now+max(my.snp$POS[my.snp$CHR == i])/2)
	bp <- c(bp, my.snp$POS[my.snp$CHR == i] + max.now)
	my.region$bp1[my.region$chr == i] = my.region$start[my.region$chr == i] + max.now
	my.region$bp2[my.region$chr == i] = my.region$end[my.region$chr == i] + max.now
}
my.snp$bp <- bp
my.snp$MRK <- 1:nrow(my.snp)
my.region$MRK <- 1:nrow(my.region)

my.plot <- my.data %>%
	filter(COVARIABLE == 2) %>% inner_join(my.snp) %>% select(CHR, bp, BF.dB.)
my.plot2 <- my.inversion %>%
	filter(COVARIABLE == 2) %>% inner_join(my.region) %>% select(chr, bp1, bp2, BF.dB.)
my.q <- my.pod %>%
	filter(COVARIABLE == 2) %>% summarize(q99=quantile(BF.dB.,probs=0.99, na.rm=T)) %>% pull()
pdf("BayPass_plot.pdf", width=8, height=3)
print(
	ggplot() + theme_classic() +
		geom_point(data=my.plot, mapping=aes(x=bp, y=BF.dB., col=CHR), pch=20, na.rm=T, show.legend=FALSE) +
		geom_segment(data=my.plot2, mapping=aes(x=bp1,xend=bp2,y=BF.dB.,yend=BF.dB.), col="red", size=2) +
		geom_hline(aes(yintercept=my.q), col="red", alpha=0.7, linetype="dotted") +
		xlab("Chromosome") + ylab(expression(BF[is]*" "*(dB))) +
		scale_colour_manual(values=rep(c("grey60","grey90"),9)) +
		scale_x_continuous(breaks=mid, labels=levels(my.plot$CHR), expand=c(0.02,0.02)) +
		scale_y_continuous(expand=c(0,0), limits=c(-12,135)) +
		theme(axis.line.x=element_blank()) +
		ggtitle(paste("MON",variable))
)
dev.off()
