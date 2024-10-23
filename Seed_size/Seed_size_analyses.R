library(tidyverse)
library(nlme)
library(gplots)
# data sets needed from Dryad: MON_seed_size_data.csv, QTL_SS_data.csv


#################################
### MON populations seed size ###
#################################

# load data
d <- read.csv("data/MON_seed_size_data.csv", stringsAsFactors = FALSE, na.strings = c("NA")) 

# fit models of seed size 
model1 <- lme(weight ~ ecotype * experiment, random = ~ 1 | pop, data=d)
model2 <- lme(weight ~ ecotype + experiment, random = ~ 1 | pop, data=d)
model3 <- lme(weight ~ experiment, random = ~ 1 | pop, data=d)
model4 <- lme(weight ~ ecotype, random = ~ 1 | pop, data=d)

# compare models
anova(model1, model2)
# no significant interaction between ecotype and experiment (F=0.039, DF=19, p-value=0.85) (LR = 2.4, DF = 1, p-value=0.12)
anova(model2, model3)
anova(model2, model4)
# both ecotype and experiment are important predictors
summary(model2)
# dune seeds 3.9 mg heavier (F=62.8, DF=20, p-value<0.001) (LR = 30, DF = 1, p-value=<0.001)
# wild collected seeds 2.4 mg heavier (F-13.4, DF=20, p-value=0.0016) (LR = 12, DF = 1, p-value=<0.001)

# calculate mean and 95% CIs for each group
group_means <- d %>% group_by(experiment, ecotype) %>% 
  summarise(mean = mean(weight), sd = sd(weight), n=n()) %>% 
  mutate(UCI = mean + 1.96*sd/sqrt(n), LCI = mean - 1.96*sd/sqrt(n))

# plot graph of seed size (Fig 3)
stripchart(weight~pop, data=d[d$ecotype=="dune",], las=2, vertical=TRUE, pch=16, method="jitter", col=alpha("darkorange1", 0.8), at=c(25,26,1:5,11,12,18), cex=1.2,
           xlim=c(0,33), ylim=c(2,16.5), ylab = "Seed weight (mg)")
nd_pops <- d[d$ecotype!="dune",]$pop %>% unique()
axis(1, at=c(6:10, 13:17,27:29), nd_pops, las=2)
stripchart(weight~pop, data=d[d$ecotype!="dune",], las=2, vertical=TRUE, pch=16, method="jitter", col=alpha("lightgreen", 0.8), at=c(6:10, 13:17,27:29), cex=1.2, add=TRUE)
axis(1, at=c(20,21,31,32), rep(c("Dune", "Non-dune"),2), las=2)
plotCI(x=c(31,32,20,21), y=group_means$mean, ui=group_means$UCI, li=group_means$LCI, err="y", pch=3, col=rep(c(alpha("darkorange1", 0.8), alpha("lightgreen", 0.8)),2), gap=0, lwd=3, sfrac=0.005, add=TRUE)
abline(v=23)
mtext("Natural Populations", 3, 0.5, at=10, adj=0.5, cex=1.2)
mtext("Common Garden", 3, 0.5, at=28.5, adj=0.5, cex=1.2)
legend(23.5,16.5, legend=c("Dune", "Non-dune"), pch=c(15,15), col=c("darkorange1", "lightgreen"), pt.cex=2, bty='n', y.intersp=2)
# save as 7.5 x 5


#########################################
### QTL mapping populations seed size ###
#########################################

d <- read.csv("data/QTL_SS_data.csv", stringsAsFactors = FALSE, na.strings = c("NA"))
table(d$type)

# determine the effect of cytoplasm on seed size for each mapping population (Table S6)
hybrid_cyto_models <- d %>% filter(generation=="F1" | generation=="F2") %>% 
  group_by(map_pop, generation) %>% do(model = lm(weight ~ cytoplasm, data = .))
for (i in 1:8) { print(anova(hybrid_cyto_models[[3]][[i]])) }

# plot seed size for each mapping population (Fig S12)
stripchart(data=d, at=c(1:6, 8:13, 15:20, 22:27), weight~type, las=2, vertical=TRUE, pch=16, method="jitter", col=alpha("black", 0.2), ylab = "Seed weight (mg)", xaxt='n', xlim=c(0,28))
axis(1, at=c(1:6, 8:13, 15:20, 22:27), rep(c("F0D", "F1(D)", "F1(N)", "F2(D)", "F2(N)", "F0N"),4), las=2)
mtext("GSD1", 1, 3.5, at=3.5, adj=0.5, cex=1.2)
mtext("GSD2", 1, 3.5, at=10.5, adj=0.5, cex=1.2)
mtext("MON1", 1, 3.5, at=17.5, adj=0.5, cex=1.2)
mtext("MON2", 1, 3.5, at=24.5, adj=0.5, cex=1.2)
# save as 7.5 x 5






