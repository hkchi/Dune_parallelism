library(tidyverse)
library(missMDA)
# data sets needed from Dryad: Habitat_data_clean.csv


# load data
d <- read.csv("data/Habitat_data_clean.csv", stringsAsFactors = FALSE, na.strings = c("NA")) %>% 
  mutate(class=paste0(location, type))
d$location <- as.factor(d$location)
d$type <- as.factor(d$type)

# run stats on each variable at each location (Table S2)
summary(lm(mean.percent.cover ~ type, data = d[d$location=="GSD",]))
summary(lm(mean.percent.cover ~ type, data = d[d$location=="MON",]))
summary(lm(percent.grass ~ type, data = d[d$location=="GSD",]))
summary(lm(percent.grass ~ type, data = d[d$location=="MON",]))
summary(lm(total.N ~ type, data = d[d$location=="GSD",]))
summary(lm(total.N ~ type, data = d[d$location=="MON",]))
summary(lm(P ~ type, data = d[d$location=="GSD",]))
summary(lm(P ~ type, data = d[d$location=="MON",]))
summary(lm(K ~ type, data = d[d$location=="GSD",]))
summary(lm(K ~ type, data = d[d$location=="MON",]))
summary(lm(Mg ~ type, data = d[d$location=="GSD",]))
summary(lm(Mg ~ type, data = d[d$location=="MON",]))
summary(lm(Ca ~ type, data = d[d$location=="GSD",]))
summary(lm(Ca ~ type, data = d[d$location=="MON",]))

# impute missing data for PCA
d2 <- d %>% select(-c(location, pop, type, class)) 
nb <- estim_ncpPCA(d2, ncp.min=0, ncp.max=5)
d2_impute2 <- imputePCA(d2, ncp=nb$ncp)

# run PCA
soil.pca <- prcomp(d2_impute2$completeObs, scale=TRUE)

# pull PCA data
d3 <- soil.pca[[2]] %>% as.data.frame()
d4 <- soil.pca[[5]] %>% as.data.frame() %>% 
  mutate(class=d$class,
         col=case_when(class == "GSDD" ~ "firebrick",
                       class == "GSDN" ~ "lightskyblue",
                       class == "MOND" ~ "darkorange1",
                       TRUE ~ "lightgreen"),
         pch=case_when(class == "GSDD" | class == "GSDN" ~ 15, TRUE ~ 16))

# make plot (Fig 1b)
plot(d4$PC1, d4$PC2, col=alpha(d4$col, 0.8), pch=d4$pch, cex=2, xlab = "PC1 (36.2%)", ylab = "PC2 (22.5%)")
for (i in 1:nrow(d4)) {
  arrows(x0=0, y0=0, x1 = 3*d3$PC1[i], y1 = 3*d3$PC2[i])
}
legend(3.4, 2, legend=c("Dune GSD", "Dune MON", "Non-dune GSD", "Non-dune MON"), pch=c(15, 16, 15, 16), col=c("firebrick", "darkorange1", "lightskyblue", "lightgreen"), bty="n", y.intersp=1.5, pt.cex=1.5)
text(labels=c("Grass", "N", "K", "Mg", "Ca", "Cover", "P"), x=c(0.75,1.3,0.75,1.2,1,1.75,1), y=c(0.5,-0.2,2,1.6,-0.5,-0.5,-1))
# save as 8x6

