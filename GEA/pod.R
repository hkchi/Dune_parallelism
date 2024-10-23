## Produce POD samples
## Kaichi Huang 2021 Jun
## modified from Kaichi Huang 2019 Sep

#install.packages("mvtnorm")
source("~/baypass/2.1/utils/baypass_utils.R")

# Read in the covariance matrix
omega=as.matrix(read.table("all_mat_omega.out"))
# Get estimates (post. mean) of both the a_pi and b_pi parameters of the Pi Beta distribution
pi.beta.coef=read.table("all_summary_beta_params.out",h=T)$Mean
# Original data to obtain total allele count
gt.data<-geno2YN("MON_GEA.baypass.txt")
# Create the POD (Files *.<suffix> saved)
pod.data <-simulate.baypass(omega.mat=omega, nsnp=1000, sample.size=gt.data$NN, 
                            beta.pi=pi.beta.coef, pi.maf=0, suffix="MON_GEA.baypass.pod")

# Read in the covariance matrix
omega=as.matrix(read.table("all_ex_mat_omega.out"))
# Get estimates (post. mean) of both the a_pi and b_pi parameters of the Pi Beta distribution
pi.beta.coef=read.table("all_ex_summary_beta_params.out",h=T)$Mean
# Original data to obtain total allele count
gt.data<-geno2YN("MON_GEA.baypass_ex.txt")
# Create the POD (Files *.<suffix> saved)
pod.data <-simulate.baypass(omega.mat=omega, nsnp=1000, sample.size=gt.data$NN, 
                            beta.pi=pi.beta.coef, pi.maf=0, suffix="MON_GEA.baypass_ex.pod")
