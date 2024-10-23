# GEA

Analyses of habitat characteristics at each location and PCA of all environmental variables:
`Habitat_stats_and_PCA.R`

Generate input files from VCF
```
zcat MON_GEA.vcf.gz | perl vcf2baypass.pl sample_info.txt MON_GEA.baypass
zcat MON_GEA.vcf.gz | perl vcf2baypass.pl sample_info.txt MON_GEA.baypass_ex ex
# MON_GEA.baypass.pop, MON_GEA.baypass.txt, MON_GEA.baypass.snp
# MON_GEA.baypass_ex.pop, MON_GEA.baypass_ex.txt, MON_GEA.baypass_ex.snp
```

Organize environmental data input for BayPass
```
cat Population_data.cleaned.txt | perl env2baypass.pl MON_GEA.baypass.pop env.baypass
cat Population_data.cleaned.txt | perl env2baypass.pl MON_GEA.baypass_ex.pop env.baypass_ex
# env.baypass.cov, env.baypass.txt
# env.baypass_ex.cov, env.baypass_ex.txt
```

Randomly choose 1k SNPs with LD<0.1 for covariance matrix
```
plink --vcf MON_GEA.vcf.gz --allow-extra-chr --double-id --recode --out MON_GEA --threads 4
awk '{print $1"\t"$1"_"$4"\t"$3"\t"$4}' MON_GEA.map > tmp; mv tmp MON_GEA.map # Give the SNPs names
plink --file MON_GEA --allow-extra-chr --double-id --indep-pairwise 10000 10 0.1 --out MON_GEA --threads 4 # Prune the SNPs based on LD
shuf MON_GEA.prune.in | head -1000 | sed 's/_/\t/g' | sort -k1,1 -k2,2n > MON_GEA.prune.in.random1k.txt # Randomly choose 1k SNPs and subset
bcftools view -R MON_GEA.prune.in.random1k.txt MON_GEA.vcf.gz -O z -o MON_GEA.random1k.vcf.gz
zcat MON_GEA.random1k.vcf.gz | perl vcf2baypass.pl sample_info.txt MON_GEA.random1k.baypass
zcat MON_GEA.random1k.vcf.gz | perl vcf2baypass.pl sample_info.txt MON_GEA.random1k.baypass_ex ex
# MON_GEA.random1k.baypass.pop, MON_GEA.random1k.baypass.txt, MON_GEA.random1k.baypass.snp
# MON_GEA.random1k.baypass_ex.pop, MON_GEA.random1k.baypass_ex.txt, MON_GEA.random1k.baypass_ex.snp
```

Run BayPass under the core model mode to generate covariance matrix
```
i_baypass -npop 18 -gfile MON_GEA.random1k.baypass.txt -outprefix all -nthreads 4
i_baypass -npop 15 -gfile MON_GEA.random1k.baypass_ex.txt -outprefix all_ex -nthreads 4
# all_mat_omega.out
# all_ex_mat_omega.out
```

Run BayPass under the standard covariate model using importance sampling (IS) estimator 
```
i_baypass -npop 18 -gfile MON_GEA.baypass.txt -efile env.baypass.txt -scalecov -omegafile all_mat_omega.out -outprefix all -nthreads 4
i_baypass -npop 15 -gfile MON_GEA.baypass_ex.txt -efile env.baypass_ex.txt -scalecov -omegafile all_ex_mat_omega.out -outprefix all_ex -nthreads 4
```

Produce POD samples with 1,000 SNPs and run BayPass
```
Rscript pod.R
i_baypass -npop 18 -gfile G.MON_GEA.baypass.pod -efile env.baypass.txt -scalecov -omegafile all_mat_omega.out -outprefix all_pod -nthreads 4
i_baypass -npop 15 -gfile G.MON_GEA.baypass_ex.pod -efile env.baypass_ex.txt -scalecov -omegafile all_ex_mat_omega.out -outprefix all_ex_pod -nthreads 4
```

Plot BayPass results
```
Rscript BayPass_plot.R
```
