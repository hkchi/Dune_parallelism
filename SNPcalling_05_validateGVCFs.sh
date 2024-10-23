#!/bin/bash
#
#SBATCH --mem=100000
#SBATCH --workdir=~/GBS/NegPG
#SBATCH --job-name=NPGvalidateGVCFs
#SBATCH --output=2019_GBS_NegPG_validateGVCFs.log
#SBATCH --error=2019_GBS_NegPG_validateGVCFs.err

while read input
do
~/Programs/jre1.8.0_131/bin/java -Xmx32G -jar ~/Programs/GenomeAnalysisTK.jar \
    -T ValidateVariants -R ~/Reference/Ha412HOv2.0-20181130.fasta \
    -V "$input".GATK.g.vcf
done < samples_gvcf.txt

