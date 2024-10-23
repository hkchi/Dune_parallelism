#!/bin/bash
#
#SBATCH --mem=50000
#SBATCH --workdir=workdir=~/GBS/NegPG
#SBATCH --job-name=NPGfiltVCF
#SBATCH --output=2019_GBS_NegPG_filtervcf.log
#SBATCH --error=2019_GBS_NegPG_filtervcf.err

~/Programs/jre1.8.0_131/bin/java -Xmx18G -jar ~/Programs/GenomeAnalysisTK.jar \
	-T SelectVariants -R ~/Reference/Ha412HOv2.0-20181130.fasta \
	-V NegPG_GBS_NGM_GATK_total.vcf -selectType SNP -o NegPG_GBS_NGM_GATK_raw_snps.vcf

~/Programs/jre1.8.0_131/bin/java -Xmx18G -jar ~/Programs/GenomeAnalysisTK.jar \
	-T SelectVariants -R ~/Reference/Ha412HOv2.0-20181130.fasta \
	-V NegPG_GBS_NGM_GATK_total.vcf -selectType INDEL -o NegPG_GBS_NGM_GATK_raw_indels.vcf

~/Programs/jre1.8.0_131/bin/java -Xmx18G -jar ~/Programs/GenomeAnalysisTK.jar \
    -T VariantFiltration -R ~/Reference/Ha412HOv2.0-20181130.fasta \
    -V NegPG_GBS_NGM_GATK__snps.vcf \
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5" \
    --filterName "GATK_snp_filter" \
    -o NegPG_GBS_NGM_GATK_snps2.vcf

~/Programs/jre1.8.0_131/bin/java -Xmx18G -jar ~/Programs/GenomeAnalysisTK.jar \
    -T SelectVariants -R ~/Reference/Ha412HOv2.0-20181130.fasta \
    -V NegPG_GBS_NGM_GATK_snps2.vcf -env -ef -o NegPG_GBS_NGM_GATK_good_snps.vcf

~/Programs/jre1.8.0_131/bin/java -Xmx18G -jar ~/Programs/GenomeAnalysisTK.jar \
    -T VariantFiltration -R ~/Reference/Ha412HOv2.0-20181130.fasta \
    -V NegPG_GBS_NGM_GATK_raw_indels.vcf \
    --filterExpression "QD < 2.0 || FS > 200.0" \
    --filterName "GATK_indel_filter"  \
    -o NegPG_GBS_NGM_GATK_indels2.vcf

~/Programs/jre1.8.0_131/bin/java -Xmx18G -jar ~/Programs/GenomeAnalysisTK.jar \
    -T SelectVariants -R ~/Reference/Ha412HOv2.0-20181130.fasta \
    -V NegPG_GBS_NGM_GATK_indels2.vcf -env -ef -o NegPG_GBS_NGM_GATK_good_indels.vcf
