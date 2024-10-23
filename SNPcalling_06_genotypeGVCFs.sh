#!/bin/bash
#
#SBATCH --mem=100000
#SBATCH --workdir=~/GBS/NegPG
#SBATCH --job-name=NPGgenoGVCFs
#SBATCH --output=2019_GBS_NegPG_genotypeGVCFs.log
#SBATCH --error=2019_GBS_NegPG_genpoypeGVCFs.err

#ls | grep 'GATK.g.vcf' | sed 's/\.GATK\.g\.vcf.*//' | uniq > samples_gvcf.txt

tmp=""
while read prefix
do
        tmp="$tmp --variant /work/klo20/GBS/NegPG/$prefix.GATK.g.vcf"
done < samples_gvcf.txt
echo $tmp
~/Programs/jre1.8.0_131/bin/java -Xmx64G -jar ~/Programs/GenomeAnalysisTK.jar \
	-nt 4 -l INFO -R ~/Reference/Ha412HOv2.0-20181130.fasta \
    --disable_auto_index_creation_and_locking_when_reading_rods \
	-log GenotypeGVCFs.log -T GenotypeGVCFs $tmp -o NegPG_GBS_NGM_GATK_total.vcf
