#!/bin/bash
#
#SBATCH --mem=40000
#SBATCH --workdir=~/GBS/NegPG
#SBATCH --job-name=NPGHaploCaller
#SBATCH --output=2019_GBS_NegPG_HaploCaller.log
#SBATCH --error=2019_GBS_NegPG_HaploCaller.err
#SBATCH --cpus-per-task=10

## ls | grep merged.bam | sed 's/\_merged.bam//' > samples_bams.txt

while read input
do
~/Programs/jre1.8.0_131/bin/java -Xmx16G -jar  ~/Programs/GenomeAnalysisTK.jar \
	-l INFO -R ~/Reference/Ha412HOv2.0-20181130.fasta \
	-log "$input".HaplotypeCaller.log -nct 10 \
	-T HaplotypeCaller -I "$input"_merged.bam --emitRefConfidence GVCF -o "$input".GATK.g.vcf
done < samples_bams.txt
exit
