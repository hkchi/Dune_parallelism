#!/bin/bash
#
#SBATCH --mem=20000
#SBATCH --job-name=NPGmergeBams
#SBATCH --workdir=~/GBS/NegPG
#SBATCH --output=2019_GBS_NegPG_mergeBams.log
#SBATCH --error=2019_GBS_NegPG_mergeBams.err

## ls | grep merged.bam | sed 's/\_KOP.*//' | uniq -d > samples_multilane.txt

while read input
do
~/Programs/samtools-1.3.1/samtools merge \
	"$input"_all_merged.bam "$input"_*_merged.bam
~/Programs/jre1.8.0_131/bin/java -jar ~/Programs/picard-tools-2.5.0/picard.jar \
    BuildBamIndex INPUT="$input"_all_merged.bam
rm "$input"_K*_merged.ba*
done < samples_multilane.txt
exit

