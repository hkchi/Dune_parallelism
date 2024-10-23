#!/bin/bash
#
#SBATCH --mem=40000
#SBATCH --job-name=NPGallbams
#SBATCH --workdir=~/GBS/NegPG
#SBATCH --output=2019_GBS_NegPG_allbam.log
#SBATCH --error=2019_GBS_NegPG_allbam.err


while read input sample library RG PU day
do
~/Programs/jre1.8.0_131/bin/java -jar ~/Programs/picard-tools-2.5.0/picard.jar FastqToSam \
	FASTQ="$input"_R1.fasta.gz FASTQ2="$input"_R2.fasta.gz \
	OUTPUT="$input"_fastqtosam.bam READ_GROUP_NAME="$RG" \
	SAMPLE_NAME="$sample" LIBRARY_NAME="$library" \
	PLATFORM_UNIT="$PU" PLATFORM=ILLUMINA SEQUENCING_CENTER=UBC \
	RUN_DATE="$day" TMP_DIR=~/temp
~/Programs/jre1.8.0_131/bin/java -jar ~/Programs/picard-tools-2.5.0/picard.jar MarkIlluminaAdapters \
	I="$input"_fastqtosam.bam O="$input"_markedadapters.bam \
	M="$input"_markedadapters_metrics.txt TMP_DIR=~/temp
~/Programs/jre1.8.0_131/bin/java -Xmx8G -jar ~/Programs/picard-tools-2.5.0/picard.jar SamToFastq I="$input"_markedadapters.bam FASTQ="$input".fastq CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true TMP_DIR=~/temp
~/Programs/NextGenMap-0.5.5/bin/ngm-0.5.5/ngm -q "$input".fastq -p -r ~/Reference/Ha412HOv2.0-20181130.fasta -o "$input".sam
~/Programs/jre1.8.0_131/bin/java -Xmx16G -jar ~/Programs/picard-tools-2.5.0/picard.jar MergeBamAlignment R=~/Reference/Ha412HOv2.0-20181130.fasta UNMAPPED_BAM="$input"_fastqtosam.bam ALIGNED_BAM="$input".sam O="$input"_merged.bam CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS TMP_DIR=~/temp
done < ~/Data/GBS/2019_GBS_NegPG_info.txt
exit
