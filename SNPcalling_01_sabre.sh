#!/bin/bash
#
#SBATCH --mem=20000
#SBATCH --job-name=GBS_sabre

~/Programs/sabre-master/sabre pe \
	-f ~/GBS/KOP1.1.fastq.gz \
	-r ~/GBS/KOP1.2.fastq.gz \
	-b ~/GBS/KOP1_barcodes.txt \
	-u ~/GBS/unknown_KOP1_R1.fastq \
	-w ~/GBS/unknown_KOP1_R2.fastq \
	> ~/GBS/KOP1_sabre_output.txt

~/Programs/sabre-master/sabre pe \
    -f ~/GBS/KOP2.1.fastq.gz \
    -r ~/GBS/KOP2.2.fastq.gz \
    -b ~/GBS/KOP2_barcodes.txt \
    -u ~/GBS/unknown_KOP2_R1.fastq \
    -w ~/GBS/unknown_KOP2_R2.fastq \
    > ~/GBS/KOP2_sabre_output.txt

~/Programs/sabre-master/sabre pe \
    -f ~/GBS/KOP6A.1.fastq.gz \
    -r ~/GBS/KOP6A.2.fastq.gz \
    -b ~/GBS/KOP6A_barcodes.txt \
    -u ~/GBS/unknown_KOP6A_R1.fastq \
    -w ~/GBS/unknown_KOP6A_R2.fastq \
    > ~/GBS/KOP6A_sabre_output.txt

~/Programs/sabre-master/sabre pe \
    -f ~/GBS/KOP6B.1.fastq.gz \
    -r ~/GBS/KOP6B.2.fastq.gz \
    -b ~/GBS/KOP6B_barcodes.txt \
    -u ~/GBS/unknown_KOP6B_R1.fastq \
    -w ~/GBS/unknown_KOP6B_R2.fastq \
    > ~/GBS/KOP6B_sabre_output.txt

~/Programs/sabre-master/sabre pe \
    -f ~/GBS/KOP7A.1.fastq.gz \
    -r ~/GBS/KOP7A.2.fastq.gz \
    -b ~/GBS/KOP7A_barcodes.txt \
    -u ~/GBS/unknown_KOP7A_R1.fastq \
    -w ~/GBS/unknown_KOP7A_R2.fastq \
    > ~/GBS/KOP7A_sabre_output.txt

~/Programs/sabre-master/sabre pe \
    -f ~/GBS/KOP7B.1.fastq.gz \
    -r ~/GBS/KOP7B.2.fastq.gz \
    -b ~/GBS/KOP7B_barcodes.txt \
    -u ~/GBS/unknown_KOP7B_R1.fastq \
    -w ~/GBS/unknown_KOP7B_R2.fastq \
    > ~/GBS/KOP7B_sabre_output.txt
