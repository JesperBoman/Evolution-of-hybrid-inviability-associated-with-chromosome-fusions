#!/bin/bash -l


#Trimming 

module load bioinfo-tools TrimGalore/0.6.1 cutadapt/3.1 FastQC/0.11.9 MultiQC/1.11


while IFS= read -r sample
do

R1=$(grep $sample <(ls fastq_links) | grep "R1")
R2=$(grep $sample <(ls fastq_links) | grep "R2")

trim_galore \
            --quality 30 \
            --paired \
            --illumina \
            --phred33 \
            --stringency 1 \
            -e 0.1 \
            --length 30 \
            --gzip \
            --output_dir trimmed_files \
            --fastqc \
            fastq_links/$R1 \
            fastq_links/$R2


done < "samples.list"

#Get MultiQC stats on trimmed fastq files
multiqc trimmed_files
