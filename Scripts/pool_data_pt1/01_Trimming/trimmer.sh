#Trimming 

module load bioinfo-tools TrimGalore/0.6.1 cutadapt/3.1 FastQC/0.11.9 MultiQC/1.11


R1=$(grep $1 <(ls fastq_links) | grep "R1")
R2=$(grep $1 <(ls fastq_links) | grep "R2")

trim_galore \
            --quality 30 \
            --paired \
            --illumina \
            --phred33 \
            --stringency 1 \
            -e 0.1 \
            --length 30 \
	    --three_prime_clip_R1 7 \
	    --three_prime_clip_R2 7 \
            --gzip \
	    --cores 4 \
            --output_dir trimmed_files \
            --fastqc \
            fastq_links/$R1 \
            fastq_links/$R2
