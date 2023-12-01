#!/bin/bash -l
 
 module load bioinfo-tools GATK/4.2.0.0
 
 sample=$1
 ref=$2
 dir=$3
 id=$4
 
 gatk --java-options "-Xmx6g" BaseRecalibrator \
   -I $dir/02_Mapping/$sample.$id.dedup.bam \
   -R $ref \
   --known-sites $dir/03_Variant_calling/preBQSR.HQSNPs.dedup.g.vcf.gz 	\
   -O recal_tables/$sample.$id.recal_data.table
