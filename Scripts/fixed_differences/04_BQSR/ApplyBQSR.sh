#!/bin/bash -l
 
 module load bioinfo-tools GATK/4.2.0.0
 
 sample=$1
 ref=$2
 dir=$3
 id=$4

gatk --java-options "-Xmx6g" ApplyBQSR \
   -R $ref \
   -I  $dir/02_Mapping/$sample.$id.dedup.bam \
   --bqsr-recal-file recal_tables/$sample.$id.recal_data.table \
   -O $sample.$id.dedup.recal.bam
