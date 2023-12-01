#!/bin/bash -l


#Initial first round of variant calling

module load bioinfo-tools GATK/4.2.0.0 bcftools/1.14

sample=$1
id=$2
scaffold=$3
ref=$4
bamdir=$5

mkdir $scaffold.dir

date
echo "$sample.$id.$scaffold, begun"


#Variant calling
gatk --java-options "-Xmx6g" HaplotypeCaller \
  -R $ref \
  --emit-ref-confidence GVCF \
  -I $bamdir/$sample.$id.dedup.bam \
  -O $scaffold.dir/$sample.$id.$scaffold.dedup.g.vcf \
  -L $scaffold


date
echo "$sample.$id.$scaffold is DONE"
