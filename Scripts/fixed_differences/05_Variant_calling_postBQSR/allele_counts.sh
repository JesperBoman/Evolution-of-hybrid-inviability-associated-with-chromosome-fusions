#!/bin/bash -l

module load bioinfo-tools vcftools/0.1.16

vcftools --gzvcf postBQSR.SNPs.dedup.g.vcf.gz --keep samples.list.Cat --counts --out postBQSR.SNPs.Cat &

vcftools --gzvcf postBQSR.SNPs.dedup.g.vcf.gz --keep samples.list.Swe --counts --out postBQSR.SNPs.Swe &

vcftools --gzvcf postBQSR.SNPs.dedup.g.vcf.gz --keep samples.list.Out --counts --out postBQSR.SNPs.Out &

wait
