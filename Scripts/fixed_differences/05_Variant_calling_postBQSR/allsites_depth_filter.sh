#!/bin/bash -l

#Load modules
module load bioinfo-tools bcftools/1.14 GATK/4.2.0.0

#Remove genotype information if depth is lower than 5 and larger than 25 reads. Following Näsvall et al. (2023)
bcftools filter -i 'FMT/DP > 5 & FMT/DP < 25' --set-GTs '.' -o postBQSR.allsites.depthFilt.dedup.g.vcf.gz postBQSR.allsites.dedup.g.vcf.gz

gatk --java-options "-Xmx6g" IndexFeatureFile \
  -I postBQSR.allsites.depthFilt.dedup.g.vcf.gz


#Ref:
#Näsvall K, Boman J, Höök L, Vila R, Wiklund C, et al. (2023) Nascent evolution of recombination rate differences as a consequence of chromosomal rearrangements. PLOS Genetics 19(8): e1010717. https://doi.org/10.1371/journal.pgen.1010717
