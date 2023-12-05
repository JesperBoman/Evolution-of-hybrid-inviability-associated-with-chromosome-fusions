#!/bin/bash -l

ml bioinfo-tools pixy/1.2.5.beta1

#See pixy documentation for details: https://pixy.readthedocs.io/en/latest/

vcf="../05_Variant_calling_postBQSR/postBQSR.allsites.depthFilt.dedup.g.vcf.gz"

wsize=10000

#Note that pixy use the 1-coordinate system for the start position, so it's not a regular bed-file. 
#I've indicated it here by adding the prefix PF for pixy format.
#See: https://pixy.readthedocs.io/en/latest/companions.html#bed-file
bed="PF_non_sr_10kb.bed"

#When you use a specific .bed-file
pixy --stats pi dxy fst \
--vcf $vcf \
--populations samples.list \
--n_cores 8 \
--fst_type 'hudson' \
--bed_file $bed \
--output_prefix non_sr_10kb_depthFilt

#For genome-wide analyses
pixy --stats pi dxy fst \
--vcf $vcf \
--populations samples.list \
--window_size $wsize \
--n_cores 8 \
--fst_type 'hudson' \
--output_prefix depthFilt
