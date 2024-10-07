#!/bin/bash -l

module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.14
module load GATK/4.2.0.0
module load picard/2.23.4

#Ref genome
ref="Lsin_DToL_CatAllele_for_fixedDiff_noDel.fa"

#BWA index
bwa index -a bwtsw $ref

# BWA is a program that takes different commands as the first argument, here index, used with the bwtsw algorithm

#Samtools and Picard indices
samtools faidx $ref
java -Xmx7g -jar $PICARD_ROOT/picard.jar CreateSequenceDictionary R=$ref O=Lsin_DToL_CatAllele_for_fixedDiff_noDel.dict
