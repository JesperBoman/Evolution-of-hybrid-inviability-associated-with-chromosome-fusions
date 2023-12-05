#!/bin/bash -l

#Change all inferred fixed differences in the genome fasta file into the allele of the Swedish (SWE) population

GENOME="../reference/Lsin_DToL.fasta"


ml bioinfo-tools BEDTools/2.29.2

SNP_list="../05_Variant_calling_postBQSR/fixeddiff_noDel.list"

awk '{{print $1 "\t" $2 "\t" $3}}' $SNP_list > $SNP_list.swe

grep -w --no-group-separator "A" $SNP_list.swe > snps_A

n=1
awk '{print $1 "\t" $2-1 "\t" $2}' snps_A > snpsA.bed


bedtools maskfasta -fi $GENOME \
 -bed snpsA.bed \
 -fo Lsin_DToL${n}.fa \
 -mc A

rm snps_A
rm snpsA.bed

o=$((n+1))
grep -w --no-group-separator "T" $SNP_list.swe > snps_T
awk '{print $1 "\t" $2-1 "\t" $2}' snps_T > snpsT.bed

bedtools maskfasta -fi Lsin_DToL${n}.fa \
 -bed snpsT.bed \
 -fo Lsin_DToL${o}.fa \
 -mc T

rm Lsin_DToL${n}.fa
rm snps_T
rm snpsT.bed

n=$o
o=$((n+1))
grep -w --no-group-separator "C" $SNP_list.swe > snps_C
awk '{print $1 "\t" $2-1 "\t" $2}' snps_C > snpsC.bed

bedtools maskfasta -fi Lsin_DToL${n}.fa \
 -bed snpsC.bed \
 -fo Lsin_DToL${o}.fa \
 -mc C

rm Lsin_DToL${n}.fa
rm snps_C
rm snpsC.bed

n=$o
grep -w --no-group-separator "G" $SNP_list.swe > snps_G
awk '{print $1 "\t" $2-1 "\t" $2}' snps_G > snpsG.bed

bedtools maskfasta -fi Lsin_DToL${n}.fa \
 -bed snpsG.bed \
 -fo Lsin_DToL_SweAllele_for_fixedDiff_noDel.fa \
 -mc G

rm Lsin_DToL${n}.fa
rm snps_G
rm snpsG.bed
