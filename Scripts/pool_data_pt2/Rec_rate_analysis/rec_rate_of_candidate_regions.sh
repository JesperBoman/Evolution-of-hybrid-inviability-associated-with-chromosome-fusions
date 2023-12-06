#!/bin/bash -l


#Reformat LASTZ output into bed.file
awk 'NR>1{print $1 "\t" $2-1 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7-1 "\t" $8 "\t" $9 "\t" $10}' LASTZ_SweM_chr_to_DToL.output > SweM_DToL.bed


#Remove singleton chromosome combination alignments
awk '{a[NR]=$0; b[NR]=$1; c[NR]=$6;  d[$1 "\t" $6]++} END{for(i in a){if( d[b[i] "\t" c[i]]>1 ){print a[i]}}}' SweM_DToL.bed | sort -k 1,1 -k 2,2n > SweM_DToL_noSingles.bed 

wc -l DToL_SweM.bed 
#16042

wc -l DToL_SweM_noSingles.bed
#15994 DToL_SweM_noSingles.bed

#This table comes from Näsvall et al. (2023, Plos Genetics)
sed 's/"//g' rec_rate_221027_update.table > rec_rate_221027_update.table_mod


#Obtaining SWE recombination rates, some modifications needed to get CAT recombination rates
awk 'BEGIN{prev=1} $2=="ls_swe"{if(prevChr == $9){print $9 "\t" prev-1 "\t" $16 "\t" (10^6 *(($18-prevGenPos)/($16-prev+1))); prev=$16; prevChr=$9; prevGenPos=$18 } else{print $9 "\t" 0 "\t" $16 "\t" "NA"; prev=$16; prevChr=$9; prevGenPos=$18 }}' <(sort -k 19,19n -k 16,16n -k 3,3n rec_rate_221027_update.table_mod)  > rec_rate_221027_update.SWE.bed


ml bioinfo-tools BEDTools/2.29.2

sr="../sign_regions_ADvDEAD_dedup_wave0.05.bed"
aln="../../genome_alignment/lastz/SweM_DToL_noSingles.bed"


bedtools intersect -a $aln -b rec_rate_221027_update.SWE.bed -wao | awk '{OFS="\t"; print $6, $7, $8, $9, $10, $1, $2, $3, $4, $5, $11, $12, $13, $14, $15}' > DToL_SweM_rec_rate_SWE.bed

bedtools intersect -b $sr -a DToL_SweM_rec_rate_SWE.bed -wao > DToL_SweM_rec_rate_SWE_srADvD.bed
awk '{if($16 != "." && $14 != "NA" && $11 != "."){numerator[$16 "\t" $17 "\t" $18 "\t" $19]+=($14*$15); denominator[$16 "\t" $17 "\t" $18  "\t" $19]+=$15 }} END{for(region in numerator){print region "\t" numerator[region]/denominator[region] }}' DToL_SweM_rec_rate_SWE_srADvD.bed > DToL_SweM_rec_rate_SWE_ONLYsrADvD.bed 

