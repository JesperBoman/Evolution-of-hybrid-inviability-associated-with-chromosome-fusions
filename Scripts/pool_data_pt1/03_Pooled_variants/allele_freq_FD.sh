#!/bin/bash -l

module load bioinfo-tools BEDTools/2.29.2 samtools/1.14 

mapgd_dir="../pool_data_pt1/MAPGD-master/bin"

sample=$1

#The following step takes some time - 37 minutes for F2RP
bedtools intersect -a ../02_Mapping/$sample.dedup.bam -b fixeddiff_noDel.bed > fixed.diff.bams/$sample.fixed.diff.bam 

samtools index fixed.diff.bams/$sample.fixed.diff.bam

samtools view -H  fixed.diff.bams/$sample.fixed.diff.bam > fixed.diff.bams/$sample.fixed.diff.header

#The following step seems to use at most 3-4 cores, took 40 minutes

#MAPGD filtering
map_qual=20
LR=6

#Maximum likelihood variant calling using MAPGD (https://github.com/LynchLab/MAPGD)
samtools mpileup -q $map_qual fixed.diff.bams/$sample.fixed.diff.bam  | $mapgd_dir/mapgd proview -H fixed.diff.bams/$sample.fixed.diff.header -o $sample.fixed.diff
$mapgd_dir/mapgd pool -a $LR -i $sample.fixed.diff.pro -o $sample.fixed.diff_noDel.mq$map_qual.af$LR


fdlist="../../fixed_differences/05_Variant_calling_postBQSR/fixeddiff_noDel.list"
fdlist_anc="../../fixed_differences/05_Variant_calling_postBQSR/fixeddiff_out_noDel_anc.list"

#Pro file format: A/C/G/T

#Since MAPGD only outputs variant sites we need to manually add markers which are monomorphic in a pool

awk -f fixed_in_pool.awk $fdlist $sample.fixed.diff.pro > $sample.tmp.fix

awk 'NR==FNR{sweAl[$1,$2]=$3; catAl[$1,$2]=$4} NR!=FNR && FNR >2{split($8, f, "/"); if(sweAl[$1,$2] == $4){print $1 "\t" $2 "\t" $6 "\t" f[1] "\t" sweAl[$1,$2] "\t" catAl[$1,$2]} else if(catAl[$1,$2] == $4){print $1 "\t" $2 "\t" $6 "\t" 1-f[1] "\t" sweAl[$1,$2] "\t" catAl[$1,$2]}}'  $fdlist $sample.fixed.diff_noDel.mq$map_qual.af$LR.pol | head -n -1 | cat - $sample.tmp.fix > $sample.fixed.diff_noDel.mq$map_qual.af$LR.SweFreq.pol

awk 'NR==FNR{ancAl[$1,$2]=$6} NR!=FNR {if(ancAl[$1,$2] == "SWE_anc"){print $0 "\t" "SWE_anc" "\t" $4} else if(ancAl[$1,$2] == "CAT_anc"){print $0 "\t" "CAT_anc" "\t" 1-$4} else{print $0 "\t" "NA" "\t" "NA"}  }' $fdlist_anc $sample.fixed.diff_noDel.mq$map_qual.af$LR.SweFreq.pol > $sample.fixed.diff_noDel.mq$map_qual.af$LR.SweFreq_AncFreq.pol

rm $sample.tmp.fix
rm $sample.fixed.diff.pro

echo "$sample is DONE"
