#!/bin/bash -l

# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ==============================================================================================================================
# A bootstrapping method to calculate the empirical p-value of base pair overlap between two .bed-files
# ==============================================================================================================================
# Jesper Boman                      1 oct 2024
# ==============================================================================================================================
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#This ver.2. handles the case of underenrichment better. 
#In particular, when there are many 0 overlaps among the shuffled samples and the observed is a 0. It takes the average position.
#Also now outputs two-sided p-values by default


mkdir Monte_Carlo_overlap


ml bioinfo-tools BEDTools/2.29.2
faidx="../reference/Lsin_DToL.fasta.fai" #A fasta index file
dir="."


#In this example we are testing the overlap between chromosomes that fissioned in the CAT lineage and the large-effect QTL seq loci
Group1="fi_cat_simple_rearrangements_fullChrs" #Reference
Group2="QTLseq_sim_CI95"


group1_len=$(awk '{sum+=($3-$2)} END{print sum}' $Group1.bed)
group2_len=$(awk '{sum+=($3-$2)} END{print sum}' $Group2.bed)
genome_len=$(awk '{sum+=$2} END{print sum}' $faidx)


obs_overlap=$(bedtools intersect -wo -a <(cut -f1-3 $dir/$Group2.bed) -b <(cut -f1-3 $dir/$Group1.bed) | awk '{sum+=$7}END{print sum}' | awk '{if($1==""){print 0} else{print $1}}')

RANDOM=$(date +%N | cut -b4-9)



resamples=1000

for i in $(seq 1 $resamples)
do


bedtools intersect -wo -a <(awk -v genome_len=$genome_len -v RANDOM=$RANDOM -f region_shuffler.awk $faidx $dir/$Group2.bed ) \
	-b <(cut -f1-3 $dir/$Group1.bed) | cut -f7 | awk '{sum+=$1}END{print sum}' >> $Group1.v.$Group2.resample_overlap

done


#This is to correct for handling when there are many ties (e.g. many 0's)
lnumb_upper=$(cat <( awk '{if($1==""){print 0} else{print $1}}' $Group1.v.$Group2.resample_overlap) <(echo $obs_overlap) | sort -n | grep -n -w "$obs_overlap" | cut -f1 -d ":" | head -n1 )
lnumb_lower=$(cat <( awk '{if($1==""){print 0} else{print $1}}' $Group1.v.$Group2.resample_overlap) <(echo $obs_overlap) | sort -n | grep -n -w "$obs_overlap" | cut -f1 -d ":" | tail -n1 )

#In case of ties, just take the average position
lnumb=$(awk -v lnumb_upper=$lnumb_upper -v lnumb_lower=$lnumb_lower 'BEGIN{print (lnumb_lower+lnumb_upper)/2 }')

pval=$(awk -v lnumb=$lnumb -v resamples=$resamples 'BEGIN{r=resamples+1-lnumb; print r/(resamples)}')


# Group2 as a fraction of genome, Group1 as a fraction of genome, Group1 and Group2 overlap, odds ratio of Group2 overlapping Group1, empirical p-value 
awk -v Group1=$Group1 -v Group2=$Group2 -v group1_len=$group1_len -v genome_len=$genome_len -v group2_len=$group2_len -v obs_overlap=$obs_overlap -v pval=$pval 'BEGIN{if(pval<0.5){pval_two_sided=pval*2} else{pval_two_sided=(1-pval)*2} print group2_len/genome_len "\t" group1_len/genome_len "\t"  obs_overlap  "\t" (obs_overlap/group2_len)/(group1_len/genome_len) "\t" pval_two_sided "\t" Group1 ".v." Group2 }' >> Monte_Carlo_overlap2/stats_$Group1.v.$Group2.overlap


rm $Group1.v.$Group2.resample_overlap
