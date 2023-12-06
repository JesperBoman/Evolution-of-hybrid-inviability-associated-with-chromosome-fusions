#!/bin/bash -l


## Sign (also known as candidate) region resampling overlap analysis ##



ml bioinfo-tools BEDTools/2.29.2
faidx="/crex/proj/uppstore2017185/b2014034_nobackup/Jesper/F2/reference/Lsin_DToL.fasta.fai"
dir="/crex/proj/uppstore2017185/b2014034_nobackup/Jesper/F2/pool_data_pt2/Rec_rate_analysis"

pop="CAT"
popsmall="Cat"

Group1="DToL_${popsmall}M_rec_rate_${pop}_srADvD" 
Group2="DToL_${popsmall}M_rec_rate_${pop}_ONLYsrADvD"

genome_len=$(awk '{sum+=$2} END{print sum}' $faidx)

#obs_rec_rate=2.43688
obs_rec_rate=3.32214


RANDOM=$(date +%N | cut -b4-9)

resamples=100000

for i in $(seq 1 $resamples)
do


bedtools intersect -wo -a <(awk -v genome_len=$genome_len -v RANDOM=$RANDOM -f region_shuffler.awk $faidx $dir/$Group2.bed) \
	-b $Group1.bed | awk '{if($17 != "NA" && $14 != "."){numerator[$1 "\t" $2 "\t" $3]+=($17*$18); denominator[$1 "\t" $2 "\t" $3]+=$18 }} END{for(region in numerator){print region "\t" numerator[region]/denominator[region] }}' | awk '{sum+=$4; i++} END{print sum/i}' >> rr.resamp.${pop}_rec.rate

done


lnumb=$(cat rr.resamp.${pop}_rec.rate <(echo $obs_rec_rate) | sort -n | grep -n -w "$obs_rec_rate" | cut -f1 -d ":" | head -n1 )

resamples=$(wc -l rr.resamp.${pop}_rec.rate | cut -f1 -d ' ')

prob=$(awk -v lnumb=$lnumb -v resamples=$resamples 'BEGIN{r=resamples+1-lnumb; print r/(resamples)}')

echo "$pop"
echo "Prob is $prob"

#Two-sided p

awk -v prob=$prob 'BEGIN{if(prob<=0.5){print prob*2} else{print (1-prob)*2}}'


#One-sided p

awk -v prob=$prob 'BEGIN{print (1-prob)}'
