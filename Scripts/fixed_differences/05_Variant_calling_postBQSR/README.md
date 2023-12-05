This folder contains the relevant scripts for variant calling after base-quality score recalibration (BQSR).

It also contains scripts for inferring fixed differences.


## Get fixed differences ##
awk -f fix_diff_2pops.awk postBQSR.SNPs.Swe.frq.count postBQSR.SNPs.Cat.frq.count > fixeddiff.list

#Removed fixed indels (since they can be problematic for PoolSeq data)
grep -v '*' fixeddiff.list > fixeddiff_noDel.list

#Transform into .bed-format
awk '{print $1 "\t" $2-1 "\t" $2 "\t" $3 "\t" $4}' fixeddiff_noDel.list > fixeddiff_noDel.bed

#Polarixation of fixed differences
awk -f fix_diff_outgroups.awk postBQSR.SNPs.Out.frq.count fixeddiff.list > fixeddiff_out.list

grep -v '*' fixeddiff_out.list > fixeddiff_out_noDel.list

awk '$5 != "NA"{if($5 == $3){print $0 "\t" "SWE_anc"} else if($5 == $4){print $0 "\t" "CAT_anc"} else{print $0 "\t" "Other"}}' fixeddiff_out_noDel.list > fixeddiff_out_noDel_anc.list
