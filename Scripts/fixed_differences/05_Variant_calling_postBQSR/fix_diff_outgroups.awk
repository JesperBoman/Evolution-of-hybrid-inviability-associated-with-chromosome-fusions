#u!/usr/bin/awk -f
#Usage: awk -f fix_diff_outgroups.awk postBQSR.SNPs.Out.frq.count fixeddiff.list > fixeddiff_out.list

#Five out of eight gene copies must have a certain allele. That gene allele is considered ancestral.
#So far only adapted to biallelic sites (in the entire vcf: SWE+CAT+Outgroup samples)

NR==FNR && $4 > 4{
split($5, ref, ":");
split($6, alt, ":");
if(ref[2] > 4){outgrp[$1 "\t" $2]=ref[1]};
if(alt[2] > 4){outgrp[$1 "\t" $2]=alt[1]};

}
NR!=FNR {
locus=$1 "\t" $2
if(locus in outgrp){print $0 "\t" outgrp[locus]}
else{ print $0 "\t" "NA" }

}
