#u!/usr/bin/awk -f

#Input: Allele counts from vcftools, for example:
#CHROM	POS	N_ALLELES	N_CHR	{ALLELE:COUNT}
#Chr_10	681	2	0	G:0	C:0

#Usage: awk -f fix_diff_2pops.awk postBQSR.SNPs.Swe.frq.count postBQSR.SNPs.Cat.frq.count > fixeddiff.list

#Only alleles with complete information for both populations are used here. 
#This is a conservative measure to increase the probability that an observed "fixed difference" in the sample is a true fixed difference.

#So far only adapted to biallelic sites (in the entire vcf: SWE+CAT+Outgroup samples)

NR==FNR && $4 == 20{
split($5, ref, ":");
split($6, alt, ":");
if(ref[2] == 20){pop1[$1 "\t" $2]=ref[1]};
if(alt[2] == 20){pop1[$1 "\t" $2]=alt[1]};

}
NR!=FNR && $4 == 20{
split($5, ref, ":");
split($6, alt, ":");
if(ref[2] == 20){pop2[$1 "\t" $2]=ref[1]};
if(alt[2] == 20){pop2[$1 "\t" $2]=alt[1]};
}

END{
for(locus in pop1){
	if(locus in pop2 && pop1[locus] != pop2[locus] ){print locus "\t" pop1[locus] "\t" pop2[locus]
	}
}
}
