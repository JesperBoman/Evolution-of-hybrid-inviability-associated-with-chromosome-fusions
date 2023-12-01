#!/bin/bash -l

#Load modules
module load bioinfo-tools GATK/4.2.0.0 bcftools/1.14

#Concatenate scaffold-separate allsites vcfs
bcftools concat -o preBQSR.allsites.dedup.g.vcf.gz Chr*.dir/*.allsites.vcf.gz

#Extract SNPs only

ref="../reference/Lsin_DToL.fasta"

gatk --java-options "-Xmx6g" IndexFeatureFile \
  -I preBQSR.allsites.dedup.g.vcf.gz

gatk --java-options "-Xmx6g" SelectVariants \
  -R $ref \
  -V preBQSR.allsites.dedup.g.vcf.gz\
  --select-type-to-include SNP \
  -O preBQSR.SNPs.dedup.g.vcf.gz

#Extract quality metrics
bcftools query preBQSR.SNPs.dedup.g.vcf.gz -f'%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\n' > preBQSR_FS.SOR.MQRS.RPRS.QD.MQ.DP.txt

awk '{for (i=1; i<=NF; i++){if($i == ".")next} print $0}'  preBQSR_FS.SOR.MQRS.RPRS.QD.MQ.DP.txt > preBQSR_FS.SOR.MQRS.RPRS.QD.MQ.DP.filt.txt

#These can be visualized and explored for example using R

#Apply filtering

bcftools filter -i 'FS<60.0 && SOR<3 && MQ>40 && MQRankSum>-12.5 && QD>2 && ReadPosRankSum>-8 && INFO/DP<746.2598' -O z -o preBQSR.HQSNPs.dedup.g.vcf.gz preBQSR.SNPs.dedup.g.vcf.gz


#Indexing HQSNPs vcf-file

gatk --java-options "-Xmx6g" IndexFeatureFile \
  -I preBQSR.HQSNPs.dedup.g.vcf.gz

