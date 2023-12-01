#!/bin/bash -l

module load bioinfo-tools GATK/4.2.0.0


id=$1
scaffold=$2



s1="37.$id.$scaffold.dedup.g.vcf"
s2="38.$id.$scaffold.dedup.g.vcf"
s3="39.$id.$scaffold.dedup.g.vcf"
s4="40.$id.$scaffold.dedup.g.vcf"
s5="41.$id.$scaffold.dedup.g.vcf"
s6="42.$id.$scaffold.dedup.g.vcf"
s7="43.$id.$scaffold.dedup.g.vcf"
s8="44.$id.$scaffold.dedup.g.vcf"
s9="46.$id.$scaffold.dedup.g.vcf"
s10="47.$id.$scaffold.dedup.g.vcf"
s11="13.$id.$scaffold.dedup.g.vcf"
s12="17.$id.$scaffold.dedup.g.vcf"
s13="29.$id.$scaffold.dedup.g.vcf"
s14="31.$id.$scaffold.dedup.g.vcf"
s15="Swe-sin-101C.$id.$scaffold.dedup.g.vcf"
s16="Swe-sin-102C.$id.$scaffold.dedup.g.vcf"
s17="Swe-sin-1C.$id.$scaffold.dedup.g.vcf"
s18="Swe-sin-2C.$id.$scaffold.dedup.g.vcf"
s19="Swe-sin-31C.$id.$scaffold.dedup.g.vcf"
s20="Swe-sin-32C.$id.$scaffold.dedup.g.vcf"
s21="Swe-sin-61C.$id.$scaffold.dedup.g.vcf"
s22="Swe-sin-62C.$id.$scaffold.dedup.g.vcf"
s23="Swe-sin-91C.$id.$scaffold.dedup.g.vcf"
s24="Swe-sin-92C.$id.$scaffold.dedup.g.vcf"


date
echo "GenomicsDB $id.$scaffold, begun"

gatk --java-options "-Xmx28g" GenomicsDBImport \
  -V $scaffold.dir/$s1 \
  -V $scaffold.dir/$s2 \
  -V $scaffold.dir/$s3 \
  -V $scaffold.dir/$s4 \
  -V $scaffold.dir/$s5 \
  -V $scaffold.dir/$s6 \
  -V $scaffold.dir/$s7 \
  -V $scaffold.dir/$s8 \
  -V $scaffold.dir/$s9 \
  -V $scaffold.dir/$s10 \
  -V $scaffold.dir/$s11 \
  -V $scaffold.dir/$s12 \
  -V $scaffold.dir/$s13 \
  -V $scaffold.dir/$s14 \
  -V $scaffold.dir/$s15 \
  -V $scaffold.dir/$s16 \
  -V $scaffold.dir/$s17 \
  -V $scaffold.dir/$s18 \
  -V $scaffold.dir/$s19 \
  -V $scaffold.dir/$s20 \
  -V $scaffold.dir/$s21 \
  -V $scaffold.dir/$s22 \
  -V $scaffold.dir/$s23 \
  -V $scaffold.dir/$s24 \
  --genomicsdb-workspace-path genomicsdb_preBQSR.$id.$scaffold \
  --reader-threads 1 \
  --tmp-dir $SNIC_TMP \
  -L $scaffold
  
date
echo "GenomicsDB $id.$scaffold, DONE"
