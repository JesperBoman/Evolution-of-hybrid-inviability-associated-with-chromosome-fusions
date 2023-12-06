#!/bin/bash -l


dir="../../reference"

cp $dir/Lsin_DToL.fasta $SNIC_TMP
cp $dir/LsinapisSpaM_chr.fasta $SNIC_TMP
cd $SNIC_TMP


TARGET=Lsin_DToL.fasta
QUERY=LsinapisSpaM_chr.fasta

module load bioinfo-tools lastz/1.04.00

lastz_32 $QUERY[multiple] $TARGET M=254 K=4500 L=3000 Y=15000 C=2 T=2 --matchcount=10000 --ambiguous=iupac --format=general:name1,start1,end1,length1,strand1,name2,start2,end2,length2,strand2
