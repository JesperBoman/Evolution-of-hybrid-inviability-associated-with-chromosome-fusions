#!/bin/bash -l

date
echo "$1 has BEGUN"

module load bioinfo-tools
module load bwa/0.7.17
module load picard/2.23.4

id=$2

#Ref genome
ref=$3

#Alignment

fqdir="../01_Trimming/trimmed_files"

R1=$(grep -E "^$1" <(ls $fqdir | grep ".fq.gz") | grep "R1")
R2=$(grep -E "^$1" <(ls $fqdir | grep ".fq.gz") | grep "R2")

#Map using bwa mem
bwa mem -t 19 $ref $fqdir/$R1 $fqdir/$R2 > $1.$id.sam

java -Xmx7g -jar $PICARD_ROOT/picard.jar AddOrReplaceReadGroups INPUT=$1.$id.sam OUTPUT=$1.$id.bam SORT_ORDER=coordinate RGID=$1-id RGLB=$1-lib RGPL=ILLUMINA RGPU=$1-01 RGSM=$1

java -Xmx7g -jar $PICARD_ROOT/picard.jar BuildBamIndex INPUT=$1.$id.bam

java -Xmx7g -jar $PICARD_ROOT/picard.jar MarkDuplicates \
  INPUT=$1.$id.bam \
  OUTPUT=$1.$id.dedup.bam \
  METRICS_FILE=$1.$id.dedup.bam.metrics\
  REMOVE_DUPLICATES=true\
  READ_NAME_REGEX=null

java -Xmx7g -jar $PICARD_ROOT/picard.jar BuildBamIndex INPUT=$1.$id.dedup.bam

#Reads were not deduplicated in our run, but duplication level was very low so shouldn't be a problem. I've added it here for completeness.

date
echo "$1 is DONE"

rm $1.$id.sam
rm $1.$id.bam
rm $1.$id.bai
