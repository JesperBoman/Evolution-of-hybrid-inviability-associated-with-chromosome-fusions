#!/bin/bash -l


#

id="DToL"
ref="../reference/Lsin_DToL.fasta"
bamdir="../02_Mapping"


while IFS= read -r sample
do

while IFS= read -r scaffold
do

#Haplotype calling
sbatch -J $sample.$id.$scaffold -o slurm/$sample.$id.$scaffold.output -e slurm/$sample.$id.$scaffold.error  -A "project ID" -t 0-14:00:00 -p core -n 1 hap_call.sh $sample $id $scaffold $ref $bamdir

done < "scaffolds.list"

done < "samples.list"
