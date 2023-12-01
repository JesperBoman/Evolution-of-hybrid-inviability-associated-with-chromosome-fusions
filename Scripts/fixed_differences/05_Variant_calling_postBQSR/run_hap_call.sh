#!/bin/bash -l


#
id="DToL"
ref="../reference/Lsin_DToL.fasta"
bamdir="../04_BQSR"


while IFS= read -r sample
do

while IFS= read -r scaffold
do

#Haplotype calling
sbatch -J $sample.$id.$scaffold.hap_call -o slurm/$sample.$id.$scaffold.hap_call.output -e slurm/$sample.$id.$scaffold.hap_call.error  -A "Project ID" -t 0-16:00:00 -p core -n 1 hap_call.sh $sample $id $scaffold $ref $bamdir

done < "scaffolds.list"

done < "samples.list"
