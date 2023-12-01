#!/bin/bash -l


ref="../reference/Lsin_DToL.fasta"

dir="../fixed_differences"

id="DToL"

mkdir recal_tables

while IFS= read -r sample
do

sbatch -J $sample.$id -o slurm/$sample.$id.ApplyBQSR.output -e slurm/$sample.$id.ApplyBQSR.error  -A "Project ID" -t 0-12:00:00 -p core -n 1 ApplyBQSR.sh $sample $ref $dir $id

done < "samples.list"
