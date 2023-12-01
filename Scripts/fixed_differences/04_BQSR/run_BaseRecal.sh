#!/bin/bash -l


ref="../reference/Lsin_DToL.fasta"

dir="../fixed_differences"

id="DToL"

mkdir recal_tables

while IFS= read -r sample
do

sbatch -J $sample.$id -o slurm/$sample.$id.BaseRecal.output -e slurm/$sample.$id.BaseRecal.error  -A "Project ID" -t 0-08:00:00 -p core -n 1 BaseRecal.sh $sample $ref $dir $id

done < "samples.list"
