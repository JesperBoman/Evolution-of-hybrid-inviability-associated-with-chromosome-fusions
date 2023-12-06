#!/bin/bash -l

id="Cat"

while IFS= read -r sample
do

sbatch -J $sample -o slurm/$sample.$id.output -e slurm/$sample.$id.error  -A "Project ID"  -t 0-03:00:00 -p core -n 3 allele_freq_FD.sh $sample.$id

done < "samples.list"
