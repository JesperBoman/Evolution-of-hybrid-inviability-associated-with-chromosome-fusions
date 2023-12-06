#!/bin/bash -l


while IFS= read -r sample
do

id="Cat"

#Map using BWA
sbatch -J $sample -o slurm/$sample.$id.output -e slurm/$sample.$id.error  -A "Project ID"  -t 1-12:00:00 -p node -n 20 map_bwa.sh $sample $id


done < "samples.list"
