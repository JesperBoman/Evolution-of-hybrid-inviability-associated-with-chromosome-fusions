#!/bin/bash -l



while IFS= read -r sample
do

sbatch -J $sample -o slurm/$sample.output -e slurm/$sample.error  -A "Project ID" -t 0-12:00:00 -p node -n 16 trimmer.sh $sample

done < "samples.list"
