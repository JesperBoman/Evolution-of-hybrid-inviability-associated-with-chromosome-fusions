#!/bin/bash -l

id="DToL"

while IFS= read -r scaffold
do

#Building genomics DB per scaffold
sbatch -J genomicsDB.$id.$scaffold -o slurm/genomicsDB.$id.$scaffold.output -e slurm/genomicsDB.$id.$scaffold.error  -A "Project ID" -t 0-03:00:00 -p core -n 1 genomicsDB.sh $id $scaffold

done < "scaffolds.list"
