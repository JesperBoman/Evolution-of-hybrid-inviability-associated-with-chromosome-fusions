#!/bin/bash -l

id="DToL"
ref="../reference/Lsin_DToL.fasta"

while IFS= read -r scaffold
do

#Genotype gVCFs
sbatch -J genotypeGVCFs.$id.$scaffold -o slurm/genotypeGVCFs.$id.$scaffold.output -e slurm/genotypeGVCFs.$id.$scaffold.error -A "Project ID"  -t 0-18:00:00 -p core -n 1 genotypeGVCFs.sh $id $scaffold $ref

done < "scaffolds.list
