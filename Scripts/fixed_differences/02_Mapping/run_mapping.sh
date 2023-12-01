#!/bin/bash -l


while IFS= read -r sample
do

id="DToL" #If you want to map reads to several reference genomes, the id suffix keeps track on which reference you've mapped to
ref="PATH/reference/Lsin_DToL.fasta" #Add path to reference genome

#Map using BWA
sbatch -J $sample -o slurm/$sample.$id.output -e slurm/$sample.$id.error  -A "project ID"  -t 03:00:00 -p node -n 20 map_bwa.sh $sample $id $ref

done < "samples.list"
