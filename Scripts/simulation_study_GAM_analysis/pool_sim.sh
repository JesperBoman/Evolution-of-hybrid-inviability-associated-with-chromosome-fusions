#!/bin/bash -l


ml R_packages/4.3.1

for pf in 0.5 0.625 0.75 0.875 1.0;
do

for loci in 100 500 1000;
do

for rr in $(seq 0 3);
do

#args[1] = constant pool fraction
#args[2] = number of loci
#args[3] = rec rate

Rscript GAMsim.R $pf $loci $rr >> sim_res.csv

done

done

done
