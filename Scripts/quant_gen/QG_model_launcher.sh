#!/bin/bash -l

ml R_packages/4.1.1

#If you want to run the models on e.g. a cluster. However, with 10^5 iterations both model are quite fast on a regular desktop.

Rscript --vanilla devtime_survival_model.R
#Rscript --vanilla additive_survival_models.R
