#!/bin/bash -l
#$ -wd .../
#$ -N LGB_olink_soma5k
#$ -l h_vmem=30G
#$ -l h_rt=7:00:00
#$ -o .../logs/olink_to_5ksoma_enet_$TASK_ID.log
#$ -e .../logs/olink_to_5ksoma_enet_$TASK_ID.log
#$ -t 1-20

source .../software/scripts/useuse
use .gcc-13.2.0
use Anaconda3
source activate R_v4.6_env

Rscript .../Codes/revision1/reviewer1/reviewer1_benchmarking_other_meths/c_train_olink_to_soma_elasticnet.R ${SGE_TASK_ID}