#!/bin/bash -l
#$ -wd .../
#$ -N LGB_soma5k_olink
#$ -l h_vmem=30G
#$ -l h_rt=3:00:00
#$ -o .../logs/5ksoma_to_olink_enet_$TASK_ID.log
#$ -e .../logs/5ksoma_to_olink_enet_$TASK_ID.log
#$ -t 1-10

source .../software/scripts/useuse
use .gcc-13.2.0
use Anaconda3
source activate R_v4.6_env

Rscript .../Codes/revision1/reviewer1/reviewer1_benchmarking_other_meths/c_train_soma_to_olink_elasticnet.R ${SGE_TASK_ID}