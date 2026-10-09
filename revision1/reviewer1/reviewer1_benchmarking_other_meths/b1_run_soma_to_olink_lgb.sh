#!/bin/bash -l
#$ -wd .../
#$ -N LGB_soma5k_olink
#$ -l h_vmem=30G
#$ -l h_rt=8:00:00
#$ -o .../logs/5ksoma_to_olink_lgb_$TASK_ID.log
#$ -e .../logs/5ksoma_to_olink_lgb_$TASK_ID.log
#$ -t 1-10

source .../software/scripts/useuse
use .gcc-7.3.0
use Anaconda3
source activate python_env

i=${SGE_TASK_ID}

papermill .../Codes/revision1/reviewer1/reviewer1_benchmarking_other_meths/b_train_log2scaled_soma_to_olink_lgb.ipynb /dev/null -p batch_number ${i}
