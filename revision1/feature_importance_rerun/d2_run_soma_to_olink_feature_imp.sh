#!/bin/bash -l
#$ -wd ...
#$ -N FeatImp_soma_olink
#$ -l h_vmem=30G
#$ -l h_rt=3:00:00
#$ -o .../logs/feature_imp_soma_olink_$TASK_ID.log
#$ -e .../logs/feature_imp_soma_olink_$TASK_ID.log
#$ -t 1


source .../software/scripts/useuse
use .gcc-13.2.0
use Anaconda3
source activate python_env

papermill .../Codes/revision1/feature_importance_rerun/d_feature_importance_SO.ipynb /dev/null
