#!/bin/bash -l
#$ -wd ...
#$ -N FeatImp_olink_soma
#$ -l h_vmem=60G
#$ -l h_rt=8:00:00
#$ -o .../logs/feature_imp_olink_soma_$TASK_ID.log
#$ -e .../logs/feature_imp_olink_soma_$TASK_ID.log
#$ -t 1

source .../software/scripts/useuse
use .gcc-13.2.0
use Anaconda3
source activate python_env

papermill .../Codes/revision1/feature_importance_rerun/d_feature_importance_OS.ipynb /dev/null
