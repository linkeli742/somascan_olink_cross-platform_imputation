#!/bin/bash -l
#$ -wd ...
#$ -N bmi_fullpan
#$ -pe smp 4 -binding linear:4 -l h_vmem=16G
#$ -l h_rt=4:00:00
#$ -o .../logs/ukbolink_chssoma_bmi.log
#$ -e .../logs/ukbolink_chssoma_bmi.log
#$ -t 1


source /broad/software/scripts/useuse
use .gcc-13.2.0
use Anaconda3
source activate R_v4.6_env

Rscript .../Codes/revision1/application1_bmi_rerun/b_everything_chssoma.R &
Rscript .../Codes/revision1/application1_bmi_rerun/a_everything_ukbolink.R &
wait