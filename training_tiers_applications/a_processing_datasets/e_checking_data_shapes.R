##### 01 - library #####
library(dplyr)
library(data.table)


gc()
rm(list = ls())


setwd("...")


##### 02 - all files #####
mesa_soma <- fread("Data/final_datasets/MESA_SOMA_5k_log2scaled_nodup_final.tsv.gz", data.table = F) # 5325 x 4974
mesa_olink <- fread("Data/final_datasets/MESA_Olink_3k_log2_nodup_rm10na_knn_final.tsv.gz", data.table = F) # 5325 x 2872


chs_soma <- fread("Data/final_datasets/CHS_SOMA_log2scaled_nodup_final.tsv.gz", data.table = F) # 3171 x 4974

identical(colnames(mesa_soma), colnames(chs_soma)) # TRUE

ukb_olink <- fread("Data/final_datasets/UKB_Olink_3k_log2_nodup_rm10na_knn_final.tsv.gz", data.table = F)

identical(colnames(mesa_olink), colnames(ukb_olink)) # TRUE

#### these dataset are ready to use in training
