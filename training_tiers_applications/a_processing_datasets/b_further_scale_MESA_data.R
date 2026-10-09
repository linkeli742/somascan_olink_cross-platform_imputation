##### 01 - library #####
library(dplyr)
library(data.table)
library(tidyr)


gc()
rm(list = ls())

setwd("...")


##### 02 - files #####
## mesa's soma raw (no log2 transformed)
mesa_soma_5k_raw <- fread("Data/mesa_olink_processed/mesa_soma_5k_full.tsv.gz", data.table = F)
mesa_soma_5k_log2 <- cbind(mesa_soma_5k_raw[,1:9],log2(mesa_soma_5k_raw[,10:ncol(mesa_soma_5k_raw)]))



##### 03 - process #####
#### chs
chs_raw <- fread(".../CHS_soma5k_bmiwaisttohipratiolipids.txt.gz", data.table = F)[,-1]
chs_proteins <- colnames(chs_raw)[11:ncol(chs_raw)]


##### 03 - NOW WE DONT NEED OVERLAPPING PROTEINS!! #####
soma_5k_ref <- fread("Data/soma_5kprot_ref", data.table = F, header = T)$idno
## note the soma5k is obtained from CHS and olink3k is obtained from ukb

length(intersect(colnames(mesa_soma_5k_log2), soma_5k_ref)) # 4972
mesa_olink_clean_nocontrol_short_filtered_wide <- fread("Data/MESA_Olink_3k_log2_nodup_final.tsv.gz", data.table = F)


##### 04 finalizing each dataset #####
# SOMA --------------------------------------------------------------------

# 4987 columns because +2 columns of id and bmi

mesa_soma_5k_log2_final <- mesa_soma_5k_log2 %>%
    filter(id %in% mesa_olink_clean_nocontrol_short_filtered_wide$id) %>%
    select(id, bmi, any_of(soma_5k_ref)) # 4972+2


write.table(mesa_soma_5k_log2_final, gzfile("Data/MESA_SOMA_5k_log2_nodup_final.tsv.gz"), sep = "\t", row.names = F, quote = F)



# CHS soma ----------------------------------------------------------------
chs_clean <- chs_raw %>%
    rename(id = idno) %>%
    select(id, bmi, all_of(soma_5k_ref))

chs_clean_log2 <- cbind(chs_clean[,1:2],log2(chs_clean[,3:ncol(chs_clean)]))


write.table(chs_clean_log2, gzfile("Data/CHS_SOMA_log2_nodup_final.tsv.gz"), sep = "\t", row.names = F, quote = F)




