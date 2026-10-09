##### 01 - library #####
library(dplyr)
library(data.table)
library(tidyr)
library(impute)

gc()
rm(list = ls())

setwd("...")


### for distance based imputation like KNN, it is suggested to standardize before imputation


##### 02 - files #####
mesa_olink <- fread("Data/MESA_Olink_3k_log2_nodup_final.tsv.gz", data.table = F)
mesa_soma_5k <- fread("Data/MESA_SOMA_5k_log2_nodup_final.tsv.gz", data.table = F)
sum(is.na(mesa_soma_5k))

length(intersect(colnames(mesa_olink)[-1], colnames(mesa_soma_5k)[-1:-2])) # 1759

## NA only present in the olink dataset

### quality checking olink dataset
row_ok <- rowMeans(is.na(mesa_olink[,-1])) <= 0.1
mesa_olink_short <- mesa_olink[row_ok,]


##### 03 - quality control #####
# isolate olink columns
qc_cols <- setdiff(names(mesa_olink_short), "id")
## keeps the 10% missingness
keep_cols <- qc_cols[colMeans(is.na(mesa_olink_short[qc_cols])) <= 0.1]

# but row(pat wise) this will be different
mesa_olink_qc_na10 <- mesa_olink_short %>%
    select(id, all_of(keep_cols)) %>%
    filter(rowMeans(is.na(across(all_of(keep_cols)))) <= 0.1)
sum(is.na(mesa_olink_qc_na10)) # 15666



## KNN imptation
mesa_olink_qc10_knn <- cbind(mesa_olink_qc_na10 %>% select(id), data.frame(impute.knn(as.matrix(mesa_olink_qc_na10 %>% select(-id)), k = 10, rng.seed = 613)$data))
sum(is.na(mesa_olink_qc10_knn))

write.table(mesa_olink_qc10_knn, gzfile("Data/MESA_Olink_3k_log2_nodup_rm10na_knn_final.tsv.gz"), sep = "\t", row.names = F, quote = F)




# scaling soma data -------------------------------------------------------
### BMI is also scaled in this case
mesa_soma_5k_scaled <- cbind(mesa_soma_5k %>% select(id), data.frame(scale(mesa_soma_5k %>% select(-id))))
mesa_soma_5k_scaled_final <- mesa_soma_5k_scaled %>%
    filter(id %in% mesa_olink_qc10_knn$id)

write.table(mesa_soma_5k_scaled_final, gzfile("Data/MESA_SOMA_5k_log2scaled_nodup_final.tsv.gz"), sep = "\t", row.names = F, quote = F)


### also scaling CHS soma data
chs_soma <- fread("Data/CHS_SOMA_log2_nodup_final.tsv.gz", data.table = F)
chs_soma_scaled <- cbind(chs_soma %>% select(id), data.frame(scale(chs_soma %>% select(-id)))) %>%
    select(id, all_of(colnames(mesa_soma_5k_scaled_final)))

write.table(chs_soma_scaled, gzfile("Data/CHS_SOMA_log2scaled_nodup_final.tsv.gz"), sep = "\t", row.names = F, quote = F)







# impute ukb missing ------------------------------------------------------
ukb_olink <- fread("Data/UKB_Olink_3k_log2_nodup_final.tsv.gz", data.table = F)
ukb_olink <- ukb_olink %>%
    rename(ERVV.1 = `ERVV-1`, HLA.A = `HLA-A`, HLA.DRA = `HLA-DRA`, HLA.E = `HLA-E`) %>%
    select(id, all_of(colnames(mesa_olink_qc10_knn)))


### imputing
# force imputation
ukb_knn_impute <- cbind(ukb_olink %>% select(id), data.frame(impute.knn(as.matrix(ukb_olink %>% select(-id)), k = 10, 
                                                                        rng.seed = 613)$data))

write.table(ukb_knn_impute, gzfile("Data/UKB_Olink_3k_log2_nodup_rm10na_knn_final.tsv.gz"), sep = "\t", row.names = F, quote = F)








