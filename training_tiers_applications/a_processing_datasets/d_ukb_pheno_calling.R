##### 01 - library #####
library(dplyr)
library(data.table)


gc()
rm(list = ls())

setwd("...")

##### 02 - files #####
ukb_id <- fread("Data/final_datasets/UKB_Olink_3k_log2_nodup_rm10na_knn_final.tsv.gz", data.table = F) %>%
    select(id)
ukb_pheno <- fread(".../NEW_3Kproteins_adiposityindex.txt.gz", drop = 1, data.table = F)
ukb_bmi3k <- ukb_pheno %>%
    select(id = eid)

## make covariates
covar <- fread(".../UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.REAL.txt", data.table = F)
covs$IDs_toRemove_SampleQCv3 = ifelse((!(covs$Submitted_Gender == covs$Inferred_Gender) |  covs$Non_Consented== 1),1,0)
covs = covs[which(covs$IDs_toRemove_SampleQCv3==0),]

covs_short <- covar %>%
    select(id, Sex, enroll_age.y, in_white_British_ancestry_subset, paste("PC",1:10,sep=""), ever_smoked, Prev_Diabetes_Type_2, BMI) %>%
    left_join(ukb_bmi3k, by = "id")

sum(is.na(covs_short))

covs_short[which(is.na(covs_short)),]

write.table(covs_short, "Data/mesa_olink_processed/UKB_CHIP_covs.tsv", sep = "\t", row.names = F)




##### 03 - statin #####
ukb_drugs <- fread(".../lipid_drug_2visits_withlipids.txt", data.table = F) %>%
    rename(id = eid) %>%
    filter(id %in% ukb_id$id)

ukb_drugs_clean <- ukb_drugs %>%
    select(id, ldl, hdl, trig, chol, statin0, statin1, ldladj, choladj, trigadj, hdladj)


write.table(ukb_drugs_clean, "Data/mesa_olink_processed/UKB_ldl_adj.tsv", sep = "\t", row.names = F)










