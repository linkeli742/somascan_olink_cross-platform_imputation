##### 01 - library #####
library(dplyr)
library(data.table)


gc()
rm(list = ls())

setwd("...")




##### 02 - files #####
no_consent <- data.frame(fread("/medpop/esp2/projects/UK_Biobank/withdrawn_samples/w7089_20250818.csv"))$V1
pheno <- fread("/medpop/esp2/lli/UKB_pheno_curation/Results/20250711_pheno_noGP_clean.tsv.gz", data.table = F)

pheno_short <- pheno %>%
    select(id, contains("Diabetes_all"), contains("Heart_Failure"), contains("Coronary_Artery_Disease"), contains("Atrial_fibrillation_or_flutter"))

pheno_final <- pheno_short %>%
    select(id, contains("Heart_Failure"))

write.table(pheno_final, "Data/bingyu_hf/pheno_heart_failure.tsv", sep = "\t", row.names = F, quote = F)

pheno_covs <- pheno_short %>%
    select(id, prevalent_disease_Diabetes_All, contains("prevalent_disease_Coronary_Artery_Disease"), prevalent_disease_Atrial_fibrillation_or_flutter)


# covs <- fread("/medpop/esp2/mzekavat/UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.REAL.txt", data.table = F)
# covs$IDs_toRemove_SampleQCv3 = ifelse((!(covs$Submitted_Gender == covs$Inferred_Gender) |  covs$Non_Consented== 1),1,0)
# covs = covs[which(covs$IDs_toRemove_SampleQCv3==0),]
# 
# covs_short <- covs %>%
#     select(id, Sex, enroll_age.y, in_white_British_ancestry_subset, paste("PC",1:10,sep=""), ever_smoked, BMI, SBP, DBP) %>%
#     filter(!id %in% no_consent)
# 
# write.table(covs_short, "Data/bingyu_hf/covs_short.tsv", sep = "\t", row.names = F, quote = F)


covs_short <- fread("Data/bingyu_hf/covs_short.tsv", data.table = F)
ukb_recruitment_site <- fread("Data/ukb_recruitment_site.tsv", data.table = F) %>%
    rename(id=eid, recruitment_site = `54-0.0`)

hypertensive_med <- fread("/medpop/esp2/smjcho/outgoing/for_mike/ukb_antihtn.csv", data.table = F) %>%
    rename(id = eid) %>%
    filter(!id %in% no_consent) %>%
    select(id, antihypertensive_meds = antihtnbase)
egfr <- fread("/medpop/esp2/smjcho/outgoing/for_aeron/ukb_egfr_baseline.csv") %>%
    rename(id = eid)

covs_all <- merge(covs_short, egfr, by = "id") %>%
    left_join(ukb_recruitment_site, by ="id") %>%
    left_join(hypertensive_med,by="id") %>%
    left_join(pheno_covs, by = "id") %>%
    mutate(hypertension = if_else((SBP >140 & DBP > 90) | antihypertensive_meds == 1,1,0))

table(covs_all$recruitment_site, useNA = "ifany")


write.table(covs_all, "Data/bingyu_hf/covs_final.tsv", sep = "\t", row.names = F, quote = F)





