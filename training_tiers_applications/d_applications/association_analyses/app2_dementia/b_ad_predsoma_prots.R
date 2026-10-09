##### 01 - library #####
library(dplyr)
library(data.table)
library(survival)
library(glue)
library(impute)
library(skimr)

gc()
rm(list = ls())


setwd("...")



##### 02 - files #####
outcomes <- fread("Data/keenan_walker/pheno_ad_diabetes.tsv", data.table = F)
exposures <- fread("Results/UKB_pred/HGB_predicted_UKB_3kolink_log2scaled_5ksoma_r.tsv.gz", data.table = F)
covs_df_raw <- fread("Data/keenan_walker/UKB_AD_covs.tsv", data.table = F)

covs_df_raw <- covs_df_raw %>%
    select(-c(SBP, DBP, rs429358_C, rs7412_C, antihypertensive_meds, paste("PC",1:10, sep = "")))


## covs
covs_knn <- impute.knn(as.matrix(covs_df_raw %>% select(-c(id, Sex, education, apoe_e4, hypertension))),
                       k = 10, rng.seed = 613)$data

covs_df <- cbind(covs_df_raw %>% select(id, Sex, education, apoe_e4, hypertension), covs_knn)


id_450k <- fread("/medpop/esp2/mesbah/datasets/CHIP/UKBB/450k/ukb450k.eid.list")$V1
no_consent <- data.frame(fread("/medpop/esp2/projects/UK_Biobank/withdrawn_samples/w7089_20250818.csv"))$V1
relate_remove <- fread("/medpop/esp2/zyu/chip_protemoics/listofremovefor450K.txt")$IID


master_df <- merge(exposures, covs_df, by = "id") %>%
    left_join(outcomes, by = "id") %>%
    filter(!id %in% no_consent) %>%
    filter(id %in% id_450k) %>%
    filter(!id %in% relate_remove) %>%
    mutate(
        BMI = scale(BMI), 
        enroll_age.y = scale(enroll_age.y),
        egfrbase = scale(egfrbase)
    )



## variable section
outcome_var <- "Alzheimers_or_Dementia"
exposure_vars <- setdiff(names(exposures), "id")
covs <- paste(c(setdiff(names(covs_df), "id"), "prevalent_disease_Diabetes_All"), collapse = "+")

temp <- master_df %>%
    filter(prevalent_disease_Alzheimers_or_Dementia == 0)


temp_skim <- temp %>%
    select(-c(id, all_of(exposure_vars))) %>%
    skim()



##### 03 - analysis #####
results <- data.frame()

for (i in 1:length(exposure_vars)){
    
    ## exposure variables
    exposure_var <- exposure_vars[i]
    
    ## log message
    print(glue("processing soma proteins: {exposure_var}"))
    
    f <- as.formula(glue("Surv(followup_{outcome_var}, incident_disease_{outcome_var}) ~ {exposure_var} + {covs}"))
    r <- coxph(f, data = temp)
    
    result <- data.frame(
        outcomes = outcome_var,
        exposures = exposure_var,
        
        hr = summary(r)$coefficients[1,2],
        lci = exp(confint(r))[1,1],
        uci = exp(confint(r))[1,2],
        pval = summary(r)$coefficients[1,5]
    )
    results <- rbind(results, result)
    
    ## save results
    write.table(results, "Results/keenan_walker/UKB_AD_somas_ad_no10pc_hypertension01_edunew.tsv", sep = "\t", row.names = F, quote = F)
}


