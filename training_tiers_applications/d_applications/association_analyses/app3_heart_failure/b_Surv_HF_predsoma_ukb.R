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
outcomes <- fread("/medpop/esp2/lli/UKB_pheno_curation/Results/20251126_CLD_v2_all_cause_dementia_HF_BYver.tsv.gz", data.table = F) %>%
    select(id, contains("Heart_Failure"))
exposures <- fread("Results/UKB_pred/HGB_predicted_UKB_3kolink_log2scaled_5ksoma_r.tsv.gz", data.table = F)
recruitment_site_coding <- fread("Data/recruitment_site_coding", data.table = F)

covs_df_raw <- fread("Data/bingyu_hf/covs_final.tsv", data.table = F)
covs_df_raw <- covs_df_raw %>%
    select(id, Sex, enroll_age.y, in_white_British_ancestry_subset, ever_smoked, BMI, egfrbase, prevalent_disease_Diabetes_All, 
           # recruitment_site,
           prevalent_disease_Coronary_Artery_Disease_SOFT, hypertension, prevalent_disease_Atrial_fibrillation_or_flutter, paste("PC",1:10, sep=""))


# ## covs
covs_knn <- impute.knn(as.matrix(covs_df_raw %>% select(-c(id, Sex, hypertension, contains("PC")))), k = 10, rng.seed = 613)$data
covs_df <- cbind(covs_df_raw %>% select(c(id, Sex, hypertension, contains("PC"))), covs_knn)


no_consent <- data.frame(fread("/medpop/esp2/projects/UK_Biobank/withdrawn_samples/w7089_20250818.csv"))$V1


master_df <- left_join(exposures, covs_df, by = "id") %>%
    left_join(outcomes, by = "id") %>%
    filter(!id %in% no_consent) %>%
    mutate(
        BMI = scale(BMI), 
        enroll_age.y = scale(enroll_age.y),
        egfrbase = scale(egfrbase)
    )


## variable section
outcome_var <- "Heart_Failure_BY"
exposure_vars <- setdiff(names(exposures), "id")
covs <- paste(c(setdiff(names(covs_df), "id")), collapse = "+")

temp <- master_df %>%
    filter(prevalent_disease_Heart_Failure_BY == 0)

temp_skim <- temp %>%
    select(-c(id, all_of(exposure_vars))) %>%
    skim()


##### 03 - analysis #####
results <- data.frame()

for (i in 1:length(exposure_vars)){
    
    ## exposure variables
    exposure_var <- exposure_vars[i]
    
    ## log message
    print(glue("processing olink proteins: {exposure_var}"))
    
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
    write.table(results, "Results/bingyu_hf/UKB_HF_BY_predsoma_impute_10pcs_cadsoft_hypertension01.tsv", sep = "\t", row.names = F, quote = F)
}







