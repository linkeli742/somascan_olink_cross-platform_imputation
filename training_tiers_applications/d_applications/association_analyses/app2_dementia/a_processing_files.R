##### 01 - library #####
library(dplyr)
library(data.table)

gc()
rm(list = ls())


setwd("...")



##### 02 - files #####
somaolink_tiers <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F)
latelife_incd_prots <- fread("Data/keenan_walker/latelife_incidence_dementia_proteins.txt", data.table = F, header = F)$V1

somaolink_tiers_short <- somaolink_tiers %>% filter(protein %in% latelife_incd_prots)


### APOE4 coded as WT or >1 allele
### education less than high school/high school, general education diploma [GED], or vocational school/college, graduate or professional school
### eGFR, BMI, SBP >140mmHg or DBP >= 90 mmhg or use of hypertensive medication
### diabetes, current smoker or not
### age, sex, race, education, apoe4 (model 1)
### model 1 + eGFR
### BMI, diabetes, hypertension, and smoking status
### bonferroni corrected two-sided P < 0.05 (0.05/4800)
### APOE： rs429358 (T>C), rs7412(T>C)
### Each of the haplotypes of the ApoE gene (ε2: T/T, ε3: C/T, ε4: C/C) was generated from combinations of these two SNPs.


### bgenix:
# .../software/bgen.tgz/build/apps/bgenix -g /broad/ukbb/imputed_v3/ukb_imp_chr19_v3.bgen \
#     -incl-rsids rs429358 rs7412 \
#     > .../UKB_e4_geno_bgen

# .../software/./plink2 \
#     --bgen .../UKB_e4_geno_bgen ref-first \
#     --sample .../UK_Biobank/linkers/app7089/ukb22828_c1_b0_v3_s487207.sample \
#     --recode A \
#     --recode-allele .../apoe_e4_alleles.txt \
#     --out .../Data/keenan_walker/UKB_e4_geno_a





### UKB covariates
id_450k <- fread("/medpop/esp2/mesbah/datasets/CHIP/UKBB/450k/ukb450k.eid.list")$V1
no_consent <- data.frame(fread("/medpop/esp2/projects/UK_Biobank/withdrawn_samples/w7089_20250818.csv"))$V1
relate_remove <- fread("/medpop/esp2/zyu/chip_protemoics/listofremovefor450K.txt")$IID


# ## make covariates
# covs <- fread("/medpop/esp2/mzekavat/UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.REAL.txt", data.table = F)
# covs$IDs_toRemove_SampleQCv3 = ifelse((!(covs$Submitted_Gender == covs$Inferred_Gender) |  covs$Non_Consented== 1),1,0)
# covs = covs[which(covs$IDs_toRemove_SampleQCv3==0),]
# 
# covs_short <- covs %>%
#     select(id, Sex, enroll_age.y, in_white_British_ancestry_subset, paste("PC",1:10,sep=""), ever_smoked, BMI, SBP, DBP) %>%
#     filter(id %in% id_450k) %>%
#     filter(!id %in% no_consent) %>%
#     filter(!id %in% relate_remove)
# 
# write.table(covs_short, "Data/keenan_walker/covs_short.tsv", sep = "\t", row.names = F, quote = F)


covs_short <- fread("Data/keenan_walker/covs_short.tsv", data.table = F)

egfr <- fread("/medpop/esp2/smjcho/outgoing/for_aeron/ukb_egfr_baseline.csv") %>%
    rename(id = eid)

hypertensive_med <- fread("/medpop/esp2/smjcho/outgoing/for_mike/ukb_antihtn.csv", data.table = F) %>%
    rename(id = eid) %>%
    filter(!id %in% no_consent) %>%
    filter(id %in% id_450k) %>%
    filter(!id %in% relate_remove) %>%
    select(id, antihypertensive_meds = antihtnbase)


## defining e4e4 genotype
apoe_e4_raw <- fread("Data/keenan_walker/UKB_e4_geno_a.raw", data.table = F) %>%
    dplyr::select(IID, rs429358_C, rs7412_C) %>%
    mutate(apoe_e4 = ifelse(rs429358_C >0 & rs7412_C >0, 1, 0)) %>%
    rename(id = IID) %>%
    filter(id %in% id_450k) %>%
    filter(!id %in% no_consent) %>%
    filter(!id %in% relate_remove)

table(apoe_e4_raw$apoe_e4, useNA = "ifany")

# zcat /medpop/esp2/projects/UK_Biobank/baskets/2008463/ukb47823.csv.gz | awk -F',' '{print $1 "," $6762}' > /medpop/esp2/lli/somalink_olink/Data/keenan_walker/ukb_education.csv 

ukb_education_raw <- fread("Data/keenan_walker/ukb_education.tsv", data.table = F)
ukb_education <- ukb_education_raw %>%
    rename(id = eid, edu_raw = `6138-0.0`) %>%
    # less than high school/high school, general education diploma [GED], or vocational school/college, graduate or professional school
    mutate(education = case_when(
        edu_raw %in% c(3,4,-7) ~ "less_than_and_high_school",
        edu_raw %in% c(1,6) ~ "graduate_professional_school",
        edu_raw %in% c(2) ~ "general_education_diploma",
        edu_raw %in% c(5) ~ "vocational_school",
        TRUE ~ NA
    )) %>%
    select(id, education)


covs_all <- merge(covs_short, apoe_e4_raw, by = "id") %>%
    left_join(ukb_education, by = "id") %>%
    left_join(egfr, by = "id")  %>%
    left_join(hypertensive_med,by="id") %>%
    mutate(hypertension = if_else((SBP >140 & DBP > 90) | antihypertensive_meds == 1,1,0))



write.table(covs_all, "Data/keenan_walker/UKB_AD_covs.tsv", row.names = F, sep = "\t", quote = F)
    




##### pheno files
pheno <- fread(".../20250711_pheno_noGP_clean.tsv.gz", data.table = F)
pheno_short <- pheno %>%
    select(id, contains("Dementia"), contains("Diabetes"))

write.table(pheno_short, "Data/keenan_walker/pheno_ad_diabetes.tsv", sep = "\t", row.names = F, quote = F)


