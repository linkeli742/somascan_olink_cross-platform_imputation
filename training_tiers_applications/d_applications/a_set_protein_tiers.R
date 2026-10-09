##### 01 - library #####
library(dplyr)
library(data.table)
library(ggplot2)
library(colorspace)

gc()
rm(list = ls())

setwd("...")


#### 02 - SOMA OLINK files #####
soma5k_olink3k_metrics_raw <- fread("Results/metrics/HGB_log2scaled_5ksoma_3kolink_results.tsv", data.table = F) %>%
    filter(has_soma == T) %>%
    select(-has_soma)


##### 03 - adding rank values #####
soma5k_olink3k_ranks <- soma5k_olink3k_metrics_raw %>%
    mutate(tiers = case_when(
        pretrain_r >= 0.7 & posttrain_r >= 0.7 ~ 1,
        pretrain_r < 0.1 & posttrain_r < 0.1 ~ 4,
        posttrain_r >= 0.4 ~ 2,
        posttrain_r >= 0.1 & posttrain_r < 0.4 ~ 3,
        TRUE ~ NA
    ))

table(soma5k_olink3k_ranks$tiers, useNA = "ifany") 

soma5k_olink3k_ranks[is.na(soma5k_olink3k_ranks$tiers),]  



### no NA categorization
soma5k_olink3k_ranks <- soma5k_olink3k_metrics_raw %>%
    mutate(tiers = case_when(
        pretrain_r >= 0.7 & posttrain_r >= 0.7 ~ 1,
        pretrain_r < 0.1 & posttrain_r < 0.1 ~ 4,
        posttrain_r >= 0.4 ~ 2,
        posttrain_r < 0.4 ~ 3,
        TRUE ~ NA
    ))   

table(soma5k_olink3k_ranks$tiers, useNA = "ifany")    




# 04 - Rank by Gain ------------------------------------------------------------
soma5k_olink3k_ranks1 <- soma5k_olink3k_ranks %>%
    mutate(gain = posttrain_r - pretrain_r) %>%
    group_by(tiers) %>%
    arrange(gain) %>%
    mutate(index = row_number(gain)) %>%
    ungroup() %>%
    arrange(tiers, desc(index)) %>%
    rename(pretrain_r2 = pretrain_rsq, posttrain_r2 = posttrain_rsq)

write.table(soma5k_olink3k_ranks1, "Results/rank_by_gain/soma5k_olink3k.tsv", sep = "\t", row.names = F, quote = F)




#### 02 - OLINK SOMA files #####
olink3k_soma5k_metrics_raw <- fread("Results/metrics/HGB_3kolink_log2scaled_5ksoma_results.tsv", data.table = F) %>%
    filter(has_soma == T) %>%
    select(-has_soma)

##### 03 - adding rank values #####
olink3k_soma5k_ranks <- olink3k_soma5k_metrics_raw %>%
    mutate(tiers = case_when(
        pretrain_r >= 0.7 & posttrain_r >= 0.7 ~ 1,
        pretrain_r < 0.1 & posttrain_r < 0.1 ~ 4,
        posttrain_r >= 0.4 ~ 2,
        posttrain_r >= 0.1 & posttrain_r < 0.4 ~ 3,
        TRUE ~ NA
    ))

table(olink3k_soma5k_ranks$tiers, useNA = "ifany")


### no NA categorization
olink3k_soma5k_ranks <- olink3k_soma5k_metrics_raw %>%
    mutate(tiers = case_when(
        pretrain_r >= 0.7 & posttrain_r >= 0.7 ~ 1,
        pretrain_r < 0.1 & posttrain_r < 0.1 ~ 4,
        posttrain_r >= 0.4 ~ 2,
        posttrain_r < 0.4 ~ 3,
        TRUE ~ NA
    ))

table(olink3k_soma5k_ranks$tiers, useNA = "ifany")    
write.table(olink3k_soma5k_ranks, "Results/rank_by_gain/olink3k_soma5k.tsv", sep = "\t", row.names = F)



