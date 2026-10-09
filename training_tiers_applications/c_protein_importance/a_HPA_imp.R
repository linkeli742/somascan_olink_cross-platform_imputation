##### 01 - library #####
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)

gc()
rm(list = ls())


setwd("...")



##### 02 - files ######
hpa <- fread("Data/HPA/proteinatlas.tsv.gz", data.table = F) %>%
    rename(protein  = Gene)
names <- data.frame(colnames(hpa))

## overlapping proteins that has ranks based on Gain
## this will be the same for soma to olink or olink to soma 
## thus only one dataset will be read as reference protein list
soma_to_olink <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F)
soma_to_olink_merged <- left_join(soma_to_olink, hpa, by = "protein")
soma_to_olink_merged <- soma_to_olink_merged %>%
    select(protein, `Disease involvement`, `Protein class`)

## I think one to confirm the scoring system is to check with malacards 
# after picking out top proteins based on our system to check

## there are 899 proteins with no associated tags
sum(soma_to_olink_merged$`Disease involvement` == "", na.rm = T)
## there are actually 5 proteins that is not in the HPA database
## correspondingly, these proteins have less associations based on MalaCards

soma_to_olink_merged[is.na(soma_to_olink_merged$`Disease involvement`),] 
# there are 5 proteins without records that is in tier 2



#### long format tags -- among proteins with tags
soma_olink_long_tags <- soma_to_olink_merged %>%
    mutate(tags = str_squish(`Disease involvement`)) %>%
    separate_rows(tags, sep = "\\s*,\\s*") %>%   # split on commas, trim spaces
    filter(!is.na(tags), tags != "")


soma_olink_long_tags_prot <- soma_to_olink_merged %>%
    mutate(protein_tags = str_squish(`Protein class`)) %>%
    separate_rows(protein_tags, sep = "\\s*,\\s*") %>%   # split on commas, trim spaces
    filter(!is.na(protein_tags), protein_tags != "")


# total count in 1737 proteins
tag_counts <- soma_olink_long_tags %>%
    count(tags, sort = TRUE)

protein_tag_counts <- soma_olink_long_tags_prot %>%
    count(protein_tags, sort = TRUE)


# 452 disease variant
# 422 cancer-related genes
# 199 FDA approved drug targets
# there are also potential drug targets tag

protein_tag_flags <- soma_olink_long_tags %>%
    mutate(flag = 1L) %>%
    pivot_wider(
        names_from = tags,
        values_from = flag,
        values_fill = 0
    ) 

protein_prot_tag_flags <- soma_olink_long_tags_prot %>%
    mutate(flag = 1L) %>%
    pivot_wider(
        names_from = protein_tags,
        values_from = flag,
        values_fill = 0
    ) 

# how many tags do a protein have?
protein_tag_flags <- protein_tag_flags %>%
    select(-`Protein class`) %>%
    rowwise() %>%
    mutate(rowsum = sum(c_across(-c(protein, `Disease involvement`)), 
                        na.rm = TRUE), .after = protein) %>%
    ungroup() %>%
    arrange(desc(rowsum))
## there are 899 missing tags, 5 missing, thus the final dataset have 833 proteins with different numbers of tags

protein_protein_tag_flags <- protein_prot_tag_flags %>%
    select(-c(`Protein class`, `Disease involvement`)) %>%
    rowwise() %>%
    mutate(rowsum = sum(c_across(-c(protein)), 
                        na.rm = TRUE), .after = protein) %>%
    ungroup() %>%
    arrange(desc(rowsum)) %>%
    relocate(`Disease related genes`, .after = rowsum) %>%
    relocate(`Human disease related genes`, .after = `Disease related genes`) %>%
    relocate(`FDA approved drug targets`, .after = `Human disease related genes`) %>%
    relocate(`Potential drug targets` , .after = `FDA approved drug targets`)


write.table(protein_protein_tag_flags, "Results/protein_importance/HPA_importance_protclasstags.tsv", 
            sep = "\t", row.names = F, quote = F)







    










