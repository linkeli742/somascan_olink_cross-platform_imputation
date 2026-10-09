##### 01 - library #####
library(dplyr)
library(data.table)

gc()
rm(list = ls())

setwd("...")

##### 02 - files #####
# 3k olink (UKB limitation)
# olink has missing value
ukb_olink_raw <- fread(".../olink_3kfull.tsv.gz", data.table = F)
ukb_olink_raw <- ukb_olink_raw %>%
    select(-ins_index)

row_ok <- rowMeans(is.na(ukb_olink_raw[-1])) <= 0.1
ukb_olink_short <- ukb_olink_raw[row_ok,]

# isolate olink columns
qc_cols <- setdiff(names(ukb_olink_short), "id")
## keeps the 10% missingness
keep_cols <- qc_cols[colMeans(is.na(ukb_olink_short[qc_cols])) <= 0.1]

# but row(pat wise) this will be different
ukb_olink_qc_na10 <- ukb_olink_short %>%
    select(id, all_of(keep_cols)) %>%
    filter(rowMeans(is.na(across(all_of(keep_cols)))) <= 0.1)
sum(is.na(ukb_olink_qc_na10)) # 1949698


## There still are lots of missingness to impute but now it had passed QC
olink_3k_ref <- setdiff(names(ukb_olink_qc_na10), "id")

write.table(olink_3k_ref, "Data/olink_3kprot_ref", sep = " ", row.names = F, quote = F)
write.table(ukb_olink_qc_na10, gzfile("Data/UKB_Olink_3k_log2_nodup_final.tsv.gz"), sep = "\t", row.names = F, quote = F)



# mesa's new olink data
mesa_olink_raw <- fread(".../data/mesa/IRC_Data_Transer_MESA_12232024/SMP_IntensityNormalized_12232024.csv")
mesa_id_sample <- fread(".../data/mesa/IRC_Data_Transer_MESA_12232024/Mapping_SMP_Plate_12232024.csv")
mesa_id_sample_clean <- mesa_id_sample %>%
    filter(Exam == 1) %>%
    select(sampleid, sidno) # 6116
mesa_soma_raw <- fread("Data/mesa_olink_processed/mesa_soma_5k_full.tsv.gz", data.table = F) # 5737


##### 03 - process #####
mesa_olink_clean <- left_join(mesa_id_sample_clean, mesa_olink_raw, by = "sampleid") %>%
    rename(id = sidno) %>%
    relocate(id, .before = 1) %>%
    filter(id %in% mesa_soma_raw$id)

# remove control batches
mesa_olink_clean_nocontrol <- mesa_olink_clean %>%
    filter(!grepl("control", Assay, ignore.case = T)) %>%
    # quality control parameters
    filter(QC_Warning == "PASS") %>%
    select(id, Assay, Panel, NPX)


olink_3k_ref <- fread("Data/olink_3kprot_ref", data.table = F)$x


# finalizing --------------------------------------------------------------
mesa_olink_clean_nocontrol_short <- mesa_olink_clean_nocontrol %>%
    filter(Assay %in% olink_3k_ref)
length(unique(mesa_olink_clean_nocontrol_short$Assay))

## check olink
duplicates <- mesa_olink_clean_nocontrol_short %>%
    group_by(id, Assay) %>%
    filter(n() > 1) %>%
    ungroup()

duplicates_clean <- duplicates %>%
    select(-id, -NPX) %>%
    filter(!duplicated(.))

duplicates_clean_assay <- duplicates_clean[duplicated(duplicates_clean$Assay) | duplicated(duplicates_clean$Assay, fromLast = T),]

dup_proteins <- unique(duplicates_clean_assay$Assay)

# olink clean
mesa_olink_clean_nocontrol_short_filtered <- mesa_olink_clean_nocontrol_short %>%
    filter(
        !(
            (Assay == "TNF"   & Panel != "Cardiometabolic") |
                (Assay == "IL6"   & Panel != "Oncology") |
                (Assay == "CXCL8" & Panel != "Oncology") |
                (Assay == "IDO1"  & Panel != "Cardiometabolic_II") |
                (Assay == "LMOD1" & Panel != "Neurology_II") |
                (Assay == "SCRIB" & Panel != "Cardiometabolic_II")
        )
    ) %>% 
    filter(!duplicated(.))

length(unique(mesa_olink_clean_nocontrol_short_filtered$Assay)) # 2920 olink proteins

dup <- mesa_olink_clean_nocontrol_short_filtered %>%
    filter(!duplicated(.)) %>%
    group_by(id, Assay) %>%
    filter(n() > 1) %>%
    ungroup()
dup_id <- unique(dup$id)   


# there are 7 identical rows for some reason
# these rows have different raw count in the original file but have the same NPX normalized value
mesa_olink_clean_nocontrol_short_filtered_wide <- mesa_olink_clean_nocontrol_short_filtered %>%
    # just keep one of the values since the NPX value is the same
    arrange(id, Assay) %>%
    distinct(id, Assay, .keep_all = TRUE) %>%
    pivot_wider(names_from = Assay, id_cols = id, values_from = NPX)


write.table(mesa_olink_clean_nocontrol_short_filtered_wide, gzfile("Data/MESA_Olink_3k_log2_nodup_final.tsv.gz"), sep = "\t", row.names = F, quote = F)








