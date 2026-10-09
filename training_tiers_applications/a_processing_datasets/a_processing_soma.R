##### 01 - library #####
library(dplyr)
library(data.table)
library(SomaDataIO)

gc()
rm(list = ls())

setwd("...")



##### 02 - files #####
# LINKER RENAME PROTEINS --------------------------------------------------
## Read linker file
linker <- fread(".../data/chs/5K_Processed_ANML_Normalized/CHS_Markers_5K_ANML-Norm_2022_04.csv",
                data.table = F)
linker <- linker %>% select(AptName, SeqId)


# Remove trailing numbers from AptName
linker$AptName <- sub("\\.\\d+\\.\\d+$", "", linker$AptName)

# Make AptName unique by appending _2, _3, etc., to duplicates
make_unique <- function(names_vector) {
    name_count <- table(names_vector)
    name_usage <- rep(0, length(name_count))
    names(name_usage) <- names(name_count)
    
    for (i in seq_along(names_vector)) {
        name <- names_vector[i]
        name_usage[name] <- name_usage[name] + 1
        if (name_usage[name] > 1) {
            names_vector[i] <- paste(name, name_usage[name], sep = "_")
        }
    }
    return(names_vector)
}

linker$AptName <- make_unique(linker$AptName)
# create a mapping between soma seq and protein Aptamer name
name_mapping <- setNames(linker$AptName, linker$SeqId)


# Rename columns in both datasets based on linker file
rename_columns <- function(dataset, mapping) {
    colnames(dataset) <- sapply(colnames(dataset), function(x) {
        if (x %in% names(mapping)) {
            mapping[x]
        } else {
            x
        }
    })
    return(dataset)
}




# processing soma needs conda environment for the installation of soma data io package
# MESA --------------------------------------------------------------------
mesa_soma_path <- ".../data/mesa/Proteomics_SomaLogic_SIDNO/SS-2212385_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.adat"

## read in the data
mesa_soma_7k_raw <- loadAdatsAsList(mesa_soma_path, collapse = TRUE)


## getting protein columns and positions
# 7596 columns, total 7629 columns
proteomics_cols <- grep("^seq\\.", colnames(mesa_soma_7k_raw), value = TRUE)


## extracting protein data
mesa_soma_7k_prots <- mesa_soma_7k_raw[, proteomics_cols, drop = FALSE]
mesa_soma_7k_attr <- attr(mesa_soma_7k_raw, "Col.Meta")
meta_cols_7k <- mesa_soma_7k_attr[, c("SeqId", "SeqIdVersion", "SomaId", 
                                 "TargetFullName", "Target", "UniProt", 
                                 "EntrezGeneID", "EntrezGeneSymbol"), drop = FALSE]
## 7k metadata
mesa_soma_7k_meta_cols <- c("SampleId", "SampleType", "PercentDilution")  
mesa_soma_7k_meta <- mesa_soma_7k_raw[, mesa_soma_7k_meta_cols, drop = FALSE]

### saving a copy of 7k meta data
write.table(meta_cols_7k, "Data/mesa_olink_processed/mesa_soma_7k_metacols.tsv", sep = "\t", quote = F, row.names = F)


### merging meta data with the proteomics data
mesa_soma_7k_clean <- cbind(mesa_soma_7k_prots, mesa_soma_7k_meta) %>%
    filter(SampleType == "Sample")


### linker part
linker <- fread(".../data/mesa/Proteomics_SomaLogic_SIDNO/MESA_A370_BridgeID_SIDNO_20240906.txt",
                data.table = F) 
linker=linker %>% select(sidno, Barcode_2D1_E1)
linker$Barcode_2D1_E1 <- as.character(linker$Barcode_2D1_E1)
mesa_soma_7k_clean_linker <- merge(linker, mesa_soma_7k_clean, by.x = "Barcode_2D1_E1", by.y = "SampleId")


## cleaning colnames
colnames(mesa_soma_7k_clean_linker) <- gsub("^seq\\.", "", colnames(mesa_soma_7k_clean_linker))
colnames(mesa_soma_7k_clean_linker) <- gsub("\\.", "-", colnames(mesa_soma_7k_clean_linker))

## rename sidno into id
mesa_soma_7k_clean_linker <- mesa_soma_7k_clean_linker %>%
    rename(id = sidno)


## other data
main <- fread(".../data/mesa/SHARe_MesaExam1Main_DS.txt", data.table = F)
#dim 6452 1238
main <- main %>% 
    select(id=subject_id, site=site1c, gender=gender1, race=race1c, 
           age=age1c, bmi=bmi1c, waist=waistcm1, hip=hipcm1) %>% 
    mutate(waisttohipratio=waist/hip)

mesa_soma_7k_clean_linker <- mesa_soma_7k_clean_linker %>% 
    select(-c(Barcode_2D1_E1, SampleType, PercentDilution))
mesa_soma_7k_full <- merge(main, mesa_soma_7k_clean_linker, by = "id")

## renaming colnames to symbols (from seq numbers)
mesa_soma_7k_full_symbol <- rename_columns(mesa_soma_7k_full, name_mapping)
head(mesa_soma_7k_full_symbol[1:10])
#### saving processed 7k data
write.table(mesa_soma_7k_full_symbol, 
            file = gzfile("Data/mesa_olink_processed/mesa_soma_7k_full.tsv.gz"), 
            sep = "\t", row.names = F, quote = F)





# MESA 5K ----------------------------------------------------------------
## lifting over to 5k using soma's software
## this didn't change the version just changed the signal space from 7k to 5k
mesa_soma_5k_raw <- lift_adat(mesa_soma_7k_raw, bridge = "7k_to_5k") # 3487 extra proteins
7629-2326 # 5270 left


# Newer versions of SomaScan typically have additional content, i.e. new reagents added to the multi-plex assay that bind to additional proteins. 
# When lifting to a previous SomaScan version, new reagents that do not exist in the "earlier" assay version assay are scaled by 1.0, and thus maintained, 
# unmodified in the returned object. Users may need to drop these columns in order to combine these data with a previous study from an earlier SomaScan version, 
# e.g. with collapseAdats().

# what column remained unchanged after lifting?
common_cols <- intersect(names(mesa_soma_7k_raw), names(mesa_soma_5k_raw))
unchanged <- purrr::keep(common_cols, ~ identical(mesa_soma_7k_raw[[.x]], mesa_soma_5k_raw[[.x]]))
unchanged_proteins <- unchanged[grep("^seq", unchanged)]


## and do everything for mesa_5k data as well
## getting protein columns and positions
proteomics_cols <- grep("^seq\\.", colnames(mesa_soma_5k_raw), value = TRUE)


## extracting protein data
mesa_soma_5k_prots <- mesa_soma_5k_raw[, proteomics_cols, drop = FALSE] %>%
    select(!all_of(unchanged_proteins))
mesa_soma_5k_attr <- attr(mesa_soma_5k_raw, "Col.Meta")
meta_cols_5k <- mesa_soma_5k_attr[, c("SeqId", "SeqIdVersion", "SomaId", 
                                      "TargetFullName", "Target", "UniProt", 
                                      "EntrezGeneID", "EntrezGeneSymbol"), drop = FALSE]
## 5k metadata
mesa_soma_5k_meta_cols <- c("SampleId", "SampleType", "PercentDilution")  
mesa_soma_5k_meta <- mesa_soma_5k_raw[, mesa_soma_5k_meta_cols, drop = FALSE]

### saving a copy of 5k meta data
write.table(meta_cols_5k, "Data/mesa_olink_processed/mesa_soma_5k_metacols.tsv", sep = "\t", quote = F, row.names = F)


### merging meta data with the proteomics data
mesa_soma_5k_clean <- cbind(mesa_soma_5k_prots, mesa_soma_5k_meta) %>%
    filter(SampleType == "Sample")


### linker part
linker <- fread(".../data/mesa/Proteomics_SomaLogic_SIDNO/MESA_A370_BridgeID_SIDNO_20240906.txt",
                data.table = F) 
linker=linker %>% select(sidno, Barcode_2D1_E1)
linker$Barcode_2D1_E1 <- as.character(linker$Barcode_2D1_E1)
mesa_soma_5k_clean_linker <- merge(linker, mesa_soma_5k_clean, by.x = "Barcode_2D1_E1", by.y = "SampleId")


## cleaning colnames
colnames(mesa_soma_5k_clean_linker) <- gsub("^seq\\.", "", colnames(mesa_soma_5k_clean_linker))
colnames(mesa_soma_5k_clean_linker) <- gsub("\\.", "-", colnames(mesa_soma_5k_clean_linker))

## rename sidno into id
mesa_soma_5k_clean_linker <- mesa_soma_5k_clean_linker %>%
    rename(id = sidno)


## other data
main <- fread(".../data/mesa/SHARe_MesaExam1Main_DS.txt", data.table = F)
#dim 6452 1238
main <- main %>% 
    select(id=subject_id, site=site1c, gender=gender1, race=race1c, 
           age=age1c, bmi=bmi1c, waist=waistcm1, hip=hipcm1) %>% 
    mutate(waisttohipratio=waist/hip)

mesa_soma_5k_clean_linker <- mesa_soma_5k_clean_linker %>% 
    select(-c(Barcode_2D1_E1, SampleType, PercentDilution))
mesa_soma_5k_full <- merge(main, mesa_soma_5k_clean_linker, by = "id")

## rename columns
mesa_soma_5k_full_symbol <- rename_columns(mesa_soma_5k_full, name_mapping)
head(mesa_soma_5k_full_symbol[1:10])

#### saving processed 5k data
write.table(mesa_soma_5k_full_symbol, 
            file = gzfile("Data/mesa_olink_processed/mesa_soma_5k_full.tsv.gz"), 
            sep = "\t", row.names = F, quote = F)




# Reference --------------------------------------------------------------
# 5k soma (CHS) limitation
# soma has no missing value so it's fine
chs_raw <- fread(".../CHS_soma5k_bmiwaisttohipratiolipids.txt.gz", data.table = F)[,-1]
chs_proteins <- chs_raw[,c(1,11:ncol(chs_raw))]
soma_5k_ref <- setdiff(colnames(chs_proteins), "id")
write.table(soma_5k_ref, "Data/soma_5kprot_ref", sep = " ", row.names = F, col.names = F, quote = F)


