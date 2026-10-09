##### 01 - library #####
library(dplyr)
library(data.table)
library(readxl)
library(stringr)
library(glue)


gc()
rm(list = ls())

setwd("...")


##### 02 - files #####

### need to separate batches for manual curation
soma_olink <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F) %>%
    select(protein)

# olink_soma <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = F)
# 
# protlist <- soma_olink$protein
# 
# batch_size = 10
# split_protlist <- split(protlist, ceiling(seq_along(protlist) / batch_size))

# separating and writing chunks
# for (i in seq_along(split_protlist)) {
#     chunk <- split_protlist[[i]]
#     chunk_string <- paste(chunk, collapse = ",")
#     # Write to tab-separated txt file
#     write.table(chunk_string, file = paste0("Data/malacards_protlist/protlist_chunk", i, ".txt"), sep = "\t", row.names = F, col.names = F, quote = FALSE)
# }

protlist <- soma_olink$protein



##### 03 - processing queried results #####
all_files <- list.files("Data/curated_malacards/")
n_files <- length(list.files("Data/curated_malacards/"))

malacards <- data.frame()

for (i in 1:n_files){
    
    cat(glue("\n----processing {i}th file----\n"))
    
    file_name <- all_files[i]
    
    df <- read_excel(glue("Data/curated_malacards/{file_name}"), sheet = "MalaCardsInferredDisorders")
    df_short <- df %>%
        select(protein = Symbol, disease = Accession, Sources, Publications, GeneScore, DisorderScore) %>%
        # number of pubs equals to the number of || + 1
        mutate(n_pubs = str_count(Publications, "\\|\\|") + 1) %>%
        mutate(n_pubs = ifelse(is.na(n_pubs),0,n_pubs)) %>%
        relocate(c(GeneScore, DisorderScore, n_pubs), .after = disease)
    
    malacards <- rbind(malacards, df_short)
    
}

malacards_imp <- left_join(soma_olink, malacards, by = "protein")

## fill in the NA of the n_pubs with 0, those proteins were missing from the database after manual checking
malacards_imp$n_pubs[is.na(malacards_imp$n_pubs)] <- 0
## save importance dataframe (without the sum)
write.table(malacards_imp, gzfile("Results/protein_importance/MALACARDS_proteins.tsv.gz"), sep = "\t", row.names = F)


## calculate sum of pubs as importance ranking of proteins within each tier
malacards_imp_sum <- malacards_imp %>%
    group_by(protein) %>%
    mutate(tot_pubs = sum(n_pubs)) %>%
    ungroup() %>%
    select(protein, tot_pubs) %>%
    filter(!duplicated(.))
write.table(malacards_imp_sum, gzfile("Results/protein_importance/MALACARDS_proteins_sumofpubs.tsv.gz"), sep = "\t", row.names = F)



