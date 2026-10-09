##### 01 - library #####
library(dplyr)
library(jsonlite)
library(glue)
library(stringr)
library(purrr)

gc()
rm(list=ls())

setwd("...")

#### 02 - files #####
all_folders <- setdiff(list.files("Results/alphafold3/"), "terms_of_use.md")
all_proteins <- toupper(str_match(all_folders, "\\d{2}_\\d{2}_\\d{2}_\\d{2}_\\d{2}_([^_]+)")[,2])

results_list <- list()

for (i in seq_along(all_folders)) {

    #folder name
    folder_name <- all_folders[i]
    #  protein name
    protein_name <- toupper(str_match(folder_name, "\\d{2}_\\d{2}_\\d{2}_\\d{2}_\\d{2}_([^_]+)")[,2])

    ## all of the files
    json_paths <- glue("Results/alphafold3/{folder_name}/fold_{folder_name}_summary_confidences_{0:4}.json")

    ## if all the files exist
    existing_paths <- json_paths[file.exists(json_paths)]

    if (length(existing_paths) > 0) {
        protein_data <- map_dfr(existing_paths, jsonlite::fromJSON, .id = "seed_id") %>%
            mutate(
                protein = protein_name,
                folder = folder_name,
                .before = 1
            )
        results_list[[i]] <- protein_data
    }
}

final_df <- bind_rows(results_list)

write.table(final_df, "Results/alphafold3/alphafold3_data_aggregated.tsv", sep="\t", row.names = F, quote = F)


final_df <- fread("Results/alphafold3/alphafold3_data_aggregated.tsv", data.table = F)

final_df_1 <- final_df %>%
    group_by(protein) %>%
    # mutate(mean_fraction_disordered = mean(fraction_disordered)) %>%
    mutate(mean_fraction_disordered = mean(fraction_disordered), sd_fraction_disordered = sd(fraction_disordered)) %>%
    ungroup() %>%
    select(protein, mean_fraction_disordered, sd_fraction_disordered) %>%
    filter(!duplicated(.))


