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

### merging and cleaning
soma_to_olink <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F)
soma_to_olink_merged <- left_join(soma_to_olink, hpa, by = "protein")
soma_to_olink_merged <- soma_to_olink_merged %>%
    select(protein, `Disease involvement`, `Protein class`)



### separate protein class
soma_olink_long_tags_prot <- soma_to_olink_merged %>%
    mutate(protein_tags = str_squish(`Protein class`)) %>%
    separate_rows(protein_tags, sep = "\\s*,\\s*") %>%   # split on commas, trim spaces
    filter(!is.na(protein_tags), protein_tags != "")

# table(soma_olink_long_tags_prot$protein_tags)

protein_class <- soma_olink_long_tags_prot %>%
    filter(!grepl("genes|drug targets", protein_tags))


table(protein_class$protein_tags)

## deriving my own class
protein_class1 <- protein_class %>%
    group_by(protein) %>%
    summarise(
        derived_class = case_when(
            any(protein_tags %in% c("Plasma proteins", "Predicted secreted proteins")) ~ "Secreted_Plasma",
            any(protein_tags == "Predicted membrane proteins") ~ "Membrane",
            any(protein_tags %in% c("Enzymes", "Predicted intracellular proteins")) ~ "Intracellular",
            TRUE ~ NA
        )
    )

membrane_and_plasma_list <- protein_class %>%
    group_by(protein) %>%
    filter(any(protein_tags == "Predicted membrane proteins") & 
               any(protein_tags == "Plasma proteins")) %>%
    ungroup()


## needing couple more proteins
setdiff(soma_to_olink$protein, protein_class1$protein) ## HPA just dont have these
table(protein_class1$derived_class)
# head(protein_class1[protein_class1$derived_class == "Membrane",])


##### 03 - appending classes ####
# LILRA3, MENT, LEG1: Secreted_Plasma
# DDX58: Intracellular
# KIR2DL2: membrane
# source: https://www.uniprot.org/uniprotkb/P43627/entry#subcellular_location (this is for KIR2DL2)


write.table(protein_class1, "Results/protein_importance/HPA_derived_protein_class.tsv", 
            sep="\t", row.names = F, quote = F)



