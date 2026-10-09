##### 01 - library #####
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggh4x)

gc()
rm(list=ls())

setwd("...")


##### 02 - files #####
soma_olink_raw <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F)
olink_soma_raw <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = F)

# LILRA3, MENT, LEG1: Secreted_Plasma
# DDX58: Intracellular
# KIR2DL2: membrane

protein_class1 <- fread("Results/protein_importance/HPA_derived_protein_class.tsv", data.table = F)
protein_class2 <- data.frame(
    protein = c("LILRA3", "MENT", "LEG1", "DDX58", "KIR2DL2"),
    derived_class = c("Secreted_Plasma","Secreted_Plasma","Secreted_Plasma","Intracellular", "Membrane")
)

protein_class <- rbind(protein_class1, protein_class2) %>%
    mutate(derived_class=ifelse(derived_class == "Secreted_Plasma", "Secreted or Plasma", derived_class))


soma_olink_class <- left_join(soma_olink_raw, protein_class, by = "protein")
olink_soma_class <- left_join(olink_soma_raw, protein_class, by = "protein")

### all the correlations
baseline_corr <- soma_olink_class %>% select(protein, pretrain_r, derived_class) %>%
    rename(corr = pretrain_r) %>%
    mutate(impute_label = "Baseline Correlation",
           direction_label = "Bi-directional Identical")

soma_olink_class_posttrain <- soma_olink_class %>% select(protein, posttrain_r, derived_class) %>%
    rename(corr = posttrain_r) %>%
    mutate(impute_label = "Imputated Correlation",
           direction_label = "SomaScan to Olink")

olink_soma_class_posttrain <- olink_soma_class %>% select(protein, posttrain_r, derived_class) %>%
    rename(corr = posttrain_r) %>%
    mutate(impute_label = "Imputated Correlation",
           direction_label = "Olink to SomaScan"
           )


## merging to long
corr_df <- rbind(baseline_corr, soma_olink_class_posttrain, olink_soma_class_posttrain) %>%
    mutate(direction_label = factor(direction_label, levels = c("Bi-directional Identical", "SomaScan to Olink", "Olink to SomaScan")))



#### 03 - plotting #####
## SO direction
ggplot(corr_df, 
       aes(x = derived_class, y = corr, fill = derived_class)) +
    facet_nested(.~ impute_label + direction_label, scales="free_x", space = "free_x") +
    geom_violin(alpha = 0.7) +
    geom_jitter(width = 0.3, alpha = 0.2, size = 0.7) +
    labs(
         y = "Pearson Correlation",
         x = "") +
    theme_bw(base_size = 12) +
    guides(fill = "none")



ggsave("Results/figures/imputed_corr_so_os.png", dpi=300, width = 10, height = 5, bg = "white")





