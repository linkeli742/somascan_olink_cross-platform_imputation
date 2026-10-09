##### 01 - library #####
library(dplyr)
library(data.table)
library(tidyverse)
library(scales)
library(RColorBrewer)
library(ggalluvial)
library(ggh4x)

gc()
rm(list = ls())

setwd("...")


##### 02 - files ####
### ukb results
files <- list.files(path = "Results/ukb_bmi_assoc_metrics", pattern = "feature_imp", full.names = TRUE)

ukb_top10_all <- data.frame()

for (i in 1:length(files)){
    file <- files[i]
    df <- fread(file) %>%
        mutate(label = file) %>%
        arrange(desc(Importance)) %>%
        head(6)
    ukb_top10_all <- rbind(ukb_top10_all, df)
}

ukb_top10_all <- ukb_top10_all %>%
    mutate(
        label = case_when(
            label == "Results/ukb_bmi_assoc_metrics/bencholink_overlapping_proteins_BMI_baseline_feature_imp.tsv" ~ "UKB OLINK",
            label == "Results/ukb_bmi_assoc_metrics/proteins_BMI_baseline_feature_imp.tsv" ~ "UKB Predidcted SOMA",
            label == "Results/ukb_bmi_assoc_metrics/tier1onlyolinksoma_actualolink_overlapping_proteins_BMI_baseline_feature_imp.tsv" ~ "UKB OLINK (Tier 1 Only)",
            label == "Results/ukb_bmi_assoc_metrics/tier1onlyolinksoma_predsoma_overlapping_proteins_BMI_baseline_feature_imp.tsv" ~ "UKB Imputed SOMA (Tier 1)"
        ),
        model = "OLINK -> SOMA"
    )



### chs results
files <- list.files(path = "Results/chs_bmi_assoc_metrics", pattern = "feature_imp", full.names = TRUE)

chs_top10_all <- data.frame()

for (i in 1:length(files)){
    file <- files[i]
    df <- fread(file) %>%
        mutate(label = file) %>%
        arrange(desc(Importance)) %>%
        head(6)
    chs_top10_all <- rbind(chs_top10_all, df)
}

chs_top10_all <- chs_top10_all %>%
    mutate(
        label = case_when(
            label == "Results/chs_bmi_assoc_metrics/proteins_bmi_feature_imp.tsv" ~ "CHS Imputed OLINK",
            label == "Results/chs_bmi_assoc_metrics/benchsoma_overlapping_proteins_bmi_feature_imp.tsv" ~ "CHS SOMA",
            label == "Results/chs_bmi_assoc_metrics/tier1onlysomaolink_actualsoma_overlapping_proteins_bmi_feature_imp.tsv" ~ "CHS SOMA (Tier 1 Only)",
            label == "Results/chs_bmi_assoc_metrics/tier1onlysomaolink_predolink_overlapping_proteins_bmi_feature_imp.tsv" ~ "CHS Imputed OLINK (Tier 1)"
        ),
        model = "SOMA -> OLINK"
    )


##### 03 - plotting sankey ####
ukb_top10_all1 <- ukb_top10_all %>%
    group_by(label) %>%
    mutate(Importance_prop = Importance / sum(Importance)) %>%
    select(label, Protein, Importance, Importance_prop)

chs_top10_all1 <- chs_top10_all %>%
    group_by(label) %>%
    mutate(Importance_prop = Importance / sum(Importance)) %>%
    select(label, Protein, Importance, Importance_prop)

order <- c("CHS SOMA", "UKB OLINK",
           "CHS Imputed OLINK (Tier 1)", "CHS SOMA (Tier 1 Only)", "UKB Imputed SOMA (Tier 1)", "UKB OLINK (Tier 1 Only)")
order_labels <- c("CHS SomaScan", "UKB Olink",
            "CHS Imputed Olink (Tier 1)", "CHS SomaScan (Tier 1)", "UKB Imputed \nSomaScan (Tier 1)", "UKB Olink (Tier 1)")

df_base <- rbind(chs_top10_all1, ukb_top10_all1) %>%
    mutate(annot = paste0(Protein," (",round(Importance,3),")")) %>%
    filter(label %in% c("CHS SOMA", "UKB OLINK",
                        "CHS Imputed OLINK (Tier 1)", "CHS SOMA (Tier 1 Only)", "UKB Imputed SOMA (Tier 1)", "UKB OLINK (Tier 1 Only)")) %>%
    mutate(
        label = factor(label, levels = order, labels = order_labels),
        panel = ifelse(label %in% c("CHS SomaScan", "UKB Olink"), "Default Situation", "Tier Membership Informed")   
    )

protein_global_importance <- df_base %>%
    group_by(Protein) %>%
    summarise(mean_importance = mean(Importance), .groups = 'drop') %>%
    arrange(mean_importance)

df_alluvium <- df_base %>%
    group_by(Protein) %>%
    mutate(Protein_ID = cur_group_id()) %>%
    ungroup() %>%
    mutate(Protein = factor(Protein, levels = protein_global_importance$Protein))


## coloring
all_proteins <- unique(df_alluvium$Protein)
my_palette <- colorRampPalette(brewer.pal(9, "Set1"))(22)
protein_colors <- setNames(my_palette, all_proteins)


ggplot(df_alluvium,
       aes(x = label, y = Importance_prop,
           alluvium = Protein_ID, 
           stratum = Protein, 
           fill = Protein, 
           label = Protein)) +
    facet_nested(.~panel, scales="free_x", space = "free_x") + 
    geom_alluvium() + 
    geom_stratum(alpha = 0.5) +
    geom_text(stat = "stratum", size = 3.5, color = "black") +
    scale_fill_manual(values = protein_colors) +
    theme_bw(base_size = 12) +
    labs(
        x = "Models",
        y = "Importance"
    ) +
    theme(
        axis.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 12, color="black"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none")


ggsave("Results/figures/BMI_important_features.png", dpi = 300, width = 15, height = 6.5, bg = "white")







##### 04 - feature importance results representation #####
feature_imp_wide <- df_base %>%
    select(label, annot) %>%
    mutate(rank = row_number()) %>%
    pivot_wider(id_cols = c("rank"),
                names_from = "label",
                values_from = "annot"
    ) %>%
    select(-rank) %>%
    select(!!order)

write.table(feature_imp_wide, "Results/bmi_feature_importances.tsv", sep = "\t", row.names = F, quote = F)



