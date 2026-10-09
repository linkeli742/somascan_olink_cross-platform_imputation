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


##### 02 - file mapping #####
## feature importance files >> label
feature_imp_map <- tribble(
    ~file, ~label,
    "Results/revision1/ukb_assoc_metrics/actualolink_overlapping_proteins_BMI_baseline_feature_imp.tsv", "UKB OLINK (Full Panel)",
    "Results/revision1/chs_assoc_metrics/somaolink_actualsoma_overlapping_proteins_bmi_feature_imp.tsv", "CHS SOMA (Full Panel)",
    
    "Results/revision1/ukb_assoc_metrics/tier1onlyolinksoma_actualolink_overlapping_proteins_BMI_baseline_feature_imp.tsv", "UKB OLINK (Tier 1)",
    "Results/revision1/ukb_assoc_metrics/tier1onlyolinksoma_predsoma_overlapping_proteins_BMI_baseline_feature_imp.tsv",    "UKB Imputed SOMA (Tier 1)",
    "Results/revision1/chs_assoc_metrics/tier1onlysomaolink_actualsoma_overlapping_proteins_bmi_feature_imp.tsv",           "CHS SOMA (Tier 1)",
    "Results/revision1/chs_assoc_metrics/tier1onlysomaolink_predolink_overlapping_proteins_bmi_feature_imp.tsv",            "CHS Imputed OLINK (Tier 1)"
)

## R2 files >> label (same label set, mapped to the corresponding _baseline_rsq.tsv file)
rsq_map <- tribble(
    ~file, ~label,
    "Results/revision1/ukb_assoc_metrics/actualolink_overlapping_proteins_BMI_baseline_baseline_rsq.tsv", "UKB OLINK (Full Panel)",
    "Results/revision1/chs_assoc_metrics/somaolink_actualsoma_overlapping_proteins_bmi_baseline_rsq.tsv", "CHS SOMA (Full Panel)",
    
    "Results/revision1/ukb_assoc_metrics/tier1onlyolinksoma_actualolink_overlapping_proteins_BMI_baseline_baseline_rsq.tsv", "UKB OLINK (Tier 1)",
    "Results/revision1/ukb_assoc_metrics/tier1onlyolinksoma_predsoma_overlapping_proteins_BMI_baseline_baseline_rsq.tsv",    "UKB Imputed SOMA (Tier 1)",
    "Results/revision1/chs_assoc_metrics/tier1onlysomaolink_actualsoma_overlapping_proteins_bmi_baseline_rsq.tsv",           "CHS SOMA (Tier 1)",
    "Results/revision1/chs_assoc_metrics/tier1onlysomaolink_predolink_overlapping_proteins_bmi_baseline_rsq.tsv",            "CHS Imputed OLINK (Tier 1)"
)

order <- c("CHS SOMA (Full Panel)", "UKB OLINK (Full Panel)",
           "CHS SOMA (Tier 1)", "CHS Imputed OLINK (Tier 1)", 
           "UKB Imputed SOMA (Tier 1)", "UKB OLINK (Tier 1)")
order_labels <- c("CHS SomaScan\n(Full Panel)", "UKB Olink\n(Full Panel)",
                  "CHS SomaScan\n(Tier 1)", "CHS Imputed Olink\n(Tier 1)", 
                  "UKB Imputed SomaScan\n(Tier 1)", "UKB Olink\n(Tier 1)")

##### 03 - load feature importance (top 6 per model) #####
load_feature_imp <- function(file, label) {
    full_df <- fread(file)
   
    full_df %>%
        arrange(desc(Importance)) %>%
        head(5) %>%
        mutate(label = label)
}

feature_imp_all <- pmap_dfr(feature_imp_map, load_feature_imp)

##### 04 - load R2 metrics #####
load_rsq <- function(file, label) {
    fread(file) %>%
        mutate(label = label)
}

rsq_all <- pmap_dfr(rsq_map, load_rsq) %>%
    mutate(
        r2_annotation = paste0("R\u00B2 = ", round(rsq, 3), "\nAdj. R\u00B2 = ", round(adjusted_rsq, 3))
    ) %>%
    select(label, rsq, adjusted_rsq, cv_r2_training_fold, r2_annotation)



##### 05 - build plotting dataframe #####
df_base <- feature_imp_all %>%
    left_join(rsq_all %>% select(label, rsq, adjusted_rsq), by = "label") %>%
    group_by(label) %>%
    mutate(
        Importance_prop = (Importance / sum(Importance))
    ) %>%
    ungroup() %>%
    mutate(
        label = factor(label, levels = order, labels = order_labels),
        panel = ifelse(label %in% c("CHS SomaScan\n(Full Panel)", "UKB Olink\n(Full Panel)"),
                       "Default Situation", "Tier Membership Informed")
    )

all_proteins_levels <- sort(unique(feature_imp_all$Protein))

df_alluvium <- df_base %>%
    mutate(Protein = factor(Protein, levels = all_proteins_levels)) %>%
    group_by(Protein) %>%
    mutate(Protein_ID = cur_group_id()) %>%
    ungroup()

rsq_annot <- rsq_all %>%
    mutate(
        label = factor(label, levels = order, labels = order_labels),
        panel = ifelse(label %in% c("CHS SomaScan\n(Full Panel)", "UKB Olink\n(Full Panel)"),
                       "Default Situation", "Tier Membership Informed")
    )


##### 06 - coloring #####
all_proteins <- unique(df_alluvium$Protein)
my_palette <- colorRampPalette(brewer.pal(9, "Set1"))(length(all_proteins))
protein_colors <- setNames(my_palette, all_proteins)

##### 07 - plot #####
max_stack_height <- df_base %>%
    group_by(label) %>%
    summarise(total = sum(Importance_prop), .groups = "drop") %>%
    pull(total) %>%
    max()

p_sankey <- ggplot(df_alluvium,
                   aes(x = label, y = Importance_prop,
                       alluvium = Protein_ID,
                       stratum = Protein,
                       fill = Protein,
                       label = Protein)) +
    facet_nested(. ~ panel, scales = "free_x", space = "free_x") +
    geom_alluvium(decreasing = TRUE) +
    geom_stratum(alpha = 0.5, decreasing = TRUE) +
    geom_text(stat = "stratum", decreasing = TRUE, size = 3.5, color = "black") +
    geom_text(
        data = rsq_annot,
        aes(x = label, y = max_stack_height * 1.05, label = r2_annotation),
            inherit.aes = FALSE, size = 3.9, 
            # fontface = "italic", 
            lineheight = 0.9
    ) +
    scale_fill_manual(values = protein_colors) +
    scale_y_continuous(limits = c(0, max_stack_height * 1.10), expand = c(0, 0)) +
    theme_bw(base_size = 12) +
    labs(x = "Models", y = "Contribution to R\u00B2 (top 5 proteins)") +
    theme(
        axis.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none"
    )

p_sankey


ggsave("Results/revision1/bmi_rerun/bmi_feature_imp_sankey.png", p_sankey, dpi=300, bg="white", width = 13, heigh=7)


##### 08 - save combined R2 table (Reviewer 1 leakage fix + Reviewer 2 adjusted R2) #####
rsq_table_out <- rsq_all %>%
    select(label, rsq, adjusted_rsq)

write.table(rsq_table_out, "Results/revision1/bmi_rerun/bmi_r2_comparison_table.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)



##### 09 - feature importance wide table (unchanged logic, updated labels) #####
feature_imp_wide <- df_base %>%
    mutate(annot = paste0(Protein, " (", round(Importance, 3), ")")) %>%
    mutate(label = str_replace_all(label, "\n", " ")) %>%
    select(label, annot) %>%
    group_by(label) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    pivot_wider(names_from = label, values_from = annot) %>%
    select(-rank)

feature_imp_wide

write.table(feature_imp_wide, "Results/revision1/bmi_rerun/bmi_feature_importances.tsv", sep = "\t", row.names = FALSE, quote = FALSE)





