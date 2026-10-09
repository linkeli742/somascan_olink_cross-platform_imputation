##### 01 - library #####
library(dplyr)
library(data.table)
library(ggplot2)
library(colorspace)
library(ggrepel)
library(ggh4x)

gc()
rm(list = ls())

setwd("...")



##### 02 - files #####
tot_pubs <- fread("Results/protein_importance/MALACARDS_proteins_sumofpubs.tsv.gz", data.table = F)
soma_olink <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F)
soma_olink <- soma_olink %>%
    select(-index) %>%
    mutate(tiers = as.character(tiers))

## importing protein tags
protein_hpa_tags <- fread("Results/protein_importance/HPA_importance_protclasstags.tsv", data.table = F) 

# separating drug tags
fda_drug_targets <- unique(protein_hpa_tags[protein_hpa_tags$`FDA approved drug targets` == 1,]$protein)
potential_drug_targets <- unique(protein_hpa_tags[protein_hpa_tags$`Potential drug targets` == 1,]$protein)

# disease gene tags
cancer_genes <- unique(protein_hpa_tags[protein_hpa_tags$`Cancer-related genes` == 1,]$protein)
disease_genes <- unique(protein_hpa_tags[protein_hpa_tags$`Disease related genes` == 1,]$protein)


## full dataset
soma_olink_full <- left_join(soma_olink, tot_pubs, by = "protein") %>%
    ## appending tags
    mutate(
        drug_group = case_when(
            protein %in% fda_drug_targets ~ "FDA Drug Targets",
            protein %in% potential_drug_targets ~ "Potential Drug Targets",
            TRUE ~ "None"
        ),
        gene_group = case_when(
            protein %in% disease_genes ~ "Disease related genes",
            protein %in% cancer_genes ~ "Cancer related genes",
            TRUE ~ "N/A"
        ),
        ## factoring
        drug_group = factor(drug_group, levels = c("FDA Drug Targets", "Potential Drug Targets", "None")),
        gene_group = factor(gene_group, levels = c("Cancer related genes", "Disease related genes", "N/A"))
    )


## ranking by gain, color by pubs
plot_data <- soma_olink_full %>%
    arrange(desc(gain)) %>%
    mutate(index = row_number(gain))


## what's the the proportion of FDA drug targets
plot_data %>%
    filter(gain > 0) %>%
    summarize(n = n(),
              n_fda = sum(drug_group=="FDA Drug Targets"),
              prop_fda = n_fda / n,
              n_pfda = sum(drug_group=="Potential Drug Targets"),
              prop_pfda = n_pfda / n,
              n_cancer = sum(gene_group == "Cancer related genes"),
              prop_cancer = n_cancer / n,
              n_disease = sum(gene_group == "Disease related genes"),
              prop_disease = n_disease / n
    )

## by tiers
summary <- plot_data %>%
    filter(gain > 0) %>%
    group_by(tiers) %>%
    summarise(n = n(),
              n_fda = sum(drug_group=="FDA Drug Targets"),
              prop_fda = n_fda / n,
              n_pfda = sum(drug_group=="Potential Drug Targets"),
              prop_pfda = n_pfda / n,
              n_cancer = sum(gene_group == "Cancer related genes"),
              prop_cancer = n_cancer / n,
              n_disease = sum(gene_group == "Disease related genes"),
              prop_disease = n_disease / n
    )
summary
# 1. Test if the proportion of 'FDA Drug Targets' (prop_fda) is the same across tiers
fda_test <- prop.test(x = summary$n_fda, n = summary$n)
print(fda_test) # 0.20

# 2. Test if the proportion of 'Potential Drug Targets' (prop_pfda) is the same across tiers
pfda_test <- prop.test(x = summary$n_pfda, n = summary$n)
print(pfda_test) # 0.06

# 3. Test if the proportion of 'Cancer related genes' (prop_cancer) is the same across tiers
cancer_test <- prop.test(x = summary$n_cancer, n = summary$n)
print(cancer_test) # 0.46

# 4. Test if the proportion of 'Disease related genes' (prop_disease) is the same across tiers
disease_test <- prop.test(x = summary$n_disease, n = summary$n)
print(disease_test)  # 0.005




##### I think bins looks better
tier_cols <- c("Tier 1"="#fc9c26", "Tier 2"="#74c365", "Tier 3"="#ff5959", "Tier 4"="#26e3fc")





# color by quadrant -------------------------------------------------------
global_max_totpub <- max(soma_olink_full$tot_pubs)
plot_data <- soma_olink_full %>%
    group_by(tiers) %>%
    arrange(desc(gain)) %>%
    mutate(
        index = row_number(gain),
        tiers = case_when(
            tiers == "1" ~ "Tier 1",
            tiers == "2" ~ "Tier 2",
            tiers == "3" ~ "Tier 3",
            tiers == "4" ~ "Tier 4",
            TRUE ~ NA # there should be any but this is my habit
        )
    )


# top 6 FDA approval drug for annotation
annot_set1 <- plot_data %>%
    group_by(tiers) %>%
    slice_max(order_by = tot_pubs, n = 6)

annot_set2 <- plot_data %>%
    filter(protein %in% c("MAPT", "TNNI3"))

annot <- rbind(annot_set1, annot_set2)


### ggplot time
ggplot(data = plot_data, aes(x = gain, y = tot_pubs, color = gene_group, shape = drug_group, alpha = gene_group)) +
    geom_point() +
    facet_wrap(~tiers,
                # scales = "free_x",
               axes = "all", 
               ncol = 1
                 ) + 
    # ggh4x::scale_x_facet(tiers %in% c("Tier 1", "Tier 2"), limits = c(-0.05,0.4)) +
    # ggh4x::scale_x_facet(tiers == "Tier 4", limits = c(-0.2, 0.2)) + 
    scale_x_continuous(limits = c(-0.2, 0.45), breaks = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4)) +
    geom_label_repel(data = annot, aes(x = gain, y = tot_pubs, label = protein, color = gene_group), inherit.aes = F,
                    size = 3, max.overlaps = Inf, fontface = "bold", box.padding = 0.4, min.segment.length = 0, show.legend = F) +
    ## shape based on drug targets
    scale_shape_manual(values = c("None" = 1, "FDA Drug Targets" = 8, "Potential Drug Targets" = 17)) +
    ## color based on gene groups
    scale_color_manual(values = c("Disease related genes" = "firebrick3", 
                                  "Cancer related genes" = "royalblue3", 
                                  "N/A" = "grey6")) +
    ## if it's grey points then it's lighter
    scale_alpha_manual(values = c("Disease related genes" = 1, 
                                  "Cancer related genes" = 1, 
                                  "N/A" = 0.5), guide = "none", ) +
    ## labels
    labs(
        shape = "Drug Targets",
        color = "Gene Group",
        x = "Gain",
        y = "Sum of Publications"
        # title = "SOMA - OLINK Protein Categories"
    ) +
    theme_minimal()


ggsave("Results/figures/soma_olink_quadrant_scatter_proteinclass.png", dpi = 300, width = 11, height = 7, bg = "white")


