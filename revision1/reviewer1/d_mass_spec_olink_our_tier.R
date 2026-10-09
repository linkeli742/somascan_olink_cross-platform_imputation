##### 01 - library #####
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggh4x)
library(gridExtra)
library(dunn.test)
library(patchwork)
library(ggsignif)


gc()
rm(list=ls())

setwd("...")


##### 02 - files #####
# 2025 paper --------------------------------------------------------------
## this their supplementary table 8 (mass spec vs Olink)
ms_olink_corr <- readxl::read_xlsx("Data/mass_spec_olink_corr/42004_2025_1753_MOESM11_ESM.xlsx")
ms_olink_corr_short <- ms_olink_corr %>%
    select(Gene.Name, Assay, Correlation_Spearman, Correlation_Pearson)


## this is their supplementary table 14 (mass spec vs olink vs eldjarn)
ms_olink_eldjarn <- readxl::read_xlsx("Data/mass_spec_olink_corr/42004_2025_1753_MOESM17_ESM.xlsx")
ms_olink_eldjarn_short <- ms_olink_eldjarn %>%
    select(Gene.Name, Correlation_Sissala, Correlation_Eldjarn, confidence.tier_Eldjarn)


# our own data ------------------------------------------------------------
so_corr_ll <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F)
os_corr_ll <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = F)



##### 03 - visualize #####
## our tier data, tagged by direction
so_tiers <- so_corr_ll %>% select(protein, posttrain_r, tiers) %>% mutate(direction = "SomaScan to Olink")
os_tiers <- os_corr_ll %>% select(protein, posttrain_r, tiers) %>% mutate(direction = "Olink to SomaScan")
our_tiers_all <- bind_rows(so_tiers, os_tiers) %>%
    mutate(tiers = factor(tiers, levels = 1:4, labels = paste("Tier", 1:4)))


## merge with Sissala et al. MS-Olink correlations, matched by gene/protein name
comparison_df <- our_tiers_all %>%
    inner_join(ms_olink_corr_short, by = c("protein" = "Gene.Name"))

## how much overlap do we actually have?
## there are 821 overlapping proteins
comparison_df %>%
    group_by(direction) %>%
    summarise(n_matched = n_distinct(protein))


# plot A ------------------------------------------------------------------
comparison_long <- comparison_df %>%
    select(protein, tiers, direction, Correlation_Spearman, Correlation_Pearson) %>%
    mutate(direction=factor(direction, levels=c("SomaScan to Olink", "Olink to SomaScan"))) %>%
    pivot_longer(cols = c(Correlation_Spearman, Correlation_Pearson),
                 names_to = "metric", values_to = "correlation") %>%
    mutate(metric = recode(metric,
                           "Correlation_Spearman" = "Spearman",
                           "Correlation_Pearson" = "Pearson"))


## compute KW p-value per direction x metric combination
kw_pvals <- comparison_long %>%
    group_by(direction, metric) %>%
    summarise(
        kw_p = kruskal.test(correlation ~ tiers)$p.value,
        .groups = "drop"
    ) %>%
    mutate(
        kw_label = ifelse(kw_p < 0.001, "Kruskal-Wallis p < 0.001",
                          paste0("Kruskal-Wallis p = ", formatC(kw_p, digits = 3, format = "f")))
    )


## small helper to format p-values consistently (scientific below 0.001)
format_p <- function(p) {
    ifelse(p < 0.001,
           formatC(p, format = "e", digits = 2),
           formatC(p, format = "f", digits = 3))
}


## helper: run dunn.test per direction x metric, extract ADJACENT tier comparisons,
## keeping BOTH raw and bonferroni-adjusted p-values
get_adjacent_dunn <- function(sub_df) {
    dunn_res <- dunn.test(sub_df$correlation, sub_df$tiers, method = "bonferroni", kw = FALSE, table = FALSE)
    dunn_df <- tibble(comparison = dunn_res$comparisons,
                      p_raw = dunn_res$P,
                      p_adj = dunn_res$P.adjusted)
    
    ## keep only adjacent pairs: Tier1-Tier2, Tier2-Tier3, Tier3-Tier4
    adjacent_pairs <- c("Tier 1 - Tier 2", "Tier 2 - Tier 3", "Tier 3 - Tier 4")
    dunn_df %>% filter(comparison %in% adjacent_pairs)
}


## build significance annotations per direction x metric facet
sig_annotations <- comparison_long %>%
    group_by(direction, metric) %>%
    group_modify(~ get_adjacent_dunn(.x)) %>%
    ungroup() %>%
    mutate(
        label = paste0("Raw P=", toupper(format_p(p_raw)), "\nBonf adj P=", toupper(format_p(p_adj))),
        group1 = sub(" - .*", "", comparison),
        group2 = sub(".* - ", "", comparison)
    )

## there's no buildin facet wrap, therefore just put 4 plots together
build_subplot <- function(dir, met) {
    sub_data <- comparison_long %>% filter(direction == dir, metric == met)
    sub_sig <- sig_annotations %>% filter(direction == dir, metric == met)
    
    ## step room for the annotation bars
    step <- (y_max - min(sub_data$correlation, na.rm = TRUE)) * 0.1

    y_max <- max(sub_data$correlation, na.rm = TRUE)
    y_min <- min(sub_data$correlation, na.rm = TRUE)
    annot_top <- y_max + step * 3.4
    y_upper   <- annot_top + step * 1.3

    
    ggplot(sub_data, aes(x = tiers, y = correlation, fill = tiers)) +
        geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
        geom_jitter(width = 0.15, alpha = 0.3, size = 0.6) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
        geom_signif(
            comparisons = list(c("Tier 1", "Tier 2"), c("Tier 2", "Tier 3"), c("Tier 3", "Tier 4")),
            annotations = sub_sig$label,
            y_position = c(y_max + step, y_max + step * 1.7, y_max + step * 2.5),
            tip_length = 0.01, textsize = 2.6, vjust = -0.1
        ) +
        coord_cartesian(ylim = c(y_min, y_upper), clip = "off") +
        scale_fill_manual(values = c("Tier 1" = "#6beeff", "Tier 2" = "#93db86",
                                     "Tier 3" = "#ffb963", "Tier 4" = "#ff7878")) +
        labs(x = NULL, y = glue::glue("Independent MS-Olink correlation ({met})\n(Sissala et al. 2025)"),
             title = dir) +
        theme_classic(base_size = 12) +
        guides(fill = "none")
}


## reorder combos so metric is the outer grouping: Pearson row first, Spearman row second
combos <- comparison_long %>%
    distinct(direction, metric) %>%
    mutate(metric = factor(metric, levels = c("Pearson", "Spearman"))) %>%
    arrange(metric, direction)

subplots <- purrr::pmap(combos, function(direction, metric) build_subplot(direction, metric))
names(subplots) <- paste(combos$metric, combos$direction, sep = "_")

## split into the two metric rows, each row = one composite subpanel
pearson_row  <- wrap_plots(subplots[grepl("^Pearson_",  names(subplots))], ncol = 2)
spearman_row <- wrap_plots(subplots[grepl("^Spearman_", names(subplots))], ncol = 2)


p_tier_vs_ms_combined <- (wrap_elements(full = pearson_row) / wrap_elements(full = spearman_row)) +
    plot_annotation(tag_levels = "A")

p_tier_vs_ms_combined


##### 04 - save #####
ggsave(
    filename = "Results/revision1/figures/our_tier_vs_msolink_sissala.png",
    plot = p_tier_vs_ms_combined,
    width = 10, height = 9, dpi = 300
)

