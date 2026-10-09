##### 01 - library #####
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggh4x)
library(purrr)
library(dunn.test)
library(ggsignif)
library(patchwork)
library(ggpubr)

gc()
rm(list = ls())
setwd("...")


##### 02 - files #####
alphafold3_full <- fread("Results/alphafold3/qiayi/batch_results_prediction_quality.tsv", data.table = F)
so_corr_ll <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F)
os_corr_ll <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = F)

##### 03 - build dual-direction tier data #####
## AF3's own "tier" column IS the SO tier -- drop it and pull fresh from so_corr_ll
## for consistency/traceability, then attach OS tier separately
af3_base <- alphafold3_full %>%
    select(-tier)

af3_so <- af3_base %>%
    inner_join(so_corr_ll %>% select(protein, tiers), by = "protein") %>%
    mutate(direction = "SomaScan to Olink")

af3_os <- af3_base %>%
    inner_join(os_corr_ll %>% select(protein, tiers), by = "protein") %>%
    mutate(direction = "Olink to SomaScan")

af3_dual <- bind_rows(af3_so, af3_os) %>%
    mutate(
        tier_label = paste("Tier", tiers),
        tier_label = factor(tier_label, levels = paste("Tier", 1:4)),
        direction = factor(direction, levels = c("SomaScan to Olink", "Olink to SomaScan"))
    )

## quick check: how many proteins have data in both directions vs. only SO
af3_dual %>% count(direction)


##### 04 - long format & pairwise comparison setup #####
continuous_metrics_cols <- c("fraction_disordered", "avg_plddt", "ptm", "ranking_score")
continuous_metrics <- c("Fraction Disordered", "Average pLDDT", "pTM", "Ranking Score")

af3_long <- af3_dual %>%
    select(protein, tier_label, direction, all_of(continuous_metrics_cols)) %>%
    pivot_longer(cols = all_of(continuous_metrics_cols), names_to = "metric", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(metric = factor(
        metric,
        levels = c("fraction_disordered", "avg_plddt", "ptm", "ranking_score"),
        labels = c("Fraction Disordered", "Average pLDDT", "pTM", "Ranking Score")
    ))

## exact desired comparison order, left-to-right on x-axis
ordered_pairs <- list(
    c("Tier 1", "Tier 2"),
    c("Tier 2", "Tier 3"),
    c("Tier 3", "Tier 4"),
    c("Tier 1", "Tier 3"),
    c("Tier 2", "Tier 4"),
    c("Tier 1", "Tier 4")
)
comparison_levels <- sapply(ordered_pairs, paste, collapse = " vs ")


##### 05 - pairwise wilcoxon tests: build & save full results table #####

compute_pairwise_table <- function(sub_data, metric_name, dir_val) {
    pvals <- sapply(ordered_pairs, function(pr) {
        g1 <- sub_data$value[sub_data$tier_label == pr[1]]
        g2 <- sub_data$value[sub_data$tier_label == pr[2]]
        if (length(g1) < 2 || length(g2) < 2) return(NA_real_)
        wilcox.test(g1, g2)$p.value
    })
    data.frame(
        metric = metric_name,
        direction = dir_val,
        comparison = comparison_levels,
        p_raw = pvals,
        p_adj = p.adjust(pvals, method = "bonferroni"),
        stringsAsFactors = FALSE
    )
}

pairwise_results <- af3_long %>%
    group_by(metric, direction) %>%
    group_split() %>%
    purrr::map_dfr(~ compute_pairwise_table(.x, as.character(unique(.x$metric)), as.character(unique(.x$direction))))

pairwise_results <- pairwise_results %>%
    mutate(
        sig_label = case_when(
            is.na(p_adj) ~ "n/a",
            p_adj < 0.001 ~ "***",
            p_adj < 0.01 ~ "**",
            p_adj < 0.05 ~ "*",
            TRUE ~ "ns"
        ),
        comparison = factor(comparison, levels = comparison_levels),
        metric = factor(metric, levels = continuous_metrics),
        direction = factor(direction, levels = c("SomaScan to Olink", "Olink to SomaScan")),
        neg_log10_p = -log10(p_adj)
    )

write.table(pairwise_results, "Results/revision1/alphafold3/af3_tier_pairwise_wilcox.tsv", sep = "\t", row.names = F, quote = F)


##### 06 - heatmap (only figure) #####
p_af3_heatmap <- ggplot(pairwise_results, aes(x = comparison, y = metric, fill = neg_log10_p)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sig_label), size = 4) +
    facet_wrap(~direction) +
    scale_fill_gradient(low = "white", high = "red", na.value = "grey90", name = "-log10(adj. p)") +
    labs(x = "Tier comparison", y = "Metric") +
    theme_classic2(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 15, hjust = 1),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold")
    )

p_af3_heatmap

ggsave("Results/revision1/figures/alphafold3_tier_pairwise_heatmap.png", p_af3_heatmap,
       dpi = 300, width = 12, height = 6, bg = "white")


