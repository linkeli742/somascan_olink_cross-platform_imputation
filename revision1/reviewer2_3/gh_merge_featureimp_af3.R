##### 01 - library #####
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggh4x)
library(purrr)
library(patchwork)
library(ggpubr)

gc()
rm(list = ls())
setwd("...")


##### 02 - PANEL A: feature importance concentration by overlap status #####
olink_prots <- names(fread("Data/final_datasets/MESA_Olink_3k_log2_nodup_rm10na_knn_final.tsv.gz", nrows = 0))
soma_prots  <- names(fread("Data/final_datasets/MESA_SOMA_5k_log2scaled_nodup_final.tsv.gz", nrows = 0))
olink_prots <- setdiff(olink_prots, "id")
soma_prots  <- setdiff(soma_prots, c("id", "bmi"))

summarize_importance_file <- function(path, target_protein) {
    imp <- fread(path, data.table = FALSE)
    if (nrow(imp) == 0) return(NULL)
    
    imp <- imp %>% arrange(desc(importance_mean))
    total <- sum(imp$importance_mean)
    if (total <= 0) return(NULL)
    
    p <- imp$importance_mean / total
    
    tibble(
        protein = target_protein,
        n_features_scored = nrow(imp),
        top1_feature = imp$feature[1],
        top1_is_self = imp$feature[1] == target_protein,
        self_in_top5 = target_protein %in% head(imp$feature, 5),
        self_in_top10 = target_protein %in% head(imp$feature, 10),
        top5_share = sum(head(p, 5)),
        top10_share = sum(head(p, 10)),
        top20_share = sum(head(p, 20))
    )
}

read_direction <- function(dir_path, suffix = "_feature_imp_shap.tsv") {
    files <- list.files(dir_path, pattern = paste0(suffix, "$"), full.names = TRUE)
    proteins <- sub(suffix, "", basename(files), fixed = TRUE)
    map2_dfr(files, proteins, summarize_importance_file)
}

imp_os <- read_direction("Results/feature_imp/olink_soma") %>%
    mutate(direction = "Olink to SomaScan", overlapping = protein %in% olink_prots)

imp_so <- read_direction("Results/feature_imp/soma_olink") %>%
    mutate(direction = "SomaScan to Olink", overlapping = protein %in% soma_prots)

imp_all <- bind_rows(imp_os, imp_so) %>%
    mutate(overlap_group = if_else(overlapping,
                                   "Overlapping (has counterpart)",
                                   "Purely imputed (no counterpart)"))

## group comparison tests (needed for the p-value annotation only; not written to disk)
run_group_tests <- function(df, direction_label) {
    d <- df %>% filter(direction == direction_label)
    wilcox_top5 <- wilcox.test(top5_share ~ overlapping, data = d)
    tibble(
        direction = direction_label,
        test = "top5_share (Wilcoxon)",
        p_value = wilcox_top5$p.value
    )
}

test_results <- bind_rows(
    run_group_tests(imp_all, "Olink to SomaScan"),
    run_group_tests(imp_all, "SomaScan to Olink")
)

label_stats <- imp_all %>%
    filter(!is.na(top5_share)) %>%
    group_by(direction, overlap_group) %>%
    summarise(
        mean_val = mean(top5_share, na.rm = TRUE),
        median_val = median(top5_share, na.rm = TRUE),
        max_val = max(top5_share, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(label = sprintf("Mean: %.3f\nMedian: %.3f", mean_val, median_val))

format_p <- function(p) {
    if (p < 0.0001) format(p, digits = 3, scientific = TRUE) else sprintf("%.4f", p)
}

p_annot <- test_results %>%
    filter(test == "top5_share (Wilcoxon)") %>%
    mutate(label = paste0("p = ", sapply(p_value, format_p)))

y_top     <- max(label_stats$max_val, na.rm = TRUE)
y_bracket <- y_top * 1.30
y_label   <- y_top * 1.40

imp_all      <- imp_all      %>% mutate(direction = factor(direction, levels = c("SomaScan to Olink", "Olink to SomaScan")))
label_stats  <- label_stats  %>% mutate(direction = factor(direction, levels = c("SomaScan to Olink", "Olink to SomaScan")))
p_annot      <- p_annot      %>% mutate(direction = factor(direction, levels = c("SomaScan to Olink", "Olink to SomaScan")))

p_concentration <- ggplot(imp_all, aes(x = overlap_group, y = top5_share, fill = overlap_group)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2.2,
                 fill = "red", color = "black") +
    geom_text(data = label_stats, aes(x = overlap_group, y = max_val * 1.12, label = label),
              inherit.aes = FALSE, size = 2.7, lineheight = 0.85, vjust = 0) +
    facet_wrap(~ direction) +
    geom_segment(data = p_annot, aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket),
                 inherit.aes = FALSE, linewidth = 0.3) +
    geom_text(data = p_annot, aes(x = 1.5, y = y_label, label = label),
              inherit.aes = FALSE, size = 3.2) +
    scale_fill_manual(values = c("Overlapping (has counterpart)" = "#6beeff",
                                 "Purely imputed (no counterpart)" = "#ff7878")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    labs(x = NULL, y = "Share of total SHAP importance\nheld by top 5 predictors") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none")


##### 03 - PANEL B: AlphaFold3 tier pairwise heatmap #####
alphafold3_full <- fread("Results/alphafold3/qiayi/batch_results_prediction_quality.tsv", data.table = F)
so_corr_ll <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F)
os_corr_ll <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = F)

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

ordered_pairs <- list(
    c("Tier 1", "Tier 2"),
    c("Tier 2", "Tier 3"),
    c("Tier 3", "Tier 4"),
    c("Tier 1", "Tier 3"),
    c("Tier 2", "Tier 4"),
    c("Tier 1", "Tier 4")
)
comparison_levels <- sapply(ordered_pairs, paste, collapse = " vs ")

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


##### 04 - merge panels: A = feature importance, B = AlphaFold3 #####
p_merged <- p_concentration / p_af3_heatmap +
    plot_annotation(tag_levels = "A") +
    plot_layout(heights = c(1, 1)) &
    theme(plot.tag = element_text(face = "bold"))

p_merged

ggsave("Results/revision1/figures/merged_feature_importance_alphafold3.png",
       p_merged, dpi = 300, width = 12, height = 11, bg = "white")