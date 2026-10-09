##### 01 - library #####
library(dplyr)
library(tidyr)
library(data.table)
library(purrr)
library(ggplot2)
library(ggh4x)

gc()
rm(list = ls())

setwd("...")


##### 02 - overlap status:  #####
olink_prots <- names(fread("Data/final_datasets/MESA_Olink_3k_log2_nodup_rm10na_knn_final.tsv.gz", nrows = 0))
soma_prots  <- names(fread("Data/final_datasets/MESA_SOMA_5k_log2scaled_nodup_final.tsv.gz", nrows = 0))
olink_prots <- setdiff(olink_prots, "id")
soma_prots  <- setdiff(soma_prots, c("id", "bmi"))



##### 03 - read per-protein SHAP importance files and summarize (one row per protein) #####
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

## target = SomaScan protein, predictors = Olink panel
imp_os <- read_direction("Results/feature_imp/olink_soma") %>%
    mutate(direction = "Olink to SomaScan", overlapping = protein %in% olink_prots)

## target = Olink protein, predictors = SomaScan panel
imp_so <- read_direction("Results/feature_imp/soma_olink") %>%
    mutate(direction = "SomaScan to Olink", overlapping = protein %in% soma_prots)

imp_all <- bind_rows(imp_os, imp_so) %>%
    mutate(overlap_group = if_else(overlapping,
                                   "Overlapping (has counterpart)",
                                   "Purely imputed (no counterpart)"))

cat("=== proteins summarized per direction ===\n")
print(table(imp_all$direction, imp_all$overlap_group))

write.table(imp_all, "Results/revision1/feature_importance_by_overlap/per_protein_importance_summary.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)


##### 04 - group-level summary: overlapping vs. purely-imputed, per direction #####
group_summary <- imp_all %>%
    group_by(direction, overlap_group) %>%
    summarise(
        n_proteins = n(),
        pct_top1_is_self = mean(top1_is_self, na.rm = TRUE) * 100,
        median_top5_share = median(top5_share, na.rm = TRUE),
        median_top10_share = median(top10_share, na.rm = TRUE),
        median_top20_share = median(top20_share, na.rm = TRUE),
        .groups = "drop"
    )

group_summary

write.table(group_summary, "Results/revision1/feature_importance_by_overlap/group_summary.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)


##### 05 - statistical tests: do the two groups actually differ? (per direction) #####
run_group_tests <- function(df, direction_label) {
    d <- df %>% filter(direction == direction_label)
    
    wilcox_top5  <- wilcox.test(top5_share ~ overlapping, data = d)
    wilcox_top10 <- wilcox.test(top10_share ~ overlapping, data = d)
    wilcox_top20 <- wilcox.test(top20_share ~ overlapping, data = d)
    fisher_top1  <- fisher.test(table(d$overlapping, d$top1_is_self))
    
    tibble(
        direction = direction_label,
        test = c("top5_share (Wilcoxon)", "top10_share (Wilcoxon)", "top20_share (Wilcoxon)", "top1_is_self (Fisher)"),
        p_value = c(wilcox_top5$p.value, wilcox_top10$p.value, wilcox_top20$p.value, fisher_top1$p.value)
    )
}

test_results <- bind_rows(
    run_group_tests(imp_all, "Olink to SomaScan"),
    run_group_tests(imp_all, "SomaScan to Olink")
)

cat("=== Overlapping vs. purely-imputed: group comparison tests ===\n")
print(test_results)

write.table(test_results, "Results/revision1/feature_importance_by_overlap/group_comparison_tests.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)



##### 06 - stratify by tier: is the overlap effect just a prediction-quality effect in disguise? #####
so_corr_ll <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = FALSE) %>%
    select(protein, tiers)
os_corr_ll <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = FALSE) %>%
    select(protein, tiers)

imp_all <- imp_all %>%
    left_join(
        bind_rows(
            os_corr_ll %>% mutate(direction = "Olink to SomaScan"),
            so_corr_ll %>% mutate(direction = "SomaScan to Olink")
        ),
        by = c("protein", "direction")
    )

cat("=== proteins missing a tier assignment (dropped from tier-stratified summary) ===\n")
print(sum(is.na(imp_all$tiers)))

tier_overlap_summary <- imp_all %>%
    filter(!is.na(tiers)) %>%
    group_by(direction, tiers, overlap_group) %>%
    summarise(
        n_proteins = n(),
        pct_top1_is_self = mean(top1_is_self, na.rm = TRUE) * 100,
        median_top5_share = median(top5_share, na.rm = TRUE),
        median_top10_share = median(top10_share, na.rm = TRUE),
        median_top20_share = median(top20_share, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(direction, tiers, overlap_group)

tier_overlap_summary

write.table(tier_overlap_summary, "Results/revision1/feature_importance_by_overlap/tier_overlap_summary.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)


##### 06c - does self-reliance differ across tiers? (overlapping proteins only) #####
## Restricted to overlapping proteins throughout: tiers comes from pretrain_r, which is only
## defined for proteins measured on both platforms, so purely-imputed proteins have no tier
## to stratify by in the first place (see tier_overlap_summary above -- it's all "Overlapping").
run_tier_chisq <- function(df, direction_label, self_col) {
    d <- df %>% filter(direction == direction_label, overlapping, !is.na(tiers))
    tbl <- table(d$tiers, d[[self_col]])
    test <- chisq.test(tbl)
    tibble(
        direction = direction_label,
        metric = self_col,
        chisq_stat = unname(test$statistic),
        df = unname(test$parameter),
        p_value = test$p.value
    )
}


tier_self_chisq_results <- bind_rows(
    run_tier_chisq(imp_all, "Olink to SomaScan", "top1_is_self"),
    run_tier_chisq(imp_all, "SomaScan to Olink", "top1_is_self"),
    run_tier_chisq(imp_all, "Olink to SomaScan", "self_in_top5"),
    run_tier_chisq(imp_all, "SomaScan to Olink", "self_in_top5"),
    run_tier_chisq(imp_all, "Olink to SomaScan", "self_in_top10"),
    run_tier_chisq(imp_all, "SomaScan to Olink", "self_in_top10")
)

cat("=== self-reliance (top1 / top5 / top10) across tiers, overlapping proteins only ===\n")
print(tier_self_chisq_results)


write.table(tier_self_chisq_results, "Results/revision1/feature_importance_by_overlap/tier_self_reliance_chisq.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)




##### 07 - figure: importance concentration by overlap group, per direction #####
## per-group mean/median labels, positioned above each box individually
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

## actual (unrounded to threshold) p-value, formatted for readability
format_p <- function(p) {
    if (p < 0.0001) format(p, digits = 3, scientific = TRUE) else sprintf("%.4f", p)
}

p_annot <- test_results %>%
    filter(test == "top5_share (Wilcoxon)") %>%
    mutate(label = paste0("p = ", sapply(p_value, format_p)))

y_top     <- max(label_stats$max_val, na.rm = TRUE)
y_bracket <- y_top * 1.30
y_label   <- y_top * 1.40

imp_all <- imp_all %>%
    mutate(direction = factor(direction, levels=c("SomaScan to Olink", "Olink to SomaScan")))
label_stats <- label_stats %>%
    mutate(direction = factor(direction, levels=c("SomaScan to Olink", "Olink to SomaScan")))
p_annot <- p_annot %>%
    mutate(direction = factor(direction, levels=c("SomaScan to Olink", "Olink to SomaScan")))

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

p_concentration

ggsave("Results/revision1/figures/feature_importance_concentration_by_overlap_status.png",
       p_concentration, dpi = 300, width = 10, height = 5, bg = "white")





##### 08 - detailed top5 predictor summary for specific proteins of interest #####
proteins_of_interest <- c("FABP4", "SVEP1", "PILRA", "SERPINA5")

get_top5_detail <- function(target_protein, dir_path, direction_label, suffix = "_feature_imp_shap.tsv") {
    path <- file.path(dir_path, paste0(target_protein, suffix))
    if (!file.exists(path)) return(NULL)
    
    imp <- fread(path, data.table = FALSE)
    if (nrow(imp) == 0) return(NULL)
    
    imp <- imp %>% arrange(desc(importance_mean))
    total <- sum(imp$importance_mean)
    if (total <= 0) return(NULL)
    
    imp %>%
        slice_head(n = 5) %>%
        mutate(
            target_protein = target_protein,
            direction = direction_label,
            rank = row_number(),
            share_of_total = importance_mean / total,
            is_self = feature == target_protein
        ) %>%
        select(target_protein, direction, rank, predictor_feature = feature,
               importance_mean, share_of_total, is_self)
}

## check both directions for each protein, keep whichever exists (should be one)
top5_detail <- map_dfr(proteins_of_interest, function(prot) {
    bind_rows(
        get_top5_detail(prot, "Results/feature_imp/olink_soma", "Olink to SomaScan"),
        get_top5_detail(prot, "Results/feature_imp/soma_olink", "SomaScan to Olink")
    )
})

cat("=== top5 predictors for proteins of interest ===\n")
print(top5_detail)

write.table(top5_detail, "Results/revision1/feature_importance_by_overlap/top5_detail_proteins_of_interest.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)


##### 08b - wide/readable version: one row per protein, top5 predictors concatenated #####
top5_wide <- top5_detail %>%
    group_by(target_protein, direction) %>%
    summarise(
        top5_predictors = paste0(rank, ". ", predictor_feature,
                                 " (", sprintf("%.1f%%", share_of_total * 100),
                                 if_else(is_self, ", self", "")) %>%
            paste(collapse = " | "),
        .groups = "drop"
    )

top5_wide

write.table(top5_wide, "Results/revision1/feature_importance_by_overlap/top5_summary_proteins_of_interest.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)