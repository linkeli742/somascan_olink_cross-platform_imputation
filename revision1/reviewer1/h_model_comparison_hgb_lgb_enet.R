##### 01 - library #####
library(dplyr)
library(tidyr)
library(data.table)
library(purrr)
library(ggplot2)
library(ggh4x)
library(irr)
library(car)
library(cocor)

gc()
rm(list = ls())

setwd("...")

##### 02 - files #####
# HGB (main model) -- established tiers
so_corr_ll <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = FALSE)
os_corr_ll <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = FALSE) %>%
    mutate(gain = posttrain_r - pretrain_r)

# LightGBM evaluation metrics (has pretrain_r/posttrain_r but no tiers of its own yet)
lgb_so <- fread("Results/metrics/LGB_log2scaled_5ksoma_3kolink_results.tsv", data.table = FALSE)
lgb_os <- fread("Results/metrics/LGB_3kolink_log2scaled_5ksoma_results.tsv", data.table = FALSE)

# Elastic net -- PLACEHOLDER, uncomment once available
enet_so <- fread("Results/metrics/ENET_log2scaled_5ksoma_3kolink_results.tsv", data.table = FALSE)
enet_os <- fread("Results/metrics/ENET_log2scaled_3kolink_5ksoma_results.tsv", data.table = FALSE)

# raw wide matrices for calibration (needed per model)
y_pred_so_hgb <- fread("Results/mesa_test_bootstrap_inputs/mesa_5ksoma_pred_3kolink_y_pred_wide.tsv.gz", data.table = FALSE)
y_test_so_hgb <- fread("Results/mesa_test_bootstrap_inputs/mesa_5ksoma_pred_3kolink_y_test_wide.tsv.gz", data.table = FALSE)
y_pred_os_hgb <- fread("Results/mesa_test_bootstrap_inputs/mesa_3kolink_pred_5ksoma_y_pred_wide.tsv.gz", data.table = FALSE)
y_test_os_hgb <- fread("Results/mesa_test_bootstrap_inputs/mesa_3kolink_pred_5ksoma_y_test_wide.tsv.gz", data.table = FALSE)

y_pred_so_lgb <- fread("Results/mesa_test_bootstrap_inputs/lgb_mesa_5ksoma_pred_3kolink_y_pred_wide.tsv.gz", data.table = FALSE)
y_test_so_lgb <- fread("Results/mesa_test_bootstrap_inputs/lgb_mesa_5ksoma_pred_3kolink_y_test_wide.tsv.gz", data.table = FALSE)
y_pred_os_lgb <- fread("Results/mesa_test_bootstrap_inputs/lgb_mesa_3kolink_pred_5ksoma_y_pred_wide.tsv.gz", data.table = FALSE)
y_test_os_lgb <- fread("Results/mesa_test_bootstrap_inputs/lgb_mesa_3kolink_pred_5ksoma_y_test_wide.tsv.gz", data.table = FALSE)

# PLACEHOLDER -- elastic net wide matrices
y_pred_so_enet <- fread("Results/mesa_test_bootstrap_inputs_enet/mesa_5ksoma_pred_3kolink_y_pred_wide.tsv.gz", data.table = FALSE)
y_test_so_enet <- fread("Results/mesa_test_bootstrap_inputs_enet/mesa_5ksoma_pred_3kolink_y_test_wide.tsv.gz", data.table = FALSE)
y_pred_os_enet <- fread("Results/mesa_test_bootstrap_inputs_enet/mesa_3kolink_pred_5ksoma_y_pred_wide.tsv.gz", data.table = FALSE)
y_test_os_enet <- fread("Results/mesa_test_bootstrap_inputs_enet/mesa_3kolink_pred_5ksoma_y_test_wide.tsv.gz", data.table = FALSE)

# external MS validation data (Sissala et al. 2025)
ms_olink_corr <- readxl::read_xlsx("Data/mass_spec_olink_corr/42004_2025_1753_MOESM11_ESM.xlsx") %>%
    select(Gene.Name, Correlation_Spearman, Correlation_Pearson)



##### 03 - re-derive tiers per method (same logic, method-specific posttrain_r) #####
assign_tier_vec <- function(pretrain_r, posttrain_r) {
    case_when(
        pretrain_r >= 0.7 & posttrain_r >= 0.7 ~ 1,
        pretrain_r < 0.1 & posttrain_r < 0.1 ~ 4,
        posttrain_r >= 0.4 ~ 2,
        posttrain_r < 0.4 ~ 3,
        TRUE ~ NA_real_
    )
}

## pretrain_r is computed from measured (not model-dependent) data, so it is
## identical across methods -- only posttrain_r changes per method
lgb_so_tiered <- lgb_so %>%
    select(protein, posttrain_r) %>%
    left_join(so_corr_ll %>% select(protein, pretrain_r), by = "protein") %>%
    mutate(tiers_lgb = assign_tier_vec(pretrain_r, posttrain_r))

lgb_os_tiered <- lgb_os %>%
    select(protein, posttrain_r) %>%
    left_join(os_corr_ll %>% select(protein, pretrain_r), by = "protein") %>%
    mutate(tiers_lgb = assign_tier_vec(pretrain_r, posttrain_r))

# elastic net tiers, uncomment once available
enet_so_tiered <- enet_so %>%
    select(protein, posttrain_r) %>%
    left_join(so_corr_ll %>% select(protein, pretrain_r), by = "protein") %>%
    mutate(tiers_enet = assign_tier_vec(pretrain_r, posttrain_r))
enet_os_tiered <- enet_os %>%
    select(protein, posttrain_r) %>%
    left_join(os_corr_ll %>% select(protein, pretrain_r), by = "protein") %>%
    mutate(tiers_enet = assign_tier_vec(pretrain_r, posttrain_r))




##### 04 - tier membership concordance (do the methods agree on tier assignment) #####
tier_concordance_so <- so_corr_ll %>%
    select(protein, tiers_hgb = tiers) %>%
    inner_join(lgb_so_tiered %>% select(protein, tiers_lgb), by = "protein") %>%
    inner_join(enet_so_tiered %>% select(protein, tiers_enet), by = "protein")

tier_concordance_os <- os_corr_ll %>%
    select(protein, tiers_hgb = tiers) %>%
    inner_join(lgb_os_tiered %>% select(protein, tiers_lgb), by = "protein") %>%
    inner_join(enet_os_tiered %>% select(protein, tiers_enet), by = "protein")



cat("=== SO direction: 3-way tier confusion matrix (HGB x LightGBM x ElasticNet) ===\n")
print(table(tier_concordance_so$tiers_hgb, tier_concordance_so$tiers_lgb, tier_concordance_so$tiers_enet))

cat("=== OS direction: 3-way tier confusion matrix ===\n")
print(table(tier_concordance_os$tiers_hgb, tier_concordance_os$tiers_lgb, tier_concordance_os$tiers_enet))

## pairwise Cohen's kappa for every method pair, both directions
run_kappa <- function(df, col_a, col_b) {
    kappa2(df %>% select(all_of(c(col_a, col_b))))
}

kappa_pairs <- list(
    list(pair = "HGB vs LightGBM",       cols = c("tiers_hgb", "tiers_lgb")),
    list(pair = "HGB vs ElasticNet",      cols = c("tiers_hgb", "tiers_enet")),
    list(pair = "LightGBM vs ElasticNet", cols = c("tiers_lgb", "tiers_enet"))
)

kappa_summary <- map_dfr(kappa_pairs, function(k) {
    k_so <- run_kappa(tier_concordance_so, k$cols[1], k$cols[2])
    k_os <- run_kappa(tier_concordance_os, k$cols[1], k$cols[2])
    cat(sprintf("=== %s: Cohen's kappa ===\n", k$pair))
    cat(sprintf("  SO: kappa = %.3f (p = %.3g, n = %d)\n", k_so$value, k_so$p.value, k_so$subjects))
    cat(sprintf("  OS: kappa = %.3f (p = %.3g, n = %d)\n", k_os$value, k_os$p.value, k_os$subjects))
    tibble(
        method_pair = k$pair,
        direction = c("SomaScan to Olink", "Olink to SomaScan"),
        kappa = c(k_so$value, k_os$value),
        p_value = c(k_so$p.value, k_os$p.value),
        n = c(k_so$subjects, k_os$subjects)
    )
})

kappa_summary

write.table(kappa_summary, "Results/revision1/model_comparison/tier_concordance_kappa.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tier_concordance_so, "Results/revision1/model_comparison/SO_tier_concordance_per_protein.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tier_concordance_os, "Results/revision1/model_comparison/OS_tier_concordance_per_protein.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

##### 05 - external MS validation: eta-squared per method (how much does each method's #####
#####      tiering explain independent, ground-truth reliability variation)            #####

compute_eta_squared <- function(ms_corr, tier) {
    df <- data.frame(ms_corr = ms_corr, tier = factor(tier)) %>% na.omit()
    fit <- aov(ms_corr ~ tier, data = df)
    ss <- summary(fit)[[1]]
    list(eta_squared = ss[1, "Sum Sq"] / sum(ss[, "Sum Sq"]), n = nrow(df))
}

ms_hgb_so <- so_corr_ll %>% select(protein, tiers, posttrain_r) %>%
    inner_join(ms_olink_corr, by = c("protein" = "Gene.Name"))
ms_lgb_so <- lgb_so_tiered %>% select(protein, tiers_lgb, posttrain_r) %>%
    inner_join(ms_olink_corr, by = c("protein" = "Gene.Name"))
ms_enet_so <- enet_so_tiered %>% select(protein, tiers_enet, posttrain_r) %>%
    inner_join(ms_olink_corr, by = c("protein" = "Gene.Name"))

ms_hgb_os <- os_corr_ll %>% select(protein, tiers, posttrain_r) %>%
    inner_join(ms_olink_corr, by = c("protein" = "Gene.Name"))
ms_lgb_os <- lgb_os_tiered %>% select(protein, tiers_lgb, posttrain_r) %>%
    inner_join(ms_olink_corr, by = c("protein" = "Gene.Name"))
ms_enet_os <- enet_os_tiered %>% select(protein, tiers_enet, posttrain_r) %>%
    inner_join(ms_olink_corr, by = c("protein" = "Gene.Name"))

eta_hgb_so  <- compute_eta_squared(ms_hgb_so$Correlation_Spearman, ms_hgb_so$tiers)
eta_lgb_so  <- compute_eta_squared(ms_lgb_so$Correlation_Spearman, ms_lgb_so$tiers_lgb)
eta_enet_so <- compute_eta_squared(ms_enet_so$Correlation_Spearman, ms_enet_so$tiers_enet)
eta_hgb_os  <- compute_eta_squared(ms_hgb_os$Correlation_Spearman, ms_hgb_os$tiers)
eta_lgb_os  <- compute_eta_squared(ms_lgb_os$Correlation_Spearman, ms_lgb_os$tiers_lgb)
eta_enet_os <- compute_eta_squared(ms_enet_os$Correlation_Spearman, ms_enet_os$tiers_enet)

eta_squared_summary <- tibble(
    method = rep(c("HGB", "LightGBM", "ElasticNet"), times = 2),
    direction = rep(c("SomaScan to Olink", "Olink to SomaScan"), each = 3),
    eta_squared = c(eta_hgb_so$eta_squared, eta_lgb_so$eta_squared, eta_enet_so$eta_squared,
                    eta_hgb_os$eta_squared, eta_lgb_os$eta_squared, eta_enet_os$eta_squared),
    n_matched = c(eta_hgb_so$n, eta_lgb_so$n, eta_enet_so$n,
                  eta_hgb_os$n, eta_lgb_os$n, eta_enet_os$n)
)

eta_squared_summary

write.table(eta_squared_summary, "Results/revision1/model_comparison/eta_squared_ms_validation.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

##### 06 - Steiger's Z: does HGB's posttrain_r agree with MS ground truth more than LightGBM's does #####
## (comparing two overlapping correlations computed on the same set of proteins)

## cocor() only compares two overlapping correlations at a time, so with 3 methods
## we run 3 pairwise comparisons (HGB-LGB, HGB-ENET, LGB-ENET) rather than one omnibus test
run_cocor_comparison <- function(ms_a_df, ms_b_df) {
    merged <- ms_a_df %>%
        select(protein, ms_corr = Correlation_Spearman, r_a = posttrain_r) %>%
        inner_join(ms_b_df %>% select(protein, r_b = posttrain_r), by = "protein") %>%
        na.omit()
    
    cocor(~ ms_corr + r_a | ms_corr + r_b, data = merged)
}

cocor_pairs <- list(
    list(pair = "HGB vs LightGBM",       so = list(ms_hgb_so, ms_lgb_so),   os = list(ms_hgb_os, ms_lgb_os)),
    list(pair = "HGB vs ElasticNet",      so = list(ms_hgb_so, ms_enet_so),  os = list(ms_hgb_os, ms_enet_os)),
    list(pair = "LightGBM vs ElasticNet", so = list(ms_lgb_so, ms_enet_so),  os = list(ms_lgb_os, ms_enet_os))
)

extract_steiger <- function(cocor_result) {
    st <- cocor_result@steiger1980
    tibble(statistic = unname(st$statistic), p_value = unname(st$p.value))
}

steiger_summary <- map_dfr(cocor_pairs, function(cp) {
    cocor_so <- run_cocor_comparison(cp$so[[1]], cp$so[[2]])
    cocor_os <- run_cocor_comparison(cp$os[[1]], cp$os[[2]])
    
    cat(sprintf("=== SO direction: Steiger's Z, %s agreement with MS ground truth ===\n", cp$pair))
    print(cocor_so)
    cat(sprintf("=== OS direction: Steiger's Z, %s ===\n", cp$pair))
    print(cocor_os)
    
    bind_rows(
        extract_steiger(cocor_so) %>% mutate(direction = "SomaScan to Olink"),
        extract_steiger(cocor_os) %>% mutate(direction = "Olink to SomaScan")
    ) %>% mutate(method_pair = cp$pair, .before = 1)
})

steiger_summary

write.table(steiger_summary, "Results/revision1/model_comparison/steiger_z_ms_validation.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)


