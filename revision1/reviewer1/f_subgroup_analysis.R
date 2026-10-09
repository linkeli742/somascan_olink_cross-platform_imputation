##### 01 - library #####
library(dplyr)
library(tidyr)
library(data.table)
library(purrr)
library(ggplot2)
library(ggh4x)
library(ggsignif)
library(patchwork)

gc()
rm(list = ls())


setwd("...")


##### 02 - files #####
mesa_demo <- fread("Data/table1_datasets/mesa_demos.tsv", data.table = FALSE)
overlapping_prots <- fread("Results/rank_by_gain/soma5k_olink3k.tsv")$protein

so_corr_ll <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = FALSE)
os_corr_ll <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = FALSE)

y_pred_so <- fread("Results/mesa_test_bootstrap_inputs/mesa_5ksoma_pred_3kolink_y_pred_wide.tsv.gz", data.table = FALSE) %>%
    select(id, all_of(overlapping_prots))
x_test_so <- fread("Results/mesa_test_bootstrap_inputs/mesa_5ksoma_pred_3kolink_x_test_protein_wide.tsv.gz", data.table = FALSE) %>%
    select(id, all_of(overlapping_prots))

y_pred_os <- fread("Results/mesa_test_bootstrap_inputs/mesa_3kolink_pred_5ksoma_y_pred_wide.tsv.gz", data.table = FALSE) %>%
    select(id, all_of(overlapping_prots))
x_test_os <- fread("Results/mesa_test_bootstrap_inputs/mesa_3kolink_pred_5ksoma_x_test_protein_wide.tsv.gz", data.table = FALSE) %>%
    select(id, all_of(overlapping_prots))



##### 03 - core functions #####
## per-protein, per-subgroup correlation
subgroup_performance <- function(y_pred, x_test, covariates, group_var) {
    ids <- y_pred$id
    strata <- covariates[[group_var]][match(ids, covariates$id)]
    proteins <- setdiff(intersect(colnames(y_pred), colnames(x_test)), "id")
    
    map_dfr(unique(strata), function(g) {
        idx <- which(strata == g)
        if (length(idx) < 10) return(NULL)
        
        r_vals <- sapply(proteins, function(p) cor(y_pred[idx, p], x_test[idx, p], use = "complete.obs"))
        tibble(protein = proteins, group = g, group_var = group_var, r = r_vals, n = length(idx))
    })
}



## paired Wilcoxon (sex, 2 groups)
run_paired_sex_test <- function(perf_df) {
    wide <- perf_df %>% 
        select(protein, group, r) %>% 
        pivot_wider(names_from = group, values_from = r)
    
    group_cols <- setdiff(colnames(wide), "protein") # strata groups
    
    wtest <- wilcox.test(wide[[group_cols[1]]], wide[[group_cols[2]]], paired = TRUE)
    tibble(test_method = "Paired Wilcoxon", p_value = wtest[["p.value"]], n_proteins = nrow(wide))
}


## Friedman test (race, >2 groups, paired/repeated-measures)
run_friedman_race_test <- function(perf_df) {
    wide <- perf_df %>% 
        select(protein, group, r) %>% 
        pivot_wider(names_from = group, values_from = r)
    
    group_cols <- setdiff(colnames(wide), "protein")
    mat <- as.matrix(wide[, group_cols])
    ftest <- friedman.test(mat)
    tibble(test_method = "Friedman", p_value = ftest[["p.value"]], n_proteins = nrow(wide))
}

## tier-stratified wrapper functions to test within tier differences
## the above functions test all the proteins, not just overlapping proteins
## tiers would make the 
run_test_by_tier <- function(perf_df, tier_lookup, test_fun) {
    perf_df %>%
        inner_join(tier_lookup, by = "protein") %>%
        group_by(tiers) %>%
        group_modify(~ test_fun(.x)) %>%
        ungroup()
}

## the ultimate wrapper: runs everything for a direction/variable, returns per-protein data +
## overall test + tier-stratified test + summary table (long, with p-values attached)
run_subgroup_analysis <- function(y_pred, x_test, covariates, group_var, tier_lookup, test_fun) {
    perf_df <- subgroup_performance(y_pred, x_test, covariates, group_var)
    
    overall_test <- test_fun(perf_df)
    tier_tests <- run_test_by_tier(perf_df, tier_lookup, test_fun)
    
    ## summary: per tier x group, median r, with tier-specific p attached
    summary_by_tier <- perf_df %>%
        inner_join(tier_lookup, by = "protein") %>%
        group_by(group, tiers) %>%
        summarise(median_r = median(r, na.rm = TRUE), n_proteins = n(), .groups = "drop") %>%
        left_join(tier_tests %>% select(tiers, tier_p_value = p_value), by = "tiers")
    
    ## overall (all proteins, tier-agnostic) summary
    summary_overall <- perf_df %>%
        group_by(group) %>%
        summarise(median_r = median(r, na.rm = TRUE), n_proteins = n(), .groups = "drop") %>%
        mutate(tiers = "Overall", overall_p_value = overall_test$p_value)
    
    list(
        per_protein = perf_df,
        overall_test = overall_test,
        tier_tests = tier_tests,
        summary_by_tier = summary_by_tier,
        summary_overall = summary_overall
    )
}



##### 04 - running for both directions ######
tier_lookup_so <- so_corr_ll %>% 
    select(protein, tiers)
tier_lookup_os <- os_corr_ll %>% 
    select(protein, tiers)

### these will give list results
sex_so  <- run_subgroup_analysis(y_pred_so, x_test_so, mesa_demo, "sex",  tier_lookup_so, run_paired_sex_test)
race_so <- run_subgroup_analysis(y_pred_so, x_test_so, mesa_demo, "race", tier_lookup_so, run_friedman_race_test)
sex_os  <- run_subgroup_analysis(y_pred_os, x_test_os, mesa_demo, "sex",  tier_lookup_os, run_paired_sex_test)
race_os <- run_subgroup_analysis(y_pred_os, x_test_os, mesa_demo, "race", tier_lookup_os, run_friedman_race_test)

##### 05 - save tables #####
save_all <- function(result, direction_label, var_label) {
    prefix <- paste0("Results/revision1/subgroup_perf/", direction_label, "_", var_label)
    write.table(result$per_protein, paste0(prefix, "_per_protein.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(result$summary_by_tier, paste0(prefix, "_summary_by_tier.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(result$summary_overall, paste0(prefix, "_summary_overall.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}

save_all(sex_so,  "SO", "sex")
save_all(race_so, "SO", "race")
save_all(sex_os,  "OS", "sex")
save_all(race_os, "OS", "race")
