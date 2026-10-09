##### 01 - library #####
library(dplyr)
library(data.table)

gc()
rm(list=ls())

setwd("...")


##### 02 - files #####
so_corr_ll <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = FALSE)
os_corr_ll <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = FALSE)

## SO direction (SomaScan -> Olink)
y_pred_so <- fread("Results/mesa_test_bootstrap_inputs/mesa_5ksoma_pred_3kolink_y_pred_wide.tsv.gz", data.table = F)
y_test_so <- fread("Results/mesa_test_bootstrap_inputs/mesa_5ksoma_pred_3kolink_y_test_wide.tsv.gz", data.table = F) 

## OS direction (Olink -> SomaScan)
y_pred_os <- fread("Results/mesa_test_bootstrap_inputs/mesa_3kolink_pred_5ksoma_y_pred_wide.tsv.gz", data.table = F)
y_test_os <- fread("Results/mesa_test_bootstrap_inputs/mesa_3kolink_pred_5ksoma_y_test_wide.tsv.gz", data.table = F)


# both true
identical(y_pred_so$id, y_test_so$id)
identical(y_pred_os$id, y_test_os$id)



##### 03 - analysis #####
build_conformal_intervals <- function(y_pred, y_test, alpha = 0.05) {
    proteins <- setdiff(intersect(colnames(y_pred), colnames(y_test)), "id")
    
    map_dfr(proteins, function(p) {
        pred <- y_pred[, p]
        actual <- y_test[, p]
        if (all(is.na(actual))) return(NULL)
        
        abs_resid <- abs(actual - pred)
        ## finite-sample-adjusted quantile (Lei et al. 2018)
        n <- sum(!is.na(abs_resid))
        q_level <- min(1, ceiling((1 - alpha) * (n + 1)) / n)
        conformal_width <- quantile(abs_resid, q_level, na.rm = TRUE)
        
        tibble(protein = p, 
               conformal_half_width = conformal_width,
               q_adj = q_level,
               n_calibration = n)
    })
}

conformal_so <- build_conformal_intervals(y_pred_so, y_test_so)
conformal_os <- build_conformal_intervals(y_pred_os, y_test_os)



##### 04- summaries ####
conformal_so_summary <- conformal_so %>%
    left_join(so_corr_ll %>% select(protein, tiers), by = "protein") %>%
    mutate(tiers = as.character(tiers)) %>%
    group_by(tiers) %>%
    summarise(
        mean_half_width   = mean(conformal_half_width, na.rm = TRUE),
        median_half_width = median(conformal_half_width, na.rm = TRUE),
        max_half_width    = max(conformal_half_width, na.rm = TRUE),
        min_half_width    = min(conformal_half_width, na.rm = TRUE),
        n_proteins        = n()
    )

conformal_so_summary

conformal_os_summary <- conformal_os %>%
    left_join(os_corr_ll %>% select(protein, tiers), by = "protein") %>%
    mutate(tiers = as.character(tiers)) %>%
    group_by(tiers) %>%
    summarise(
        mean_half_width   = mean(conformal_half_width, na.rm = TRUE),
        median_half_width = median(conformal_half_width, na.rm = TRUE),
        max_half_width    = max(conformal_half_width, na.rm = TRUE),
        min_half_width    = min(conformal_half_width, na.rm = TRUE),
        n_proteins        = n()
    )

conformal_os_summary

##### 04b - overall summary (all proteins, no tier stratification) ####
conformal_so_overall <- conformal_so %>%
    summarise(
        tiers             = "All",
        mean_half_width   = mean(conformal_half_width, na.rm = TRUE),
        median_half_width = median(conformal_half_width, na.rm = TRUE),
        max_half_width    = max(conformal_half_width, na.rm = TRUE),
        min_half_width    = min(conformal_half_width, na.rm = TRUE),
        n_proteins        = n()
    )

conformal_os_overall <- conformal_os %>%
    summarise(
        tiers             = "All",
        mean_half_width   = mean(conformal_half_width, na.rm = TRUE),
        median_half_width = median(conformal_half_width, na.rm = TRUE),
        max_half_width    = max(conformal_half_width, na.rm = TRUE),
        min_half_width    = min(conformal_half_width, na.rm = TRUE),
        n_proteins        = n()
    )

## append to per-tier summaries before pivoting wider
conformal_so_summary <- bind_rows(conformal_so_summary, conformal_so_overall)
conformal_os_summary <- bind_rows(conformal_os_summary, conformal_os_overall)


## combining and pivot wider
conformal_combined <- bind_rows(
    conformal_so_summary %>% mutate(direction = "SomaScan to Olink"),
    conformal_os_summary %>% mutate(direction = "Olink to SomaScan")
)

conformal_combined_wide <- conformal_combined %>%
    pivot_wider(
        id_cols = tiers,
        names_from = direction,
        values_from = c(mean_half_width, median_half_width, max_half_width, min_half_width, n_proteins),
        names_glue = "{.value}_{direction}"
    ) %>%
    relocate(ends_with("SomaScan to Olink"), .after=tiers)

conformal_combined_wide



#### saving results
# summary level
write.table(conformal_combined_wide, "Results/revision1/uncertainty_estimate/conformal_combined_wide.tsv", 
            sep = "\t", row.names = F, quote = F)


# individual protein level
write.table(conformal_so, 
            "Results/revision1/uncertainty_estimate/SO_conformal_per_protein.tsv", 
            sep = "\t", row.names = F, quote = F)
write.table(conformal_os, 
            "Results/revision1/uncertainty_estimate/OS_conformal_per_protein.tsv", 
            sep = "\t", row.names = F, quote = F)








