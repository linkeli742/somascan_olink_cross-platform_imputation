##### 01 - library #####
library(dplyr)
library(tidyr)
library(data.table)
library(purrr)


gc()
rm(list=ls())

setwd("...")



##### 02 - files #####
# raw files ---------------------------------------------------------------
so_corr_ll <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F)
os_corr_ll <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = F)

## SO direction
y_pred_so <- fread("Results/mesa_test_bootstrap_inputs/mesa_5ksoma_pred_3kolink_y_pred_wide.tsv.gz", data.table = F) 
y_test_so <- fread("Results/mesa_test_bootstrap_inputs/mesa_5ksoma_pred_3kolink_y_test_wide.tsv.gz", data.table = F) 


## OS direction
y_pred_os <- fread("Results/mesa_test_bootstrap_inputs/mesa_3kolink_pred_5ksoma_y_pred_wide.tsv.gz", 
                   data.table = F)
y_test_os <- fread("Results/mesa_test_bootstrap_inputs/mesa_3kolink_pred_5ksoma_y_test_wide.tsv.gz", 
                   data.table = F)


y_pred_so_mat <- y_pred_so %>% select(-id) %>% as.matrix()
y_test_so_mat <- y_test_so %>% select(-id) %>% as.matrix()


y_pred_os_mat <- y_pred_os %>% select(-id) %>% as.matrix()
y_test_os_mat <- y_test_os %>% select(-id) %>% as.matrix()






##### 03 - analyses #####
# calibration slope, intercept, and agreement (RMSE / CCC) ----------------
# CCC (Lin's concordance correlation coefficient) is computed manually
# CCC = 2*cov(pred,actual) / (var(pred) + var(actual) + (mean(pred)-mean(actual))^2)

calibration_check <- function(y_pred_mat, y_actual_mat) {
    proteins <- intersect(colnames(y_pred_mat), colnames(y_actual_mat))
    map_dfr(proteins, function(p) {
        pred <- y_pred_mat[, p]
        actual <- y_actual_mat[, p]
        
        fit <- lm(actual ~ pred)
        
        # CCC on complete cases for this protein
        cc <- complete.cases(pred, actual)
        pred_c <- pred[cc]
        actual_c <- actual[cc]
        n <- length(pred_c)
        
        if (n > 1) {
            mp <- mean(pred_c)
            ma <- mean(actual_c)
            vp <- var(pred_c) * (n - 1) / n   # population variance
            va <- var(actual_c) * (n - 1) / n
            cv <- cov(pred_c, actual_c) * (n - 1) / n
            ccc_val <- 2 * cv / (vp + va + (mp - ma)^2)
        } else {
            ccc_val <- NA_real_
        }
        
        tibble(
            protein = p,
            calib_slope = coef(fit)[2],
            calib_intercept = coef(fit)[1],
            rmse = sqrt(mean((actual - pred)^2, na.rm = TRUE)),
            mae = mean(abs(actual - pred), na.rm = TRUE),
            sd_ratio = sd(pred, na.rm = TRUE) / sd(actual, na.rm = TRUE),  # variance preservation
            ccc = ccc_val
        )
    })
}


calib_so <- calibration_check(y_pred_so_mat, y_test_so_mat)
calib_os <- calibration_check(y_pred_os_mat, y_test_os_mat)



calib_so_table <- calib_so %>%
    left_join(so_corr_ll %>% select(protein, tiers), by = "protein") %>%
    group_by(tiers) %>%
    summarise(
        n_proteins= n(),
        median_slope = median(calib_slope, na.rm = TRUE),
        median_intercept = median(calib_intercept, na.rm = TRUE),
        median_rmse = median(rmse, na.rm = TRUE),
        median_sd_ratio = median(sd_ratio, na.rm = TRUE),
        median_ccc = median(ccc, na.rm = TRUE)
    )

summary(calib_so$calib_intercept)


write.table(calib_so_table, "Results/revision1/calibration/SO_cali_slope_int_rmse_sdratio.tsv",
            sep="\t", row.names = F, quote = F)

calib_so_per_prot <- calib_so %>%
    left_join(so_corr_ll %>% select(protein, tiers), by = "protein") %>%
    relocate(tiers, .after=protein)


write.table(calib_so_per_prot, 
            "Results/revision1/calibration/SO_cali_per_protein.tsv",
            sep="\t", row.names = F, quote = F)




### os direction
calib_os_table <- calib_os %>%
    left_join(os_corr_ll %>% select(protein, tiers), by = "protein") %>%
    group_by(tiers) %>%
    summarise(
        n_proteins= n(),
        median_slope = median(calib_slope, na.rm = TRUE),
        median_intercept = median(calib_intercept, na.rm = TRUE),
        median_rmse = median(rmse, na.rm = TRUE),
        median_sd_ratio = median(sd_ratio, na.rm = TRUE),
        median_ccc = median(ccc, na.rm = TRUE)
    )


write.table(calib_os_table, "Results/revision1/calibration/OS_cali_slope_int_rmse_sdratio.tsv",
            sep="\t", row.names = F, quote = F)

calib_os_per_prot <- calib_os %>%
    left_join(os_corr_ll %>% select(protein, tiers), by = "protein") %>%
    relocate(tiers, .after=protein)

write.table(calib_os_per_prot, 
            "Results/revision1/calibration/OS_cali_per_protein.tsv",
            sep="\t", row.names = F, quote = F)



# all proteins summary ----------------------------------------------------
calib_so_table_allprot <- calib_so %>%
    summarise(
        tiers="All",
        n_proteins= n(),
        median_slope = median(calib_slope, na.rm = TRUE),
        median_intercept = median(calib_intercept, na.rm = TRUE),
        median_rmse = median(rmse, na.rm = TRUE),
        median_sd_ratio = median(sd_ratio, na.rm = TRUE),
        median_ccc = median(ccc, na.rm = TRUE)
    )


calib_os_table_allprot <- calib_os %>%
    summarise(
        tiers="All",
        n_proteins= n(),
        median_slope = median(calib_slope, na.rm = TRUE),
        median_intercept = median(calib_intercept, na.rm = TRUE),
        median_rmse = median(rmse, na.rm = TRUE),
        median_sd_ratio = median(sd_ratio, na.rm = TRUE),
        median_ccc = median(ccc, na.rm = TRUE)
    )




write.table(calib_so_table_allprot, "Results/revision1/calibration/SO_allprots_cali_slope_int_rmse_sdratio.tsv",
            sep="\t", row.names = F, quote = F)
write.table(calib_os_table_allprot, "Results/revision1/calibration/OS_allprots_cali_slope_int_rmse_sdratio.tsv",
            sep="\t", row.names = F, quote = F)




# pivoting wider ----------------------------------------------------------
calib_table_full <- rbind(
    rbind(calib_so_table, calib_so_table_allprot) %>% mutate(direction = "SomaScan to Olink"),
    rbind(calib_os_table, calib_os_table_allprot) %>% mutate(direction = "Olink to SomaScan")
)

calib_table_full_wide <- calib_table_full %>%
    pivot_wider(
        id_cols = tiers,
        names_from = direction,
        values_from = c(n_proteins, median_slope, median_intercept, median_rmse, median_sd_ratio, median_ccc),
        names_glue = "{.value}_{direction}"
    ) %>%
    relocate(ends_with("SomaScan to Olink"), .after=tiers)

calib_table_full_wide

write.table(calib_table_full_wide, "Results/revision1/calibration/ALLdirection_allprots_cali_slope_int_rmse_sdratio_wide.tsv",
            sep="\t", row.names = F, quote = F)


