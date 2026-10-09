##### 01 - library #####
library(data.table)
library(dplyr)
library(glue)
library(purrr)

gc()
rm(list=ls())

setwd("...")


##### 02 - files #####
overlapping_proteins <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F)$protein
# raw files ---------------------------------------------------------------
so_corr_ll <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F)
os_corr_ll <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = F)



## SO direction
y_pred_so <- fread("Results/mesa_test_bootstrap_inputs/mesa_5ksoma_pred_3kolink_y_pred_wide.tsv.gz", data.table = F) %>%
    select(id, all_of(overlapping_proteins))
y_test_so <- fread("Results/mesa_test_bootstrap_inputs/mesa_5ksoma_pred_3kolink_y_test_wide.tsv.gz", data.table = F) %>%
    select(id, all_of(overlapping_proteins))
x_test_so <- fread("Results/mesa_test_bootstrap_inputs/mesa_5ksoma_pred_3kolink_x_test_protein_wide.tsv.gz", data.table = F) %>%
    select(id, all_of(overlapping_proteins))


## OS direction
y_pred_os <- fread("Results/mesa_test_bootstrap_inputs/mesa_3kolink_pred_5ksoma_y_pred_wide.tsv.gz", 
                   data.table = F) %>%
    select(id, all_of(overlapping_proteins))
y_test_os <- fread("Results/mesa_test_bootstrap_inputs/mesa_3kolink_pred_5ksoma_y_test_wide.tsv.gz", 
                   data.table = F) %>%
    select(id, all_of(overlapping_proteins))
x_test_os <- fread("Results/mesa_test_bootstrap_inputs/mesa_3kolink_pred_5ksoma_x_test_protein_wide.tsv.gz", 
                   data.table = F) %>%
    select(id, all_of(overlapping_proteins))

y_pred_so_mat <- y_pred_so %>% select(-id) %>% as.matrix()
y_test_so_mat <- y_test_so %>% select(-id) %>% as.matrix()
x_test_so_mat <- x_test_so %>% select(-id) %>% as.matrix()

y_pred_os_mat <- y_pred_os %>% select(-id) %>% as.matrix()
y_test_os_mat <- y_test_os %>% select(-id) %>% as.matrix()
x_test_os_mat <- x_test_os %>% select(-id) %>% as.matrix()



##### 03 - test #####
## rank-based inverse normal transformation
inverse_normal_transform <- function(x) {
    qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

## apply INT to your existing y_pred/y_test wide matrices, recompute correlation
compute_int_correlation <- function(y_pred_mat, x_test_mat) {
    proteins <- intersect(colnames(y_pred_mat), colnames(x_test_mat))
    map_dfr(proteins, function(p) {
        pred_int <- inverse_normal_transform(y_pred_mat[, p])
        actual_int <- inverse_normal_transform(x_test_mat[, p])
        tibble(protein = p, posttrain_r_int = cor(pred_int, actual_int, use = "complete.obs"))
    })
}

int_comparison_so <- compute_int_correlation(y_pred_so_mat, x_test_so_mat) %>%
    left_join(so_corr_ll %>% select(protein, tiers, posttrain_r), by = "protein")
int_comparison_os <- compute_int_correlation(y_pred_os_mat, x_test_os_mat) %>%
    left_join(os_corr_ll %>% select(protein, tiers, posttrain_r), by = "protein")

## does INT meaningfully change the correlation or tier ordering?
cor(int_comparison_so$posttrain_r, int_comparison_so$posttrain_r_int, use = "complete.obs")
cor(int_comparison_os$posttrain_r, int_comparison_os$posttrain_r_int, use = "complete.obs")

int_comparison_so %>%
    group_by(tiers) %>%
    summarise(median_original = median(posttrain_r, na.rm = TRUE),
              median_int = median(posttrain_r_int, na.rm = TRUE))


int_comparison_os %>%
    group_by(tiers) %>%
    summarise(median_original = median(posttrain_r, na.rm = TRUE),
              median_int = median(posttrain_r_int, na.rm = TRUE))















