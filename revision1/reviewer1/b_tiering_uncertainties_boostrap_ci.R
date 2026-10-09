##### 01 - library #####
library(data.table)
library(dplyr)
library(purrr)
library(tibble)

gc()
rm(list=ls())


setwd("/medpop/esp2/lli/somalink_olink/")


##### 02 - files #####
overlapping_prots <- fread("Results/rank_by_gain/soma5k_olink3k.tsv")$protein
## SO direction (SomaScan -> Olink)
y_pred_so <- fread("Results/mesa_test_bootstrap_inputs/mesa_5ksoma_pred_3kolink_y_pred_wide.tsv.gz", data.table = F) %>%
    select(id, all_of(overlapping_prots))
y_test_so <- fread("Results/mesa_test_bootstrap_inputs/mesa_5ksoma_pred_3kolink_y_test_wide.tsv.gz", data.table = F) %>%
    select(id, all_of(overlapping_prots))
x_test_so <- fread("Results/mesa_test_bootstrap_inputs/mesa_5ksoma_pred_3kolink_x_test_protein_wide.tsv.gz", data.table = F) %>%
    select(id, all_of(overlapping_prots))

## OS direction (Olink -> SomaScan)
y_pred_os <- fread("Results/mesa_test_bootstrap_inputs/mesa_3kolink_pred_5ksoma_y_pred_wide.tsv.gz", 
                   data.table = F) %>%
    select(id, all_of(overlapping_prots))
y_test_os <- fread("Results/mesa_test_bootstrap_inputs/mesa_3kolink_pred_5ksoma_y_test_wide.tsv.gz", 
                   data.table = F) %>%
    select(id, all_of(overlapping_prots))
x_test_os <- fread("Results/mesa_test_bootstrap_inputs/mesa_3kolink_pred_5ksoma_x_test_protein_wide.tsv.gz", 
                   data.table = F) %>%
    select(id, all_of(overlapping_prots))


y_pred_so_mat <- y_pred_so %>% select(-id) %>% as.matrix()
y_test_so_mat <- y_test_so %>% select(-id) %>% as.matrix()
x_test_so_mat <- x_test_so %>% select(-id) %>% as.matrix()

y_pred_os_mat <- y_pred_os %>% select(-id) %>% as.matrix()
y_test_os_mat <- y_test_os %>% select(-id) %>% as.matrix()
x_test_os_mat <- x_test_os %>% select(-id) %>% as.matrix()



#### vectorized pearson corr
colwise_cor <- function(A, B) {
    n <- nrow(A) # assuming A and B are equal length (which they are and should be)
    # formula is sigma(ZxZy) / (n-1)
    A_z <- scale(A)
    B_z <- scale(B)
    colSums(A_z * B_z, na.rm = FALSE) / (n - 1)
}


set.seed(613)
n_boot=1000

## running bootstrap
run_bootstrap <- function(y_pred, y_test, x_test, n_boot = 1000) {
    
    n <- nrow(y_pred)
    
    boot_posttrain_r <- matrix(NA_real_, nrow = n_boot, ncol = length(overlapping_prots),
                               dimnames = list(NULL, overlapping_prots))
    boot_testset_r <- matrix(NA_real_, nrow = n_boot, ncol = length(overlapping_prots),
                             dimnames = list(NULL, overlapping_prots))
    
    for (b in seq_len(n_boot)) {
        
        idx <- sample.int(n, size = n, replace = TRUE)
        boot_posttrain_r[b, ] <- colwise_cor(y_pred[idx, ], x_test[idx, ])
        boot_testset_r[b, ]   <- colwise_cor(y_test[idx, ], x_test[idx, ])
        if (b %% 250 == 0) cat("bootstrap iter:", b, "/", n_boot, "\n")
        
    }
    
    list(posttrain_r = boot_posttrain_r, testset_r = boot_testset_r, proteins = overlapping_prots)
}


boot_so <- run_bootstrap(y_pred_so_mat, y_test_so_mat, x_test_so_mat, n_boot = n_boot)
boot_os <- run_bootstrap(y_pred_os_mat, y_test_os_mat, x_test_os_mat, n_boot = n_boot)



##### 03 - combine with the full set to see the distribution of the correlation #####
soma_olink_orig <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = FALSE)
olink_soma_orig <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = FALSE)

pretrain_r_so <- setNames(soma_olink_orig$pretrain_r, soma_olink_orig$protein)
pretrain_r_os <- setNames(olink_soma_orig$pretrain_r, olink_soma_orig$protein)

assign_tier_vec <- function(pretrain_r, posttrain_r) {
    case_when(
        pretrain_r >= 0.7 & posttrain_r >= 0.7 ~ 1,
        pretrain_r < 0.1 & posttrain_r < 0.1 ~ 4,
        posttrain_r >= 0.4 ~ 2,
        posttrain_r < 0.4 ~ 3,
        TRUE ~ NA_real_
    )
}

compute_tier_probs <- function(boot_result, pretrain_r_lookup, direction_label) {
    map_dfr(boot_result$proteins, function(p) {
        pr <- pretrain_r_lookup[[p]]
        if (is.null(pr) || is.na(pr)) return(NULL)
        
        # this is feeding n_boot times baseline r to match with the bootstrapping number
        boot_tiers <- assign_tier_vec(pr, boot_result$posttrain_r[, p])
        
        tibble(tier = boot_tiers) %>%
            count(tier) %>% # for each tier
            mutate(
                protein = p,
                # how many times did protein landed in this tier?
                probability = n / sum(n),
                # this is just the training direction
                direction = direction_label,
                pretrain_r_fixed = pr,
                posttrain_r_ci_low  = quantile(boot_result$posttrain_r[, p], 0.025, na.rm = TRUE),
                posttrain_r_ci_high = quantile(boot_result$posttrain_r[, p], 0.975, na.rm = TRUE)
            )
    })
}

tier_probs_so <- compute_tier_probs(boot_so, pretrain_r_so, "SomaScan to Olink")
tier_probs_os <- compute_tier_probs(boot_os, pretrain_r_os, "Olink to SomaScan")


write.table(tier_probs_so, glue::glue("Results/mesa_test_bootstrap_tier_probs/so_{n_boot}bootstraps.tsv"), 
            sep = "\t", row.names = F, quote = F)
write.table(tier_probs_os, glue::glue("Results/mesa_test_bootstrap_tier_probs/os_{n_boot}bootstraps.tsv"), 
            sep = "\t", row.names = F, quote = F)


