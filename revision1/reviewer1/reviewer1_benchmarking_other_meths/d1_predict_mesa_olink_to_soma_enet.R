##### 01 - library #####
library(glmnet)
library(data.table)
library(dplyr)
library(tidyr)
library(glue)

gc()
rm(list=ls())

setwd("...")

##### 02 - files #####
MESA_Olink <- fread("Data/final_datasets/MESA_Olink_3k_log2_nodup_rm10na_knn_final.tsv.gz", data.table = FALSE)
MESA_Soma  <- fread("Data/final_datasets/MESA_SOMA_5k_log2scaled_nodup_final.tsv.gz", data.table = FALSE) %>%
  select(-bmi)

overlap_ids <- intersect(unique(MESA_Olink$id), unique(MESA_Soma$id))

MESA_Olink <- MESA_Olink %>% 
    filter(id %in% overlap_ids) %>% 
    arrange(id)
MESA_Soma  <- MESA_Soma  %>% 
    filter(id %in% overlap_ids) %>% 
    arrange(id)

stopifnot(all(MESA_Olink$id == MESA_Soma$id))

rownames(MESA_Olink) <- MESA_Olink$id
rownames(MESA_Soma)  <- MESA_Soma$id

MESA_Olink <- MESA_Olink %>%
    select(-id)
MESA_Soma <- MESA_Soma %>%
    select(-id)

train_index <- fread("Data/helper_data/train_index.csv", data.table = FALSE)$id
test_index  <- fread("Data/helper_data/test_index.csv", data.table = FALSE)$id

## OS-direction target list (SomaScan proteins)
protlist <- colnames(MESA_Soma)
olink_cols <- colnames(MESA_Olink)

##### 03 - predict on MESA test set using saved elastic net models #####
X_test <- as.matrix(MESA_Olink[as.character(test_index), , drop = FALSE])

all_test_values <- list()

for (protein_name in protlist) {
    model_path <- file.path(
        "Data/soma_trained_models/olink_soma/3kolink_5k_log2scaled_soma_enet/",
        glue("{protein_name}_enet.rds")
    )
    
    if (!file.exists(model_path)) {
        cat("skip", protein_name, "- no saved elastic net model found\n")
        next
    }
    
    fit <- readRDS(model_path)
    cat("predicting", protein_name, "\n")
    
    y_test <- MESA_Soma[as.character(test_index), protein_name]
    y_pred <- as.numeric(predict(fit, newx = X_test, s = "lambda.min"))
    
    if (protein_name %in% olink_cols) {
        x_test_protein <- MESA_Olink[as.character(test_index), protein_name]
    } else {
        x_test_protein <- NA
    }
    
    all_test_values[[protein_name]] <- data.frame(
        id = test_index,
        protein = protein_name,
        y_test = y_test,
        y_pred = y_pred,
        X_test_protein = x_test_protein
    )
}

raw_test_values_os_df <- bind_rows(all_test_values)


##### 03b - compute pre/post-training correlation metrics (mirrors LGB output) #####
results_list <- list()

for (protein_name in unique(raw_test_values_os_df$protein)) {
    sub <- raw_test_values_os_df %>% filter(protein == protein_name)
    y_test <- sub$y_test
    y_pred <- sub$y_pred
    x_test_protein <- sub$X_test_protein
    
    model_r  <- cor(y_test, y_pred, use = "complete.obs")
    model_r2 <- model_r^2
    
    has_olink <- protein_name %in% olink_cols
    
    if (has_olink) {
        y_train         <- MESA_Soma[as.character(train_index), protein_name]
        X_train_protein <- MESA_Olink[as.character(train_index), protein_name]
        
        pretrain_r  <- cor(y_train, X_train_protein, use = "complete.obs")
        pretrain_r2 <- pretrain_r^2
        
        posttrain_r  <- cor(y_pred, x_test_protein, use = "complete.obs")
        posttrain_r2 <- posttrain_r^2
        
        testset_r  <- cor(y_test, x_test_protein, use = "complete.obs")
        testset_r2 <- testset_r^2
    } else {
        pretrain_r <- pretrain_r2 <- posttrain_r <- posttrain_r2 <- testset_r <- testset_r2 <- NA
    }
    
    results_list[[protein_name]] <- data.frame(
        protein = protein_name,
        has_olink = has_olink,
        model_r = model_r, model_r2 = model_r2,
        pretrain_r = pretrain_r, pretrain_rsq = pretrain_r2,
        posttrain_r = posttrain_r, posttrain_rsq = posttrain_r2,
        testset_r = testset_r, testset_r2 = testset_r2
    )
    cat("done evaluating", protein_name, "\n")
}

results_df <- bind_rows(results_list)



write.table(results_df, "Results/metrics/ENET_log2scaled_3kolink_5ksoma_results.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)




##### 04 - pivot to wide (n x p matrices) and save #####
y_test_wide <- raw_test_values_os_df %>%
  select(id, protein, y_test) %>%
  pivot_wider(names_from = protein, values_from = y_test)

y_pred_wide <- raw_test_values_os_df %>%
  select(id, protein, y_pred) %>%
  pivot_wider(names_from = protein, values_from = y_pred)

x_test_protein_wide <- raw_test_values_os_df %>%
  select(id, protein, X_test_protein) %>%
  pivot_wider(names_from = protein, values_from = X_test_protein)

cat(dim(y_test_wide), "\n")
cat(dim(y_pred_wide), "\n")
cat(dim(x_test_protein_wide), "\n")

write.table(y_test_wide, gzfile("Results/mesa_test_bootstrap_inputs_enet/mesa_3kolink_pred_5ksoma_y_test_wide.tsv.gz"),
       sep = "\t", quote = F, row.names = F)
write.table(y_pred_wide, gzfile("Results/mesa_test_bootstrap_inputs_enet/mesa_3kolink_pred_5ksoma_y_pred_wide.tsv.gz"),
       sep = "\t", row.names = F)
write.table(x_test_protein_wide, gzfile("Results/mesa_test_bootstrap_inputs_enet/mesa_3kolink_pred_5ksoma_x_test_protein_wide.tsv.gz"),
       sep = "\t", row.names = F)
