##### 01 - library #####
library(rsample)
library(glmnet)
library(dplyr)
library(data.table)


gc()
rm(list = ls())

setwd("...")


##### 02 - UKB phenotype #####
no_consent <- data.frame(fread("/medpop/esp2/projects/UK_Biobank/withdrawn_samples/w7089_20250818.csv"))$V1
overlapping_prots <- fread("Results/rank_by_gain/olink3k_soma5k.tsv")$protein
ukb_raw <- fread("Results/UKB_pred/HGB_predicted_UKB_3kolink_log2scaled_5ksoma_r.tsv.gz", data.table = F) %>%
    select(id, all_of(overlapping_prots))
ukb_covs <- fread("Data/mesa_olink_processed/UKB_BMI_LDL_covs.tsv", data.table = F)
ukb_drugs <- fread("Data/mesa_olink_processed/UKB_ldl_adj.tsv", data.table = F)

## set variable names
unpenalized_covs <- c("enroll_age.y", "Sex_num", "Race_num")
protein_names <- setdiff(colnames(ukb_raw), c("id"))


## getting the final analytic dataset
ukb_analysis <- merge(ukb_raw, ukb_covs, by = "id") %>%
    merge(ukb_drugs, by = "id") %>%
    mutate(Race_num = ifelse(Race == "White",1,0),
           Sex_num = ifelse(Sex == "Female",0,1)) %>%
    select(-c(Race, Sex, LDL.direct, C.reactive.protein_NoOutlierBySex, ldl, hdl, trig, chol, 
              statin0, statin1, ldladj, choladj, trigadj, hdladj)) %>%
    filter(!id %in% no_consent) %>%
    filter(complete.cases(.)) %>%
    mutate(BMI_baseline = scale(BMI_baseline))

phenotypes <- c("BMI_baseline")


### analysis
for (pheno in phenotypes) {
    cat("\n============================\n")
    cat("Running analysis for:", pheno, "\n")
    cat("============================\n")
    
    ### extract X and y
    y_full <- ukb_analysis[[pheno]]
    x_full <- as.matrix(ukb_analysis %>% select(id, all_of(protein_names), all_of(unpenalized_covs)))
    
    ## penalties
    penalty_factors <- rep(1, ncol(x_full))
    penalty_factors[colnames(x_full) %in% unpenalized_covs] <- 0
    
    ## determining alpha
    alpha_list <- seq(0, 1, by = 0.1)
    cv_models <- list()
    set.seed(613)
    
    for (i in seq_along(alpha_list)) {
        alpha_val <- alpha_list[i]
        cat(paste("Running CV for alpha =", alpha_val, "\n"))
        
        cv_fit <- cv.glmnet(
            x_full, y_full,
            alpha = alpha_val,
            penalty.factor = penalty_factors,
            standardize = TRUE
        )
        cv_models[[i]] <- cv_fit
    }
    
    # Select best alpha & lambda
    cv_errors <- sapply(cv_models, function(m) min(m$cvm))
    best_alpha <- alpha_list[which.min(cv_errors)]
    final_model <- cv_models[[which.min(cv_errors)]]
    optimal_lambda <- final_model$lambda.min
    
    cat("Best alpha:", best_alpha, "\n")
    cat("Optimal lambda:", optimal_lambda, "\n")
    
    
    ### train test split using rsample
    set.seed(613)
    data_split <- initial_split(ukb_analysis, prop = 0.8)
    train_data <- training(data_split)
    test_data <- testing(data_split)
    
    y_train <- train_data[[pheno]]
    x_train <- as.matrix(train_data %>% select(id, all_of(protein_names), all_of(unpenalized_covs)))
    
    y_test <- test_data[[pheno]]
    x_test <- as.matrix(test_data %>% select(id, all_of(protein_names), all_of(unpenalized_covs)))
    
    ## final model
    final_refit <- glmnet(
        x_train, y_train,
        alpha = best_alpha,
        penalty.factor = penalty_factors
    )
    
    ## baseline performance
    baseline_predictions <- predict(final_refit, s = optimal_lambda, newx = x_test)
    baseline_rsq <- cor(y_test, baseline_predictions)^2
    
    baseline_rsq
    ##save rsq
    write.table(data.frame(pheno = pheno, rsq = baseline_rsq), 
                paste0("Results/ukb_assoc_metrics/proteins_",pheno,"_baseline_rsq.tsv"), sep = "\t", row.names = F, quote = F
    )
    
    ## permutation feature importance
    protein_importance <- list()
    for (prot in protein_names) {
        x_test_shuffled <- x_test
        x_test_shuffled[, prot] <- sample(x_test_shuffled[, prot])
        
        shuffled_predictions <- predict(final_refit, s = optimal_lambda, newx = x_test_shuffled)
        shuffled_rsq <- cor(y_test, shuffled_predictions)^2
        
        protein_importance[[prot]] <- baseline_rsq - shuffled_rsq
    }
    
    importance_df <- data.frame(
        Protein = names(protein_importance),
        Importance = unlist(protein_importance)
    ) %>%
        arrange(desc(Importance))
    
    ## save feature importance
    out_file <- paste0(
        "Results/ukb_assoc_metrics/proteins_", pheno, "_feature_imp.tsv"
    )
    
    write.table(
        importance_df,
        out_file,
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    
    cat("Saved results to:", out_file, "\n")
}
