##### 01 - library #####
library(rsample)
library(glmnet)
library(dplyr)
library(data.table)


gc()
rm(list = ls())

setwd("...")


##### 02 - chs phenotype #####
tier1_prots <- fread("Results/rank_by_gain/soma5k_olink3k.tsv") %>%
    filter(tiers == 1) %>%
    pull(protein)

chs_raw <- fread("Data/final_datasets/CHS_SOMA_log2scaled_nodup_final.tsv.gz", data.table = F) %>%
    select(id, all_of(tier1_prots))

chs_covs <- fread("/medpop/esp2/zyu/proteomics_pathway/CHS_soma5k_bmiwaisttohipratiolipids.txt.gz", data.table = F, 
                  select=c("idno", "age", "sex", "race", "bmi")) %>%
    rename(id = idno)


protein_names <- setdiff(colnames(chs_raw), c("id"))


# BMI -------------------------------------------------------------
chs_analysis <- merge(chs_raw, chs_covs, by = "id") %>%
    filter(complete.cases(.)) %>%
    mutate(bmi = scale(bmi),
           race = factor(race, levels = c("1","2","3","4","5")))

phenotypes <- c("bmi")
unpenalized_covs <- c("age", "sex", "race")

# Loop over phenotypes
for (pheno in phenotypes) {
    cat("\n============================\n")
    cat("Running analysis for:", pheno, "\n")
    cat("============================\n")
    
    ### extract X and y
    y_full <- chs_analysis[[pheno]]
    x_full <- as.matrix(chs_analysis %>% select(id, all_of(protein_names), all_of(unpenalized_covs)))
    
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
    data_split <- initial_split(chs_analysis, prop = 0.8)
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
                paste0("Results/chs_assoc_metrics/tier1onlysomaolink_actualsoma_overlapping_proteins_",pheno,"_baseline_rsq.tsv"), sep = "\t", row.names = F, quote = F
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
        "Results/chs_assoc_metrics/tier1onlysomaolink_actualsoma_overlapping_proteins_", pheno, "_feature_imp.tsv"
    )
    
    write.table(
        importance_df,
        out_file,
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    
    cat("Saved results to:", out_file, "\n")
}


















