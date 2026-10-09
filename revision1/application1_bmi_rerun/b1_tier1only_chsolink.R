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

chs_raw <- fread("Results/CHS_pred/HGB_predicted_CHS_log2scaled_5ksoma_3kolink.tsv.gz", data.table = F) %>%
    select(id, all_of(tier1_prots))

chs_covs <- fread("/medpop/esp2/zyu/proteomics_pathway/CHS_soma5k_bmiwaisttohipratiolipids.txt.gz", data.table = F,
                  select=c("idno", "age", "sex", "race", "bmi")) %>%
    rename(id = idno)

protein_names <- setdiff(colnames(chs_raw), c("id"))


# CRP and BMI -------------------------------------------------------------
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

    ### FIX (Reviewer 1 - leakage): split BEFORE any model fitting or hyperparameter search
    set.seed(613)
    data_split <- initial_split(chs_analysis, prop = 0.8)
    train_data <- training(data_split)
    test_data <- testing(data_split)

    y_train <- train_data[[pheno]]
    x_train <- as.matrix(train_data %>% select(all_of(protein_names), all_of(unpenalized_covs)))

    y_test <- test_data[[pheno]]
    x_test <- as.matrix(test_data %>% select(all_of(protein_names), all_of(unpenalized_covs)))

    n_train <- nrow(x_train)
    n_test <- nrow(x_test)
    n_predictors <- ncol(x_train)

    ## penalties
    penalty_factors <- rep(1, ncol(x_train))
    penalty_factors[colnames(x_train) %in% unpenalized_covs] <- 0

    ## determining alpha -- FIX: now uses x_train/y_train only, not x_full/y_full
    alpha_list <- seq(0, 1, by = 0.1)
    cv_models <- list()
    set.seed(613)

    for (i in seq_along(alpha_list)) {
        alpha_val <- alpha_list[i]
        cat(paste("Running CV for alpha =", alpha_val, "\n"))

        cv_fit <- cv.glmnet(
            x_train, y_train,
            alpha = alpha_val,
            penalty.factor = penalty_factors,
            standardize = TRUE
        )
        cv_models[[i]] <- cv_fit
    }

    # Select best alpha & lambda
    cv_errors <- sapply(cv_models, function(m) min(m$cvm))
    best_alpha_idx <- which.min(cv_errors)
    best_alpha <- alpha_list[best_alpha_idx]
    best_cv_model <- cv_models[[best_alpha_idx]]
    optimal_lambda <- best_cv_model$lambda.min

    cat("Best alpha:", best_alpha, "\n")
    cat("Optimal lambda:", optimal_lambda, "\n")

    ## ADD: leakage-free cross-validated R2, computed entirely within the training fold
    cv_mse_at_best <- min(best_cv_model$cvm)
    cv_r2 <- 1 - cv_mse_at_best / var(y_train)
    cat("Training-fold CV R2:", cv_r2, "\n")

    ## testing
    # best_alpha = 1
    # optimal_lambda = 0.001092252

    ## final model -- refit on training data only
    final_refit <- glmnet(
        x_train, y_train,
        alpha = best_alpha,
        penalty.factor = penalty_factors,
        lambda = optimal_lambda
    )

    preds <- predict(final_refit, newx = x_test, s = optimal_lambda)
    
    # R2
    sst <- sum((y_test - mean(y_test))^2)
    sse <- sum((y_test - preds)^2)
    r2 <- 1 - sse / sst
    
    n <- nrow(x_test)
    # number of non-zero coefficients at this lambda (excluding intercept)
    p <- sum(coef(final_refit, s = optimal_lambda) != 0) - 1
    
    adjusted_rsq <- 1 - (1 - r2) * (n - 1) / (n - p - 1)
    
    
    ##save rsq -- now includes n_train, n_test, n_predictors, adjusted_rsq, cv_r2
    write.table(data.frame(pheno = pheno, n_train = n_train, n_test = n_test,
                            n_predictors = n_predictors, best_alpha = best_alpha,
                            optimal_lambda = optimal_lambda, rsq = r2,
                            adjusted_rsq = adjusted_rsq, cv_r2_training_fold = cv_r2),
                paste0("Results/revision1/chs_assoc_metrics/tier1onlysomaolink_predolink_overlapping_proteins_",pheno,"_baseline_rsq.tsv"),
                sep = "\t", row.names = F, quote = F
    )

    ## permutation feature importance -- MUST use the same R2 formula as baseline
    protein_importance <- list()
    for (prot in protein_names) {
        x_test_shuffled <- x_test
        x_test_shuffled[, prot] <- sample(x_test_shuffled[, prot])
        
        shuffled_predictions <- as.numeric(predict(final_refit, s = optimal_lambda, newx = x_test_shuffled))
        
        ## same formula as baseline_rsq -- 1 - SSE/SST
        sse_shuffled <- sum((y_test - shuffled_predictions)^2)
        shuffled_rsq <- 1 - sse_shuffled / sst   # SST is unchanged since y_test never changes
        
        protein_importance[[prot]] <- r2 - shuffled_rsq
    }

    importance_df <- data.frame(
        Protein = names(protein_importance),
        Importance = unlist(protein_importance)
    ) %>%
        arrange(desc(Importance))

    ## save feature importance
    out_file <- paste0(
        "Results/revision1/chs_assoc_metrics/tier1onlysomaolink_predolink_overlapping_proteins_", pheno, "_feature_imp.tsv"
    )

    write.table(
        importance_df,
        out_file,
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    cat("Saved results to:", out_file, "\n")
}
