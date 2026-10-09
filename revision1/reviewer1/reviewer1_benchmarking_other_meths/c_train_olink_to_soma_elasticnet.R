##### 01 - library #####
library(glmnet)
library(data.table)
library(dplyr)
library(glue)

gc()
rm(list=ls())

setwd("...")

##### 02 - parameters #####
args <- commandArgs(trailingOnly = TRUE)
batch_number <- if (length(args) >= 1) as.integer(args[1]) else 0



##### 03 - files #####
MESA_Olink <- fread("Data/final_datasets/MESA_Olink_3k_log2_nodup_rm10na_knn_final.tsv.gz", data.table = FALSE)
MESA_Soma  <- fread("Data/final_datasets/MESA_SOMA_5k_log2scaled_nodup_final.tsv.gz", data.table = FALSE) %>%
  select(-bmi)


## overlapping subject ids
overlap_ids <- intersect(unique(MESA_Olink$id), unique(MESA_Soma$id))

MESA_Olink <- MESA_Olink %>% filter(id %in% overlap_ids) %>% 
    arrange(id)
MESA_Soma  <- MESA_Soma  %>% filter(id %in% overlap_ids) %>% 
    arrange(id)

stopifnot(all(MESA_Olink$id == MESA_Soma$id))  # fools check to see if the IDs are in the same sequence


### removing index similar to python one
rownames(MESA_Olink) <- MESA_Olink$id
rownames(MESA_Soma)  <- MESA_Soma$id

MESA_Olink <- MESA_Olink %>%
    select(-id)
MESA_Soma <- MESA_Soma %>%
    select(-id)

## SAME pickled splits used for HGB/LightGBM, exported once via z_export_pickles_for_R.py,
## so results are directly comparable across all three models
train_index <- fread("Data/helper_data/train_index.csv", data.table = FALSE)$id
test_index  <- fread("Data/helper_data/test_index.csv", data.table = FALSE)$id

## SAME chunking scheme as HGB/LightGBM batches
protlist <- fread(glue("Data/helper_data/os_protlist_chunk{batch_number}.csv"), data.table = FALSE)$protein

### testing section for 10 proteins
# protlist <- protlist[1:10]

cat("Batch", batch_number, "-", length(protlist), "proteins\n")

##### 04 - training loop (elastic net, Olink -> SomaScan) #####
## X = Olink (feature set), target = SomaScan protein
X_all <- as.matrix(MESA_Olink)


train_soma_from_olink_glmnet <- function(protein_name) {
    y <- MESA_Soma[[protein_name]]
    names(y) <- rownames(MESA_Soma)
    
    ## fools proof, but still there wont be any NA in Soma or Olink
    if (any(is.na(y[as.character(train_index)]))) {
        cat("skip protein", protein_name, "due to NA in SOMA data\n")
        return(invisible(NULL))
    }
    
    X_train <- X_all[as.character(train_index), ,drop = FALSE]
    y_train <- y[as.character(train_index)]
    
    ## elastic net (alpha = 0.5), lambda chosen via 10-fold CV -- default settings,
    ## no extensive tuning, matching the "reasonable defaults" framing used for HGB/LightGBM
    fit <- cv.glmnet(X_train, y_train, alpha = 0.5, nfolds = 5)
    
    out_dir <- "Data/soma_trained_models/olink_soma/3kolink_5k_log2scaled_soma_enet"
    model_path <- file.path(out_dir, glue("{protein_name}_enet.rds"))
    saveRDS(fit, model_path)
    cat("Saved model for", protein_name, "to", model_path, "\n")
}

##### 05 - run batch #####
t0 <- Sys.time()
for (p in protlist) {
    train_soma_from_olink_glmnet(p)
}

cat("Batch", batch_number, "elapsed:", format(Sys.time() - t0), "\n")





