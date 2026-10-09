##### 01 - library #####
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(gridExtra)

gc()
rm(list=ls())

setwd("...")


##### 02 - files #####
soma_olink <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F)
olink_soma <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = F) %>%
    mutate(gain = posttrain_r - pretrain_r)

cross_dir <- soma_olink %>%
    select(protein, r_SO = posttrain_r) %>%
    inner_join(
        olink_soma %>% select(protein, r_OS = posttrain_r),
        by = "protein"
    )

# correlation between the two directions' post-imputation performance, per protein
cor_val_p <- cor(cross_dir$r_SO, cross_dir$r_OS, use = "complete.obs", method = "pearson")
cor_val_s <- cor(cross_dir$r_SO, cross_dir$r_OS, use = "complete.obs", method = "spearman")

p <- ggplot(cross_dir, aes(x = r_SO, y = r_OS)) +
    geom_point(alpha = 0.4, size = 1.4) +
    geom_smooth(method = "lm", se = TRUE, color = "#d62728") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    annotate("text", x = min(cross_dir$r_SO, na.rm = TRUE), y = max(cross_dir$r_OS, na.rm = TRUE),
             label = paste0("pearson r = ", round(cor_val_p, 2)), hjust = 0, vjust = 1, size = 4) +
    labs(x = "Post-imputation r, SomaScan → Olink (SO)",
         y = "Post-imputation r, Olink → SomaScan (OS)") +
    theme_minimal(base_size = 12)

p





