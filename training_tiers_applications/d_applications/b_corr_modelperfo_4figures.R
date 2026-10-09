##### 01 - library #####
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(ggpubr)

gc()
rm(list = ls())

setwd("...")


##### 02 - files #####
# raw files for both direction
soma_olink_raw <- fread("Results/metrics/HGB_log2scaled_5ksoma_3kolink_results.tsv", data.table = FALSE)
soma_olink_raw <- soma_olink_raw %>%
    filter(has_soma) %>%
    filter(complete.cases(.)) %>%
    mutate(idx = row_number(pretrain_rsq)) %>%
    mutate(model = "SomaScan to Olink") %>%
    mutate(tiers = case_when(
        pretrain_r >= 0.7 & posttrain_r >= 0.7 ~ "Tier 1",
        pretrain_r < 0.1 & posttrain_r < 0.1 ~ "Tier 4",
        posttrain_r >= 0.4 ~ "Tier 2",
        posttrain_r < 0.4 ~ "Tier 3",
        TRUE ~ NA
    ))


olink_soma_raw <- fread("Results/metrics/HGB_3kolink_log2scaled_5ksoma_results.tsv", data.table = F)
olink_soma_raw <- olink_soma_raw %>%
    filter(has_soma) %>%
    filter(complete.cases(.)) %>%
    mutate(idx = row_number(pretrain_rsq)) %>%
    mutate(model = "Olink to SomaScan") %>%
    mutate(tiers = case_when(
        pretrain_r >= 0.7 & posttrain_r >= 0.7 ~ "Tier 1",
        pretrain_r < 0.1 & posttrain_r < 0.1 ~ "Tier 4",
        posttrain_r >= 0.4 ~ "Tier 2",
        posttrain_r < 0.4 ~ "Tier 3",
        TRUE ~ NA
    ))

tier_cols <- c("Tier 1"="#26e3fc", "Tier 2"="#74c365", "Tier 3"="#fc9c26", "Tier 4"="#ff5959")




##### 03 - plots #####
### baseline and post-imputation plot
so_prepost_r <- ggplot(soma_olink_raw, aes(x = pretrain_r, y = posttrain_r, color = tiers)) +
    geom_point(alpha = 0.7, size = 1.4) +
    scale_color_manual(values = tier_cols)+
    labs(
        x = "Baseline Pearson r between SomaScan & Olink",
        y = "Post-imputation Pearson r between \nSomaScan & SomaScan imputed-Olink",
        color = "Tiers"
    ) +
    geom_abline(intercept = 0, slope = 1, color = "black", size = 1, linetype = "dotted") +
    theme_classic() +
    theme(
        plot.title = element_text(size=9),
        plot.caption = element_text(hjust = 0, size = 7),
        axis.text = element_text(size=12, colour = "black"),
        axis.title = element_text(size = 13, face="bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13), 
        legend.key.width = unit(1.2, "cm")
    ) +
    guides(color = guide_legend(
        override.aes = list(size = 4)
    ))

so_prepost_r


os_prepost_r <- ggplot(olink_soma_raw, aes(x = pretrain_r, y = posttrain_r, color = tiers)) +
    geom_point(alpha = 0.7, size = 1.4) +
    # scale_y_continuous(breaks = c(-0.4, -0.2, seq(0,1,0.1))) +
    labs(
        x = "Baseline Pearson r between Olink & SomaScan",
        y = "Post-imputation Pearson r between \nOlink & Olink imputed-SomaScan",
        color = "Tiers"
    ) +
    scale_color_manual(values = tier_cols)+
    geom_abline(intercept = 0, slope = 1, color = "black", size = 1, linetype = "dotted") +
    theme_classic() +
    theme(
        plot.title = element_text(size=9),
        plot.caption = element_text(hjust = 0, size = 7),
        axis.text = element_text(size=12, colour = "black"),
        axis.title = element_text(size = 13, face="bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13), 
        legend.key.width = unit(1.2, "cm")
    ) +
    guides(color = guide_legend(
        override.aes = list(size = 4)
    ))

os_prepost_r



### model performance plots
so_modelper_r <- ggplot(soma_olink_raw, aes(x = idx, y = model_r)) +
    geom_point(alpha = 0.6, size = 1.5, aes(color = tiers)) +
    scale_color_manual(values = tier_cols)+
    scale_y_continuous(breaks = c(-0.4, -0.2, seq(0,1,0.1))) +
    labs(
        x = "Proteins (Ranked Index)",
        y = "Model Performance (Pearson r)") +
    geom_line(aes(x = idx, y = pretrain_r), color="black", size = 0.5, show.legend = F) +
    theme_classic() + 
    theme(
        axis.text = element_text(size=12, colour = "black"),
        axis.title = element_text(size = 13, face="bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13), 
        legend.key.width = unit(1.2, "cm")
    ) +
    guides(color = guide_legend(
        override.aes = list(size = 4)
    ))

so_modelper_r


os_modelper_r <- ggplot(olink_soma_raw, aes(x = idx, y = model_r)) +
    geom_point(alpha = 0.6, size = 1.5, aes(color = tiers)) +
    scale_color_manual(values = tier_cols)+
    scale_y_continuous(breaks = c(-0.4, -0.2, seq(0,1,0.1))) +
    labs(
        x = "Proteins (Ranked Index)",
        y = "Model Performance (Pearson r)") +
    geom_line(aes(x = idx, y = pretrain_r), color="black", size = 0.5, show.legend = FALSE) +
    theme_classic() +
    theme(
        plot.title = element_text(size=9),
        plot.caption = element_text(hjust = 0, size = 7),
        axis.text = element_text(size=12, colour = "black"),
        axis.title = element_text(size = 13, face="bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13), 
        legend.key.width = unit(1.2, "cm")
    ) +
    guides(color = guide_legend(
        override.aes = list(size = 4)
    ))

os_modelper_r


##### 04 - combine and output #####
combined_plot <- ggarrange(so_prepost_r, os_prepost_r,
                           so_modelper_r, os_modelper_r, 
                           labels = "AUTO",
                           common.legend = T, legend = "bottom", align="hv", nrow=2, ncol=2)
combined_plot

ggsave("Results/figures/so_os_preppost_r_modelpre_r.png", combined_plot, dpi = 300, width = 13, height = 14, bg="white")



