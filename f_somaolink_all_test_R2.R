##### 01 - library #####
library(dplyr)
library(data.table)
library(ggplot2)
library(colorspace)
library(ggrepel)
library(ggh4x)
library(glue)

gc()
rm(list = ls())

setwd("...")



##### 02 - files #####
somaolink_test_r2 <- fread("Results/metrics/HGB_log2scaled_5ksoma_3kolink_results.tsv", data.table = F)
olinksoma_test_r2 <- fread("Results/metrics/HGB_3kolink_log2scaled_5ksoma_results.tsv", data.table = F)



### processing
tier_cols <- c("\u22650.7"="#26e3fc", "[0.4-0.7)"="#74c365", "[0.1-0.4)"="#fc9c26", "<0.1"="#ff5959")
binning_somaolink <- somaolink_test_r2 %>%
    mutate(
        r_bins = case_when(
            model_r >= 0.7 ~ "\u22650.7",
            model_r < 0.1 ~ "<0.1",
            model_r >= 0.1 & model_r <0.4 ~ "[0.1-0.4)",
            model_r >= 0.4 & model_r <0.7 ~ "[0.4-0.7)",
        ),
        r_bins = factor(r_bins, levels = c("\u22650.7", "[0.4-0.7)", "[0.1-0.4)", "<0.1")),
        has_soma = ifelse(has_soma, "Overlapping", "No Overlapping")
    ) %>%
    count(has_soma, r_bins, name = "n") %>%
    group_by(has_soma) %>%
    mutate(prop = n / sum(n),
           labels = "Soma\u2192Olink") %>%
    ungroup()


binning_olinksoma <- olinksoma_test_r2 %>%
    mutate(
        r_bins = case_when(
            model_r >= 0.7 ~ "\u22650.7",
            model_r < 0.1 ~ "<0.1",
            model_r >= 0.1 & model_r <0.4 ~ "[0.1-0.4)",
            model_r >= 0.4 & model_r <0.7 ~ "[0.4-0.7)",
        ),
        r_bins = factor(r_bins, levels = c("\u22650.7", "[0.4-0.7)", "[0.1-0.4)", "<0.1")),
        has_soma = ifelse(has_soma, "Overlapping", "No Overlapping")
    ) %>%
    count(has_soma, r_bins, name = "n") %>%
    group_by(has_soma) %>%
    mutate(prop = n / sum(n),
           labels = "Olink\u2192Soma") %>%
    ungroup()



binning_full <- rbind(binning_somaolink, binning_olinksoma) %>%
    mutate(labels =factor(labels, levels = c("Soma\u2192Olink", "Olink\u2192Soma")))




ggplot(binning_full, aes(x = has_soma, y = prop, fill = r_bins)) +
    geom_col(position = "fill") +
    geom_text(
        aes(label = glue("{scales::comma(n)} ({scales::percent(prop, accuracy = 0.05)})")),
        position = position_fill(vjust = 0.5),
        color = "black", 
        size = 3
    ) +
    facet_nested(.~labels)+
    scale_fill_manual(values = tier_cols) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, 
         y = "Proportion", 
         fill = "Pearson \nCorrelation"
         # title = "Distribution of Correlation (Pearson R) in Overlapping & Non-Overlapping Proteins"
         ) +
    theme_classic() +
    theme(
        title = element_text(size = 8)
    )



ggsave("Results/figures/testset_correlation_proteins.png", dpi = 300, width = 7, height = 5, bg = "white")








# density plots -----------------------------------------------------------
r2_full <- rbind(somaolink_test_r2 %>% mutate(labels = "Soma\u2192Olink"), 
                 olinksoma_test_r2 %>% mutate(labels = "Olink\u2192Soma")) %>%
    mutate(labels =factor(labels, levels = c("Soma\u2192Olink", "Olink\u2192Soma")),
           has_soma = ifelse(has_soma, "Overlapping", "No Overlapping"))


ggplot(r2_full, aes(x = model_r)) +
    geom_density(fill = "skyblue", alpha = 0.7) +
    geom_vline(aes(xintercept = 0.4), 
               linetype = "dashed", color = "red") +
    facet_nested(.~labels + has_soma) + 
    labs(
        # title = "Distribution of Pearson R in Test set",
         x = "Pearson R Value",
         y = "Density") +
    theme_classic()


ggsave("Results/figures/somaolink_olinksoma_testset_pearsonr_density.png", dpi = 300, width = 8, height = 4, bg = "white")


# #### hex plots
# ggplot(r2_full %>% filter(has_soma == "Overlapping"), aes(x = r, y = correlation)) +
#     geom_hex(bins = 50) +
#     facet_nested_wrap(labels ~ ., axes = "x", ncol = 1)+
#     scale_y_continuous(labels = scales::percent) +
#     scale_fill_continuous(type = "viridis") + # Use a good color scale for count data
#     labs(
#         x = "Pearson R in Test set (Pred-Actual)",
#         y = "Correlation between Soma Olink in Full set", 
#         fill = "Count") +
#     theme_minimal()
# 
# 
# 
# ggsave("Results/figures/somaolink_olinksoma_testset_pearsonr_overallcorr_hex.png", dpi = 300, width = 5, height = 9, bg = "white")





