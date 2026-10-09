##### 01 - library #####
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(patchwork)

gc()
rm(list=ls())

setwd("...")


##### 02 - Panel A: pre- vs post-imputation correlation #####
soma_olink <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F)
olink_soma <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = F) %>%
    mutate(gain = posttrain_r - pretrain_r)

soma_olink2 <- soma_olink %>% mutate(direction = "SomaScan to Olink")
olink_soma2 <- olink_soma %>% mutate(direction = "Olink to SomaScan")

combined <- bind_rows(
    soma_olink2 %>% select(protein, direction, tiers, baseline = pretrain_r, post_imputation = posttrain_r),
    olink_soma2 %>% select(protein, direction, tiers, baseline = pretrain_r, post_imputation = posttrain_r)
) %>%
    mutate(direction = factor(direction, levels = c("SomaScan to Olink", "Olink to SomaScan")),
           tiers = factor(tiers, levels = c(1,2,3,4), labels = c("Tier 1","Tier 2","Tier 3","Tier 4")))


long_df <- combined %>%
    pivot_longer(cols = c(baseline, post_imputation), names_to = "stage", values_to = "r") %>%
    mutate(stage = factor(stage, levels = c("baseline","post_imputation"),
                          labels = c("Baseline","Post-imputation")))


p_A <- ggplot(long_df, aes(x = r, fill = stage)) +
    geom_density(alpha = 0.55, position = "identity", color = NA) +
    facet_wrap(~direction, ncol = 1) +
    scale_fill_manual(values = c("Baseline" = "#8c8c8c", "Post-imputation" = "#2166ac")) +
    labs(x = "Cross-platform correlation (Pearson R)", y = "Number of proteins", fill = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top", strip.text = element_text(face = "bold"))


p_A



##### 03 - Panel B: model test-set performance #####
somaolink_test_r2 <- fread("Results/metrics/HGB_log2scaled_5ksoma_3kolink_results.tsv", data.table = F)
olinksoma_test_r2 <- fread("Results/metrics/HGB_3kolink_log2scaled_5ksoma_results.tsv", data.table = F)

##### merge tier assignments into model performance data #####
soma_olink_tiers <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F) %>%
    select(protein, tiers) %>%
    mutate(labels = "SomaScan to Olink")
olink_soma_tiers <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = F) %>%
    select(protein, tiers) %>%
    mutate(labels = "Olink to SomaScan")

tier_lookup <- bind_rows(soma_olink_tiers, olink_soma_tiers) %>%
    mutate(tiers = factor(tiers, levels = c(1,2,3,4),
                          labels = c("Tier 1","Tier 2","Tier 3","Tier 4")))


r2_full <- rbind(
    somaolink_test_r2 %>% mutate(labels = "SomaScan to Olink"),
    olinksoma_test_r2 %>% mutate(labels = "Olink to SomaScan")
) %>%
    mutate(labels = factor(labels, levels = c("SomaScan to Olink", "Olink to SomaScan")),
           has_soma = ifelse(has_soma, "Overlapping", "No Overlapping"))

r2_full_tiers <- r2_full %>%
    left_join(tier_lookup, by = c("protein", "labels"))
# proteins with has_soma == "No Overlapping" will have tiers == NA, as expected

##### Panel B: split by overlap status, tier-colored where applicable #####
tier_cols <- c("Tier 1" = "#6beeff", "Tier 2" = "#93db86",
               "Tier 3" = "#ffb963", "Tier 4" = "#ff7878")


p_B <- ggplot() +
    # No Overlapping: single density, matches Panel A's blue
    geom_density(data = r2_full_tiers %>% filter(has_soma == "No Overlapping"),
                 aes(x = model_r), fill = "#2166ac", color = NA, alpha = 0.7) +
    # Overlapping: colored by tier
    geom_density(data = r2_full_tiers %>% filter(has_soma == "Overlapping"),
                 aes(x = model_r, fill = tiers), color=NA, alpha = 0.7) +
    scale_fill_manual(values = tier_cols, name = "Tier", na.translate = FALSE) +
    scale_color_manual(values = tier_cols, guide = "none") +
    facet_nested(. ~ labels + has_soma) +
    labs(x = "Model Performance (Pearson R)", y = "Density") +
    theme_classic(base_size = 12) +
    theme(legend.position = "top")

p_B



##### 05 - Panel C: SO vs OS post-imputation correlation (overlapping proteins) #####
scatter_df <- inner_join(
    soma_olink %>% select(protein, tiers, SO_r = posttrain_r),
    olink_soma %>% select(protein, tiers, OS_r = posttrain_r),
    by = "protein",
    suffix = c("_SO", "_OS")
) %>%
    mutate(tiers_SO = factor(tiers_SO, levels = c(1,2,3,4), labels = c("Tier 1","Tier 2","Tier 3","Tier 4")))

# concordance annotation: Pearson r between the two directions' post-imputation r's
# yes, a correlation of correlations :(
concordance <- cor.test(scatter_df$SO_r, scatter_df$OS_r)
concordance_label <- sprintf("r = %.2f (n = %d)", concordance$estimate, nrow(scatter_df))


p_C <- ggplot(scatter_df, aes(x = SO_r, y = OS_r, color = tiers_SO)) +
    geom_point(alpha = 0.7, size = 1.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey30", linewidth=1) +
    geom_smooth(mapping = aes(x = SO_r, y = OS_r), data = scatter_df,
                method = "lm", se = FALSE, color = "red",
                linewidth = 0.9, inherit.aes = FALSE) +
    annotate("text", x = 0.02, y = 0.98, label = concordance_label,
             hjust = 0, vjust = 1, size = 3.8, color = "black") +
    scale_color_manual(values = tier_cols, name = "Tier", guide = "none") +
    labs(x = "Post-imputation correlation, SomaScan to Olink",
         y = "Post-imputation correlation, Olink to SomaScan") +
    coord_equal(xlim = c(0,1), ylim = c(0,1)) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top")

p_C



##### 06 - merge #####
ab_col <- p_A / p_B +
    plot_layout(heights = c(1, 0.8))

merged_fig <- (ab_col | p_C) +
    plot_layout(widths = c(1, 0.7), guides = "collect") +
    plot_annotation(tag_levels = 'A') &
    theme(legend.position = "bottom", plot.tag = element_text(face = "bold"))

merged_fig

ggsave("Results/revision1/figures/pre_post_model_perf_density_plot.png",
       merged_fig, width = 13, height = 9, dpi = 300, bg = "white")



