##### 01 - library #####
library(data.table)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(patchwork)

gc()
rm(list=ls())


setwd("/medpop/esp2/lli/somalink_olink/")


##### 02 - files #####
soma_olink_orig <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = FALSE)
olink_soma_orig <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = FALSE)

tier_probs_so <- fread("Results/mesa_test_bootstrap_tier_probs/so_1000bootstraps.tsv", data.table = F)
tier_probs_os <- fread("Results/mesa_test_bootstrap_tier_probs/os_1000bootstraps.tsv", data.table = F)

tier_probs_all <- bind_rows(tier_probs_so, tier_probs_os)


## build the original tier lookup (as before)
original_tiers <- bind_rows(
    soma_olink_orig %>% select(protein, tiers) %>% mutate(direction = "SomaScan to Olink"),
    olink_soma_orig %>% select(protein, tiers) %>% mutate(direction = "Olink to SomaScan")
) %>%
    rename(original_tier = tiers)


## for each protein, find the highest probability among NON-original tiers
stability_check <- tier_probs_all %>%
    inner_join(original_tiers, by = c("protein", "direction")) %>%
    group_by(protein, direction, original_tier) %>%
    summarise(
        max_alt_tier_prob = max(probability[tier != original_tier], default = 0),
        .groups = "drop"
    ) %>%
    mutate(
        stability = ifelse(max_alt_tier_prob >= 0.05, "Unstable (≥5% alt. tier)", "Stable"),
        original_tier_label = factor(original_tier, levels = 1:4, labels = paste("Tier", 1:4))
    )

# SOME ARE WILDLY OFF OMG





##### 03 - plotting #####
stability_summary <- stability_check %>%
    count(direction, original_tier_label, stability) %>%
    group_by(direction, original_tier_label) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    mutate(direction = factor(direction, levels = c("SomaScan to Olink", "Olink to SomaScan")))


# first plot: stability of the tier assignment ----------------------------
stability_labels <- stability_summary %>%
    arrange(direction, original_tier_label, desc(stability)) %>%
    group_by(direction, original_tier_label) %>%
    mutate(
        cum_top = cumsum(prop),
        cum_bottom = cum_top - prop,
        y_mid = (cum_top + cum_bottom) / 2,
        ## push label inward: never closer than 5% to either edge
        y_mid_adj = pmin(pmax(y_mid, 0.05), 0.95)
    ) %>%
    ungroup()

p_stability <- ggplot(stability_summary, aes(x = original_tier_label, y = prop, fill = stability)) +
    geom_col(position = "stack", width = 0.6) +
    geom_text(
        data = stability_labels,
        aes(x = original_tier_label, y = y_mid_adj,
            label = paste0(n, "\n(", scales::percent(prop, accuracy = 1), ")")),
        inherit.aes = FALSE, size = 2.8, color = "white", fontface = "bold"
    ) +
    facet_wrap(~direction) +
    scale_fill_manual(values = c("Stable" = "#2166ac", "Unstable (≥5% alt. tier)" = "#d62728")) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Proportion of proteins", fill = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom")

p_stability









# second plot: the fluctuation in pearson correlation ---------------------
## combine point estimate + CI + original tier, one row per protein/direction
ci_data <- tier_probs_all %>%
    distinct(protein, direction, posttrain_r_ci_low, posttrain_r_ci_high) %>%
    inner_join(original_tiers, by = c("protein", "direction")) %>%
    ## merging the original correlation for reference
    left_join(
        bind_rows(
            soma_olink_orig %>% select(protein, posttrain_r) %>% mutate(direction = "SomaScan to Olink"),
            olink_soma_orig %>% select(protein, posttrain_r) %>% mutate(direction = "Olink to SomaScan")
        ),
        by = c("protein", "direction")
    ) %>%
    mutate(original_tier_label = factor(original_tier, levels = 1:4, labels = paste("Tier", 1:4))) %>%
    group_by(direction) %>%
    arrange(posttrain_r) %>%
    mutate(rank_order = row_number()) %>%
    ungroup() %>%
    mutate(direction = factor(direction, levels = c("SomaScan to Olink", "Olink to SomaScan")))

p_ci_ribbon <- ggplot(ci_data, aes(x = rank_order, y = posttrain_r, color = original_tier_label,
                                   fill = original_tier_label)) +
    geom_ribbon(aes(ymin = posttrain_r_ci_low, ymax = posttrain_r_ci_high), alpha = 0.25, color = NA) +
    geom_line(linewidth = 0.5) +
    geom_hline(yintercept = c(0.1, 0.4, 0.7), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    facet_wrap(~direction) +
    scale_color_manual(values = c("Tier 1" = "#6beeff", "Tier 2" = "#93db86",
                                  "Tier 3" = "#ffb963", "Tier 4" = "#ff7878")) +
    scale_fill_manual(values = c("Tier 1" = "#6beeff", "Tier 2" = "#93db86",
                                 "Tier 3" = "#ffb963", "Tier 4" = "#ff7878")) +
    labs(x = "Protein (ranked by point-estimate r)", y = "Post-imputation Pearson r",
         color = "Original tier", fill = "Original tier") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom")

p_ci_ribbon



##### 05 - save plot objects for downstream combined display item #####
# dir.create("Results/revision1/figures/rds", recursive = TRUE, showWarnings = FALSE)
saveRDS(p_stability, "Results/revision1/figures/rds/p_stability.rds")
saveRDS(p_ci_ribbon, "Results/revision1/figures/rds/p_ci_ribbon.rds")




##### 04 - merging two plots #####
combined_stability_plot <- p_stability / p_ci_ribbon +
    plot_layout(heights = c(1, 0.9), guides = "collect") + plot_annotation(tag_levels = 'A') &
    theme(legend.position = "bottom")

combined_stability_plot

ggsave("Results/revision1/figures/tier_stability_combined.png", combined_stability_plot,
       dpi = 300, width = 10, height = 11, bg = "white")




