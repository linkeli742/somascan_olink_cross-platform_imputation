##### 01 - library #####
library(dplyr)
library(purrr)
library(data.table)

gc()
rm(list=ls())

setwd("...")




##### 02 - files #####
overlapping_prots <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F)$protein

so_corr_ll <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = FALSE)
os_corr_ll <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = FALSE)

y_pred_so <- fread("Results/mesa_test_bootstrap_inputs/mesa_5ksoma_pred_3kolink_y_pred_wide.tsv.gz", data.table = FALSE) %>%
    select(id, all_of(overlapping_prots))
x_test_so <- fread("Results/mesa_test_bootstrap_inputs/mesa_5ksoma_pred_3kolink_x_test_protein_wide.tsv.gz", data.table = FALSE) %>%
    select(id, all_of(overlapping_prots))
y_pred_os <- fread("Results/mesa_test_bootstrap_inputs/mesa_3kolink_pred_5ksoma_y_pred_wide.tsv.gz", data.table = FALSE) %>%
    select(id, all_of(overlapping_prots))
x_test_os <- fread("Results/mesa_test_bootstrap_inputs/mesa_3kolink_pred_5ksoma_x_test_protein_wide.tsv.gz", data.table = FALSE) %>%
    select(id, all_of(overlapping_prots))


##### 03 - analysis #####
compute_spearman <- function(y_pred, x_test) {
    proteins <- setdiff(intersect(colnames(y_pred), colnames(x_test)), "id")
    map_dfr(proteins, function(p) {
        tibble(
            protein = p,
            posttrain_rho = cor(y_pred[[p]], x_test[[p]], method = "spearman", use = "complete.obs")
        )
    })
}

spearman_so <- compute_spearman(y_pred_so, x_test_so)
spearman_os <- compute_spearman(y_pred_os, x_test_os)


## merge with existing Pearson-based tiers for direct comparison
so_comparison <- so_corr_ll %>%
    select(protein, tiers, posttrain_r) %>%
    inner_join(spearman_so, by = "protein")

os_comparison <- os_corr_ll %>%
    select(protein, tiers, posttrain_r) %>%
    inner_join(spearman_os, by = "protein")

## tier-level summary: does Spearman tell the same story as Pearson?
so_comparison_summary <- so_comparison %>%
    group_by(tiers) %>%
    summarise(median_pearson = median(posttrain_r, na.rm = TRUE),
              median_spearman = median(posttrain_rho, na.rm = TRUE),
              n = n())

os_comparison_summary <- os_comparison %>%
    group_by(tiers) %>%
    summarise(median_pearson = median(posttrain_r, na.rm = TRUE),
              median_spearman = median(posttrain_rho, na.rm = TRUE),
              n = n())

## overall agreement between the two metrics
cor(so_comparison$posttrain_r, so_comparison$posttrain_rho, use = "complete.obs")
cor(os_comparison$posttrain_r, os_comparison$posttrain_rho, use = "complete.obs")


##### 04 - plot and more outputs #####
##### combined tier-level summary, pivoted wide (Pearson/Spearman side by side per direction) #####
combined_summary <- bind_rows(
    so_comparison_summary %>% mutate(direction = "SO"),
    os_comparison_summary %>% mutate(direction = "OS")
)

combined_summary_wide <- combined_summary %>%
    pivot_wider(
        id_cols = tiers,
        names_from = direction,
        values_from = c(median_pearson, median_spearman, n),
        names_glue = "{.value}_{direction}"
    ) %>%
    relocate(ends_with("_SO"), .after = tiers)

combined_summary_wide

write.table(combined_summary_wide, "Results/revision1/spearman/pearson_spearman_tier_summary_wide.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

##### overall (protein-level) concordance between Pearson and Spearman #####
overall_concordance <- tibble(
    direction = c("SomaScan to Olink", "Olink to SomaScan"),
    pearson_spearman_concordance = c(
        cor(so_comparison$posttrain_r, so_comparison$posttrain_rho, use = "complete.obs"),
        cor(os_comparison$posttrain_r, os_comparison$posttrain_rho, use = "complete.obs")
    )
)

overall_concordance

write.table(overall_concordance, "Results/revision1/spearman/pearson_spearman_overall_concordance.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)


##### visualization: scatter of Pearson vs Spearman, colored by tier, faceted by direction #####
plot_data <- bind_rows(
    so_comparison %>% mutate(direction = "SomaScan to Olink"),
    os_comparison %>% mutate(direction = "Olink to SomaScan")
) %>%
    mutate(
        tiers = factor(tiers, levels = 1:4, labels = paste("Tier", 1:4)),
        direction = factor(direction, levels = c("SomaScan to Olink", "Olink to SomaScan"))
    )

concordance_label <- overall_concordance %>%
    mutate(
        direction = factor(direction, levels = c("SomaScan to Olink", "Olink to SomaScan")),
        label = paste0("r = ", round(pearson_spearman_concordance, 3))
    )

p_pearson_vs_spearman <- ggplot(plot_data, aes(x = posttrain_r, y = posttrain_rho, color = tiers)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_text(data = concordance_label, aes(x = -0.3, y = 0.9, label = label),
              inherit.aes = FALSE, size = 3.5) +
    facet_wrap(~direction) +
    scale_color_manual(values = c("Tier 1" = "#6beeff", "Tier 2" = "#93db86",
                                  "Tier 3" = "#ffb963", "Tier 4" = "#ff7878")) +
    labs(x = "Pearson r", y = "Spearman rho", color = "Tier",
         title = "Concordance between Pearson and Spearman correlation") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom")

p_pearson_vs_spearman

ggsave("Results/revision1/figures/pearson_vs_spearman.png", p_pearson_vs_spearman,
       dpi = 300, width = 10, height = 6, bg = "white")


##### 05 - save #####
# dir.create("Results/revision1/spearman", showWarnings = FALSE, recursive = TRUE)
write.table(so_comparison, "Results/revision1/spearman/SO_pearson_spearman_per_protein.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(os_comparison, "Results/revision1/spearman/OS_pearson_spearman_per_protein.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
















