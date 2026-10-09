##### 01 - library #####
library(dplyr)
library(data.table)
library(ggh4x)
library(ggplot2)

gc()
rm(list=ls())

setwd("...")


##### 02 - files #####
soma_olink <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F)
olink_soma <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = F)

## original tiering cutoffs
orig_cuts <- list(tier1_cut = 0.7, tier4_cut = 0.1, tier2_cut = 0.4)

# function to cut tier the proteins for sets of proteins
assign_tiers <- function(df, tier1_cut, tier4_cut, tier2_cut) {
    df %>%
        dplyr::mutate(tiers = case_when(
            pretrain_r >= tier1_cut & posttrain_r >= tier1_cut ~ 1,
            pretrain_r < tier4_cut & posttrain_r < tier4_cut ~ 4,
            posttrain_r >= tier2_cut ~ 2,
            posttrain_r < tier2_cut ~ 3,
            TRUE ~ NA_real_
        ))
}


## original tiering, kept as the reference for reclassification comparisons
original_so <- assign_tiers(soma_olink, orig_cuts$tier1_cut, orig_cuts$tier4_cut, orig_cuts$tier2_cut) %>%
    select(protein, tiers_original = tiers)
original_os <- assign_tiers(olink_soma, orig_cuts$tier1_cut, orig_cuts$tier4_cut, orig_cuts$tier2_cut) %>%
    select(protein, tiers_original = tiers)

## denser one-at-a-time perturbation grid: vary each cutoff across a fine range,
## holding the other two fixed at their original value
delta_seq <- seq(-0.1, 0.1, by = 0.02)

build_one_at_a_time <- function(df, original, direction_label) {
    purrr::map_dfr(c("tier1_cut", "tier4_cut", "tier2_cut"), function(param) {
        purrr::map_dfr(delta_seq, function(delta) {
            cuts <- orig_cuts
            cuts[[param]] <- orig_cuts[[param]] + delta
            
            tiered <- assign_tiers(df, cuts$tier1_cut, cuts$tier4_cut, cuts$tier2_cut) %>%
                select(protein, tiers_alt = tiers)
            
            merged <- original %>% inner_join(tiered, by = "protein")
            
            tier_props <- merged %>%
                count(tiers_alt) %>%
                mutate(prop = n / sum(n))
            
            pct_changed <- mean(merged$tiers_original != merged$tiers_alt, na.rm = TRUE)
            
            tibble(
                direction = direction_label,
                param_varied = param,
                delta = delta,
                threshold_value = cuts[[param]],
                pct_reclassified = pct_changed
            )
        })
    })
}

sensitivity_curve <- bind_rows(
    build_one_at_a_time(soma_olink, original_so, "SomaScan to Olink"),
    build_one_at_a_time(olink_soma, original_os, "Olink to SomaScan")
)

## preparing for plotting
sensitivity_curve <- sensitivity_curve %>%
    mutate(direction = factor(direction, levels = c("SomaScan to Olink", "Olink to SomaScan")))


##### 03 - plotting the change #####
param_labels <- c(tier1_cut = "Tier 1 orig cutoff (0.7)", tier4_cut = "Tier 4 orig cutoff (0.1)", tier2_cut = "Tier 2 orig cutoff (0.4)")

p_sensitivity <- ggplot(sensitivity_curve, aes(x = delta, y = pct_reclassified, color = param_varied)) +
    geom_line(linewidth = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~direction, ncol = 1) +
    scale_color_manual(values = c(tier1_cut = "#d62728", tier4_cut = "#1f77b4", tier2_cut = "#2ca02c"),
                       labels = param_labels) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(
        breaks = seq(-0.15, 0.15, by = 0.05)
    ) +
    labs(x = "Threshold shift (± from original cutoff)", y = "% of proteins reclassified",
         color = "Threshold varied",
         title = "Sensitivity of tier assignment to threshold choice") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom")

p_sensitivity




##### 03 - another plot showing differences! #####
build_tier_proportions <- function(df, original, direction_label) {
    purrr::map_dfr(c("tier1_cut", "tier4_cut", "tier2_cut"), function(param) {
        purrr::map_dfr(delta_seq, function(delta) {
            cuts <- orig_cuts
            cuts[[param]] <- orig_cuts[[param]] + delta
            
            tiered <- assign_tiers(df, cuts$tier1_cut, cuts$tier4_cut, cuts$tier2_cut)
            
            tiered %>%
                count(tiers) %>%
                mutate(prop = n / sum(n),
                       direction = direction_label,
                       param_varied = param,
                       delta = delta)
        })
    })
}

tier_proportions <- bind_rows(
    build_tier_proportions(soma_olink, original_so, "SomaScan to Olink"),
    build_tier_proportions(olink_soma, original_os, "Olink to SomaScan")
    ) %>%
    filter(!is.na(tiers)) %>%
    # mutate(tiers = factor(tiers, levels = c(1, 2, 3, 4), labels = c("Tier 1", "Tier 2", "Tier 3", "Tier 4"))) %>%
    mutate(tiers = factor(tiers, levels = c(1, 3, 2, 4), labels = c("Tier 1", "Tier 3", "Tier 2", "Tier 4"))) %>%
    mutate(direction = factor(direction, levels = c("SomaScan to Olink", "Olink to SomaScan"))) %>%
    filter(
        (param_varied == "tier1_cut" & tiers == "Tier 1") | 
        (param_varied == "tier2_cut" & tiers %in% c("Tier 2", "Tier 3")) |
        (param_varied == "tier4_cut" & tiers == "Tier 4")
    )

## pick a handful of checkpoints to label, rather than every delta value
label_deltas <- c(-0.10, -0.05, 0, 0.05, 0.10)

## compute cumulative (stacked) y-position for each label, matching the
## area plot's stacking order exactl

tier_labels <- tier_proportions %>%
    filter(delta %in% label_deltas) %>%
    arrange(param_varied, direction, delta, desc(tiers)) %>%
    group_by(param_varied, direction, delta) %>%
    mutate(
        cum_top = cumsum(prop),
        cum_bottom = cum_top - prop,
        y_mid = (cum_top + cum_bottom) / 2,
        ## push edge labels inward so they don't clip past the panel boundary
        label_hjust = case_when(
            delta == min(label_deltas) ~ 0,
            delta == max(label_deltas) ~ 1,
            TRUE ~ 0.5
        )
    ) %>%
    ungroup() %>%
    mutate(tiers = factor(tiers, levels = c("Tier 1", "Tier 3", "Tier 2", "Tier 4")))

p_tier_proportions <- ggplot(tier_proportions, aes(x = delta, y = prop, fill = tiers)) +
    geom_area(position = position_stack(reverse = F)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "white") +
    geom_label(
        data = tier_labels,
        aes(x = delta, y = y_mid, label = scales::percent(prop, accuracy = 1), hjust = label_hjust),
        inherit.aes = FALSE, size = 2.4, color = "black",
        fill = "white", alpha = 0.75, label.size = 0, label.padding = unit(0.12, "lines")
    ) +
    facet_nested(. ~ param_varied + direction, labeller = labeller(param_varied = param_labels),
                 axes = "all") +
    scale_fill_manual(values = c("Tier 1" = "#6beeff", "Tier 2" = "#93db86", "Tier 3" = "#ffb963", "Tier 4" = "#ff7878")) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(expand = expansion(mult = 0.08)) +
    labs(x = "Threshold shift (± from original cutoff)", y = "Proportion changed", fill = "Tier") +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom")

p_tier_proportions



##### 04 - save plot objects for downstream combined display item #####
saveRDS(p_sensitivity, "Results/revision1/figures/rds/p_sensitivity.rds")
saveRDS(p_tier_proportions, "Results/revision1/figures/rds/p_tier_proportions.rds")


