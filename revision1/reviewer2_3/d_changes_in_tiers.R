##### 01 - library #####
library(dplyr)
library(data.table)
library(ggalluvial)
library(ggplot2)
library(patchwork)
library(ggpubr)

gc()
rm(list=ls())

setwd("...")



##### 02 - files ####
soma_olink <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F) %>%
    mutate(tiers = factor(tiers, levels=c(1,2,3,4), labels = c("Tier 1", "Tier 2", "Tier 3", "Tier 4")))

olink_soma <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = F) %>%
    mutate(gain = posttrain_r - pretrain_r) %>%
    mutate(tiers = factor(tiers, levels=c(1,2,3,4), labels = c("Tier 1", "Tier 2", "Tier 3", "Tier 4")))

tier_compare <- soma_olink %>%
    select(protein, tier_SO = tiers) %>%
    inner_join(
        olink_soma %>% select(protein, tier_OS = tiers),
        by = "protein"
    )

tier_pal <- c("Tier 1" = "#6beeff", "Tier 2" = "#93db86", "Tier 3" = "#ffb963", "Tier 4" = "#ff7878")

flow_counts <- tier_compare %>%
    count(tier_SO, tier_OS, name = "n")

# --- Plot A: alluvial diagram, SO tier -> OS tier, colored by tier (original) ---
pA <- ggplot(flow_counts,
             aes(axis1 = tier_SO, axis2 = tier_OS, y = n)) +
    geom_alluvium(aes(fill = tier_SO), width = 0.25, alpha = 0.8) +
    geom_stratum(width = 0.25, fill = "grey85", color = "grey30") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
    scale_x_discrete(limits = c("SomaScan → Olink tier", "Olink → SomaScan tier"),
                     expand = c(0.15, 0.15)) +
    scale_fill_manual(values = tier_pal) +
    labs(y = "Number of proteins", fill = "SO tier") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")

# --- Plot B: SO tier x OS tier confusion matrix, % of row total ---
tier_matrix <- tier_compare %>%
    mutate(tier_SO = factor(tier_SO, levels = c("Tier 4", "Tier 3", "Tier 2", "Tier 1")),
           tier_OS = factor(tier_OS, levels = c("Tier 4", "Tier 3", "Tier 2", "Tier 1"))
    ) %>%
    count(tier_SO, tier_OS) %>%
    group_by(tier_SO) %>%
    mutate(pct = n / sum(n)) %>%
    ungroup()

pB <- ggplot(tier_matrix, aes(x = tier_OS, y = tier_SO, fill = pct)) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste0(n, "\n(", scales::percent(pct, accuracy = 1), ")")),
              size = 3.2) +
    scale_fill_gradient(low = "#f7fbff", high = "#2166ac", labels = scales::percent) +
    labs(x = "Tier, Olink → SomaScan direction", y = "Tier, SomaScan → Olink direction",
         fill = "% of SO tier") +
    theme_minimal(base_size = 12) +
    coord_fixed()

##### 03 - enrichment of protein location within specific transitions #####
prot_loc <- fread("Results/protein_importance/HPA_derived_protein_class_full.tsv", data.table = F)

tier_change <- tier_compare %>%
    inner_join(prot_loc %>% select(protein, derived_class) %>% distinct(), by = "protein") %>%
    mutate(changed = tier_SO != tier_OS)

# background = composition across ALL proteins in the tier comparison (same for every transition)
background <- tier_change %>%
    count(derived_class, name = "n_background") %>%
    mutate(prop_background = n_background / sum(n_background))

compute_enrichment <- function(so, os) {
    observed <- tier_change %>%
        filter(tier_SO == so, tier_OS == os) %>%
        count(derived_class, name = "n_observed") %>%
        mutate(prop_observed = n_observed / sum(n_observed))
    
    enrichment <- observed %>%
        full_join(background, by = "derived_class") %>%
        mutate(
            n_observed = tidyr::replace_na(n_observed, 0),
            prop_observed = tidyr::replace_na(prop_observed, 0),
            enrichment_ratio = round(prop_observed / prop_background, 2)
        )
    
    chisq_result <- chisq.test(x = enrichment$n_observed, p = enrichment$prop_background)
    
    list(table = enrichment, chisq = chisq_result, label = paste0(so, " → ", os))
}

transitions <- list(c("Tier 3", "Tier 2"), c("Tier 3", "Tier 4"), c("Tier 4", "Tier 3"))
enrichment_results <- lapply(transitions, function(x) compute_enrichment(x[1], x[2]))

# --- one table-plot per transition, chi-sq p-value in the title ---
make_table_plot <- function(res) {
    disp <- res$table %>%
        transmute(
            `Location` = derived_class,
            `Obs n` = n_observed,
            `Obs %` = scales::percent(prop_observed, accuracy = 1),
            `Bkgd %` = scales::percent(prop_background, accuracy = 1),
            `Enr. ratio` = enrichment_ratio
        )
    
    p_label <- paste0(res$label, "  (χ² p = ", formatC(res$chisq$p.value, format = "e", digits = 2), ")")
    
    ggtexttable(disp, rows = NULL, theme = ttheme("mOrange", base_size = 9)) %>%
        tab_add_title(text = p_label, face = "italic", size = 10, padding = unit(2, "mm"))
}

table_plots <- lapply(enrichment_results, make_table_plot)

pC <- wrap_plots(table_plots, ncol = 3)


##### 04 - combine #####
left  <- (pB / pA) + plot_layout(heights = c(1, 1.3))
right <- wrap_plots(table_plots, ncol = 1)

combined <- (left | right) +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 14, face = "bold"))

combined

ggsave("Results/revision1/figures/tier_change_combined.png", combined, width = 13, height = 10, dpi = 300)
