##### 01 - library #####
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(patchwork)


gc()
rm(list=ls())

setwd("...")


##### 02 - files #####
soma_olink <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F)
olink_soma <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = F) %>%
    mutate(gain = posttrain_r - pretrain_r)


# --- add direction labels ---
soma_olink2 <- soma_olink %>% mutate(direction = "SomaScan → Olink (SO)")
olink_soma2 <- olink_soma %>% mutate(direction = "Olink → SomaScan (OS)")

prot_loc <- fread("Results/protein_importance/HPA_derived_protein_class_full.tsv", data.table = F)


# --- combine, keep just what we need for the histogram ---
combined <- bind_rows(
    soma_olink2 %>% select(protein, direction, tiers, baseline = pretrain_r, post_imputation = posttrain_r),
    olink_soma2 %>% select(protein, direction, tiers, baseline = pretrain_r, post_imputation = posttrain_r)
) %>%
    mutate(
        tiers = factor(tiers, levels=c(1,2,3,4), labels = c("Tier 1", "Tier 2", "Tier 3", "Tier 4")),
        direction = factor(direction, levels = c("SomaScan → Olink (SO)", "Olink → SomaScan (OS)"))
    ) %>%
    left_join(prot_loc, by = "protein")

table(combined$derived_class, useNA = "ifany") # there's no NA, checked out



# creating plotting dataframe ---------------------------------------------
soma_olink_neggain_loc <- soma_olink %>%
    filter(gain < 0) %>%
    left_join(prot_loc, by = "protein")


olink_soma_neggain_loc <- olink_soma %>%
    filter(gain < 0) %>%
    left_join(prot_loc, by = "protein")




##### 03 - plots #####
gain_threshold <- -0.05
neg_gain <- combined %>%
    mutate(gain = post_imputation - baseline,
           neg_flag = gain < gain_threshold)

##### 03b - enrichment: class composition of negative-gain proteins vs. background #####
# --- background composition (fixed, from full protein universe, n = 1737) ---
bg_counts <- c(Intracellular = 357, Membrane = 294, `Secreted or Plasma` = 1086)
bg_total  <- sum(bg_counts)          # 1737
bg_prop   <- bg_counts / bg_total    # expected proportions under null

gain_threshold <- -0.05  # definition of meaningful negative gain

# --- omnibus goodness-of-fit test per direction (exact via simulation, small-n safe) ---
omnibus_results <- neg_gain %>%
    filter(gain < gain_threshold) %>%
    group_by(direction) %>%
    group_modify(~ {
        obs_tab <- table(factor(.x$derived_class, levels = names(bg_counts)))
        ct <- chisq.test(x = as.vector(obs_tab), p = bg_prop,
                         simulate.p.value = TRUE, B = 1e5)
        tibble(
            n_neg_total = sum(obs_tab),
            p_omnibus   = ct$p.value
        )
    }) %>%
    ungroup()

print(omnibus_results)

# --- per-class post-hoc: exact one-sample binomial test vs. background proportion ---
class_summary <- neg_gain %>%
    filter(gain < gain_threshold) %>%
    group_by(direction) %>%
    group_modify(~ {
        d <- .x
        n_neg_total <- nrow(d)
        purrr::map_dfr(names(bg_counts), function(cls) {
            n_neg_in_class <- sum(d$derived_class == cls, na.rm = TRUE)
            bt <- binom.test(n_neg_in_class, n_neg_total, p = bg_prop[[cls]])
            tibble(
                derived_class    = cls,
                n_neg_in_class   = n_neg_in_class,
                n_neg_total      = n_neg_total,
                observed_prop    = n_neg_in_class / n_neg_total,
                expected_prop    = bg_prop[[cls]],
                fold_enrichment  = observed_prop / bg_prop[[cls]],
                ci_low           = bt$conf.int[1] / bg_prop[[cls]],
                ci_high          = bt$conf.int[2] / bg_prop[[cls]],
                p_value          = bt$p.value
            )
        })
    }) %>%
    ungroup() %>%
    group_by(direction) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(
        sig_label = case_when(
            p_adj < 0.001 ~ "***",
            p_adj < 0.01  ~ "**",
            p_adj < 0.05  ~ "*",
            TRUE          ~ ""
        )
    )


##### 03c - bar plot: fold enrichment with exact CI #####
pB_top <- ggplot(class_summary, aes(x = derived_class, y = fold_enrichment, fill = derived_class)) +
    geom_col(width = 0.65) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, color = "grey30") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    geom_text(aes(label = paste0(n_neg_in_class, "/", n_neg_total)),
              y = -Inf, vjust = -0.8, size = 3) +
    geom_text(aes(label = sig_label, y = ci_high), vjust = -0.8, size = 4.5, fontface = "bold") +
    facet_wrap(~direction, axes = "all") +
    labs(x = NULL, y = "Relative enrichment vs. background composition",
         fill = "Localization") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold"),
          plot.caption = element_text(size = 8, color = "grey40"))


##### 03d - Panel B: compact stats table underneath the figure #####
table_df <- class_summary %>%
    mutate(
        `95% CI`   = paste0(sprintf("%.2f", ci_low), " - ", sprintf("%.2f", ci_high)),
        `Relative Enrichment` = sprintf("%.2f", fold_enrichment),
        `P adj.`   = sprintf("%.3f", p_adj)
    ) %>%
    select(Direction = direction, Class = derived_class,
           `N Negative / N Total` = n_neg_in_class, `Relative Enrichment`, `95% CI`, `P adj.`, `Sig ?` = sig_label)

# collapse n_neg_in_class / n_neg_total into one column for compactness
table_df$`N Negative / N Total` <- paste0(class_summary$n_neg_in_class, " / ", class_summary$n_neg_total)

tbl_grob <- ggtexttable(table_df, rows = NULL,
                        theme = ttheme(base_size = 8, padding = unit(c(2, 2), "mm")))

##### 03e - combine plot (Panel A) + table (Panel B) into one figure #####
pB_combined <- pB_top / tbl_grob +
    plot_layout(heights = c(1.5, 1)) +
    plot_annotation(tag_levels = "A") &  
    theme(plot.tag = element_text(face = "bold"))

pB_combined



ggsave("Results/revision1/figures/neg_gain_class_enrichment.png", pB_combined, dpi=300,
       width=10, height=7, bg="white")





