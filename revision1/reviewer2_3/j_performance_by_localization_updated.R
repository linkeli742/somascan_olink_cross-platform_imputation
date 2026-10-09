##### 01 - library #####
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(patchwork)
library(dunn.test)
library(patchwork)

gc()
rm(list=ls())

setwd("...")


##### 02 - files #####
soma_olink_raw <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F)
olink_soma_raw <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = F)

# LILRA3, MENT, LEG1: Secreted_Plasma
# DDX58: Intracellular
# KIR2DL2: membrane

protein_class1 <- fread("Results/protein_importance/HPA_derived_protein_class.tsv", data.table = F)
protein_class2 <- data.frame(
    protein = c("LILRA3", "MENT", "LEG1", "DDX58", "KIR2DL2"),
    derived_class = c("Secreted_Plasma","Secreted_Plasma","Secreted_Plasma","Intracellular", "Membrane")
)

protein_class <- rbind(protein_class1, protein_class2) %>%
    mutate(derived_class=ifelse(derived_class == "Secreted_Plasma", "Secreted or Plasma", derived_class))


write.table(protein_class, "Results/protein_importance/HPA_derived_protein_class_full.tsv",
            sep="\t", row.names = F, quote = F)

soma_olink_class <- left_join(soma_olink_raw, protein_class, by = "protein")
olink_soma_class <- left_join(olink_soma_raw, protein_class, by = "protein")

### all the correlations
baseline_corr <- soma_olink_class %>% select(protein, pretrain_r, derived_class) %>%
    rename(corr = pretrain_r) %>%
    mutate(impute_label = "Baseline Correlation",
           direction_label = "Bi-directional Identical")

soma_olink_class_posttrain <- soma_olink_class %>% select(protein, posttrain_r, derived_class) %>%
    rename(corr = posttrain_r) %>%
    mutate(impute_label = "Imputated Correlation",
           direction_label = "SomaScan to Olink")

olink_soma_class_posttrain <- olink_soma_class %>% select(protein, posttrain_r, derived_class) %>%
    rename(corr = posttrain_r) %>%
    mutate(impute_label = "Imputated Correlation",
           direction_label = "Olink to SomaScan"
    )


## merging to long
corr_df <- rbind(baseline_corr, soma_olink_class_posttrain, olink_soma_class_posttrain) %>%
    mutate(direction_label = factor(direction_label, levels = c("Bi-directional Identical", "SomaScan to Olink", "Olink to SomaScan")))




#### 03a - existing correlation plot (unchanged) #####
summary_stats <- corr_df %>%
    group_by(derived_class, impute_label, direction_label) %>%
    summarise(
        mean_corr = mean(corr, na.rm = TRUE),
        median_corr = median(corr, na.rm = TRUE),
        n = n(),
        .groups = "drop"
    )

p_corr <- ggplot(corr_df, aes(x = derived_class, y = corr, fill = derived_class)) +
    facet_nested(. ~ impute_label + direction_label, scales = "free_x", space = "free_x") +
    geom_violin(alpha = 0.7) +
    geom_jitter(width = 0.3, alpha = 0.2, size = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2.5, fill = "white", color = "black") +
    stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
    geom_text(
        data = summary_stats,
        aes(x = derived_class, y = max(corr_df$corr, na.rm = TRUE) * 1.08,
            label = paste0("mean=", round(mean_corr, 2), "\nmedian=", round(median_corr, 2))),
        inherit.aes = FALSE, size = 2.6, lineheight = 0.85
    ) +
    labs(y = "Pearson Correlation", x = "") +
    theme_bw(base_size = 12) +
    guides(fill = "none")




##### 03b - gain density with KW test annotation #####
## build gain_df: same derived_class join, but gain = posttrain_r - pretrain_r per direction
soma_olink_gain <- soma_olink_class %>%
    mutate(gain = posttrain_r - pretrain_r) %>%
    select(protein, gain, derived_class) %>%
    mutate(direction_label = "SomaScan to Olink")

olink_soma_gain <- olink_soma_class %>%
    mutate(gain = posttrain_r - pretrain_r) %>%
    select(protein, gain, derived_class) %>%
    mutate(direction_label = "Olink to SomaScan")

gain_df <- rbind(soma_olink_gain, olink_soma_gain) %>%
    mutate(direction_label = factor(direction_label, levels = c("SomaScan to Olink", "Olink to SomaScan")))

## run KW + Dunn's per direction, capture results
kw_results <- list()
dunn_results <- list()


for (dir in unique(gain_df$direction_label)) {
    sub <- gain_df %>% filter(direction_label == dir)
    kw <- kruskal.test(gain ~ derived_class, data = sub)
    kw_results[[dir]] <- kw
    
    dunn <- dunn.test(sub$gain, sub$derived_class, method = "bonferroni", kw = FALSE, table = FALSE)
    dunn_results[[dir]] <- dunn
}


## build annotation text per direction: KW p-value + which pairs differ (Dunn's, adj p < 0.05)
annot_df <- lapply(names(kw_results), function(dir) {
    kw_p <- kw_results[[dir]]$p.value
    kw_label <- ifelse(kw_p < 0.001, "p < 0.001", paste0("p = ", formatC(kw_p, digits = 3, format = "f")))
    
    dunn <- dunn_results[[dir]]
    sig_pairs <- dunn$comparisons[dunn$P.adjusted < 0.05]
    sig_label <- if (length(sig_pairs) > 0) {
        paste0("Significant pairs:\n", paste(sig_pairs, collapse = "\n"))
    } else {
        "No pairs significant"
    }
    
    data.frame(
        direction_label = dir,
        label = paste0("Kruskal-Wallis ", kw_label, "\n", sig_label)
    )
}) %>% 
    bind_rows() %>%
    mutate(direction_label = factor(direction_label, levels = c("SomaScan to Olink", "Olink to SomaScan")))


## density plot with KW/Dunn annotation
p_gain_density <- ggplot(gain_df, aes(x = gain, fill = derived_class, color = derived_class)) +
    facet_nested(.~direction_label) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_density(alpha = 0.4, linewidth = 0.6) +
    geom_vline(
        data = gain_df %>% group_by(derived_class, direction_label) %>%
            summarise(mean_gain = mean(gain, na.rm = TRUE), .groups = "drop"),
        aes(xintercept = mean_gain, color = derived_class),
        linetype = "solid", linewidth = 0.9, show.legend = FALSE
    ) +
    geom_text(
        data = annot_df,
        aes(x = Inf, y = Inf, label = label),
        inherit.aes = FALSE, hjust = 1.05, vjust = 1.1, size = 2.8, lineheight = 0.95
    ) +
    labs(x = "Gain (post - baseline r)", y = "Density",
         fill = "Localization", color = "Localization") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

p_gain_density



##### 04 - merging 2 together ##### 
combined_plot <- p_corr / p_gain_density +
    plot_layout(heights = c(1.5, 1.1)) +
    plot_annotation(tag_levels = "A")

combined_plot

ggsave("Results/revision1/figures/imputed_corr_and_gain_density_combined.png", combined_plot,
       dpi = 300, width = 10, height = 8, bg = "white")


