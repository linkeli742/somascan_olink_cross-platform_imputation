##### 01 - library #####
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(glue)

gc()
rm(list = ls())


setwd("...")


##### 02 - files #####
bingyu_results <- fread("Data/bingyu_hf/bingyu_visit3_results", data.table = F)
olink_raw <- fread("Results/bingyu_hf/UKB_HF_BY_olinks_impute_10pcs_cadsoft_hypertension01.tsv", data.table = F) %>%
    mutate(pval_adj = p.adjust(pval, method = "bonferroni"))
soma_raw <- fread("Results/bingyu_hf/UKB_HF_BY_predsoma_impute_10pcs_cadsoft_hypertension01.tsv", data.table = F) %>%
    mutate(pval_adj = p.adjust(pval, method = "bonferroni"))

# OS-direction (Olink -> SomaScan) tiers from the primary HGB model
# same source table used in h_model_comparison_hgb_lgb_enet.R
os_tiers <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = FALSE) %>%
    select(protein, tiers, pretrain_r, posttrain_r)


##### 03 - plotting to compare (17 highlighted proteins from Shah et al.) #####
imp_17_prots <- fread("/medpop/esp2/lli/somalink_olink/Data/bingyu_hf/visit3_17_proteins", data.table = F, header = F)$V1
soma_17 <- soma_raw %>%
    filter(exposures %in% c(imp_17_prots,"NPPB_3")) %>%
    filter(!exposures %in% "NPPB")

olink_17 <- olink_raw %>%
    filter(exposures %in% imp_17_prots)

yb_effect_sizes <- bingyu_results %>%
    filter(gene %in% c(imp_17_prots,"NPPB_3")) %>%
    ## duplicate NPPB2
    filter(!gene %in% c("NPPB_2")) %>%
    ### double MSR1 which we dont have 
    filter(!term %in% c("SeqId_11207_3"))

names(soma_17)
names(yb_effect_sizes)

## renaming and getting ready for merging
soma_17 <- soma_17 %>%
    rename(protein = exposures) %>%
    select(-outcomes) %>%
    mutate(label = "UKB Imputed SomaScan") 

olink_17 <- olink_17 %>%
    rename(protein = exposures) %>%
    select(-outcomes) %>%
    mutate(label = "UKB Olink") 

yb_effect_sizes <- yb_effect_sizes %>%
    select(protein = gene, hr = estimate, lci = conf.low, uci = conf.high, pval = p.value, pval_adj = adjPValueBonf) %>%
    mutate(label = "Shah et al., Nature Communications 2024")

## fixing molecule names based on seq id
yb_effect_sizes[yb_effect_sizes$protein=="NPPB",]$protein <- "NPPB_2"
yb_effect_sizes[yb_effect_sizes$protein=="NPPB_3",]$protein <- "NPPB"

soma_17[soma_17$protein=="NPPB_2",]$protein <- "NPPB"
soma_17[soma_17$protein=="NPPB_3",]$protein <- "NPPB_2"

plot_data <- rbind(soma_17, olink_17, yb_effect_sizes) %>%
    mutate(sig = ifelse(pval_adj < 0.05, "p<0.05" , "ns"))

order <- yb_effect_sizes %>%
    arrange(desc(pval_adj)) %>%
    pull(protein)

plot_data <- plot_data %>%
    mutate(protein = factor(protein, levels = order))

colors <- c("Shah et al., Nature Communications 2024" = "#F8766D", "UKB Olink" = "#00BA38", "UKB Imputed SomaScan" = "#619cfF")
n <- nlevels(plot_data$protein)
sep_df <- data.frame(y = seq(1.5, n - 0.5, by = 1))

ggplot(plot_data, aes(x = hr, y = protein, color = label)) + 
    geom_hline(data = sep_df, aes(yintercept = y),
               inherit.aes = FALSE, color = "gray85", linewidth = 0.4) +
    geom_point(aes(group = label, shape = sig, color = label),
               position = position_dodge(width = 0.5), size = 1.4) +
    geom_errorbar(aes(group = label, xmin = lci, xmax = uci, color = label),
                  position = position_dodge(width = 0.5), width = 0.4) +
    scale_shape_manual(values = c("ns" = 1, "p<0.05" = 19), guide = "none") +
    scale_color_manual(values = colors) + 
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    labs(
        x = "Hazard Ratio",
        y = "Protein",
        color = "Analysis"
    ) + 
    theme_linedraw() +
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        plot.title = element_text(size=9),
        plot.caption = element_text(hjust = 0, size = 7),
        axis.text = element_text(size=10, colour = "black"),
        axis.title = element_text(size = 13, face="bold"),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 11), 
        legend.key.width = unit(1.2, "cm")
    )

ggsave("Results/figures/forest_plot_hf_replication.png", dpi = 300, width = 9, height = 6, bg = "white")


#### completing tables (17 highlighted proteins) ####
plot_data <- rbind(soma_17, olink_17,yb_effect_sizes)
hf_protein_list <- plot_data[plot_data$label == "Shah et al., Nature Communications 2024",]$protein
hf_results_long <- plot_data %>%
    complete(label, protein = hf_protein_list) %>%
    select(label, protein, hr, lci, uci, pval) %>%
    mutate(hr = round(hr, 2),
           lci = round(lci, 2),
           uci = round(uci, 2),
    ) %>%
    mutate(
        across(c(hr, lci, uci), ~ ifelse(is.na(.), "-", format(., nsmall = 2))),
        hr_95ci = ifelse(hr == "-", "-", glue("{hr} ({lci}-{uci})"))
        
    )

order <- yb_effect_sizes %>%
    arrange(pval_adj) %>%
    pull(protein)

hf_results_wide <- hf_results_long %>%
    pivot_wider(id_cols="protein", 
                names_from = "label",
                values_from = c("hr_95ci", "pval"),
    ) %>%
    select(protein, contains("Shah et al., Nature Communications 2024"), contains("UKB Imputed SomaScan"), contains("UKB Olink")) %>%
    mutate(protein = factor(protein, levels = order)) %>%
    arrange(protein)

write.table(hf_results_wide, "Results/bingyu_hf/hf_results_17prots_wide.tsv", sep = "\t", row.names = F, quote = F)


##### 04 - systematic effect-size concordance across ALL overlapping proteins #####
olink_full <- olink_raw %>%
    rename(protein = exposures) %>%
    select(protein, hr_olink = hr, lci_olink = lci, uci_olink = uci,
           pval_olink = pval, pval_adj_olink = pval_adj)

soma_full <- soma_raw %>%
    rename(protein = exposures) %>%
    select(protein, hr_soma = hr, lci_soma = lci, uci_soma = uci,
           pval_soma = pval, pval_adj_soma = pval_adj)

concordance_17 <- olink_full %>%
    inner_join(soma_full, by = "protein") %>%
    inner_join(os_tiers, by = "protein") %>%
    filter(!is.na(hr_olink), !is.na(hr_soma), protein %in% imp_17_prots) %>%
    mutate(
        log_hr_olink = log(hr_olink),
        log_hr_soma  = log(hr_soma),
        sig_olink = pval_adj_olink < 0.05,
        sig_soma  = pval_adj_soma  < 0.05,
        same_direction = sign(log_hr_olink) == sign(log_hr_soma)
    )

cat("=== N proteins in full overlapping-panel concordance set ===\n")
print(nrow(concordance_17))

## overall concordance: correlation of log(HR) (HR is multiplicative -> log scale
## is the natural one for comparing effect sizes, same logic as a forest plot's x-axis)
overall_cor <- cor.test(concordance_17$log_hr_olink, concordance_17$log_hr_soma,
                        method = "pearson")

cat("=== Overall log(HR) concordance, measured Olink vs imputed SomaScan ===\n")
print(overall_cor)

cat("=== Sign (direction) concordance, all proteins ===\n")
print(mean(concordance_17$same_direction))

## tier-stratified concordance (OS direction tiers) -- directly informs Reviewer 3
## comment #4's skepticism: do tier 2/3 proteins replicate meaningfully worse than
## tier 1, or is direction/effect-size concordance roughly preserved across tiers?
tier_concordance_summary <- concordance_17 %>%
    group_by(tiers) %>%
    summarise(
        n_proteins = n(),
        log_hr_cor = cor(log_hr_olink, log_hr_soma, use = "complete.obs"),
        pct_same_direction = mean(same_direction),
        pct_both_significant = mean(sig_olink & sig_soma),
        .groups = "drop"
    )

cat("=== Tier-stratified (OS direction) log(HR) concordance ===\n")
print(tier_concordance_summary)

write.table(concordance_all, "Results/bingyu_hf/hf_concordance_all_overlapping_proteins.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tier_concordance_summary, "Results/bingyu_hf/hf_concordance_by_tier.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

##### 05 - scatter plot: measured Olink log(HR) vs imputed SomaScan log(HR), by tier #####
p_concordance <- ggplot(concordance_all, aes(x = log_hr_olink, y = log_hr_soma, color = factor(tiers))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(alpha = 0.6, size = 1.6) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
    labs(
        x = "log(HR), measured UKB Olink",
        y = "log(HR), UKB imputed SomaScan",
        color = "Tier (OS direction)",
        caption = glue("Pearson r = {round(overall_cor$estimate, 3)}, ",
                       "p = {format.pval(overall_cor$p.value, digits = 2)}, ",
                       "n = {nrow(concordance_all)} overlapping proteins")
    ) +
    theme_linedraw() +
    theme(
        legend.position = "right",
        plot.caption = element_text(hjust = 0, size = 9),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12, face = "bold")
    )

p_concordance

ggsave("Results/figures/hf_concordance_all_overlapping_proteins.png", p_concordance,
       dpi = 300, width = 7, height = 6, bg = "white")