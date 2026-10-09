##### 01 - library #####
library(dplyr)
library(data.table)
library(tidyr)
library(glue)
library(ggplot2)
library(ggrepel)

gc()
rm(list = ls())


setwd("...")



##### 02 - keenan files #####
## overlapping
overlapping_proteins <- fread("Results/rank_by_gain/soma5k_olink3k.tsv")$protein
keenan_results <- fread("Data/keenan_walker/midlife_dementia_proteins", data.table = F, header = F)
colnames(keenan_results) <- c("gene_symbol", "somascan_id", "hr", "lci", "uci", "pval")
keenan_results <- keenan_results %>%
    filter(pval < 0.0013) %>%
    select(-somascan_id) %>%
    mutate(label = "Walker et al., Nature Aging 2021")

soma_results <- fread("Results/keenan_walker/UKB_AD_somas_ad_no10pc_hypertension01_edunew.tsv", data.table = F)
soma_short <- soma_results %>%
    filter(exposures %in% keenan_results$gene_symbol) %>%
    select(-outcomes) %>%
    rename(gene_symbol = exposures) %>%
    mutate(label = "UKB Imputed SomaScan")

# write.table(soma_short, "Results/keenan_walker/UKB_AD_somas_short.tsv", sep = "\t", row.names = F, quote = F)
merge <- rbind(keenan_results, soma_short) %>%
    mutate(sig = ifelse(pval < 0.05, "p<0.05" , "ns"))

order <- keenan_results %>% 
    arrange(desc(pval)) %>%
    pull(gene_symbol)

merge <- merge %>%
    mutate(label = factor(label)) %>%
    filter(!gene_symbol %in% olink_short$gene_symbol) %>%
    mutate(annote = ifelse(gene_symbol %in% c("SVEP1", "SVEP1_2"), glue("{round(hr,2)} ({round(lci,2)}-{round(uci,2)}); P-value={sprintf('%.2E', pval)}"), NA)) %>%
    mutate(gene_symbol = factor(gene_symbol, levels = order)) 

keenan_colors <- c("Walker et al., Nature Aging 2021" = "#F8766D", "UKB Imputed SomaScan" = "#619cfF")
n <- nlevels(merge$gene_symbol)
sep_df <- data.frame(y = seq(1.5, n - 0.5, by = 1))

ad_plot <- ggplot(merge, aes(x = hr, y = gene_symbol, color = label)) + 
    geom_hline(data = sep_df, aes(yintercept = y),
               inherit.aes = FALSE, color = "gray85", linewidth = 0.4) +
    geom_text(
        aes(x = uci, label = annote, color = label), 
        position = position_dodge(width = 0.4), 
        hjust = -0.05,
        fontface = "bold",
        size = 3.2,
        show.legend = FALSE
    ) +
    geom_point(aes(group = label, shape = sig),
               position = position_dodge(width = 0.4), size = 1.4) + 
    geom_errorbar(aes(group = label, xmin = lci, xmax = uci),
                  position = position_dodge(width = 0.4), width = 0.4) +
    scale_shape_manual(values = c("ns" = 1, "p<0.05" = 19), guide = "none") +
    scale_color_manual(values = keenan_colors) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    labs(
        x = "Hazard Ratio",
        y = "Protein",
        color = "Analysis"
    ) + 
    theme_linedraw() +
    scale_y_discrete(expand = expansion(add = c(0.6, 0.6))) + 
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




##### 03 - BY files #####
bingyu_results <- fread("Data/bingyu_hf/bingyu_visit3_results", data.table = F)
olink_raw <- fread("Results/bingyu_hf/UKB_HF_BY_olinks_impute_10pcs_cadsoft_hypertension01.tsv", data.table = F) %>%
    mutate(pval_adj = p.adjust(pval, method = "bonferroni"))
soma_raw <- fread("Results/bingyu_hf/UKB_HF_BY_predsoma_impute_10pcs_cadsoft_hypertension01.tsv", data.table = F) %>%
    mutate(pval_adj = p.adjust(pval, method = "bonferroni"))


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

hf_figure <- ggplot(plot_data, aes(x = hr, y = protein, color = label)) + 
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


hf_figure


##### 04 - combine figures #####
combined_plot <- ggarrange(ad_plot, hf_figure, labels = "AUTO", legend = "bottom", align="hv", nrow=2)
combined_plot


ggsave("Results/figures/ad_hf_forest_stacked.png", combined_plot, dpi = 300, width = 9, height = 12, bg = "white")

