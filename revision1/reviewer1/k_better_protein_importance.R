##### 01 - library #####
library(dplyr)
library(data.table)
library(ggplot2)
library(ggnewscale)
library(colorspace)
library(wesanderson)
library(ggrepel)
library(ggh4x)
library(ggpubr)
library(grid)

gc()
rm(list = ls())

setwd("...")


##### 02 - shared inputs #####
tot_pubs <- fread("Results/protein_importance/MALACARDS_proteins_sumofpubs.tsv.gz", data.table = F)
protein_hpa_tags <- fread("Results/protein_importance/HPA_importance_protclasstags.tsv", data.table = F)

fda_drug_targets      <- unique(protein_hpa_tags[protein_hpa_tags$`FDA approved drug targets` == 1, ]$protein)
potential_drug_targets<- unique(protein_hpa_tags[protein_hpa_tags$`Potential drug targets` == 1, ]$protein)
cancer_genes  <- unique(protein_hpa_tags[protein_hpa_tags$`Cancer-related genes` == 1, ]$protein)
disease_genes <- unique(protein_hpa_tags[protein_hpa_tags$`Disease related genes` == 1, ]$protein)

tier_levels <- c("1", "2", "3", "4")
tier_cols   <- c("1" = "#6beeff", "2" = "#93db86", "3" = "#ffb963", "4" = "#ff7878")
nbins <- 40


##### 03 - Panel A: bin plot (rank by gain, colored by pubs) #####
build_bin_plot <- function(rank_file) {
    dat <- fread(rank_file, data.table = F) %>%
        mutate(tiers = as.character(tiers), gain = posttrain_r - pretrain_r) %>%
        left_join(tot_pubs, by = "protein") %>%
        arrange(desc(gain)) %>%
        mutate(index = row_number(gain))
    
    pd <- dat %>%
        mutate(tiers = factor(tiers, levels = tier_levels)) %>%
        arrange(index) %>%
        mutate(xbin = cut_number(index, n = nbins))
    
    slab_base <- pd %>%
        group_by(xbin) %>%
        summarise(xmin = min(index), xmax = max(index), height = max(gain), .groups = "drop")
    
    tier_stats <- pd %>%
        group_by(xbin, tiers) %>%
        summarise(n = n(), heat = sum(tot_pubs), .groups = "drop")
    
    slabs_tier <- slab_base %>%
        left_join(tier_stats, by = "xbin") %>%
        group_by(xbin) %>%
        mutate(prop = n / sum(n)) %>%
        arrange(tiers, .by_group = TRUE) %>%
        mutate(cumprop = cumsum(prop), ymin = height * (cumprop - prop), ymax = height * cumprop) %>%
        ungroup()
    
    global_max_pubs <- max(tier_stats$heat, na.rm = T)
    
    p <- ggplot()
    for (t in tier_levels) {
        slab_t <- dplyr::filter(slabs_tier, tiers == t)
        if (nrow(slab_t) > 0) {
            p <- p +
                geom_rect(data = slab_t,
                          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = heat),
                          color = NA, alpha = 1, inherit.aes = FALSE) +
                scale_fill_gradient(
                    low = tier_cols[[t]], high = colorspace::darken(tier_cols[[t]], 0.6),
                    limits = c(0, global_max_pubs),
                    name = paste0("Sum of Publications \n(Tier ", t, ")"),
                    guide = guide_colorbar(order = match(t, tier_levels), direction = "vertical", title.position = "top")
                ) +
                ggnewscale::new_scale_fill()
        }
    }
    
    p +
        geom_line(data = dat, aes(x = index, y = gain), linewidth = 0.8, color = "black") +
        labs(x = "Index", y = "Gain") +
        scale_y_continuous(limits = c(-0.2, 0.5), breaks = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
        theme_minimal() +
        theme(legend.direction = "vertical", legend.position = "right",
              legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 10),
              axis.title = element_text(size = 12, face = "bold"),
              axis.text = element_text(size = 11, color = "black"))
}

soma_olink_bin <- build_bin_plot("Results/rank_by_gain/soma5k_olink3k.tsv")
olink_soma_bin <- build_bin_plot("Results/rank_by_gain/olink3k_soma5k.tsv")

# combined bin-plot construct (this whole block becomes outer Panel A)
bin_combo <- ggarrange(soma_olink_bin, olink_soma_bin,
                       common.legend = T, legend = "right", align = "hv", nrow = 1)


##### 04 - Panel B: quadrant scatter (gain vs. log(N Pubs)) #####
build_quadrant_plot <- function(rank_file) {
    dat <- fread(rank_file, data.table = F) %>%
        mutate(tiers = as.character(tiers), gain = posttrain_r - pretrain_r) %>%
        left_join(tot_pubs, by = "protein") %>%
        mutate(
            drug_group = case_when(
                protein %in% fda_drug_targets ~ "FDA Drug Targets",
                protein %in% potential_drug_targets ~ "Potential Drug Targets",
                TRUE ~ "None"
            ),
            gene_group = case_when(
                protein %in% disease_genes ~ "Disease related genes",
                protein %in% cancer_genes ~ "Cancer related genes",
                TRUE ~ "N/A"
            ),
            drug_group = factor(drug_group, levels = c("FDA Drug Targets", "Potential Drug Targets", "None")),
            gene_group = factor(gene_group, levels = c("Cancer related genes", "Disease related genes", "N/A")),
            # explicit log10 column per Reviewer 1's Fig 4 comment -----------------
            log_n_pubs = log10(tot_pubs)
        ) %>%
        group_by(tiers) %>%
        arrange(desc(gain)) %>%
        mutate(
            index = row_number(gain),
            tiers = case_when(
                tiers == "1" ~ "Tier 1", tiers == "2" ~ "Tier 2",
                tiers == "3" ~ "Tier 3", tiers == "4" ~ "Tier 4",
                TRUE ~ NA_character_
            )
        ) %>%
        ungroup()
    
    annot_set1 <- dat %>% group_by(tiers) %>% slice_max(order_by = tot_pubs, n = 6)
    annot_set2 <- dat %>% filter(protein %in% c("MAPT", "TNNI3"))
    annot <- rbind(annot_set1, annot_set2)
    
    ggplot(data = dat, aes(x = gain, y = log_n_pubs, color = gene_group, shape = drug_group, alpha = gene_group)) +
        geom_point() +
        facet_wrap(~tiers, axes = "all", ncol = 1, scales = "free_x") +
        geom_label_repel(data = annot, aes(x = gain, y = log_n_pubs, label = protein, color = gene_group),
                         inherit.aes = F, size = 3, max.overlaps = Inf, fontface = "bold",
                         box.padding = 0.4, min.segment.length = 0, show.legend = F) +
        scale_shape_manual(values = c("None" = 1, "FDA Drug Targets" = 8, "Potential Drug Targets" = 17),
                           breaks = c("FDA Drug Targets", "Potential Drug Targets")) +
        scale_color_manual(values = c("Disease related genes" = "firebrick3",
                                      "Cancer related genes" = "royalblue3",
                                      "N/A" = "grey6"),
                           breaks = c("Disease related genes", "Cancer related genes")) +
        scale_alpha_manual(values = c("Disease related genes" = 1, "Cancer related genes" = 1, "N/A" = 0.5),
                           guide = "none") +
        labs(shape = "Drug Targets", color = "Gene Group", x = "Gain", y = "log(N Pubs)") +
        theme_minimal() +
        theme(legend.direction = "horizontal", legend.position = "bottom",
              
              legend.box = "vertical", legend.box.just = "center",
              legend.margin = margin(t = 0, b = 0),
              axis.title = element_text(size = 12, face = "bold"),
              axis.text = element_text(size = 11, color = "black"),
              legend.title = element_text(size = 12, face = "bold"),
              legend.text = element_text(size = 12)) +
        guides(color = guide_legend(order = 1, override.aes = list(size = 4)),
               shape = guide_legend(order = 2, override.aes = list(size = 4)))
}

soma_olink_quad <- build_quadrant_plot("Results/rank_by_gain/soma5k_olink3k.tsv")
olink_soma_quad <- build_quadrant_plot("Results/rank_by_gain/olink3k_soma5k.tsv")

# combined quadrant-scatter construct (this whole block becomes outer Panel B)
quadrant_combo <- ggarrange(soma_olink_quad, olink_soma_quad + rremove("ylab"),
                            common.legend = T, legend = "bottom")


##### 05 - final merged display item #####
merged_fig <- ggarrange(bin_combo, quadrant_combo,
                        labels = c("A", "B"),
                        ncol = 2, heights = c(1.3, 1.3)) # panel B is taller (faceted by tier)

merged_fig

ggsave("Results/revision1/figures/merged_bin_quadrant_display_item.png", merged_fig,
       dpi = 300, width = 14, height = 8, bg = "white")
