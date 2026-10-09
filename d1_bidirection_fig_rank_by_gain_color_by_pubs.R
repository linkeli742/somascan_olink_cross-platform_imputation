##### 01 - library #####
library(dplyr)
library(data.table)
library(ggplot2)
library(ggnewscale)
library(colorspace)
library(wesanderson)
library(ggpubr)

gc()
rm(list = ls())

setwd("...")



##### 02 - files #####
tot_pubs <- fread("Results/protein_importance/MALACARDS_proteins_sumofpubs.tsv.gz", data.table = F)
soma_olink <- fread("Results/rank_by_gain/soma5k_olink3k.tsv", data.table = F)
soma_olink <- soma_olink %>%
    mutate(tiers = as.character(tiers)) %>%
    mutate(gain = posttrain_r - pretrain_r)


## full dataset
soma_olink_full <- left_join(soma_olink, tot_pubs, by = "protein")


## ranking by gain, color by pubs
plot_data <- soma_olink_full %>%
    arrange(desc(gain)) %>%
    mutate(index = row_number(gain))



##### I think bins looks better
nbins <- 220
tier_levels <- c("1","2","3","4")
tier_cols <- c("1"="#26e3fc", "2"="#74c365", "3"="#fc9c26", "4"="#ff5959")


# based on index, cut into bins to represent game (so bars have areas)
pd <- plot_data %>%
    mutate(tiers = factor(tiers, levels = tier_levels)) %>%
    arrange(index) %>%
    mutate(xbin = cut_number(index, n = nbins))

# define the bar width and height
slab_base <- pd %>%
    group_by(xbin) %>%
    # the start and end of each bar/bin
    summarise(xmin = min(index), xmax = max(index),
              height = max(gain), .groups = "drop")

## for each bar, calculate the mean number of publications
tier_stats <- pd %>%
    group_by(xbin, tiers) %>%
    summarise(n = n(),
              heat = sum(tot_pubs),
              .groups = "drop")

# combine, turn proportions into stacked ymin/ymax per tier within each slab
slabs_tier <- slab_base %>%
    left_join(tier_stats, by = "xbin") %>%
    # for each beam, calculate the proportion of proteins in a certain tier compared to the total number of proteins in the bin
    group_by(xbin) %>%
    mutate(prop = n / sum(n)) %>%
    arrange(tiers, .by_group = TRUE) %>%
    # calculate the height of the each tier, the higher the bar the more proteins in the bin 
    mutate(cumprop = cumsum(prop),
           ymin = height*(cumprop - prop),
           ymax = height*cumprop) %>%
    ungroup()


global_max_pubs <- max(tier_stats$heat, na.rm = T)

### create a new ggplot object
p <- ggplot()
for (t in tier_levels) {
    # for each tier
    dat <- dplyr::filter(slabs_tier, tiers == t)
    if (nrow(dat) > 0) {
        p <- p +
            geom_rect(data = dat,
                      aes(xmin = xmin, xmax = xmax, 
                          ymin = ymin, ymax = ymax, fill = heat),
                      color = NA, alpha = 1, inherit.aes = FALSE) +
            
            scale_fill_gradient(
                # reference to the color list and darken the color if the number increase
                # the gradient is specific to tiers
                low  = tier_cols[[t]],
                high = colorspace::darken(tier_cols[[t]], 0.6),
                limits = c(0, global_max_pubs),
                name = paste0("Sum of Publications \n(Tier ", t, ")"),
                guide = guide_colorbar(order = match(t, tier_levels), direction = "vertical", 
                                       title.position = "top")
                # guide = guide_colorbar(order = match(t, tier_levels), direction = "horizontal", 
                #                        title.position = "top", barwidth = unit(4, "cm"))
            ) +
            ggnewscale::new_scale_fill()
    }
}


soma_olink_plot <- p +
    geom_line(data = plot_data, aes(x = index, y = gain),
              linewidth = 0.8, color = "black") +
    labs(x = "Index", y = "Gain") +
    scale_y_continuous(limits = c(-0.2,0.5), breaks = c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5)) +
    theme_minimal() +
    theme(legend.direction = "vertical", legend.position = "right",
          legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12),
          axis.title = element_text(size = 13, face = "bold"),
          axis.text = element_text(size = 12, color="black")
          )




# olink to soma -----------------------------------------------------------
tot_pubs <- fread("Results/protein_importance/MALACARDS_proteins_sumofpubs.tsv.gz", data.table = F)
olink_soma <- fread("Results/rank_by_gain/olink3k_soma5k.tsv", data.table = F)
olink_soma <- olink_soma %>%
    mutate(tiers = as.character(tiers)) %>%
    mutate(gain = posttrain_r - pretrain_r)


## full dataset
olink_soma_full <- left_join(olink_soma, tot_pubs, by = "protein")


## ranking by gain, color by pubs
plot_data <- olink_soma_full %>%
    arrange(desc(gain)) %>%
    mutate(index = row_number(gain))



##### I think bins looks better
nbins <- 220
tier_levels <- c("1","2","3","4")
tier_cols <- c("1"="#26e3fc", "2"="#74c365", "3"="#fc9c26", "4"="#ff5959")

# based on index, cut into bins to represent game (so bars have areas)
pd <- plot_data %>%
    mutate(tiers = factor(tiers, levels = tier_levels)) %>%
    arrange(index) %>%
    mutate(xbin = cut_number(index, n = nbins))

# define the bar width and height
slab_base <- pd %>%
    group_by(xbin) %>%
    # the start and end of each bar/bin
    summarise(xmin = min(index), xmax = max(index),
              height = max(gain), .groups = "drop")

## for each bar, calculate the mean number of publications
tier_stats <- pd %>%
    group_by(xbin, tiers) %>%
    summarise(n = n(),
              heat = sum(tot_pubs),
              .groups = "drop")

# combine, turn proportions into stacked ymin/ymax per tier within each slab
slabs_tier <- slab_base %>%
    left_join(tier_stats, by = "xbin") %>%
    # for each beam, calculate the proportion of proteins in a certain tier compared to the total number of proteins in the bin
    group_by(xbin) %>%
    mutate(prop = n / sum(n)) %>%
    arrange(tiers, .by_group = TRUE) %>%
    # calculate the height of the each tier, the higher the bar the more proteins in the bin 
    mutate(cumprop = cumsum(prop),
           ymin = height*(cumprop - prop),
           ymax = height*cumprop) %>%
    ungroup()

global_max_pubs <- max(tier_stats$heat, na.rm = T)

### create a new ggplot object
p <- ggplot()

for (t in tier_levels) {
    # for each tier
    dat <- dplyr::filter(slabs_tier, tiers == t)
    if (nrow(dat) > 0) {
        p <- p +
            geom_rect(data = dat,
                      aes(xmin = xmin, xmax = xmax, 
                          ymin = ymin, ymax = ymax, fill = heat),
                      color = NA, alpha = 1, inherit.aes = FALSE) +
            
            scale_fill_gradient(
                # reference to the color list and darken the color if the number increase
                # the gradient is specific to tiers
                low  = tier_cols[[t]],
                high = colorspace::darken(tier_cols[[t]], 0.6),
                limits = c(0, global_max_pubs),
                name = paste0("Sum of Publications \n(Tier ", t, ")"),
                guide = guide_colorbar(order = match(t, tier_levels), direction = "vertical", 
                                       title.position = "top")
                # guide = guide_colorbar(order = match(t, tier_levels), direction = "vertical", 
                #                        title.position = "top", barwidth = unit(4, "cm"))
            ) +
            ggnewscale::new_scale_fill()
    }
}

olink_soma_plot <- p +
    geom_line(data = plot_data, aes(x = index, y = gain),
              linewidth = 0.8, color = "black") +
    scale_y_continuous(limits = c(-0.2,0.5), breaks = c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5)) +
    labs(x = "Index", y = "Gain") +
    theme_minimal() +
    theme(legend.direction = "vertical", legend.position = "right",
          legend.title = element_text(size = 13, face = "bold"),
          axis.title = element_text(size = 13, face = "bold"),
          axis.text = element_text(size = 12, color="black"),
          )
olink_soma_plot



# combining ---------------------------------------------------------------
so_os_stacked <- ggarrange(soma_olink_plot, olink_soma_plot, labels = "AUTO",
                           common.legend = T, legend = "right", align="hv", nrow=2)
so_os_stacked   

so_os_stacked_forpre <- ggarrange(soma_olink_plot, olink_soma_plot, labels = "AUTO",
                           common.legend = T, legend = "right", align="hv", nrow=1) + 
    theme(legend.direction = "vertical", legend.position = "right",
          legend.title = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 9, color="black"),
    )
so_os_stacked_forpre  


# ggsave("Results/figures/soma_olink_rank_by_gain_color_by_pubs.png", soma_olink_plot, dpi = 250, width = 8, height = 5.5, bg = "white")
# ggsave("Results/figures/olink_soma_rank_by_gain_color_by_pubs.png", olink_soma_plot, dpi = 250, width = 8, height = 5.5, bg = "white")

ggsave("Results/figures/so_os_rank_by_gain_color_by_pubs_stacked.png", so_os_stacked, dpi = 250, width = 8, height = 11, bg = "white")
ggsave("Results/figures/forpre_so_os_rank_by_gain_color_by_pubs_stacked.png", so_os_stacked_forpre, dpi = 250, width = 11.5, height = 9, bg = "white")


