##### 01 - library #####
library(dplyr)
library(data.table)
library(ggplot2)
library(ggnewscale)
library(colorspace)
library(wesanderson)

gc()
rm(list = ls())

setwd("...")



##### 02 - files #####
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
                name = paste0("Sum of Pubs (Tier ", t, " Proteins)"),
                guide = guide_colorbar(order = match(t, tier_levels), direction = "vertical", title.position = "top")
            ) +
            ggnewscale::new_scale_fill()
    }
}

p +
    geom_line(data = plot_data, aes(x = index, y = gain),
              linewidth = 0.8, color = "black") +
    scale_y_continuous(limits = c(-0.2,0.5), breaks = c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5)) +
    labs(x = "Index", y = "Gain") +
    theme_minimal() +
    theme(legend.position = "right",
          legend.box = "vertical")



ggsave("Results/figures/olink_soma_rank_by_gain_color_by_pubs.png", dpi = 300, width = 13, height = 8, bg = "white")







# # archived versions -------------------------------------------------------
# global_max_totpub <- max(olink_soma_full$tot_pubs)
# ggplot(data = plot_data, aes(x = index, y = gain)) +
#     geom_col()+
#     geom_col(data = subset(plot_data, tiers == "1"),
#              aes(fill = tot_pubs), width = 3) +
#     scale_fill_gradient(
#         low = "#fc9c26",
#         high = darken("#fc9c26", 0.8),
#         limits = c(0, global_max_totpub),
#         name = "# pubs (Tier 1 Protein)",
#         guide = guide_colorbar(order = 1, direction = "vertical")
#     ) +
#     new_scale_fill() +
#     geom_col(data = subset(plot_data, tiers == "2"),
#              aes(fill = tot_pubs), width = 3) +
#     scale_fill_gradient(
#         low = "#74c365",
#         high = darken("#74c365", 0.8),
#         limits = c(0, global_max_totpub),
#         name = "# pubs (Tier 2 Protein)",
#         guide = guide_colorbar(order = 2, direction = "vertical")
#     ) +
#     new_scale_fill() +
#     geom_col(data = subset(plot_data, tiers == "3"),
#              aes(fill = tot_pubs), width = 3) +
#     scale_fill_gradient(
#         low = "#ff5959",
#         high = darken("#ff5959", 0.8),
#         limits = c(0, global_max_totpub),
#         name = "# pubs (Tier 3 Protein)",
#         guide = guide_colorbar(order = 3, direction = "vertical")
#     ) +
#     new_scale_fill() +
#     geom_col(data = subset(plot_data, tiers == "4"),
#              aes(fill = tot_pubs), width = 3) +
#     scale_fill_gradient(
#         low = "#26e3fc",
#         high = darken("#26e3fc", 0.8),
#         limits = c(0, global_max_totpub),
#         name = "# pubs (Tier 4 Protein)",
#         guide = guide_colorbar(order = 4, direction = "vertical")
#     ) +
#     theme_minimal()
# 
# ggsave("Results/figures/olink_soma_rank_by_gain_color_by_pubs_nobins.png", dpi = 300, width = 11, height = 7, bg = "white")









# # individual bar but separated by quadrants -------------------------------
# ## need to reindex within each tier for the individual
# plot_data <- olink_soma_full %>%
#     group_by(tiers) %>%
#     arrange(desc(gain)) %>%
#     mutate(
#         index = row_number(gain),
#         tiers = case_when(
#             tiers == "1" ~ "Tier 1",
#             tiers == "2" ~ "Tier 2",
#             tiers == "3" ~ "Tier 3",
#             tiers == "4" ~ "Tier 4",
#             TRUE ~ NA # there should be any but this is my habit
#         )
#     )
# 
# 
# ggplot(data = plot_data, aes(x = index, y = gain, fill = tot_pubs)) +
#     facet_wrap(~tiers, nrow = 4, scales = "free_x") +
#     geom_col()+
#     geom_col(data = subset(plot_data, tiers == "Tier 1"),
#              aes(fill = tot_pubs), width = 1) +
#     scale_fill_gradient(
#         low = "#fc9c26",
#         high = darken("#fc9c26", 0.8),
#         limits = c(0, global_max_totpub),
#         name = "# pubs (Tier 1 Protein)",
#         guide = guide_colorbar(order = 1, direction = "vertical")
#     ) +
#     new_scale_fill() +
#     geom_col(data = subset(plot_data, tiers == "Tier 2"),
#              aes(fill = tot_pubs), width = 1) +
#     scale_fill_gradient(
#         low = "#74c365",
#         high = darken("#74c365", 0.8),
#         limits = c(0, global_max_totpub),
#         name = "# pubs (Tier 2 Protein)",
#         guide = guide_colorbar(order = 2, direction = "vertical")
#     ) +
#     new_scale_fill() +
#     geom_col(data = subset(plot_data, tiers == "Tier 3"),
#              aes(fill = tot_pubs), width = 1) +
#     scale_fill_gradient(
#         low = "#ff5959",
#         high = darken("#ff5959", 0.8),
#         limits = c(0, global_max_totpub),
#         name = "# pubs (Tier 3 Protein)",
#         guide = guide_colorbar(order = 3, direction = "vertical")
#     ) +
#     new_scale_fill() +
#     geom_col(data = subset(plot_data, tiers == "Tier 4"),
#              aes(fill = tot_pubs), width = 1) +
#     scale_fill_gradient(
#         low = "#26e3fc",
#         high = darken("#26e3fc", 0.8),
#         limits = c(0, global_max_totpub),
#         name = "# pubs (Tier 4 Protein)",
#         guide = guide_colorbar(order = 4, direction = "vertical")
#     ) +
#     theme_minimal()
# 
# 
# 
# ggsave("Results/figures/olink_soma_rank_by_gain_color_by_pubs_nobins_quadrants.png", dpi = 300, width = 11, height = 7, bg = "white")



