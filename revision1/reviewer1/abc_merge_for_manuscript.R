##### 01 - library #####
library(ggplot2)
library(patchwork)

gc()
rm(list = ls())

setwd("...")


##### 02 - load plot objects saved by a_ and c_ #####
rds_dir <- "Results/revision1/figures/rds"

p_sensitivity      <- readRDS(file.path(rds_dir, "p_sensitivity.rds"))
p_tier_proportions <- readRDS(file.path(rds_dir, "p_tier_proportions.rds"))
p_stability        <- readRDS(file.path(rds_dir, "p_stability.rds"))
p_ci_ribbon        <- readRDS(file.path(rds_dir, "p_ci_ribbon.rds"))


##### 03 - enforce shared display-item styling #####
strip_and_style <- function(p, base_size = 11) {
    p +
        labs(title = NULL) +
        theme_classic(base_size = base_size) +
        theme(legend.position = "bottom")
}

p_sensitivity      <- strip_and_style(p_sensitivity)
p_tier_proportions <- strip_and_style(p_tier_proportions, base_size = 10)
p_stability        <- strip_and_style(p_stability)
p_ci_ribbon        <- strip_and_style(p_ci_ribbon)


##### 04 - combine figure #####
## top row: threshold sensitivity analysis (a_)
## bottom row: bootstrap uncertainty / CI analysis (c_)
combined_robustness_plot <- (p_tier_proportions) /
    (p_stability | p_ci_ribbon) +
    plot_annotation(tag_levels = "A") &
    theme(
        plot.tag = element_text(face = "bold", size = 13),
        legend.position = "bottom"
    )

combined_robustness_plot


##### 05 - save #####
ggsave(
    "Results/revision1/figures/tiering_robustness_combined_ABCD.png",
    combined_robustness_plot,
    dpi = 300, width = 14, height = 12, bg = "white"
)
