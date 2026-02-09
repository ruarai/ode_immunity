plot_annotations <- list(
  geom_segment(
    # aes(x = r_0, y = 0.0, xend = r_0 + 0.001, yend = -0.01),
    aes(x = -0.003, y = r_0, xend = -0.01, yend = r_0),
    plot_data_year_marks
  ),
  
  annotate("linerange", x = -0.0065, ymin = bifur_zero, ymax = 0.03),
  annotate(
    "segment",
    x = -0.003,
    y = bifur_zero,
    xend = -0.01,
    yend = bifur_zero
  ),
  geom_text(
    aes(x = -0.08, y = r_0 + 0.0002, label = year_label),
    hjust = 0,
    plot_data_year_marks,
    size = 4.5
  ),
  annotate(
    "text",
    x = -0.08,
    y = 0.0275,
    label = "Fixed\npoint",
    hjust = 0,
    size = 4.5
  )
)

period_cols <- viridis::inferno(n = 5, direction = -1, begin = 0.1)[1:4]
ggplot() +
  
  # geom_tile(aes(x = eta, y = r, fill = period),
  #           plot_data_periodic %>% filter(period != 0)) +
  # geom_tile(aes(x = eta, y = r, fill = factor(5)),
  #           plot_data_quasiperiodic) +
  # geom_tile(aes(x = eta, y = r, fill = factor(6)),
  #           plot_data_chaotic) +
  geom_tile(aes(x = eta, y = r, fill = factor(8)),
            plot_data_empty) +
  geom_tile(aes(x = eta, y = r, fill = factor(f)),
            plot_data_low) +
  
  plot_annotations +
  
  scale_fill_manual(
    name = "Period",
    values = c(period_cols, "#0076BC", "#9ED7F3", "grey90", "white") %>% 
      `names<-`(1:8),
    
    labels = c("1 year", str_c(2:3, " years"), "≥4 years", "Quasiperiodic", "Chaotic", "Low minimum prevalence", "Unclassified") %>%
      `names<-`(1:8),
    
    breaks = c(1:2, 5, 6, 3:4, 7, 8)
  ) +
  
  
  coord_fixed(ratio = 16.66, ylim = c(0.000, 0.03)) +
  xlab("Seasonal forcing strength <i>η</i>") + ylab("Effective antibody decay rate <i>r</i>") +
  guides(fill = guide_legend(nrow = 3, ncol = 4),
         colour = guide_none()) +
  
  plot_theme_paper +
  theme(legend.position = "none", legend.byrow = TRUE,
        legend.key = element_rect(colour = "grey80", linewidth = 0.5))



ggsave(
  "results/grid.png",
  device = png, #bg = "white",
  width = 7, height = 7
)

color_tbl <- read_delim("data/lipari.txt", col_names = c("R", "G", "B")) %>%
  mutate(c = rgb(R, G, B))

min_cols <- color_tbl$c[seq(1, nrow(color_tbl), length.out = 128)]

p_min <- ggplot()  +
  geom_tile(aes(x = eta, y = r, fill = pmax(-14, log10(inf_min))),
            plot_data) +
  
  plot_annotations +
  scale_colour_manual(
    values = rep("white", 11),
    breaks = c(0:4, 4.5, 5:9)
  ) + 
  
  scale_fill_stepsn(
    # colours = rev(colorspace::diverging_hcl(n = 20, h = c(240, 15), c = c(60, 80), l = c(75, 5), power = c(1.2, 1.5))),
    colours = min_cols,
    name = "Minimum\ninfection\nprevalence (log10)",
    limits = c(-15, -1),
    breaks = seq(-15, -1, 0.5),
    labels = c("<-15", "", "", "", "-13", "", "", "", "-11", "", "", "", "-9", 
               "", "", "", "-7", "", "", "", "-5", "", "", "", "-3", "", "", 
               "", "-1")
  ) +
  
  coord_fixed(ratio = 16.66, ylim = c(0, 0.03)) +
  xlab("Seasonal forcing strength <i>η</i>") + ylab("Effective antibody decay rate <i>r</i>")  +
  
  plot_theme_paper +
  guides(fill = guide_colourbar(barwidth = 15),
         colour = guide_none()) +
  theme(legend.position = "none",
        legend.title = element_text(margin = margin(r = 0.7, unit = "cm")))

p_min

ggsave(
  "results/grid_2.png",
  device = png, bg = "white",
  width = 7, height = 7
)
