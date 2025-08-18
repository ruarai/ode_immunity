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
  )#,
  # geom_point(
  #   aes(x = eta, y = r),
  #   plot_data_example_points,
  #   colour = "black",
  #   size = 1.4,
  #   stroke = 1
  # ),
  # geom_point(
  #   aes(x = eta, y = r),
  #   plot_data_example_points,
  #   colour = "white",
  #   size = 0.7,
  #   stroke = 0.5
  # ),
  # geom_label(
  #   aes(x = eta + 0.01, y = r + 0.0015, label = label),
  #   plot_data_example_points,
  #   label.r = unit(0.1, "cm"),
  #   label.size = 0,
  #   fill = shades::opacity("white", 0.8)
  # )
)

period_cols <- viridis::inferno(n = 8, direction = -1, begin = 0.1)

ggplot() +
  
  geom_tile(aes(x = eta, y = r, fill = period),
            plot_data_periodic) +
  geom_tile(aes(x = eta, y = r, fill = factor(9)),
            plot_data_quasiperiodic) +
  geom_tile(aes(x = eta, y = r, fill = factor(10)),
            plot_data_chaotic) +
  # geom_tile(aes(x = eta, y = r, fill = factor(11)),
  #           plot_data_empty) +
  
  plot_annotations +
  
  scale_fill_manual(
    name = "Period",
    values = c(period_cols, "#0076BC", "#9ED7F3", "white") %>% 
      `names<-`(1:11),
    
    labels = c("1 year", str_c(2:7, " years"), "≥8 years", "Quasiperiodic", "Chaotic", "Unclassified") %>%
      `names<-`(1:11),
    
    breaks = c(9, 1:4, 10, 5:8, 11)
  ) +
  
  
  coord_fixed(ratio = 16.66, ylim = c(0, 0.03)) +
  scale_x_continuous(breaks = seq(0, 0.5, by = 0.1)) +
  xlab("Seasonal forcing strength <i>η</i>") + ylab("Antibody decay rate <i>r</i> (days<sup>-1</sup>)") +
  guides(fill = guide_legend(nrow = 3, ncol = 5),
         colour = guide_none()) +
  
  plot_theme_paper +
  theme(legend.position = "none", legend.byrow = TRUE,
        legend.key = element_rect(colour = "grey80", linewidth = 0.5))


ggsave(
  "results/grid.png",
  width = 7, height = 7
)
