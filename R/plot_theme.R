

plot_theme <- list(
  theme_minimal(),
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(),
        text = element_text(family = "Helvetica", colour = "black", size = 20),
        line = element_line(linewidth = 0.7),
        axis.text.x = element_markdown(margin = margin(t = 0.3, unit = "cm"), colour = "black"),
        axis.text.y = element_markdown(margin = margin(r = 0.3, unit = "cm"), colour = "black"),
        plot.margin = margin(t = 0.5, l = 0.5, b = 0.5, r = 0.5, unit = "cm"),
        axis.ticks.length=unit(-0.1, "cm"),
        axis.title.y = element_markdown(),
        axis.title.x = element_markdown())
)
