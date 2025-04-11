library(tidyverse)
library(rhdf5)
library(patchwork)

source("R/plot_theme.R")
source("R/read_seasonality_data.R")




orbit_labels <- expand_grid(
  period = 1:3,
  n_peaks = 1:7
) %>%
  mutate(label = str_c(n_peaks, "/", period),
         label_fct = factor(label))

plot_data_peaks <- plot_data %>% 
  filter(periodic, round(period / 365) > 0) %>%
  mutate(label = case_when(
    round(period / 365) > 3 ~ "Period > 3",
    TRUE ~ str_c(round(peak_density * (period / 365)), "/", round(period/365))
  )) %>%
  # filter(label %in% orbit_labels$label) %>%
  mutate(label_fct = factor(label, levels = orbit_labels$label))




orbit_colours <- c(
  RColorBrewer::brewer.pal(7, "Blues"),
  RColorBrewer::brewer.pal(7, "Greens"),
  RColorBrewer::brewer.pal(7, "Oranges")
) %>%
  `names<-`(orbit_labels$label)

p_orbits <- ggplot() +
  
  geom_tile(aes(x = -1, y = -1, fill = label), orbit_labels) +
  
  geom_tile(aes(x = eta, y = r, fill = label),
            plot_data_peaks) +
  
  plot_annotations +
  
  scale_fill_manual(
    name = NULL,
    values = orbit_colours,
    labels = orbit_labels
  ) +
  
  coord_fixed(ratio = 16.66, ylim = c(0, 0.03), xlim = c(-0.07, 0.5)) +
  xlab("Seasonality strength <i>η</i>") + ylab("Antibody decay rate <i>r</i>") +
  
  plot_theme_paper +
  
  guides(fill = guide_legend(direction = "horizontal", reverse = FALSE, byrow = FALSE, ncol = 8)) +
  theme(legend.position = "bottom") +
  
  ggtitle(NULL, "<b>A</b> — Number of peaks per period (period ≤ 3yr)")

p_orbits

p_peak_density <- ggplot() +
  
  geom_tile(aes(x = eta, y = r, fill = peak_density),
            plot_data) +
  
  plot_annotations +
  
  scale_fill_viridis_c(option = 5, name = "Peaks per year") +
  
  coord_fixed(ratio = 16.66, ylim = c(0, 0.03), xlim = c(-0.07, 0.5)) +
  xlab("Seasonality strength <i>η</i>") + ylab("Antibody decay rate <i>r</i>") +
  
  plot_theme_paper +
  
  guides(fill = guide_colourbar(barwidth = 15)) +
  theme(legend.position = "bottom") +
  
  ggtitle(NULL, "<b>B</b> — Average number of peaks per year")


p_full <- p_orbits | p_peak_density




ggsave(
  "results/results_supp_seasonality_density.png",
  p_full,
  device = png,
  width = 14, height = 8,
  bg = "white"
)
