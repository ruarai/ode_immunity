
library(tidyverse)
library(rhdf5)
library(patchwork)


source("R/plot_theme.R")

t_seq <- h5read("data/ex_period.jld2", "t_seq")
y <- h5read("data/ex_period.jld2", "y_seasonal")

y_tbl <- y %>%
  reshape2::melt(varnames = c("ix", "t"), value.name = "prevalence") %>%
  as_tibble() %>% 
  mutate(t = t_seq[t]) %>%
  filter(ix == 34, t > 40000) %>%
  mutate(t_mod = t %% 365)


xy_to_days <- function(x, y) {
  theta <- atan2(y, x)
  if_else(theta < 0, theta + 2 * pi, theta) / (2 * pi) * 365
}


tibble(
  t = 0:1000
) %>%
  mutate(
    theta = t * pi * 2 / 365,
    x = cos(theta), y = sin(theta),
    theta = atan2(y, x),
    theta = if_else(theta < 0, theta + 2 * pi, theta),
    t2 = xy_to_days(x, y)
  ) %>% 
  ggplot() +
  geom_line(aes(x = t, y = y))



y_tbl %>%
  mutate(theta = t_mod / 365 * 2 * pi,
         x = cos(theta), y = sin(theta)) %>%
  summarise(mean_x = weighted.mean(x, prevalence), mean_y = weighted.mean(y, prevalence)) %>%
  mutate(mean_day = xy_to_days(mean_x, mean_y),
         mean_mag = sqrt(mean_x^2 + mean_y^2), mean_var = 1 - mean_mag)

ggplot() +
  geom_path(aes(x = t_mod, y = prevalence),
            y_tbl) +
  coord_radial(expand = FALSE, inner.radius = 0.7, start = -pi / 2, direction = -1) +
  
  plot_theme_paper +
  
  scale_x_continuous(breaks = scales::breaks_extended(n = 10)) +
  
  xlab("Time of year (days)") +
  ylab("Infection prevalence") +
  
  theme(axis.line.x = element_blank())


ggsave(
  "results/radial_mean_template.pdf",
  device = pdf,
  width = 7, height = 7,
  bg = "white"
)


hex_cols <- twilight_palette %>% 
  shifter(n = 240) %>%
  `[`(30:500) %>%
  shifter(n = floor(471 / 2))


ggplot() +
  geom_tile(aes(x = t, y = y, fill = t, colour = t),
            expand_grid(t = 1:365, y = seq(0, 1, by = 0.01)), linewidth = 1) +
  coord_radial(expand = FALSE, inner.radius = 0) +
  
  plot_theme_paper +
  
  scale_x_continuous(breaks = c(0, 1)) +
  
  scale_fill_gradientn(colours = twilight_palette) +
  scale_colour_gradientn(colours = twilight_palette) +
  
  xlab("Time of year (days)") +
  ylab("Infection prevalence") +
  
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank())


ggplot() +
  geom_blank(aes(x = t, y = 1.1),
            tibble(t = 1:365)) +
  geom_tile(aes(x = t, y = y, fill = y, colour = y),
            expand_grid(t = 1:365, y = seq(0, 1, by = 0.01)), linewidth = 1) +
  coord_radial(expand = FALSE, inner.radius = 0) +
  
  plot_theme_paper +
  
  scale_x_continuous(breaks = c(0, 1)) +
  
  scale_fill_stepsn(
    name = NULL,
    colours = RColorBrewer::brewer.pal(9, "Blues"),
    breaks = seq(0, 1, length.out = 11),
    labels = c("0", "", "", "", "", "0.5", "", "", "", "", "1.0"),
    limits = c(0, 1)
  ) +
  
  scale_colour_stepsn(
    name = NULL,
    colours = RColorBrewer::brewer.pal(9, "Blues"),
    breaks = seq(0, 1, length.out = 11),
    labels = c("0", "", "", "", "", "0.5", "", "", "", "", "1.0"),
    limits = c(0, 1)
  ) +

  
  xlab("Time of year (days)") +
  ylab("Infection prevalence") +
  
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank())

