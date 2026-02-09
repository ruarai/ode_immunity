
library(tidyverse)
library(rhdf5)
library(patchwork)

source("R/plot_theme.R")

source("R/read_seasonality_data.R")

xy_to_days <- function(x, y) {
  theta <- atan2(y, x)
  if_else(theta < 0, theta + 2 * pi, theta) / (2 * pi) * 365
}

plot_data <- read_seasonality_data("data/period_over_grid.jld2") %>%
  mutate(season_day = xy_to_days(season_x, season_y),
         season_mag = sqrt(season_x ^ 2 + season_y ^ 2),
         season_var = 1 - season_mag)


shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

# hex_cols <- read_table("data/vikO25.txt", skip = 2, col_names = c("R", "G", "B", "name", "hex")) %>%
#   pull(hex) %>%
#   `[`(seq(1, 24, by = 2)) %>%
#   rev() %>% 
#   shifter(n = 10)


hex_cols <- twilight_palette %>% 
  shifter(n = 240) %>%
  `[`(30:500) %>%
  shifter(n = floor(471 / 2))


example_points <- tribble(
  ~eta, ~r, ~label,
  0.25, 0.01, "iv",
  0.25, 0.013, "iii",
  0.25, 0.0165, "ii",
  0.05, 0.0165, "i"
)

plot_annotations <- list(
  geom_point(aes(x = eta, y = r), example_points,
             colour = "black"),
  geom_point(aes(x = eta, y = r), example_points,
             colour = "white", size = 0.5),
  
  geom_label(aes(x = eta + 0.015, y = r - 0.0015, label = label), example_points,
             label.r = unit(0.1, "cm"), label.size = 0, fill = shades::opacity("white", 0.8))
)


div_factor <- 365 / 365

p_bias <- plot_data %>% 
  filter(eta > 0) %>% 
  mutate(season_day = floor(season_day / div_factor) * div_factor) %>% 
  ggplot()  +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.005, fill = "white")+
  geom_tile(aes(x = eta, y = r, fill = season_day)) +
  
  coord_fixed(ratio = 16.66, ylim = c(0, 0.03)) +
  xlab("Seasonal forcing strength <i>η</i>") + ylab("Effective antibody decay rate <i>r</i>") +
  
  scale_fill_gradientn(colours = hex_cols,
                       limits = c(0, 365),
                       breaks = 365 * c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("Jan\nWinter", "Apr\nSpring",
                                  "Jul\nSummer", "Oct\nAutumn",
                                  "Jan\nWinter"),
                       name = "Mean\n(time of year)") +
  
  
  plot_theme_paper +
  
  guides(fill = guide_colourbar(barwidth = 15)) +
  theme(legend.position = "bottom",
        legend.title = element_text(hjust = 0.5, margin = margin(r = 1.5, unit = "cm")),
        plot.title = element_markdown(size = 17)) +
  
  ggtitle("Seasonal bias in infection incidence", "<b>A</b> — Circular mean (time of year)")

p_bias

p_variance <- ggplot()  +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.005, fill = "white")+
  geom_tile(aes(x = eta, y = r, fill = season_var),
            plot_data) +
  
  coord_fixed(ratio = 16.66, ylim = c(0, 0.03)) +
  xlab("Seasonal forcing strength <i>η</i>") + ylab("Effective antibody decay rate <i>r</i>") +
  
  scale_fill_distiller(name = "Variance",
                       limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.25),
                       direction = 1) +
  
  plot_theme_paper +
  
  guides(fill = guide_colourbar(barwidth = 15)) +
  theme(legend.position = "bottom",
        legend.title = element_markdown(hjust = 0.5, margin = margin(r = 1, unit = "cm")),
        plot.title = element_markdown(size = 14)) +
  
  ggtitle("", "<b>B</b> — Circular variance")

p_variance

p_bias | p_variance




p <- p_bias | p_variance

ggsave(
  "results/results_supp_grid_seasonal_bias.png",
  p,
  device = png,
  width = 13, height = 8,
  bg = "white"
)


