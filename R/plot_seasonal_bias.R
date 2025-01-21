
library(tidyverse)
library(rhdf5)
library(patchwork)

source("R/plot_theme.R")

source("R/read_seasonality_data.R")

season_to_days <- function(y, x){
  t <- -(atan2(y, x) - pi / 2)
  if_else(t < 0, t + 2 * pi, t) / (2 * pi) * 365
}

tibble(
  d = 0:365
) %>%
  mutate(theta = d / 365 * 2 * pi,
         x = sin(theta), y = cos(theta),
         t = season_to_days(y, x)) %>% 
  ggplot() +
  geom_line(aes(x = theta, y = t))

plot_data <- read_seasonality_data("data/paper/period_over_grid.jld2") %>%
  mutate(season_day = season_to_days(season_y, season_x),
         season_mag = sqrt(season_x ^ 2 + season_y ^ 2),
         season_var = 1 - season_mag)


plot_data$season_theta


m <- 2
b <- 1
rescale_density <- function(x) {
  if_else(x > 0, plogis(m * qlogis(x) + b), -plogis(m * qlogis(-x) + b))
}


p_heatmap <- ggplot()  +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.005, fill = "white")+
  geom_tile(aes(x = eta, y = r, fill = rescale_density(season_y)),
            plot_data) +
  
  coord_fixed(ratio = 5, ylim = c(0, 0.1)) +
  xlab("Seasonality constant <i>η</i>") + ylab("Antibody decay rate <i>r</i>") +
  
  scale_fill_distiller(name = "Seasonal bias\nin infection\nincidence",
                       type = "div",
                       palette = "RdBu",
                       direction = 1,
                       limits = c(-1, 1),
                       breaks = c(-1, 0, 1),
                       labels = c("-1\nSummer", "0", "1\nWinter")) +
  
  # scale_fill_stepsn(
  #   colours = RColorBrewer::brewer.pal(11, name = "RdBu"),
  #   name = "Annual infection incidence",
  #   limits = c(-1, 1),
  #   breaks = seq(-1, 1, by = 0.25)
  # ) +

  plot_theme_paper +
  guides(fill = guide_colourbar(barwidth = 15)) +
  theme(legend.position = "none")

p_legend <- tibble(
  x = seq(-1.2, 1.2, by = 0.01)
) %>%
  mutate(density = rescale_density(x)) %>% 
  fill(density, .direction = "updown") %>% 
  ggplot() +
  geom_tile(aes(x = 1, y = x, fill = density)) +
  
  geom_hline(aes(yintercept = y), tibble(y = seq(-1, 1, by = 0.5)),
             colour = "white") +
  
  scale_fill_distiller(name = "Seasonal bias\nin infection\nincidence",
                       type = "div",
                       palette = "RdBu",
                       direction = 1) +
  
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5)) +
  
  plot_theme_paper + xlab(NULL) + ylab("Seasonal<br>Bias") +
  
  theme(legend.position = "none",
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_markdown(angle = 0, vjust = 0.5, size = 14))

(p_heatmap | p_legend) +
  plot_layout(widths = c(8, 1))



data_rgb <- read_table("data/vikO.txt", col_names = c("R", "G", "B"))


shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

hex_cols <- as.matrix(data_rgb) %>%
  shades::shade() %>%
  shifter(n = 40) %>%
  rev()

hex_cols <- read_table("data/vikO25.txt", skip = 2, col_names = c("R", "G", "B", "name", "hex")) %>%
  pull(hex) %>%
  `[`(seq(1, 24, by = 2)) %>%
  rev() %>% 
  shifter(n = 10)


hex_cols <- pals::kovesi.cyclic_wrwbw_40_90_c42_s25(n = 12)
shades::swatch(hex_cols)

hex_cols <- twilight_palette %>% 
  shifter(n = 240) %>%
  `[`(30:500) %>%
  shifter(n = floor(471 / 2))


div_factor <- 365 / 365

p_bias <- plot_data %>% 
  filter(eta > 0) %>% 
  mutate(season_day = floor(season_day / div_factor) * div_factor) %>% 
  ggplot()  +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.005, fill = "white")+
  geom_tile(aes(x = eta, y = r, fill = season_day)) +
  
  coord_fixed(ratio = 5, ylim = c(0, 0.1)) +
  xlab("Seasonality constant <i>η</i>") + ylab("Antibody decay rate <i>r</i>") +
  
  scale_fill_gradientn(colours = hex_cols,
                       limits = c(0, 365),
                       breaks = 365 * c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("Jan", "Apr", "Jul", "Oct", "Jan"),
                       name = "Mean\n(time of year)  ") +
  
  
  plot_theme_paper +
  guides(fill = guide_colourbar(barwidth = 15)) +
  theme(legend.position = "bottom",
        plot.title = element_markdown(size = 17)) +
  
  ggtitle("Seasonal bias in infection incidence", "<b>A</b> — Circular mean (time of year)")

p_variance <- ggplot()  +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.005, fill = "white")+
  geom_tile(aes(x = eta, y = r, fill = season_var),
            plot_data) +
  
  coord_fixed(ratio = 5, ylim = c(0, 0.1)) +
  xlab("Seasonality constant <i>η</i>") + ylab("Antibody decay rate <i>r</i>") +
  
  # scale_fill_distiller(limits = c(0, 1), direction = 1) +
  
  scale_fill_stepsn(
    name = "Variance  ",
    colours = RColorBrewer::brewer.pal(9, "Blues"),
    breaks = seq(0, 1, length.out = 11),
    labels = c("0", "", "", "", "", "0.5", "", "", "", "", "1.0"),
    limits = c(0, 1)
  ) +
  
  plot_theme_paper +
  guides(fill = guide_colourbar(barwidth = 15)) +
  theme(legend.position = "bottom",
        plot.title = element_markdown(size = 14)) +
  
  ggtitle("", "<b>B</b> — Circular variance")


p_bias | p_variance

p_legend <- tibble(
  x = seq(0, 365, by = div_factor)
) %>% 
  mutate(x_2 = if_else(x == 365, 0, x)) %>%
  ggplot() +
  geom_tile(aes(x = 1, y = x, fill = x_2)) +
  
  # geom_hline(aes(yintercept = y), tibble(y = seq(0, 365, by = 365 / 4)),
  #            colour = "white") +
  
  scale_fill_gradientn(colours = hex_cols,
                       limits = c(0, 365)) +
  
  scale_y_continuous(breaks = c(0, 365 * 0.25, 365 * 0.5, 365 * 0.75, 365),
                     labels = c("Jan", "Apr", "Jul", "Oct", "Jan")) +
  
  plot_theme_paper + xlab(NULL) + ylab(NULL) +
  
  theme(legend.position = "none",
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_markdown(angle = 0, vjust = 0.5, size = 14))

(p_heatmap | p_legend) +
  plot_layout(widths = c(15, 1))


p_bias +
  annotate("point", x = 0.015, y = 0.015)

