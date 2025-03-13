
library(tidyverse)
library(rhdf5)
library(patchwork)

source("R/plot_theme.R")

source("R/read_seasonality_data.R")

xy_to_days <- function(x, y) {
  theta <- atan2(y, x)
  if_else(theta < 0, theta + 2 * pi, theta) / (2 * pi) * 365
}

plot_data <- read_seasonality_data("data/paper/period_over_grid.jld2") %>%
  mutate(season_day = xy_to_days(season_x, season_y),
         season_mag = sqrt(season_x ^ 2 + season_y ^ 2),
         season_var = 1 - season_mag)


shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

hex_cols <- read_table("data/vikO25.txt", skip = 2, col_names = c("R", "G", "B", "name", "hex")) %>%
  pull(hex) %>%
  `[`(seq(1, 24, by = 2)) %>%
  rev() %>% 
  shifter(n = 10)


hex_cols <- twilight_palette %>% 
  shifter(n = 240) %>%
  `[`(30:500) %>%
  shifter(n = floor(471 / 2))


## TODO change these.
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
  xlab("Seasonality constant <i>η</i>") + ylab("Antibody decay rate <i>r</i>") +
  
  scale_fill_gradientn(colours = hex_cols,
                       limits = c(0, 365),
                       breaks = 365 * c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("Jan\nWinter", "Apr\nSpring",
                                  "Jul\nSummer", "Oct\nAutumn",
                                  "Jan\nWinter"),
                       name = "Mean\n(time of year)") +
  
  
  plot_theme_paper +
  
  plot_annotations +
  
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
  xlab("Seasonality constant <i>η</i>") + ylab("Antibody decay rate <i>r</i>") +
  
  scale_fill_distiller(name = "Variance",
                       limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.25),
                       direction = 1) +
  
  plot_theme_paper +
  
  plot_annotations + 
  
  guides(fill = guide_colourbar(barwidth = 15)) +
  theme(legend.position = "bottom",
        legend.title = element_markdown(hjust = 0.5, margin = margin(r = 1, unit = "cm")),
        plot.title = element_markdown(size = 14)) +
  
  ggtitle("", "<b>B</b> — Circular variance")

p_variance

p_bias | p_variance



x_eta <- h5read("data/paper/seasonality_bias_examples.jld2", "x_eta")
x_r <- h5read("data/paper/seasonality_bias_examples.jld2", "x_r")
y_inc <- h5read("data/paper/seasonality_bias_examples.jld2", "y_inc")
y_seasonality <- h5read("data/paper/seasonality_bias_examples.jld2", "y_seasonality")

x_labels <- str_c(
  "<b>", c("iv", "iii", "ii", "i"), ".</b> ",
  " <i>η</i> = ", scales::label_comma(accuracy = 0.01)(x_eta),
  ", <i>r</i> = ", scales::label_comma(accuracy = 0.001)(x_r), ""
)

c_levels <- 10 ^ seq(0, 8, by = 8 / 32)

t_ex_start <- 365 * 100
t_ex_end <- 365 * (100 + 8)
t_ex_yearly_end <- 365 * (100 + 60)


plot_data_ex_inc_year <- y_inc %>%
  reshape2::melt(c("i", "t"), value.name = "incidence") %>%
  mutate(eta = x_eta[i], r = x_r[i], label = x_labels[i]) %>%
  filter(t >= t_ex_start, t < t_ex_yearly_end) %>% 
  group_by(label, eta, r, t) %>%
  summarise(incidence = sum(incidence)) %>%
  
  mutate(year = floor(t / 365),
         t_year = t %% 365)

plot_data_ex_means <- tibble(
  x = y_seasonality[,1], y = y_seasonality[,2],
  eta = x_eta[1,], r = x_r[1,],
  label = x_labels
)  %>%
  mutate(season_day = xy_to_days(x, y),
         season_mag = sqrt(x ^ 2 + y ^ 2),
         season_var = 1 - season_mag) %>%
  mutate(season_var_label = scales::label_comma(accuracy = 0.01)(season_var),
         season_var_label = str_c("Variance = ", season_var_label))

p_ex_yearly_inc <- ggplot() +
  geom_line(aes(x = t_year, y = incidence, group = year),
            linewidth = 0.3, alpha = 0.5,
            plot_data_ex_inc_year) +
  
  geom_vline(aes(xintercept = season_day),
             linewidth = 1, alpha = 0.5, colour = colour_C,
             plot_data_ex_means) +
  
  geom_label(aes(x = 365 / 2, y = 0.03, label = season_var_label),
             label.size = 0, fill = "#FAEFF5",
             plot_data_ex_means) +
  
  facet_wrap(~label, ncol = 4) +
  
  xlab("Time of year") + ylab("Infection incidence") +
  
  scale_x_continuous(breaks = c(0, 90, 180, 270, 365),
                     # labels = c("0 (Jan)", "", "180 (Jul)", "", "365 (Jan)"),
                     labels = c("Jan", "",
                                "Jul", "", 
                                "Jan")) +
  
  scale_y_continuous(breaks = c(0, 0.02, 0.04)) +
  
  coord_cartesian(ylim = c(0, 0.04)) +
  
  plot_theme_paper +
  theme(panel.grid.major.x = element_line(colour = "grey50", linetype = "28", linewidth = 0.5),
        strip.text = element_markdown(),
        plot.title = element_markdown(size = 14),
        panel.spacing.x = unit(1.5, "cm")) +
  
  ggtitle(NULL, "<b>C</b> — Exemplar infection incidence by time of year")

p_ex_yearly_inc


p <- ((p_bias | p_variance) / p_ex_yearly_inc) +
  plot_layout(heights = c(3, 1))

ggsave(
  "results/results_grid_seasonal_bias.png",
  p,
  device = png,
  width = 13, height = 11,
  bg = "white"
)


