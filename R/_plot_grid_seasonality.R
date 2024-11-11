
library(tidyverse)

library(rhdf5)


x_vals <- h5read("data/paper/period_over_grid.jld2", "x_vals")
y_inf_maxima <- h5read("data/paper/period_over_grid.jld2", "y_inf_maxima")
y_period <- h5read("data/paper/period_over_grid.jld2", "y_period")

plot_data <- tibble(
  eta = x_vals[1, ], rho = x_vals[2, ],
  min = y_inf_maxima[, 1], max = y_inf_maxima[, 2],
  period = y_period[,1], period_sd = y_period[,2], period_n = y_period[,3]
) %>%
  mutate(diff = max - min)




plot_data_periodic <- plot_data %>%
  filter(eta > 0) %>% 
  filter(period_sd < 1, period_n >= 5) %>% 
  mutate(period = period / 365,
         period = pmin(period, 8),
         period = factor(round(period)))

plot_data_quasiperiodic <- plot_data %>%
  filter(period_sd >= 1, period_n > 5)

plot_data_min <- plot_data %>%
  filter(min < 1e-7)

plot_data_eta_zero <- plot_data %>% filter(eta == 0)

bifur_zero <- plot_data_eta_zero %>% filter(diff < 1e-3) %>% pull(rho) %>% head(1)

year_stops <- c(1/2, 2/3, 1, 3/2, 2, 3, 4)
year_marks <- approxfun(plot_data_eta_zero$period, plot_data_eta_zero$rho)(365 * year_stops)
plot_data_year_marks <- tibble(rho_0 = year_marks, year = year_stops) %>%
  mutate(year_label = str_c(scales::label_comma()(year), " yr"))

plot_data_example_points <- tribble(
  ~eta, ~rho,
  0.3, 0.0028,
  0.02, 0.0035,
  0.15, 0.0018,
  0.32, 0.0035
)

plot_annotations <- list(
  geom_point(aes(x = -0.01, y = rho_0), plot_data_year_marks, pch = "-", size = 6),
  
  annotate("linerange", x = -0.01, ymin = bifur_zero, ymax = 0.005),
  annotate("point", x = -0.01, y = bifur_zero, pch = "-", size = 6),
  geom_point(aes(x = eta, y = rho), plot_data_example_points, colour = "black", size = 1.4, stroke = 1),
  geom_point(aes(x = eta, y = rho), plot_data_example_points, colour = "white", size = 0.7, stroke = 0.5)
)

p_period <- ggplot() +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.005, fill = "white") +
  geom_tile(aes(x = eta, y = rho, fill = period),
            plot_data_periodic) +
  geom_tile(aes(x = eta, y = rho),
            plot_data_quasiperiodic,
            fill = "#4884CC") +
  
  plot_annotations +
  
  geom_text(aes(x = -0.07, y = rho_0, label = year_label), hjust = 0, plot_data_year_marks) +
  annotate("text", x = -0.07, y = 0.00465, label = "Fixed\npoint", hjust = 0) +
  
  scale_fill_manual(name = "Period",
                    values = viridis::inferno(n = 8, direction = -1, begin = 0.1),
                    labels = c(1:7, "≥8")) +

  # scale_fill_viridis_b(option = "inferno", direction = -1, breaks = 1:8, labels = c(1:7, "8")) +

  
  coord_fixed(ratio = 100) +
  xlab("Seasonality constant η") + ylab("Waning constant ρ") +
  
  plot_theme_paper +
  theme(legend.position = "bottom", legend.byrow = TRUE)

p_period

p_max <- plot_data %>%
  # filter(period_sd < 1, period_n >= 5, period / 365 < 10) %>%
  ggplot()  +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.005, fill = "white")+
  geom_tile(aes(x = eta, y = rho, fill = max)) +
  
  plot_annotations + 
  
  scale_fill_viridis_c(name = "Peak\ninfection\nprevalence", option = "mako") +
  coord_fixed(ratio = 100) +
  xlab("Seasonality constant η") + ylab("Waning constant ρ") +
  
  plot_theme_paper +
  guides(fill = guide_colourbar(barwidth = 15)) +
  theme(legend.position = "bottom")

p_period | p_max


ggsave(
  "results/results_grid_seasonality.png",
  device = png,
  width = 14, height = 7,
  bg = "white"
)

plot_data %>% 
  ggplot() +
  geom_tile(aes(x = eta, y = rho, fill = min < 1e-7))



plot_data %>%
  mutate(log_min = log10(min),
         log_min = floor(pmax(log_min, -16) / 1) * 1,
         ) %>% 
  ggplot() +
  geom_tile(aes(x = eta, y = rho, fill = log_min)) +
  
  plot_annotations +
  
  scale_fill_viridis_c(name = "Minimum\ninfection\nprevalence",
                       limits = c(-10, 0), breaks = seq(-10, 0, by = 2),
                       # labels = str_c("10<sup>", seq(-16, 0, by = 2), "</sup>"),
                       labels = parse(text = str_c("10^'", seq(-10, 0, by = 2), "'")),
                       na.value = 10^-10,
                       ) +
  
  coord_fixed(ratio = 100) +
  xlab("Seasonality constant η") + ylab("Waning constant ρ")  +
  
  plot_theme_paper +
  guides(fill = guide_colorbar(barwidth = 20)) +
  theme(legend.position = "bottom")









