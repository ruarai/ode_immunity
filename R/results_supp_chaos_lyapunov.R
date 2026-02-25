



library(tidyverse)
library(rhdf5)
library(patchwork)

source("R/plot_theme.R")
source("R/read_seasonality_data.R")

plot_data <- read_seasonality_data("data/period_over_grid.jld2")

inf_threshold <- 1e-6

plot_data_periodic <- plot_data %>%
  filter(eta > 0) %>% 
  filter(periodic) %>% 
  mutate(period = period / 365,
         period = pmin(period, 4),
         period = factor(round(period)))

plot_data_chaotic <- plot_data %>% filter(chaotic)

plot_data_empty <- plot_data %>% filter(!periodic, !quasiperiodic, !chaotic)

plot_data_eta_zero <- plot_data %>% filter(eta == 0)
plot_data_eta_zero_periodic <- plot_data %>% filter(eta == 0, r < 0.065)

bifur_zero <- plot_data_eta_zero %>% filter(inf_diff < 1e-3) %>% pull(r) %>% head(1)

plot_data_quasiperiodic <- plot_data %>% filter(quasiperiodic, r < bifur_zero)

grid_factor <- 110
r_ratio <- 16.66
plot_data_low <- plot_data %>% filter(inf_min <= inf_threshold) %>% 
  mutate(f = if_else(
    xor(floor(eta * grid_factor) %% 2 == 0, floor(r * grid_factor * r_ratio) %% 2 == 0), factor(8), factor(7)
  ))


hdf5_file_lya <- "data/chaos_lyapunov.jld2"

lya_tbl <- tibble(
  eta = h5read(hdf5_file_lya, "x_eta_full"),
  r = h5read(hdf5_file_lya, "x_r_full"),
  lya = h5read(hdf5_file_lya, "y_lyapunov")
)

plot_data_chaos_lya <- lya_tbl %>% filter(lya > 1e-4) %>% 
  left_join(plot_data_periodic)



year_stops <- c(1/2, 2/3, 1)
year_marks <- approxfun(plot_data_eta_zero_periodic$period, plot_data_eta_zero_periodic$r)(365 * year_stops)
plot_data_year_marks <- tibble(r_0 = year_marks, year = year_stops) %>%
  mutate(year_label = str_c(scales::label_comma()(year), " yr"))


plot_data_example_points <- tribble(
  ~eta, ~r, ~label, ~colour,
  0.0, 0.018, "i", factor(0),
  0.07, 0.018, "ii", factor(4.5),
  0.2, 0.018, "iii", factor(1),
  0.27, 0.018, "iv", factor(2),
  0.37, 0.018, "v", factor(9)
)
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


p_period <- ggplot() +
  geom_tile(aes(x = eta, y = r, fill = factor(5)),
            plot_data_quasiperiodic) +
  geom_tile(aes(x = eta, y = r, fill = factor(6)),
            plot_data_chaotic) +
  
  geom_tile(aes(x = eta, y = r, fill = period),
            plot_data_periodic %>% filter(period != 0)) +
  geom_tile(aes(x = eta, y = r, fill = factor(8)),
            plot_data_empty) +
  # geom_tile(aes(x = eta, y = r, fill = factor(f)),
  #           plot_data_low) +
  
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
  theme(legend.position = "bottom", legend.byrow = TRUE,
        legend.key = element_rect(colour = "grey80", linewidth = 0.5))


p_period_lya <- ggplot() +
  geom_tile(aes(x = eta, y = r, fill = factor(5)),
            plot_data_quasiperiodic) +
  geom_tile(aes(x = eta, y = r, fill = factor(6)),
            plot_data_chaos_lya) +
  
  geom_tile(aes(x = eta, y = r, fill = period),
            plot_data_periodic %>% filter(period != 0)) +
  geom_tile(aes(x = eta, y = r, fill = factor(8)),
            plot_data_empty) +
  # geom_tile(aes(x = eta, y = r, fill = factor(f)),
  #           plot_data_low) +
  
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.01, fill = "grey50", alpha = 0.3) +
  
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
  theme(legend.position = "bottom", legend.byrow = TRUE,
        legend.key = element_rect(colour = "grey80", linewidth = 0.5))


p_period + ggtitle(NULL, "<b>A</b> — Dynamics (using Melbourne 0-1 test)") | 
  p_period_lya + ggtitle(NULL, "<b>B</b> — Dynamics (using maximum Lyapunov exponent)")


ggsave(
  "results/results_supp_chaos_lyapunov.png",
  device = png,
  width = 13, height = 9,
  bg = "white"
)



