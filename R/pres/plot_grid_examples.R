library(tidyverse)
library(rhdf5)
library(patchwork)

source("R/plot_theme.R")
source("R/read_seasonality_data.R")



x_eta <- h5read("data/period_over_grid_examples.jld2", "x_eta")
x_r <- h5read("data/period_over_grid_examples.jld2", "x_r")
y_inf <- h5read("data/period_over_grid_examples.jld2", "y_inf")
y_sus <- h5read("data/period_over_grid_examples.jld2", "y_sus")

x_labels <- str_c(
  # "<b>",
  # c("i", "ii", "iii", "iv", "v"),
  # ".</b> ",
  c(
    "Zero seasonality",
    "Quasiperiodic",
    "Periodic (1 year)",
    "Periodic (2 years)",
    "Chaotic"
  ),
  " (<i>η</i> = ",
  scales::label_comma(accuracy = 0.01)(x_eta),
  ", <i>r</i> = ",
  scales::label_comma(accuracy = 0.001)(x_r),
  ")"
)

c_levels <- 2 ^ seq(0, 8, by = 8 / 32)

t_ex_start <- 365 * 100
t_ex_end <- 365 * (100 + 6)
t_ex_yearly_end <- 365 * (100 + 60)

y_panel_spacing <- 0.3

plot_data_ex_inf <- y_inf %>%
  reshape2::melt(c("i", "t"), value.name = "prevalence") %>%
  mutate(eta = x_eta[i], r = x_r[i], label = x_labels[i]) %>%
  filter(t >= t_ex_start, t < t_ex_end) %>%
  group_by(label, eta, r, t) %>%
  summarise(prevalence = sum(prevalence))


plot_data_ex_inf_year <- y_inf %>%
  reshape2::melt(c("i", "t"), value.name = "prevalence") %>%
  mutate(eta = x_eta[i], r = x_r[i], label = x_labels[i]) %>%
  filter(t >= t_ex_start, t < t_ex_yearly_end) %>%
  group_by(label, eta, r, t) %>%
  summarise(prevalence = sum(prevalence)) %>%
  
  mutate(year = floor(t / 365), t_year = t %% 365)


x_eta_plot <- x_eta[[5]]

ggplot() +
  geom_line(aes(x = t, y = prevalence), plot_data_ex_inf %>% filter(eta == x_eta_plot)) +
  
  facet_wrap(~label, ncol = 1) +
  
  xlab("Time (years following burn-in)") +
  ylab("Infection<br>prevalence") +
  
  scale_x_continuous(breaks = scales::breaks_width(365), labels = function(x) {
    (x - min(x, na.rm = TRUE)) / 365
  }) +
  
  # scale_y_continuous(breaks = c(0, 0.1, 0.2)) +
  
  # coord_cartesian(ylim = c(0, 0.11)) +
  
  plot_theme_paper +
  theme(
    panel.grid.major.x = element_line(
      colour = "grey50",
      linetype = "28",
      linewidth = 0.5
    ),
    strip.text.x = element_markdown(),
    panel.spacing.y = unit(y_panel_spacing, "cm")
  )

ggplot() +
  geom_line(
    aes(x = t_year, y = prevalence, group = year),
    linewidth = 0.3,
    alpha = 0.5,
    plot_data_ex_inf_year %>% filter(eta == x_eta_plot)
  ) +
  
  facet_wrap(~label, ncol = 1) +
  
  xlab("Time of year") +
  ylab("Infection prevalence") +
  
  scale_x_continuous(
    breaks = c(0, 90, 180, 270, 365),
    # labels = c("0 (Jan)", "", "180 (Jul)", "", "365 (Jan)"),
    labels = c("Jan", "", "Jul", "", "Jan")
  ) +
  
  scale_y_continuous(breaks = c(0, 0.1, 0.2)) +
  
  coord_cartesian(ylim = c(0, 0.115)) +
  
  plot_theme_paper +
  theme(
    panel.grid.major.x = element_line(
      colour = "grey50",
      linetype = "28",
      linewidth = 0.5
    ),
    strip.text.x = element_blank(),
    panel.spacing.y = unit(y_panel_spacing, "cm")
  )

x_eta_plot <- x_eta[[1]]

plot_data <- y_inf %>%
  reshape2::melt(c("i", "t"), value.name = "prevalence") %>%
  as_tibble() %>% 
  mutate(eta = x_eta[i], r = x_r[i]) %>% 
  
  filter(eta == x_eta_plot) %>%
  
  filter(t >= t_ex_start) %>%
  rename(frac_I = prevalence) %>%
  left_join(
    y_sus %>%
      reshape2::melt(c("i", "t", "ix"), value.name = "prevalence") %>%
      as_tibble() %>% 
      filter(ix == 1) %>% 
      mutate(eta = x_eta[i], r = x_r[i]) %>%
      group_by(eta, r, t) %>%
      summarise(prevalence = sum(prevalence)) %>% 
      
      filter(eta == x_eta_plot) %>%
      
      filter(t >= t_ex_start) %>%
      
      rename(frac_S = prevalence)
  )
  
plot_data %>%   
  mutate(t_year = t %% 365) %>% 
  filter(t_year < 1) %>% 
  
  ggplot() +
  geom_point(aes(x = frac_I, y = frac_S)) +
  
  scale_y_log10() +
  
  plot_theme_paper


plot_data_param <- plot_data %>%
  mutate(t_year = t %% 365,
         x = cos(t_year / 365 * 2 * pi) * (1 + 20 * frac_I), y = sin(t_year / 365 * 2 * pi) * (1 + 20 * frac_I),
         z = frac_S)


plotly::plot_ly(
  plot_data_param,
  x = ~x, y = ~y, z = ~z, 
  type = "scatter3d", mode = "lines",
  line = list(sizemode = "diameter", size = 0.3, color = "black"),
  opacity = 0.1
)
  

