library(tidyverse)
library(rhdf5)
library(patchwork)

source("R/plot_theme.R")
source("R/read_seasonality_data.R")

plot_data <- read_seasonality_data("data/period_over_grid.jld2")


plot_data_periodic <- plot_data %>%
  filter(eta > 0) %>% 
  filter(periodic) %>% 
  mutate(period = period / 365,
         period = pmin(period, 8),
         period = factor(round(period)))

plot_data_eta_zero <- plot_data %>% filter(eta == 0)
plot_data_eta_zero_periodic <- plot_data %>% filter(eta == 0, r < 0.065)

bifur_zero <- plot_data_eta_zero %>% filter(inf_diff < 1e-3) %>% pull(r) %>% head(1)

year_stops <- c(1/2, 2/3, 1, 3/2, 2, 3)
year_marks <- approxfun(plot_data_eta_zero_periodic$period, plot_data_eta_zero_periodic$r)(365 * year_stops)
plot_data_year_marks <- tibble(r_0 = year_marks, year = year_stops) %>%
  mutate(year_label = str_c(scales::label_comma()(year), " yr"))



plot_annotations <- list(
  geom_segment(
    # aes(x = r_0, y = 0.0, xend = r_0 + 0.001, yend = -0.01),
    aes(x = -0.003, y = r_0, xend = -0.01, yend = r_0),
    plot_data_year_marks
  ),
  
  annotate("linerange", x = -0.0065, ymin = bifur_zero, ymax = 0.03),
  annotate("segment", x = -0.003, y = bifur_zero, xend = -0.01, yend = bifur_zero),
  geom_text(aes(x = -0.07, y = r_0 + 0.0002, label = year_label), hjust = 0, plot_data_year_marks),
  annotate("text", x = -0.07, y = 0.0275, label = "Fixed\npoint", hjust = 0)
)

plot_data_peaks <- plot_data %>% 
  filter(periodic, round(period / 365) > 0) %>%
  mutate(label = case_when(
    round(period / 365) > 2 ~ "Period > 2",
    TRUE ~ str_c(round(peak_density * (period / 365)), "/", round(period/365))
  ))



p_orbits <- ggplot() +
  
  geom_tile(aes(x = eta, y = r, fill = label),
            plot_data_peaks) +
  
  plot_annotations +
  
  ggokabeito::scale_fill_okabe_ito(name = "[Peaks]/[Period]") +
  
  coord_fixed(ratio = 16.66, ylim = c(0, 0.03), xlim = c(-0.07, 0.5)) +
  xlab("Seasonal forcing strength <i>η</i>") + ylab("Effective antibody decay rate <i>r</i>") +
  
  plot_theme_paper +
  
  guides(fill = guide_legend(direction = "horizontal", reverse = FALSE, byrow = FALSE, ncol = 8)) +
  theme(legend.position = "bottom") +
  
  ggtitle(NULL, "<b>A</b> — Number of peaks per period (period ≤ 2yr)")

p_orbits

p_peak_density <- ggplot() +
  
  geom_tile(aes(x = eta, y = r, fill = peak_density),
            plot_data) +
  
  plot_annotations +
  
  scale_fill_viridis_c(option = 5, name = "Peaks per year") +
  
  coord_fixed(ratio = 16.66, ylim = c(0, 0.03), xlim = c(-0.07, 0.5)) +
  xlab("Seasonal forcing strength <i>η</i>") + ylab("Effective antibody decay rate <i>r</i>") +
  
  plot_theme_paper +
  
  guides(fill = guide_colourbar(barwidth = 15)) +
  theme(legend.position = "bottom") +
  
  ggtitle(NULL, "<b>B</b> — Average number of peaks per year")


p_max <- ggplot()  +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.005, fill = "white")+
  geom_tile(aes(x = eta, y = r, fill = inf_max),
            plot_data) +
  
  plot_annotations + 
  
  scale_fill_viridis_b(
    name = "Peak\ninfection\nprevalence", option = "mako",
    limits = c(0, 0.25),
    breaks = seq(0, 0.25, 0.025),
    labels = c("0", "", "0.05", "", "0.1", "", "0.15", "", "0.2", "", "0.25")
  ) +
  
  coord_fixed(ratio = 16.66, ylim = c(0, 0.03)) +
  xlab("Seasonal forcing strength <i>η</i>") + ylab("Effective antibody decay rate <i>r</i>")  +
  
  plot_theme_paper +
  guides(fill = guide_colourbar(barwidth = 15)) +
  theme(legend.position = "bottom") +
  
  ggtitle(NULL, "<b>C</b> — Peak infection prevalence")


p_attack_rate <- ggplot()  +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.005, fill = "white")+
  geom_tile(aes(x = eta, y = r, fill = inc_mean * 365),
            plot_data) +
  
  plot_annotations + 
  
  scale_fill_stepsn(
    colours = colorspace::sequential_hcl(n = 20, h = c(300, 75), c = c(40, NA, 95), l = c(15, 90), power = c(1, 1.1)),
    name = "Annual infection incidence",
    limits = c(0, 2.0),
    breaks = seq(0, 2.0, 0.2),
    labels = c("0", "", "0.4", "", "0.8", "", "1.2", "", "1.6", "", "2.0")
  ) +
  
  coord_fixed(ratio = 16.66, ylim = c(0, 0.03)) +
  xlab("Seasonal forcing strength <i>η</i>") + ylab("Effective antibody decay rate <i>r</i>")  +
  
  plot_theme_paper +
  guides(fill = guide_colourbar(barwidth = 15)) +
  theme(legend.position = "bottom") +
  
  ggtitle(NULL, "<b>D</b> — Annual average infection incidence")


p_full <- (p_orbits | p_peak_density) / (p_max | p_attack_rate)


ggsave(
  "results/results_supp_seasonality.png",
  p_full,
  device = png,
  width = 14, height = 14,
  bg = "white"
)




