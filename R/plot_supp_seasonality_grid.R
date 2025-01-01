
library(tidyverse)
library(rhdf5)
library(patchwork)

source("R/plot_theme.R")

plot_data <- read_seasonality_data("data/paper/period_over_grid.jld2")

plot_data_periodic <- plot_data %>%
  filter(eta > 0) %>% 
  filter(periodic) %>% 
  mutate(period = period / 365,
         period = pmin(period, 8),
         period = factor(round(period)))

plot_data_quasiperiodic <- plot_data %>% filter(quasiperiodic)

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
  
  annotate("linerange", x = -0.0065, ymin = bifur_zero, ymax = 0.1),
  annotate("segment", x = -0.003, y = bifur_zero, xend = -0.01, yend = bifur_zero),
  geom_text(aes(x = -0.07, y = r_0 + 0.0002, label = year_label), hjust = 0, plot_data_year_marks),
  annotate("text", x = -0.07, y = 0.085, label = "Fixed\npoint", hjust = 0)
)


period_cols <- viridis::inferno(n = 8, direction = -1, begin = 0.1)
p_period <- ggplot() +
  
  geom_tile(aes(x = eta, y = r, fill = period),
            plot_data_periodic) +
  
  geom_tile(aes(x = eta, y = r, fill = factor(4.5)),
            plot_data_quasiperiodic) +
  geom_tile(aes(x = eta, y = r, fill = factor(9)),
            plot_data_chaotic) +
  
  plot_annotations +
  
  scale_fill_manual(
    name = "Period",
    values = c(period_cols[1:4], "#2260BE", period_cols[5:8], "#BDE6F4") %>% `names<-`(c(1:4, "4.5", 5:9)),
    
    labels = c("1 yr", str_c(2:4, "yrs"), "Quasiperiodic", str_c(5:7, "yrs"), "≥8 yrs", "Chaotic") %>% `names<-`(c(1:4, "4.5", 5:9)),
    breaks = c(1:4, 4.5, 5:9)
  ) +
  
  
  coord_fixed(ratio = 5, ylim = c(0, 0.1)) +
  xlab("Seasonality strength <i>η</i>") + ylab("Antibody decay rate <i>r</i>") +
  guides(fill = guide_legend(nrow = 2, ncol = 5),
         colour = guide_none()) +
  
  plot_theme_paper +
  theme(legend.position = "bottom", legend.byrow = TRUE)


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
  
  coord_fixed(ratio = 5) +
  xlab("Seasonality constant <i>η</i>") + ylab("Antibody decay rate <i>r</i>")  +
  
  plot_theme_paper +
  guides(fill = guide_colourbar(barwidth = 15)) +
  theme(legend.position = "bottom")

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
  coord_fixed(ratio = 5) +
  xlab("Seasonality constant <i>η</i>") + ylab("Antibody decay rate <i>r</i>")  +
  
  plot_theme_paper +
  guides(fill = guide_colourbar(barwidth = 15)) +
  theme(legend.position = "bottom")

# p_max


p_min <- ggplot()  +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.005, fill = "white")+
  geom_tile(aes(x = eta, y = r, fill = pmax(-14, log10(inf_min))),
            plot_data) +
  
  plot_annotations + 
  
  scale_fill_stepsn(
    colours = rev(colorspace::diverging_hcl(n = 20, h = c(240, 15), c = c(60, 80), l = c(75, 5), power = c(1.2, 1.5))),
    name = "Minimum\ninfection\nprevalence (log10)",
    limits = c(-15, -1),
    breaks = seq(-15, -1, 0.5),
    labels = c("<-15", "", "", "", "-13", "", "", "", "-11", "", "", "", "-9", 
               "", "", "", "-7", "", "", "", "-5", "", "", "", "-3", "", "", 
               "", "-1")
  ) +
  
  coord_fixed(ratio = 5) +
  xlab("Seasonality constant <i>η</i>") + ylab("Antibody decay rate <i>r</i>")  +
  
  plot_theme_paper +
  guides(fill = guide_colourbar(barwidth = 15)) +
  theme(legend.position = "bottom")

# p_min



p_full <- (p_period | p_min) / (p_attack_rate | p_max) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))




ggsave(
  "results/results_supp_seasonality.png",
  p_full,
  device = png,
  width = 14, height = 14,
  bg = "white"
)

