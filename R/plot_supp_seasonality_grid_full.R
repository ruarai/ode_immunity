
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
         period = pmin(period, 8),
         period = factor(round(period)))

plot_data_chaotic <- plot_data %>% filter(chaotic)

plot_data_empty <- plot_data %>% filter(!periodic, !quasiperiodic, !chaotic)

plot_data_eta_zero <- plot_data %>% filter(eta == 0)
plot_data_eta_zero_periodic <- plot_data %>% filter(eta == 0, r < 0.065)

bifur_zero <- plot_data_eta_zero %>% filter(inf_diff < 1e-3) %>% pull(r) %>% head(1)

plot_data_quasiperiodic <- plot_data %>% filter(quasiperiodic, r < bifur_zero)

plot_data_quasiperiodic <- plot_data %>% filter(!chaotic, !periodic, toroidal_filling > 5)

grid_factor <- 110
r_ratio <- 16.66
plot_data_low <- plot_data %>% filter(inf_min <= inf_threshold) %>% 
  mutate(f = if_else(
    xor(floor(eta * grid_factor) %% 2 == 0, floor(r * grid_factor * r_ratio) %% 2 == 0), factor(8), factor(7)
  ))


year_stops <- c(1/2, 2/3, 1, 3/2, 2, 3)
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

period_cols <- viridis::inferno(n = 8, direction = -1, begin = 0.1)

p_period_full <- ggplot() +
  
  geom_tile(aes(x = eta, y = r, fill = period),
            plot_data_periodic) +
  geom_tile(aes(x = eta, y = r, fill = factor(9)),
            plot_data_quasiperiodic) +
  geom_tile(aes(x = eta, y = r, fill = factor(10)),
            plot_data_chaotic) +
  geom_tile(aes(x = eta, y = r, fill = factor(11)),
            plot_data_empty) +
  
  plot_annotations +
  
  scale_fill_manual(
    name = "Period",
    values = c(period_cols, "#0076BC", "#9ED7F3", "white") %>% 
      `names<-`(1:11),
    
    labels = c("1 year", str_c(2:7, " years"), "≥8 years", "Quasiperiodic", "Chaotic", "Unclassified") %>%
      `names<-`(1:11),
    
    breaks = c(9, 1:4, 10, 5:8, 11)
  ) +
  
  
  coord_fixed(ratio = 16.66, ylim = c(0.000, 0.03)) +
  xlab("Seasonal forcing strength <i>η</i>") + ylab("Effective antibody decay rate <i>r</i>") +
  guides(fill = guide_legend(nrow = 3, ncol = 5),
         colour = guide_none()) +
  
  plot_theme_paper +
  theme(legend.position = "bottom", legend.byrow = TRUE,
        legend.key = element_rect(colour = "grey80", linewidth = 0.5))

p_period_full

color_tbl <- read_delim("data/lipari.txt", col_names = c("R", "G", "B")) %>%
  mutate(c = rgb(R, G, B))

min_cols <- color_tbl$c[seq(1, nrow(color_tbl), length.out = 128)]

p_min <- ggplot()  +
  geom_tile(aes(x = eta, y = r, fill = pmax(-14, log10(inf_min))),
            plot_data) +
  
  plot_annotations +
  scale_colour_manual(
    values = rep("white", 11),
    breaks = c(0:4, 4.5, 5:9)
  ) + 
  
  scale_fill_stepsn(
    # colours = rev(colorspace::diverging_hcl(n = 20, h = c(240, 15), c = c(60, 80), l = c(75, 5), power = c(1.2, 1.5))),
    colours = min_cols,
    name = "Minimum\ninfection\nprevalence (log10)",
    limits = c(-15, -1),
    breaks = seq(-15, -1, 0.5),
    labels = c("<-15", "", "", "", "-13", "", "", "", "-11", "", "", "", "-9", 
               "", "", "", "-7", "", "", "", "-5", "", "", "", "-3", "", "", 
               "", "-1")
  ) +
  
  coord_fixed(ratio = 16.66, ylim = c(0, 0.03)) +
  xlab("Seasonal forcing strength <i>η</i>") + ylab("Effective antibody decay rate <i>r</i>")  +
  
  plot_theme_paper +
  guides(fill = guide_colourbar(barwidth = 15),
         colour = guide_none()) +
  theme(legend.position = "bottom",
        legend.title = element_text(margin = margin(r = 0.7, unit = "cm")))


(
  p_min + ggtitle(NULL, "<b>A</b> — Minimum infection prevalence") |
    p_period_full + ggtitle(NULL, "<b>B</b> — Dynamics (full)")
) 

ggsave(
  "results/results_supp_grid_seasonality_2.png",
  device = png,
  width = 13, height = 7,
  bg = "white"
)
