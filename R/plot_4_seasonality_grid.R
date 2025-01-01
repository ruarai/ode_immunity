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
plot_data_chaotic <- plot_data %>% filter(chaotic)

plot_data_eta_zero <- plot_data %>% filter(eta == 0)
plot_data_eta_zero_periodic <- plot_data %>% filter(eta == 0, r < 0.065)

bifur_zero <- plot_data_eta_zero %>% filter(inf_diff < 1e-3) %>% pull(r) %>% head(1)

year_stops <- c(1/2, 2/3, 1, 3/2, 2, 3)
year_marks <- approxfun(plot_data_eta_zero_periodic$period, plot_data_eta_zero_periodic$r)(365 * year_stops)
plot_data_year_marks <- tibble(r_0 = year_marks, year = year_stops) %>%
  mutate(year_label = str_c(scales::label_comma()(year), " yr"))


plot_data_example_points <- tribble(
  ~eta, ~r, ~label, ~colour,
  0.0, 0.06, "i", factor(0),
  0.07, 0.06, "ii", factor(4.5),
  0.2, 0.06, "iii", factor(1),
  0.27, 0.06, "iv", factor(2),
  0.37, 0.06, "v", factor(9)
)
plot_annotations <- list(
  geom_segment(
    # aes(x = r_0, y = 0.0, xend = r_0 + 0.001, yend = -0.01),
    aes(x = -0.003, y = r_0, xend = -0.01, yend = r_0),
    plot_data_year_marks
  ),
  
  annotate("linerange", x = -0.0065, ymin = bifur_zero, ymax = 0.1),
  annotate("segment", x = -0.003, y = bifur_zero, xend = -0.01, yend = bifur_zero),
  geom_text(aes(x = -0.07, y = r_0 + 0.0002, label = year_label), hjust = 0, plot_data_year_marks),
  annotate("text", x = -0.07, y = 0.085, label = "Fixed\npoint", hjust = 0),
  geom_point(aes(x = eta, y = r), plot_data_example_points, colour = "black", size = 1.4, stroke = 1),
  geom_point(aes(x = eta, y = r), plot_data_example_points, colour = "white", size = 0.7, stroke = 0.5),
  geom_label(aes(x = eta + 0.01, y = r - 0.005, label = label), plot_data_example_points,
             label.r = unit(0.1, "cm"), label.size = 0, fill = shades::opacity("white", 0.8))
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

p_period


p_min <- ggplot()  +
  geom_tile(aes(x = eta, y = r, fill = pmax(-14, log10(inf_min))),
            plot_data) +
  
  plot_annotations +
  scale_colour_manual(
    values = rep("white", 11),
    breaks = c(0:4, 4.5, 5:9)
  ) + 
  
  scale_fill_stepsn(
    colours = rev(colorspace::diverging_hcl(n = 20, h = c(240, 15), c = c(60, 80), l = c(75, 5), power = c(1.2, 1.5))),
    name = "Minimum\ninfection\nprevalence (log10)",
    limits = c(-15, -1),
    breaks = seq(-15, -1, 0.5),
    labels = c("<-15", "", "", "", "-13", "", "", "", "-11", "", "", "", "-9", 
               "", "", "", "-7", "", "", "", "-5", "", "", "", "-3", "", "", 
               "", "-1")
  ) +
  
  coord_fixed(ratio = 5, ylim = c(0, 0.1)) +
  xlab("Seasonality strength <i>η</i>") + ylab("Antibody decay rate <i>r</i>")  +
  
  plot_theme_paper +
  guides(fill = guide_colourbar(barwidth = 15),
         colour = guide_none()) +
  theme(legend.position = "bottom",
        legend.title = element_text(margin = margin(r = 0.7, unit = "cm")))

p_period | p_min


p_period + ggtitle(NULL, "<b>A</b> — Qualitative dynamics") | 
  p_min + ggtitle(NULL, "<b>B</b> — Minimum infection prevalence")

ggsave(
  "results/results_grid_seasonality.png",
  device = png,
  width = 13, height = 7,
  bg = "white"
)


# (p_period + ggtitle(NULL, "<b>A</b> — Dynamics") | 
#     p_min + ggtitle(NULL, "<b>B</b> — Minimum infection prevalence")) / (
#       ((p_ex_inf | p_ex_yearly_inf | p_ex_antibody | p_ex_yearly_antibody) +
#          plot_layout(widths = c(2, 1, 2, 1)))
#     )
# 
# 
# 
# ggsave(
#   "results/results_grid_seasonality_2.png",
#   device = png,
#   width = 13, height = 13,
#   bg = "white"
# )
