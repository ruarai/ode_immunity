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

plot_data_chaotic <- plot_data %>% filter(chaotic)

plot_data_empty <- plot_data %>% filter(!periodic, !quasiperiodic, !chaotic)

plot_data_eta_zero <- plot_data %>% filter(eta == 0)
plot_data_eta_zero_periodic <- plot_data %>% filter(eta == 0, r < 0.065)

bifur_zero <- plot_data_eta_zero %>% filter(inf_diff < 1e-3) %>% pull(r) %>% head(1)

plot_data_quasiperiodic <- plot_data %>% filter(quasiperiodic, r < bifur_zero)

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
  ),
  geom_point(
    aes(x = eta, y = r),
    plot_data_example_points,
    colour = "black",
    size = 1.4,
    stroke = 1
  ),
  geom_point(
    aes(x = eta, y = r),
    plot_data_example_points,
    colour = "white",
    size = 0.7,
    stroke = 0.5
  ),
  geom_label(
    aes(x = eta + 0.01, y = r + 0.0015, label = label),
    plot_data_example_points,
    label.r = unit(0.1, "cm"),
    label.size = 0,
    fill = shades::opacity("white", 0.8)
  )
)

period_cols <- viridis::inferno(n = 8, direction = -1, begin = 0.1)

p_period <- ggplot() +
  
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
  
  
  coord_fixed(ratio = 16.66, ylim = c(0, 0.03)) +
  xlab("Seasonality strength <i>η</i>") + ylab("Antibody decay rate <i>r</i>") +
  guides(fill = guide_legend(nrow = 3, ncol = 5),
         colour = guide_none()) +
  
  plot_theme_paper +
  theme(legend.position = "bottom", legend.byrow = TRUE,
        legend.key = element_rect(colour = "grey80", linewidth = 0.5))


p_period


color_tbl <- read_delim("~/Downloads/ScientificColourMaps8/lipari/lipari.txt", col_names = c("R", "G", "B")) %>%
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
  xlab("Seasonality strength <i>η</i>") + ylab("Antibody decay rate <i>r</i>")  +
  
  plot_theme_paper +
  guides(fill = guide_colourbar(barwidth = 15),
         colour = guide_none()) +
  theme(legend.position = "bottom",
        legend.title = element_text(margin = margin(r = 0.7, unit = "cm")))

p_min

p_period | p_min



x_eta <- h5read("data/period_over_grid_examples.jld2", "x_eta")
x_r <- h5read("data/period_over_grid_examples.jld2", "x_r")
y_inf <- h5read("data/period_over_grid_examples.jld2", "y_inf")
y_sus <- h5read("data/period_over_grid_examples.jld2", "y_sus")

x_labels <- str_c(
  "<b>",
  c("i", "ii", "iii", "iv", "v"),
  ".</b> ",
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

plot_data_labels <- plot_data_ex_inf %>%
  distinct(label)


p_ex_inf <- ggplot() +
  geom_line(aes(x = t, y = prevalence), plot_data_ex_inf) +
  
  geom_richtext(
    aes(x = t_ex_start, y = 0.1, label = label),
    hjust = 0,
    label.r = unit(0, "cm"),
    label.size = 0,
    size = 4.5,
    fill = "grey97",
    plot_data_labels
  ) +
  
  facet_wrap(~label, ncol = 1) +
  
  xlab("Time (years following burn-in)") +
  ylab("Infection prevalence") +
  
  scale_x_continuous(breaks = scales::breaks_width(365), labels = function(x) {
    (x - min(x, na.rm = TRUE)) / 365
  }) +
  
  scale_y_continuous(breaks = c(0, 0.1, 0.2)) +
  
  coord_cartesian(ylim = c(0, 0.11)) +
  
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

p_ex_inf


plot_data_ex_inf_year <- y_inf %>%
  reshape2::melt(c("i", "t"), value.name = "prevalence") %>%
  mutate(eta = x_eta[i], r = x_r[i], label = x_labels[i]) %>%
  filter(t >= t_ex_start, t < t_ex_yearly_end) %>%
  group_by(label, eta, r, t) %>%
  summarise(prevalence = sum(prevalence)) %>%
  
  mutate(year = floor(t / 365), t_year = t %% 365)

p_ex_yearly_inf <- ggplot() +
  geom_line(
    aes(x = t_year, y = prevalence, group = year),
    linewidth = 0.3,
    alpha = 0.5,
    plot_data_ex_inf_year
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

p_ex_yearly_inf

(p_ex_inf | p_ex_yearly_inf) +
  plot_layout(widths = c(2, 1))


plot_data_ex_mean <- y_sus %>%
  reshape2::melt(c("i", "t", "ix"), value.name = "prevalence") %>%
  mutate(eta = x_eta[i], r = x_r[i], label = x_labels[i], c = c_levels[ix]) %>%
  filter(t >= t_ex_start, t < t_ex_end) %>%
  group_by(label, eta, r, t) %>%
  summarise(mean = sum(prevalence * c))

p_ex_antibody <- ggplot() +
  geom_line(aes(x = t, y = mean), plot_data_ex_mean) +
  
  facet_wrap(~label, ncol = 1) +
  
  xlab("Time (years)") +
  ylab("Mean antibody concentration") +
  
  scale_x_continuous(breaks = scales::breaks_width(365), labels = function(x) {
    (x - min(x, na.rm = TRUE)) / 365
  }) +
  
  scale_y_continuous(
    trans = "log2",
    labels = scales::label_log(base = 2),
    breaks = 2^c(1, 3, 5)
  ) +
  
  coord_cartesian(ylim = 10^c(2.5, 6.5)) +
  
  plot_theme_paper +
  theme(
    panel.grid.major.x = element_line(
      colour = "grey50",
      linetype = "28",
      linewidth = 0.5
    ),
    strip.text = element_markdown(size = 12)
  )


plot_data_ex_mean_year <- y_sus %>%
  reshape2::melt(c("i", "t", "ix"), value.name = "prevalence") %>%
  mutate(eta = x_eta[i], r = x_r[i], label = x_labels[i], c = c_levels[ix]) %>%
  filter(t >= t_ex_start, t < t_ex_yearly_end) %>%
  group_by(label, eta, r, t) %>%
  summarise(mean = sum(prevalence * c)) %>%
  
  mutate(year = floor(t / 365), t_year = t %% 365)


p_ex_yearly_antibody <- ggplot() +
  geom_line(
    aes(x = t_year, y = mean, group = year),
    linewidth = 0.3,
    alpha = 0.5,
    plot_data_ex_mean_year
  ) +
  
  facet_wrap(~label, ncol = 1) +
  
  xlab("Time of year") +
  ylab("Mean antibody concentration") +
  
  scale_x_continuous(
    breaks = c(0, 90, 180, 270, 365),
    labels = c("Jan", "", "Jul", "", "Jan")
  ) +
  
  scale_y_continuous(
    trans = "log2",
    labels = scales::label_log(base = 2),
    breaks = 2^c(1, 3, 5)
  ) +
  
  coord_cartesian(ylim = 2^c(0, 6.5)) +
  
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


p_examples <- ((p_ex_inf +
                  ggtitle(NULL, "<b>C</b> — Infection prevalence") |
                  p_ex_yearly_inf +
                  ggtitle(NULL, "<b>D</b> — Infection prevalence<br>by time of year") |
                  p_ex_yearly_antibody +
                  ggtitle(NULL, "<b>E</b> — Mean antibody<br>concentration")) +
                 plot_layout(widths = c(2, 1, 1)))

p_examples




p <- (
  p_period + ggtitle(NULL, "<b>A</b> — Dynamics") |
  p_min + ggtitle(NULL, "<b>B</b> — Minimum infection prevalence")
) / p_examples




ggsave(
  "results/results_grid_seasonality_2.png",
  p,
  device = png,
  width = 13, height = 13,
  bg = "white"
)

# p_legend <- cowplot::plot_grid(x[[17]])
# 
# 
# ggsave(
#   "results/results_grid_seasonality_legend.pdf",
#   p_legend,
#   device = cairo_pdf,
#   width = 13/2, height = 2,
#   bg = "white"
# )


