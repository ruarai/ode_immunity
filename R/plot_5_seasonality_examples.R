library(tidyverse)
library(rhdf5)
library(patchwork)


source("R/plot_theme.R")


x_eta <- h5read("data/paper/period_over_grid_examples.jld2", "x_eta")
x_r <- h5read("data/paper/period_over_grid_examples.jld2", "x_r")
y_inf <- h5read("data/paper/period_over_grid_examples.jld2", "y_inf")
y_sus <- h5read("data/paper/period_over_grid_examples.jld2", "y_sus")

x_labels <- str_c(
  "<b>", c("i", "ii", "iii", "iv", "v"), ".</b> ",
  c("Zero seasonality", "Quasiperiodic", "Periodic (1 year)", "Periodic (2 years)", "Chaotic"),
  " (<i>η</i> = ", scales::label_comma(accuracy = 0.01)(x_eta),
  ", <i>r</i> = ", scales::label_comma(accuracy = 0.01)(x_r), ")"
)

c_levels <- 10 ^ seq(0, 8, by = 8 / 32)

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
  geom_line(aes(x = t, y = prevalence),
            plot_data_ex_inf) +
  
  geom_richtext(aes(x = t_ex_start, y = 0.1, label = label),
                hjust = 0, label.r = unit(0, "cm"), label.size = 0,
                fill = "grey97", 
                plot_data_labels) +
  
  facet_wrap(~label, ncol = 1) +
  
  xlab("Time (years)") + ylab("Infection prevalence") +
  
  scale_x_continuous(breaks = scales::breaks_width(365),
                     labels = function(x) { (x - min(x, na.rm = TRUE)) / 365 }) +
  
  scale_y_continuous(breaks = c(0, 0.1, 0.2)) +
  
  coord_cartesian(ylim = c(0, 0.11)) +
  
  plot_theme_paper +
  theme(panel.grid.major.x = element_line(colour = "grey50", linetype = "28", linewidth = 0.5),
        strip.text.x = element_blank(),
        panel.spacing.y = unit(y_panel_spacing, "cm"))

p_ex_inf


plot_data_ex_inf_year <- y_inf %>%
  reshape2::melt(c("i", "t"), value.name = "prevalence") %>%
  mutate(eta = x_eta[i], r = x_r[i], label = x_labels[i]) %>%
  filter(t >= t_ex_start, t < t_ex_yearly_end) %>% 
  group_by(label, eta, r, t) %>%
  summarise(prevalence = sum(prevalence)) %>%
  
  mutate(year = floor(t / 365),
         t_year = t %% 365)

p_ex_yearly_inf <- ggplot() +
  geom_line(aes(x = t_year, y = prevalence, group = year),
            linewidth = 0.3, alpha = 0.5,
            plot_data_ex_inf_year) +
  
  facet_wrap(~label, ncol = 1) +
  
  xlab("Time of year") + ylab("Infection prevalence") +
  
  scale_x_continuous(breaks = c(0, 90, 180, 270, 365),
                     # labels = c("0 (Jan)", "", "180 (Jul)", "", "365 (Jan)"),
                     labels = c("Jan", "", "Jul", "", "Jan")) +
  
  scale_y_continuous(breaks = c(0, 0.1, 0.2)) +
  
  coord_cartesian(ylim = c(0, 0.115)) +
  
  plot_theme_paper +
  theme(panel.grid.major.x = element_line(colour = "grey50", linetype = "28", linewidth = 0.5),
        strip.text.x = element_blank(),
        panel.spacing.y = unit(y_panel_spacing, "cm"))

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
  geom_line(aes(x = t, y = mean),
            plot_data_ex_mean) +
  
  facet_wrap(~label, ncol = 1) +
  
  xlab("Time (years)") + ylab("Mean antibody concentration") +
  
  scale_x_continuous(breaks = scales::breaks_width(365),
                     labels = function(x) { (x - min(x, na.rm = TRUE)) / 365 }) +
  
  scale_y_continuous(trans = "log10", labels = scales::label_log(base = 10),
                     breaks = 10^c(1, 3, 5)) +
  
  coord_cartesian(ylim = 10^c(2.5, 6.5)) +
  
  plot_theme_paper +
  theme(panel.grid.major.x = element_line(colour = "grey50", linetype = "28", linewidth = 0.5),
        strip.text = element_markdown(size = 12))


plot_data_ex_mean_year <- y_sus %>%
  reshape2::melt(c("i", "t", "ix"), value.name = "prevalence") %>%
  mutate(eta = x_eta[i], r = x_r[i], label = x_labels[i], c = c_levels[ix]) %>%
  filter(t >= t_ex_start, t < t_ex_yearly_end) %>% 
  group_by(label, eta, r, t) %>% 
  summarise(mean = sum(prevalence * c)) %>%
  
  mutate(year = floor(t / 365),
         t_year = t %% 365)


p_ex_yearly_antibody <- ggplot() +
  geom_line(aes(x = t_year, y = mean, group = year),
            linewidth = 0.3, alpha = 0.5,
            plot_data_ex_mean_year) +
  
  facet_wrap(~label, ncol = 1) +
  
  xlab("Time of year") + ylab("Mean antibody concentration") +
  
  scale_x_continuous(breaks = c(0, 90, 180, 270, 365),
                     labels = c("Jan", "", "Jul", "", "Jan")) +
  
  scale_y_continuous(trans = "log10", labels = scales::label_log(base = 10),
                     breaks = 10^c(1, 3, 5)) +
  
  coord_cartesian(ylim = 10^c(2.5, 6.5)) +
  
  plot_theme_paper +
  theme(panel.grid.major.x = element_line(colour = "grey50", linetype = "28", linewidth = 0.5),
        strip.text.x = element_blank(),
        panel.spacing.y = unit(y_panel_spacing, "cm"))


p_examples <- (
  (
    p_ex_inf + ggtitle(NULL, "<b>C</b> — Infection prevalence") | 
      p_ex_yearly_inf + ggtitle(NULL, "<b>D</b> — Infection prevalence<br>by time of year") |
      p_ex_yearly_antibody + ggtitle(NULL, "<b>E</b> — Mean antibody<br>concentration") 
  ) +
    plot_layout(widths = c( 2, 1, 1))
)

p_examples

ggsave(
  "results/results_examples.png",
  p_examples,
  device = png,
  width = 10, height = 6,
  bg = "white"
)


