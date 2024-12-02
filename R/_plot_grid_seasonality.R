library(tidyverse)
library(rhdf5)
library(patchwork)


source("R/plot_theme.R")


x_vals <- h5read("data/paper/period_over_grid.jld2", "x_vals")
y_inf_summary <- h5read("data/paper/period_over_grid.jld2", "y_inf_summary")
y_period <- h5read("data/paper/period_over_grid.jld2", "y_period")

plot_data <- tibble(
  eta = x_vals[1, ], r = x_vals[2, ],
  inf_min = y_inf_summary[, 1], inf_max = y_inf_summary[, 2],  inf_mean = y_inf_summary[, 3],
  inc_min = y_inf_summary[, 4], inc_max = y_inf_summary[, 5],  inc_mean = y_inf_summary[, 6],
  period = y_period[,1], period_sd = y_period[,2], period_n = y_period[,3]
) %>%
  mutate(inf_diff = inf_max - inf_min,
         periodic = (period_sd < 1) & (period_n >= 5),
         quasiperiodic = (period_sd >= 1) & (period_n >= 5))




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


plot_data_example_points <- tribble(
  ~eta, ~r, ~label,
  0.3, 0.06, "i",
  0.3, 0.05, "ii",
  0.02, 0.06, "iii",
  0.3, 0.04, "iv",
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
  geom_label(aes(x = eta + 0.02, y = r - 0.0015, label = label), plot_data_example_points,
             label.r = unit(0.1, "cm"), label.size = 0, fill = shades::opacity("white", 0.8))
)


period_cols <- viridis::inferno(n = 8, direction = -1, begin = 0.1)
p_period <- ggplot() +
  geom_tile(aes(x = eta, y = r, fill = period),
            plot_data_periodic) +
  
  geom_tile(aes(x = eta, y = r, fill = factor(4.5)),
            plot_data_quasiperiodic) +
  geom_tile(aes(x = eta, y = r, fill = factor(9)),
            plot_data %>% filter(!periodic, !quasiperiodic)) +
  
  plot_annotations +
  
  
  scale_fill_manual(name = "Period",
                    values = c(period_cols[1:4], "#2260BE", period_cols[5:8], "lightblue") %>% `names<-`(c(1:4, "4.5", 5:9)),
                    
                    labels = c("1 yr", str_c(2:4, "yrs"), "Quasiperiodic", str_c(5:7, "yrs"), "≥8 yrs", "?") %>% `names<-`(c(1:4, "4.5", 5:9)),
                    breaks = c(1:4, 4.5, 5:9)
  ) +
  
  
  coord_fixed(ratio = 5) +
  xlab("Seasonality constant <i>η</i>") + ylab("Mean antibody decay rate <i>r</i>") +
  guides(fill = guide_legend(nrow = 2, ncol = 5)) +
  
  plot_theme_paper +
  theme(legend.position = "bottom", legend.byrow = TRUE)

p_period

x_eta <- h5read("data/paper/period_over_grid_examples.jld2", "x_eta")
x_r <- h5read("data/paper/period_over_grid_examples.jld2", "x_r")
y_sol <- h5read("data/paper/period_over_grid_examples.jld2", "sol_t")

x_labels <- str_c(
  "<b>", c("i", "ii", "iiii", "iv"), ".</b>",
  " <i>η</i> = ", scales::label_comma(accuracy = 0.01)(x_eta),
  ", <i>r</i> = ", scales::label_comma(accuracy = 0.01)(x_r), " — ",
  c("periodic (1 year)", "periodic (2 years)", "quasiperiodic", "chaotic [?]")
)

c_levels <- 10 ^ seq(0, 8, by = 8 / 32)

plot_data_ex <- y_sol %>%
  reshape2::melt(c("i", "t", "class", "ix"), value.name = "prevalence") %>%
  mutate(eta = x_eta[i], r = x_r[i], label = x_labels[i]) %>%
  as_tibble()


plot_data_ex_inf <- plot_data_ex %>%
  filter(class == 2, t >= 365 * 166, t < 365 * (166 + 8)) %>% 
  group_by(label, eta, r, t) %>%
  summarise(prevalence = sum(prevalence))

p_ex_inf <- ggplot() +
  geom_line(aes(x = t, y = prevalence),
            plot_data_ex_inf) +
  
  facet_wrap(~label, ncol = 1, scales = "free_x") +
  
  xlab("Time (years)") + ylab("Infection prevalence") +
  
  scale_x_continuous(breaks = scales::breaks_width(365),
                     labels = function(x) { (x - min(x, na.rm = TRUE)) / 365 }) +
  
  scale_y_continuous(breaks = c(0, 0.1, 0.2)) +
  
  coord_cartesian(ylim = c(0, 0.2)) +
  
  plot_theme_paper +
  theme(panel.grid.major.x = element_line(colour = "grey50", linetype = "28", linewidth = 0.5),
        strip.text = element_markdown(size = 12))



plot_data_ex_mean <- plot_data_ex %>%
  mutate(c = c_levels[ix]) %>%
  filter(class == 1, t >= 365 * 166, t < 365 * (166 + 8)) %>% 
  group_by(label, eta, r, t) %>% 
  summarise(mean = sum(prevalence * c))

plot_data_ex_mean_year <- plot_data_ex %>%
  mutate(c = c_levels[ix]) %>%
  filter(class == 1, t >= 365 * 200, t < 365 * 220) %>% 
  group_by(label, eta, r, t) %>% 
  summarise(mean = sum(prevalence * c)) %>%
  
  mutate(year = floor(t / 365),
         t_year = t %% 365)


p_ex_yearly_antibody <- ggplot() +
  geom_line(aes(x = t_year, y = mean, group = year),
            linewidth = 0.3, alpha = 0.3,
            plot_data_ex_mean_year) +
  
  facet_wrap(~label, ncol = 1, scales = "free_x") +
  
  xlab("Time of year (days)") + ylab("Mean antibody concentration") +
  
  scale_x_continuous(breaks = c(0, 90, 180, 270, 365)) +
  
  scale_y_continuous(trans = "log10", labels = scales::label_log(base = 10),
                     breaks = 10^c(1, 3, 5)) +
  
  plot_theme_paper +
  theme(panel.grid.major.x = element_line(colour = "grey50", linetype = "28", linewidth = 0.5),
        strip.text = element_markdown(colour = "white"))


# p_ex_yearly_antibody

(p_period | p_ex_inf | p_ex_yearly_antibody) +
  plot_layout(widths = c(3, 2, 1)) +
  plot_annotation(tag_levels = list(c("A", "B", ""))) &
  theme(plot.tag = element_text(face = "bold", size = 15))


ggsave(
  "results/results_grid_seasonality.png",
  device = png,
  width = 14, height = 6.75,
  bg = "white"
)



