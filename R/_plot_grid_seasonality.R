
library(tidyverse)

library(rhdf5)


source("../ode_immunity_multi/R/plot_theme.R")


x_vals <- h5read("data/paper/period_over_grid.jld2", "x_vals")
y_inf_summary <- h5read("data/paper/period_over_grid.jld2", "y_inf_summary")
y_period <- h5read("data/paper/period_over_grid.jld2", "y_period")
y_attack_rate <- h5read("data/paper/period_over_grid.jld2", "y_attack_rate")

plot_data <- tibble(
  eta = x_vals[1, ], rho = x_vals[2, ],
  inf_min = y_inf_summary[, 1], inf_max = y_inf_summary[, 2],  inf_mean = y_inf_summary[, 3],
  inc_min = y_inf_summary[, 4], inc_max = y_inf_summary[, 5],  inc_mean = y_inf_summary[, 6],
  period = y_period[,1], period_sd = y_period[,2], period_n = y_period[,3],
  attack_rate = y_attack_rate
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
plot_data_quasiperiodic <- plot_data %>% filter(!quasiperiodic, !periodic)

plot_data_eta_zero <- plot_data %>% filter(eta == 0)

bifur_zero <- plot_data_eta_zero %>% filter(inf_diff < 1e-3) %>% pull(rho) %>% head(1)

year_stops <- c(1/2, 2/3, 1, 3/2, 2, 3, 4)
year_marks <- approxfun(plot_data_eta_zero$period, plot_data_eta_zero$rho)(365 * year_stops)
plot_data_year_marks <- tibble(rho_0 = year_marks, year = year_stops) %>%
  mutate(year_label = str_c(scales::label_comma()(year), " yr"))

plot_data_example_points <- tribble(
  ~eta, ~rho, ~label,
  0.3, 0.0028, "iv",
  0.02, 0.0028, "iii",
  0.3, 0.0023, "i",
  0.3, 0.0033, "ii"
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
  geom_tile(aes(x = eta, y = rho, fill = factor(9)),
            plot_data_quasiperiodic) +
  
  plot_annotations +
  
  geom_label(aes(x = eta + 0.015, y = rho - 0.00015, label = label), plot_data_example_points,
             label.r = unit(0.1, "cm"), label.size = 0, fill = shades::opacity("white", 0.8),
             family = "bold") +
  
  geom_text(aes(x = -0.07, y = rho_0, label = year_label), hjust = 0, plot_data_year_marks) +
  annotate("text", x = -0.07, y = 0.00465, label = "Fixed\npoint", hjust = 0) +
  
  scale_fill_manual(name = "Period",
                    values = c(viridis::inferno(n = 8, direction = -1, begin = 0.1), "#2260BE"),
                    labels = c(1:7, "≥8", "Quasiperiodic")) +

  # scale_fill_viridis_b(option = "inferno", direction = -1, breaks = 1:8, labels = c(1:7, "8")) +

  
  coord_fixed(ratio = 100) +
  xlab("Seasonality constant <i>η</i>") + ylab("Waning constant <i>ρ</i>") +
  guides(fill = guide_legend(nrow = 1)) +
  
  plot_theme_paper +
  theme(legend.position = "bottom", legend.byrow = TRUE)

p_period

p_max <- ggplot()  +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.005, fill = "white")+
  geom_tile(aes(x = eta, y = rho, fill = inf_max),
            plot_data) +
  
  plot_annotations + 
  
  scale_fill_viridis_c(name = "Peak\ninfection\nprevalence", option = "mako") +
  coord_fixed(ratio = 100) +
  xlab("Seasonality constant <i>η</i>") + ylab("Waning constant <i>ρ</i>")  +
  
  plot_theme_paper +
  guides(fill = guide_colourbar(barwidth = 15)) +
  theme(legend.position = "bottom")

p_period | p_max



x_eta <- h5read("data/paper/period_over_grid_examples.jld2", "x_eta")
x_rho <- h5read("data/paper/period_over_grid_examples.jld2", "x_rho")
y_sol <- h5read("data/paper/period_over_grid_examples.jld2", "sol_t")

x_labels <- str_c(
  "<b>", c("iv", "iii", "i", "ii"), ".</b>",
  " <i>η</i> = ", scales::label_comma(accuracy = 0.01)(x_eta),
  ", <i>ρ</i> = ", scales::label_comma(accuracy = 0.0001)(x_rho), " — ",
  c("chaotic?", "quasiperiodic", "periodic (1 year)", "periodic (2 years)")
)

c_levels <- 2^seq(0, 8, by = 8 / 32)

plot_data_ex <- y_sol %>%
  reshape2::melt(c("i", "t", "class", "ix"), value.name = "prevalence") %>%
  mutate(eta = x_eta[i], rho = x_rho[i], label = x_labels[i]) %>%
  as_tibble()


plot_data_ex_inf <- plot_data_ex %>%
  filter(class == 2, t >= 365 * 166, t < 365 * (166 + 8)) %>% 
  group_by(label, eta, rho, t) %>%
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
  group_by(label, eta, rho, t) %>% 
  summarise(mean = sum(prevalence * c))

plot_data_ex_mean_year <- plot_data_ex %>%
  mutate(c = c_levels[ix]) %>%
  filter(class == 1, t >= 365 * 200, t < 365 * 230) %>% 
  group_by(label, eta, rho, t) %>% 
  summarise(mean = sum(prevalence * c)) %>%
  
  mutate(year = floor(t / 365),
         t_year = t %% 365)


p_ex_yearly_antibody <- ggplot() +
  geom_line(aes(x = t_year, y = mean, group = year),
            linewidth = 0.3, alpha = 0.3,
            plot_data_ex_mean_year) +
  
  facet_wrap(~label, ncol = 1, scales = "free_x") +
  
  xlab("Time (days)") + ylab("Mean antibody concentration") +
  
  scale_x_continuous(breaks = c(0, 90, 180, 270, 365)) +
  
  scale_y_continuous(trans = "log2", labels = scales::label_log(base = 2),
                     breaks = 2^c(1, 3, 5)) +
  
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



plot_data_attack_rate <- plot_data %>% 
  filter(rho > 0.0003) %>% 
  group_by(rho) %>% 
  filter(any(eta == 0)) %>% 
  mutate(inc_mean_eta_zero = inc_mean[eta == 0],
         inc_mean_diff = inc_mean / inc_mean_eta_zero,
         log_diff = pmax(log2(inc_mean_diff), -0.4),
         period_year = approxfun(plot_data_eta_zero$rho, plot_data_eta_zero$period)(rho) / 365 ) %>% 
  filter(eta %in% c(0, 0.1, 0.3, 0.5), rho < bifur_zero) %>% 
  mutate(eta_label = str_c("<i>η</i> = ", eta))

plot_data_year_marks_ar <- plot_data_year_marks %>% filter(year %in% c(0.5, 1, 1.5, 2, 3, 4))


p_attack_rate <- ggplot() +
  
  geom_line(aes(x = 1 / period_year, y = inc_mean * 365),
            colour = colour_C, linetype = "82",
            plot_data_attack_rate %>% filter(eta == 0) %>% ungroup() %>% select(-eta_label)) +
  
  geom_line(aes(x = 1 / period_year, y = inc_mean * 365),
            linewidth = 0.7,
            plot_data_attack_rate %>% filter(eta > 0)) +
  
  facet_wrap(~eta_label, ncol = 1) +
  
  coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 0.4)) +
  
  xlab("Natural frequency (years<sup>-1</sup>)") + ylab("Attack rate") +
  
  plot_theme_paper +
  theme(strip.text = element_markdown(),
        axis.text.x.top = element_text(margin = margin(b = 0.25, unit = "cm")),
        panel.grid.major.x = element_line(colour = "grey50", linetype = "28", linewidth = 0.5)) +
  
  ggtitle(NULL, "Yearly infection attack rate")

p_attack_rate


p_attack_rate_diff <- ggplot() +
  geom_hline(yintercept = 0, colour = colour_C, linetype = "82") +
  
  geom_line(aes(x = 1 / period_year, y = log_diff),
            linewidth = 0.7,
            plot_data_attack_rate %>% filter(eta > 0)) +
  
  facet_wrap(~eta_label, ncol = 1) +
  
  scale_y_continuous(breaks = log2(c(0.7, 0.8, 1, 1.25)),  labels = function(x) {scales::label_comma(accuracy = 0.01)(2^x)}) +
  
  coord_cartesian(ylim = c(-0.4, 0.4),
                  xlim = c(0, 2.5)) +
  
  xlab("Natural frequency (years<sup>-1</sup>)") + ylab("Proportional change") +
  
  plot_theme_paper +
  theme(strip.text = element_markdown(colour = "white"),
        axis.text.x.top = element_text(margin = margin(b = 0.25, unit = "cm")),
        panel.grid.major.x = element_line(colour = "grey50", linetype = "28", linewidth = 0.5)) +
  
  ggtitle(NULL, "Proportional change in\nyearly infection attack rate")


rho_breaks <- c(0.001, 0.002, 0.003)

rho_breaks_freq <- approxfun(plot_data_eta_zero$rho, 365 / plot_data_eta_zero$period)(rho_breaks)

p_axes_rho <- ggplot() +
  geom_blank(aes(x = 0)) +
  plot_theme_paper +
  
  scale_x_continuous(limits = c(0, 2.5),
                     breaks = c(0, rho_breaks_freq, 2.456),
                     labels = c(0, rho_breaks, 0.004)) +
  
  xlab("Waning constant <i>ρ</i>")

(p_attack_rate / p_axes_rho + plot_layout(heights = c(20, 1))) |
  (p_attack_rate_diff / p_axes_rho + plot_layout(heights = c(20, 1)))


ggsave(
  "results/results_seasonality_attack_rate.pdf",
  device = cairo_pdf,
  width = 9, height = 6.75,
  bg = "white"
)

