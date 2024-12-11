library(tidyverse)
library(rhdf5)
library(patchwork)


source("R/plot_theme.R")


x_vals <- h5read("data/paper/period_over_grid.jld2", "x_vals")
y_inf_summary <- h5read("data/paper/period_over_grid.jld2", "y_inf_summary")
y_period <- h5read("data/paper/period_over_grid.jld2", "y_period")

plot_data <- tibble(
  eta = x_vals[1, ], r = x_vals[2, ],
  inf_min = y_inf_summary[, 1], inf_max = y_inf_summary[, 2],
  inf_mean = y_inf_summary[, 3], inf_chaos = y_inf_summary[, 4],
  
  inc_min = y_inf_summary[, 5], inc_max = y_inf_summary[, 6],  
  inc_mean = y_inf_summary[, 7], inc_chaos = y_inf_summary[ , 8],
  lyapunov = y_inf_summary[, 9],
  
  period = y_period[,1], period_sd = y_period[,2], period_n = y_period[,3]
) %>%
  mutate(inf_diff = inf_max - inf_min,
         periodic = (period_sd < 1) & (period_n >= 5),
         quasiperiodic = (period_sd >= 1) & (period_n >= 5))


plot_data_eta_zero_periodic <- plot_data %>% filter(eta == 0, r < 0.065)
bifur_zero <- plot_data %>% filter(eta == 0, inf_diff < 1e-3) %>% pull(r) %>% head(1)


plot_data_incidence <- plot_data %>% 
  filter(r > 0.003) %>% 
  group_by(r) %>% 
  filter(any(eta == 0)) %>% 
  mutate(inc_mean_eta_zero = inc_mean[eta == 0],
         inc_mean_diff = inc_mean / inc_mean_eta_zero,
         log_diff = pmax(log2(inc_mean_diff), -0.4),
         period_year = approxfun(plot_data_eta_zero_periodic$r, plot_data_eta_zero_periodic$period)(r) / 365 ) %>% 
  filter(eta %in% c(0, 0.1, 0.3, 0.5)) %>% 
  mutate(eta_label = str_c("Seasonality strength <i>η</i> = ", eta))


freq_breaks <- seq(0.5, 2.5, by = 0.5)
freq_breaks_r <- approxfun(365 / plot_data_eta_zero_periodic$period, plot_data_eta_zero_periodic$r)(freq_breaks)

p_incidence <- ggplot() +
  geom_vline(aes(xintercept = r), tibble(r = c(0, freq_breaks_r)),
             linetype = "44", colour = "grey70") +
  
  annotate("rect", xmin = bifur_zero, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.2) +
  
  geom_line(aes(x = r, y = inc_mean * 365),
            linewidth = 0.7,
            plot_data_incidence %>% filter(eta > 0)) +
  
  geom_line(aes(x = r, y = inc_mean * 365),
            colour = colour_C, linetype = "82",
            plot_data_incidence %>% filter(eta == 0) %>% ungroup() %>% select(-eta_label)) +
  
  facet_wrap(~eta_label, ncol = 1) +
  
  coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 2.0)) +
  
  xlab("Mean antibody decay rate <i>r</i>") + ylab("Yearly infection incidence") +
  
  plot_theme_paper +
  theme(strip.text = element_markdown(),
        axis.text.x.top = element_text(margin = margin(b = 0.25, unit = "cm"))) +
  
  ggtitle(NULL, "Yearly infection incidence")

p_incidence

p_incidence_diff <- ggplot() +
  geom_vline(aes(xintercept = r), tibble(r = c(0, freq_breaks_r)),
             linetype = "44", colour = "grey70") +
  
  annotate("rect", xmin = bifur_zero, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.2) +
  
  geom_hline(yintercept = 1, colour = colour_C, linetype = "82") +
  
  geom_line(aes(x = r, y = inc_mean_diff),
            linewidth = 0.7,
            plot_data_incidence %>% filter(eta > 0)) +
  
  facet_wrap(~eta_label, ncol = 1) +
  
  scale_y_continuous(breaks = c(0.75, 1, 1.25)) +
  coord_cartesian(xlim = c(0, 0.1)) +
  
  xlab("Mean antibody decay rate <i>r</i>") + ylab("Proportional difference") +
  
  plot_theme_paper +
  theme(strip.text = element_markdown(colour = "white"),
        axis.text.x.top = element_text(margin = margin(b = 0.25, unit = "cm"))) +
  
  ggtitle(NULL, "Proportional difference in<br>yearly infection incidence")

p_incidence_diff 

p_axes_freq <- ggplot() +
  geom_blank(aes(x = 0)) +
  annotate("linerange", xmin = -Inf, xmax = bifur_zero, y = 0.0, linewidth = 0.7)  +
  annotate("linerange", xmin = bifur_zero, xmax = Inf, y = 0.0, linewidth = 0.7, linetype = "44") +
  plot_theme_paper +
  
  geom_linerange(aes(ymin = 0, ymax = 0.02, x = x), tibble(x = c(0.0, freq_breaks_r)),
                 linewidth = 0.7) +
  
  annotate("point", x = bifur_zero, y = 0, size = 2.5) +
  annotate("point", x = bifur_zero, y = 0, size = 1.25, colour = "white") +
  
  
  scale_x_continuous(limits = c(0, 0.1),
                     breaks = c(0.0, freq_breaks_r, bifur_zero),
                     labels = c(0.0, freq_breaks, "")) +
  
  coord_cartesian(ylim = c(0, 0)) +
  
  xlab("Natural frequency without seasonality (years<sup>-1</sup>)") + ylab(NULL) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

p_axes_freq


p_incidence / p_axes_freq + plot_layout(heights = c(20, 1))

((p_incidence / p_axes_freq + plot_layout(heights = c(20, 1))) |
    (p_incidence_diff / p_axes_freq + plot_layout(heights = c(20, 1))))


ggsave(
  "results/results_seasonality_resonance.pdf",
  device = cairo_pdf,
  width = 11, height = 6.75,
  bg = "white"
)
