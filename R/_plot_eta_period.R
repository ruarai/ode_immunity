

library(tidyverse)
library(rhdf5)


source("../ode_immunity_multi/R/plot_theme.R")

t_steps <- h5read("data/paper/period_over_eta.jld2", "t_seq")
x_eta <- h5read("data/paper/period_over_eta.jld2", "x_eta")

y_inf <- h5read("data/paper/period_over_eta.jld2", "y_inf")
y_period <- h5read("data/paper/period_over_eta.jld2", "y_period")


plot_data_inf <- y_inf %>%
  reshape2::melt(varnames = c("eta", "t"), value.name = "prevalence") %>%
  mutate(eta = x_eta[eta],
         t = t_steps[t])

plot_data_period <- tibble(
  eta = x_eta,
  mean = y_period[,1], sd = y_period[,2], n = y_period[,3]
) %>%
  mutate(stable = (n > 4) & (sd < 0.1),
         period = if_else(stable, mean, NA),
         period_years = period / 365) 

ggplot(plot_data_period %>% filter(stable)) +
  geom_point(aes(x = period_years, y = sd))

plot_data_inf %>%
  filter(eta == eta[1]) %>%
  ggplot() +
  geom_line(aes(x = t, y = prevalence))

eta_width_half <- x_eta[2] / 2


p_period <- ggplot() +
  
  geom_rect(aes(xmin = eta - eta_width_half, xmax = eta + eta_width_half, ymin = -Inf, ymax = Inf),
           linewidth = 0.7, alpha = 0.1,
           fill = colour_C,
           plot_data_period %>% filter(is.na(period))) +
  
  geom_step(aes(x = eta, y = period_years),
            linewidth = 0.2,
            plot_data_period) +
  
  geom_linerange(aes(xmin = eta - eta_width_half, xmax = eta + eta_width_half, y = period_years, label = eta),
                 linewidth = 1.0,
                 plot_data_period %>% drop_na(period_years)) +
  
  scale_y_continuous(minor_breaks = 0:10, breaks = seq(0, 10, by = 2)) +
  
  coord_cartesian(ylim = c(0, 10), xlim = c(0, 0.5)) +
  
  ylab("Period (years)") +
  xlab("Seasonality factor η") +
  
  
  plot_theme_paper +
  theme(panel.grid.major.y = element_gridline)

p_period

eta_examples <- c(0.121, 0.2, 0.3, 0.4)

label_eta <- function(eta, period) {
  str_c("η = ", eta, ", period = ", period, " years")
}

plot_data_period_sub <- plot_data_period %>%
  filter(eta %in% eta_examples) %>%
  expand_grid(mult = seq(0, 40, by = 2)) %>%
  mutate(t = mult * period) %>%
  filter(t < 365 * 19.5) %>%
  mutate(label = label_eta(eta, round(period / 365)))

plot_data_sub <- plot_data_inf %>%
  filter(eta %in% eta_examples, t > 30000, t < 30000 + 365 * 20) %>%
  left_join(plot_data_period_sub %>% select(eta, label) %>% distinct(eta, label)) %>%
  mutate(t = t - 30000)


p_examples <- ggplot() +
  geom_line(aes(x = t, y = prevalence),
            plot_data_sub) +
  
  geom_rect(aes(xmin = t, xmax = pmin(t + period, 365 * 20), ymin = -Inf, ymax = Inf),
            alpha = 0.1, fill = colour_A,
            plot_data_period_sub) +
  
  scale_x_continuous(labels = scales::label_comma(scale = 1/365),
                     breaks = seq(0, 365 * 20, by = 365 * 2),) +
  
  facet_wrap(~label, ncol = 2, scales = "free_x") +
  
  coord_cartesian(xlim = c(0, 365 * 20)) +
  
  xlab("Time (years)") +
  ylab("Infection prevalence") +
  
  plot_theme_paper +
  theme(panel.grid.major.x = element_gridline,
        panel.grid.minor.x = element_gridline,
        panel.spacing.x = unit(1, "cm"))


p_period_narrow <- (plot_spacer() | p_period | plot_spacer()) +
  plot_layout(widths = c(1, 5, 1))

(p_period / p_examples) +
  plot_layout(heights = c(1, 2)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))


