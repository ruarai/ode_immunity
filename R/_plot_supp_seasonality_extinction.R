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
  mutate(eta_label = str_c("Seasonality constant <i>η</i> = ", eta))


ggplot() +
  geom_line(aes(x = r, y = inf_min),
            linewidth = 0.7,
            plot_data %>% filter(eta == 0) %>% rename(eta_label_null = eta_label),
            colour = colour_C) +
  geom_line(aes(x = r, y = inf_min),
            linewidth = 0.7,
            plot_data %>% filter(eta %in% c(0.1, 0.3, 0.5))) +
  
  scale_y_log10(labels = scales::label_log(),
                breaks = 10^c(-12, -8, -4, 0)) +
  
  facet_wrap(~eta_label, ncol = 1) +
  
  coord_cartesian(ylim = 10^c(-12, 0)) +
  
  ylab("Minimum infection prevalence") +
  xlab("Antibody decay rate <i>r</i>") +
  
  plot_theme_paper +
  
  theme(strip.text = element_markdown())


ggsave(
  "results/results_supp_seasonality_extinction.pdf",
  device = cairo_pdf,
  width = 7, height = 6.75,
  bg = "white"
)
