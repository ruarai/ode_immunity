library(tidyverse)
library(rhdf5)
library(patchwork)


source("R/plot_theme.R")

plot_data <- read_seasonality_data("data/period_over_grid.jld2") %>%
  mutate(eta_label = str_c("Seasonal forcing strength <i>η</i> = ", eta))


p_min <- ggplot() +
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



# p_max <- ggplot() +
#   geom_line(aes(x = r, y = inf_max),
#             linewidth = 0.7,
#             plot_data %>% filter(eta == 0) %>% rename(eta_label_null = eta_label),
#             colour = colour_C) +
#   geom_line(aes(x = r, y = inf_max),
#             linewidth = 0.7,
#             plot_data %>% filter(eta %in% c(0.1, 0.3, 0.5))) +
#   
#   facet_wrap(~eta_label, ncol = 1) +
#   
#   coord_cartesian(ylim = c(0, NA)) +
#   
#   ylab("Maximum infection prevalence") +
#   xlab("Antibody decay rate <i>r</i>") +
#   
#   plot_theme_paper +
#   
#   theme(strip.text = element_markdown())
# 
# p_min | p_max


p_min

ggsave(
  "results/results_supp_seasonality_extinction.pdf",
  device = cairo_pdf,
  width = 7, height = 6.75,
  bg = "white"
)
