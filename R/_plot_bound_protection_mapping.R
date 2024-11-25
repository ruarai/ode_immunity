
library(tidyverse)
library(rhdf5)
library(patchwork)

source("R/plot_theme.R")


plot_data <- expand_grid(
  p_bound = seq(0.0, 1, by = 0.001),
  #log_a = seq(-2, 2, by = 2),
  prob_half = c(0.25, 0.5, 0.75),
  h = c(0.5, 1, 2)
) %>%
  mutate(logit_half = qlogis(prob_half),
         p_neut = plogis(h * qlogis(p_bound) - h * logit_half),
         h_label = str_c("<i>h</i> = ", h),
         prob_half_label = str_c("<i>p</i><sub>half</sub> = ", prob_half))


p_logit <- ggplot() +
  geom_line(aes(x = p_bound, y = p_neut, colour = factor(prob_half), group = prob_half),
            plot_data,
            linewidth = 1) +
  
  facet_wrap(~h_label, scales = "free") +
  
  xlab("*p*~bound~") +
  ylab("*p*~protection~") +
  
  ggokabeito::scale_colour_okabe_ito(order = c(5, 9, 3),
                                     name = "<i>p</i><sub>half</sub>") +
  
  coord_cartesian(xlim = c(0.05, 0.95), ylim = c(0.05, 0.95)) +
  scale_x_continuous(trans = "logit", breaks = c(0.05, 0.2, 0.5, 0.8, 0.95)) +
  scale_y_continuous(trans = "logit", breaks = c(0.05, 0.2, 0.5, 0.8, 0.95)) +
  
  plot_theme_paper +
  theme(panel.grid.major = element_gridline,
        legend.title = element_markdown(),
        strip.text = element_markdown(),
        panel.spacing.x = unit(1, "cm"),
        legend.position = "none")
p_logit

  
p_natural <- ggplot() +
  geom_line(aes(x = p_bound, y = p_neut, colour = factor(prob_half), group = prob_half),
            plot_data,
            linewidth = 1) +
  
  geom_point(aes(x = p_bound, y = p_neut, colour = factor(prob_half)),
             plot_data %>% filter(p_neut == 0.5)) +
  
  facet_wrap(~h_label, scales = "free") +
  
  xlab("*p*~bound~") +
  ylab("*p*~protection~") +
  
  ggokabeito::scale_colour_okabe_ito(order = c(5, 9, 3),
                                     name = "<i>p</i><sub>half</sub>") +
  
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  
  plot_theme_paper +
  
  guides(colour = guide_legend(override.aes = list(size = 0))) +
  
  theme(panel.grid.major = element_gridline,
        legend.title = element_markdown(),
        strip.text = element_markdown(),
        panel.spacing.x = unit(1, "cm"),
        legend.position = "bottom")

p_natural

(p_logit / p_natural) +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))


ggsave(
  "results/bound_protection_mapping.pdf",
  width = 10, height = 7.8,
  device = cairo_pdf,
  bg = "white"
)

ggplot() +
  geom_line(aes(x = p_bound, y = p_neut, colour = factor(prob_half), group = prob_half),
            plot_data,
            linewidth = 1) +
  
  facet_wrap(~h, scales = "free") +
  
  xlab("*p*~bound~") +
  ylab("*p*~protection~") +
  
  ggokabeito::scale_colour_okabe_ito(order = c(5, 9, 3),
                                     name = "h") +
  
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  
  plot_theme_paper +
  theme(panel.grid.major = element_gridline,
        legend.title = element_markdown(),
        strip.text = element_markdown(),
        panel.spacing.x = unit(1, "cm"),
        legend.position = "bottom")


tibble(
  c = 10^seq(-5,5, by = 0.1),
  h = 1, b = 2
) %>%
  mutate(y1 = c ^ b / (h^b + c^b),
         y2 = 1 / (1 + exp(-b * (log(c) - log(h))))
         ) %>%
  
  ggplot() +
  geom_line(aes(x = c, y = y1)) +
  geom_line(aes(x = c, y = y2),
            colour = "red", linetype = "44") +
  
  scale_x_log10()



