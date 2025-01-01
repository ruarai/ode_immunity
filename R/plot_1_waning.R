

library(tidyverse)
library(patchwork)

source("R/plot_theme.R")

C <- 8
k <- 32
rho <- 0.005
t_max <- 150

num_data <- expand_grid(
  t = seq(0, t_max, by = 0.1),
  i = 0:k
) %>%
  mutate(
    p_t = dpois(i, rho * k * t),
    log10_c = -(C * i / k),
    c = 10^log10_c,
    
    t_group = factor(t) %>% fct_rev()
  )


closed_data <- expand_grid(
  t = 0:t_max
) %>%
  mutate(mean = -C * rho * log(10) * t,
         log_var = C^2 * rho * log(10)^2 * t / k,
         
         mean_natural = rho * k * t * (10^(-C / k) - 1))

col_heatmap <- colorspace::sequential_hcl(n = 128, h = c(262, 136), c = c(39, 72, 0), l = c(13, 98), power = c(1.65, 1.1))

rescale_density <- function(x) {
  plogis(0.6 * qlogis(x) + 2)
}

trans_density <- scales::trans_new(
  "density", 
  function(x) plogis(0.6 * qlogis(x) + 2),
  function(x) plogis((qlogis(x) - 2) / 0.6),
  domain = c(0, 1)
)

log_to_i_scale <- -k / (C * log(10))

p_prob <- ggplot() +
  geom_tile(aes(x = t, y = i, fill = p_t, colour = p_t),
            num_data) +
  
  geom_line(aes(x = t, y = mean * log_to_i_scale),
            colour = "white", linewidth = 0.7,
            closed_data) +
  
  geom_line(aes(x = t, y = (mean + sqrt(log_var)) * log_to_i_scale),
            colour = "white", linetype = "24", linewidth = 0.7,
            closed_data) +
  
  geom_line(aes(x = t, y = (mean - sqrt(log_var)) * log_to_i_scale),
            colour = "white", linetype = "24", linewidth = 0.7,
            closed_data) +
  
  scale_fill_viridis_c(option = 8, name = "Probability",
                       breaks = c(0, 0.01, 0.1, 1.0),
                       transform = trans_density) +
  
  scale_colour_viridis_c(option = 8, name = "Probability",
                         breaks = c(0, 0.01, 0.1, 1.0),
                         transform = trans_density) +
  
  scale_y_continuous(
    transform = "reverse",
    breaks = c(0, 8, 16, 24, 32),
    sec.axis = sec_axis(
      transform = identity,
      breaks = c(0, 8, 16, 24, 32),
      labels = function(x) scales::label_log()(10^-(C * x / k)),
      name = "Decay <i>D</i>(<i>t</i>)"
    )
  ) +
  
  coord_cartesian(ylim = c(k, 0)) +
  
  xlab("Time (days)") + ylab("Transitions <i>I</i>(<i>t</i>)") +
  
  plot_theme_paper +
  theme(legend.position = "none",
        legend.title = element_text(margin = margin(r = 0.7, unit = "cm")))


p_legend <- tibble(
  x = seq(-0.05, 1.05, by = 0.001)
) %>%
  mutate(density = rescale_density(x)) %>% 
  fill(density, .direction = "updown") %>% 
  ggplot() +
  geom_tile(aes(x = x, y = 1, fill = density)) +
  
  geom_vline(aes(xintercept = x), tibble(x = seq(0, 1, by = 0.25)),
             colour = "white") +
  
  scale_fill_viridis_c(option = 8,
                       breaks = rescale_density(c(1e-8, 0.01, 0.05, 0.2, 1 - 1e-3)),
                       labels = c(0, 0.01, 0.05, 0.2, 1.0)) +
  
  scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
  
  plot_theme_paper + ylab(NULL) + xlab("Probability") +
  
  theme(legend.position = "none",
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), axis.line.x = element_blank(),
        axis.ticks.x = element_blank())


p_legend_space <- (plot_spacer() | p_legend | plot_spacer()) +
  plot_layout(widths = c(0.5, 5, 0.5))


(p_prob / p_legend_space) +
  plot_layout(heights = c(15, 1))


ggsave(
  "results/results_waning_probability.pdf",
  device = cairo_pdf(),
  width = 5, height = 5,
  bg = "white"
)



p_prob_2 <- ggplot() +
  geom_tile(aes(x = t, y = i, fill = p_t, colour = p_t),
            num_data) +
  
  geom_line(aes(x = t, y = mean * log_to_i_scale),
            colour = "white", linewidth = 0.7,
            closed_data) +
  
  geom_line(aes(x = t, y = mean_natural * log_to_i_scale),
            colour = "white", linetype = "44", linewidth = 0.7,
            closed_data) +
  
  scale_fill_viridis_c(option = 8, name = "Probability",
                       breaks = c(0, 0.01, 0.1, 1.0),
                       transform = trans_density) +
  
  scale_colour_viridis_c(option = 8, name = "Probability",
                         breaks = c(0, 0.01, 0.1, 1.0),
                         transform = trans_density) +
  
  scale_y_continuous(
    transform = "reverse",
    breaks = c(0, 8, 16, 24, 32),
    sec.axis = sec_axis(
      transform = identity,
      breaks = c(0, 8, 16, 24, 32),
      labels = function(x) scales::label_log()(10^-(C * x / k)),
      name = "Decay <i>D</i>(<i>t</i>)"
    )
  ) +
  
  coord_cartesian(ylim = c(k, 0)) +
  
  xlab("Time (days)") + ylab("Transitions <i>I</i>(<i>t</i>)") +
  
  plot_theme_paper +
  theme(legend.position = "none",
        legend.title = element_text(margin = margin(r = 0.7, unit = "cm")))


(p_prob_2 / p_legend_space) +
  plot_layout(heights = c(15, 1))



ggsave(
  "results/results_supp_waning_probability.png",
  device = png,
  width = 7, height = 7,
  bg = "white"
)


