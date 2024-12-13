

library(tidyverse)
library(patchwork)

source("R/plot_theme.R")

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

#\frac{C^2\rho\log_e(10)^2t}{k}
closed_data <- expand_grid(
  t = 0:t_max
) %>%
  mutate(mean = -C * rho * log(10) * t,
         log_var = C^2 * rho * log(10)^2 * t / k )

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

log_to_i_scale <- -32 / (C * log(10))

ggplot() +
  geom_tile(aes(x = t, y = i, fill = p_t),
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
  
  scale_y_continuous(
    transform = "reverse",
    breaks = c(0, 8, 16, 24, 32),
    sec.axis = sec_axis(
      transform = identity,
      breaks = c(0, 8, 16, 24, 32),
      labels = function(x) scales::label_log()(10^-(C * x / k)),
      name = "Decay factor <i>D</i>(<i>t</i>)"
    )
  ) +
  
  xlab("Time (days)") + ylab("Compartment <i>I</i>(<i>t</i>)") +
  
  plot_theme_paper +
  theme(legend.position = "bottom",
        legend.title = element_text(margin = margin(r = 0.7, unit = "cm")))

ggsave(
  "results/results_waning_probability.png",
  device = png,
  width = 5, height = 5,
  bg = "white"
)





