

library(tidyverse)
library(patchwork)

source("R/plot_theme.R")




k <- 32
lambda <- 0.01
C <- 8
t_max <- 150


c_levels <- 10^(C * ((0:k) / k))
log_c_levels <- log10(c_levels)
plot_data <- expand_grid(
  i = 0:k,
  t = 0:t_max
) %>%
  mutate(p = ((lambda * k * t)^i / factorial(i)) * exp(-lambda * k * t),
         c_loss = 1 - 10^(-C * i / k),
         log_c = log_c_levels[k - i + 1],
         c = c_levels[k - i + 1]) %>%
  
  group_by(t) %>%
  mutate(p = if_else(i == k, 1 - sum(p) + p, p)) %>%
  ungroup()

plot_data_avg <- plot_data %>% 
  summarise(mean = sum(p * c_levels[k - i + 1]),
            var = sum(p * c_levels[k - i + 1] ^ 2) - mean ^ 2,
            .by = "t")

plot_data_closed <- tibble(t = 0:t_max) %>%
  mutate(
    mean = 10^8 * exp((-1 + 10^(-C / k)) * k * lambda * t),
    var = 10^(2 * C) * exp(lambda * k * t * (10^(-2 * C / k) - 1)) - 10^(2 * C) * exp(2 * lambda * k * t * (10^(-C / k) - 1))
  )
  
p_heatmap <- ggplot() +
  geom_tile(aes(x = t, y = c, fill = pmin(0.1, p)),
            plot_data) +
  
  geom_line(aes(x = t, y = mean, colour = "Compartmental model"), plot_data_avg,
            linewidth = 1.5) +
  
  scale_y_continuous(trans = "log10", breaks = 10^c(0, 2, 4, 6, 8),
                     labels = scales::label_log(base = 10),
                     sec.axis = sec_axis(
                       transform = "identity", 
                       labels = function(x) log10(x) * 4 + 0,
                       breaks = 10^c(0, 2, 4, 6, 8),
                       name = "Strata"
                     )) +
  
  scale_colour_manual(
    name = "Mean",
    values = c("Compartmental model" = ggokabeito::palette_okabe_ito(3))
  ) +
  
  ylab("Serum antibody<br>concentration c(t)") +
  
  xlab("Time <i>t</i> (days)") +
  
  scale_fill_viridis_c(option = "B", name = "Proportion",
                       breaks = c(0.001, 0.05, 0.1), labels = c("0.00", "0.05", "≥0.1")) +
  plot_theme_paper

p_heatmap

p_mean <- ggplot() +
  
  geom_line(aes(x = t, y = mean, colour = "Compartmental model"), plot_data_avg,
            linewidth = 1.5) +
  
  geom_line(aes(x = t, y = mean, colour = "Closed form approximation"), plot_data_closed,
            linewidth = 1.0) +
  
  scale_y_continuous(trans = "log10",
                     breaks = 10^c(0, 2, 4, 6, 8),
                     labels = scales::label_log(base = 10)) +
  
  
  scale_colour_manual(
    name = "Mean",
    values = c("Compartmental model" = ggokabeito::palette_okabe_ito(3),
               "Closed form approximation" = ggokabeito::palette_okabe_ito(6)),
    breaks = c("Compartmental model", "Closed form approximation")
  ) +
  
  ylab("Serum antibody<br>concentration c(t)") +
  
  xlab("Time <i>t</i> (days)") +
  
  plot_theme_paper


p_diff <- plot_data_avg %>%
  select(t, mean_exp = mean) %>%
  left_join(plot_data_closed, by = "t") %>%
  
  mutate(mean_diff = (mean - mean_exp)) %>%
  
  ggplot() +
  geom_line(aes(x = t, y = mean_diff)) +
  
  xlab("Time <i>t</i> (days)") +
  ylab("Bias") +
  
  plot_theme_paper

(p_heatmap / p_mean ) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))


ggsave(
  "results/results_supp_waning_model.png",
  device = png,
  width = 12, height = 7,
  bg = "white"
)

plot_data_closed_k <- expand_grid(t = 0:150, k = c(8, 32, 128, 512, 2048)) %>%
  mutate(
    mean = 10^8 * exp((-1 + 10^(-C / k)) * k * lambda * t),
    var = exp(lambda * k * t * (10^(-2 * C / k) - 1)) - exp(2 * lambda * k * t * (10^(-C / k) - 1))
  )

p_sd <- ggplot() +
  geom_line(aes(x = t, y = sqrt(var), colour = factor(k), group = k ),
            linewidth = 1.0,
            plot_data_closed_k) +
  
  xlab("Time <i>t</i> (days)") +
  
  ylab("sd") +
  
  scale_y_log10(labels = scales::label_log()) +
  scale_colour_manual(name = "<i>k</i>",
                      values = RColorBrewer::brewer.pal(8, "YlGn")[4:8]) +
  
  plot_theme_paper +
  theme(legend.title = element_markdown())

p_sd

p_sd_ratio <- ggplot() +
  geom_line(aes(x = t, y = sqrt(var) / mean, colour = factor(k), group = k ),
            linewidth = 1.0,
            plot_data_closed_k) +
  
  xlab("Time <i>t</i> (days)") +
  
  ylab("sd/mean") +
  
  scale_y_log10(labels = scales::label_log())  +
  scale_colour_manual(name = "<i>k</i>",
                      values = RColorBrewer::brewer.pal(8, "YlGn")[4:8]) +
  
  plot_theme_paper +
  theme(legend.title = element_markdown())

(p_sd / p_sd_ratio) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))


ggsave(
  "results/results_supp_waning_model_k.png",
  device = png,
  width = 10, height = 7,
  bg = "white"
)




