

library(tidyverse)
library(rhdf5)
library(patchwork)


source("../ode_immunity_multi/R/plot_theme.R")

t_steps <- h5read("data/paper/ex_period.jld2", "t")
y <- h5read("data/paper/ex_period.jld2", "y")

y_tbl <- y %>%
  reshape2::melt(varnames = c("ix", "t"), value.name = "prevalence") %>%
  mutate(t = t_steps[t])

plot_data <- y_tbl %>%
  
  filter(t > 365 * 35.5, t < 365 * 45) %>% 
  
  group_by(ix) %>%
  mutate(prev_norm = prevalence - prevalence[1],
         dist = prev_norm ^ 2,
         class = if_else(ix > 33, "I", "S"))


plot_data_dist <- plot_data %>%
  group_by(t) %>%
  summarise(dist = sum(dist))


p_usual <- ggplot() +
  geom_line(aes(x = t, y = prevalence, group = ix, colour = class),
            plot_data) +
  
  scale_x_continuous(labels = scales::label_comma()) +
  scale_y_continuous() +
  
  ggokabeito::scale_colour_okabe_ito(order = c(5, 6)) +
  
  xlab("Time *t*") +
  ylab("__x__(*t*~0~) - __x__(*t*~i~)") +
  
  plot_theme_paper


p_usual

p_diff <- ggplot() +
  geom_line(aes(x = t, y = prev_norm, group = ix, colour = class),
            plot_data) +
  
  scale_x_continuous(labels = scales::label_comma()) +
  scale_y_continuous() +
  
  ggokabeito::scale_colour_okabe_ito(order = c(5, 6)) +
  
  xlab("Time *t*") +
  ylab("__x__(*t*~0~) - __x__(*t*~i~)") +
  
  plot_theme_paper


p_diff_norm <- ggplot() +
  geom_line(aes(x = t, y = dist),
            plot_data_dist) +
  
  geom_hline(yintercept = 10^-6, linetype = "44") +
  
  scale_x_continuous(labels = scales::label_comma()) +
  scale_y_log10(labels = scales::label_log()) +
  
  xlab("Time *t*") +
  ylab("||__x__(*t*~0~) - __x__(*t*~i~)||") +
  
  plot_theme_paper

p_usual / p_diff / p_diff_norm


ggsave(
  "results/results_ex_period_alg.png",
  device = png,
  width = 12, height = 12,
  bg = "white"
)


plot_sub <- y_tbl %>%
  filter(t > 365*36) %>%
  mutate(t_year = (t - min(t)) / 365)


p_heatmap <- plot_sub %>%
  filter(ix <= 33) %>% 
  mutate(prevalence = pmin(prevalence, 0.05)) %>% 
  ggplot() +
  geom_tile(aes(x = t_year, y = ix, fill = prevalence)) +
  
  scale_x_continuous(breaks = 0:5) +
  
  scale_fill_viridis_c(option = "B") +
  
  plot_theme_paper +
  theme(legend.position = "none")

