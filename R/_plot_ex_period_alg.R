

library(tidyverse)
library(rhdf5)
library(patchwork)


source("R/plot_theme.R")

t_seq <- h5read("data/paper/ex_period.jld2", "t_seq")
y <- h5read("data/paper/ex_period.jld2", "y")

y_tbl <- y %>%
  reshape2::melt(varnames = c("ix", "t"), value.name = "prevalence") %>%
  mutate(t = t_seq[t])

t_start <- 365 * 50
t_end <- 365 * (50 + 5)

plot_data <- y_tbl %>%
  
  filter(t > t_start, t < t_end) %>% 
  
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
  
  ggokabeito::scale_colour_okabe_ito(name = "Class", order = c(9, 8)) +
  
  xlab("Time *t*") +
  ylab("__x__(*t*~i~)") +
  
  plot_theme_paper

p_usual

p_diff <- ggplot() +
  geom_line(aes(x = t, y = prev_norm, group = ix, colour = class),
            plot_data) +
  
  scale_x_continuous(labels = scales::label_comma()) +
  scale_y_continuous() +
  
  ggokabeito::scale_colour_okabe_ito(name = "Class", order = c(9, 8)) +
  
  xlab("Time *t*") +
  ylab("__x__(*t*~i~) - __x__(*t*~0~)") +
  
  plot_theme_paper


p_diff_norm <- ggplot() +
  geom_line(aes(x = t, y = dist),
            plot_data_dist) +
  
  geom_hline(yintercept = 10^-6, linetype = "44") +
  
  scale_x_continuous(labels = scales::label_comma()) +
  scale_y_log10(labels = scales::label_log()) +
  
  xlab("Time *t*") +
  ylab("||__x__(*t*~i~) - __x__(*t*~0~)||") +
  
  plot_theme_paper

p_usual / p_diff / p_diff_norm


ggsave(
  "results/results_ex_period_alg.png",
  device = png,
  width = 12, height = 12,
  bg = "white"
)
