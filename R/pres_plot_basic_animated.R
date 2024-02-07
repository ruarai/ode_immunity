

library(tidyverse)
library(rhdf5)
library(ggtext)
library(gganimate)

source("R/plot_theme.R")



sol_S <- h5read("data/anziam2024/basic.jld2", "sol_S")
sol_I <- h5read("data/anziam2024/basic.jld2", "sol_I")


c_levels <- 1:ncol(sol_I) / 2


sol <- data.table::data.table(t = 1:nrow(sol_I), c = rep(c_levels, each = nrow(sol_I)), I = c(sol_I), S = c(sol_S)) %>%
  as_tibble() %>%
  filter(t <= 1100)



p_summ <- sol %>%
  summarise(I = sum(I), .by = t) %>% 
  ggplot() +
  geom_line(aes(x = t, y = I),
            linewidth = 1.0) +
  
  coord_cartesian(xlim = c(0, 1100)) +
  
  scale_x_continuous(breaks = seq(0, 365 * 4, by = 365),
                     labels = 0:4) +
  
  xlab("Time _t_ (years)") + ylab("Total Infected I") +
  
  plot_theme

p_summ

animate(
  p_summ + transition_reveal(t),
  duration = 16,
  renderer = av_renderer(file = "results/summ.webm"),
  width = 680,
  height = 400
)

p_heatmap <- sol %>%
  
  mutate(S = pmin(S, 0.1)) %>% 
  
  ggplot() +
  geom_tile(aes(x = t, y = c, fill = S)) +
  
  scale_fill_viridis_c(option = "B", name = "Proportion",
                       breaks = c(0.001, 0.05, 0.1), labels = c("0.00", "0.05", "≥0.10")) +
  
  scale_x_continuous(breaks = seq(0, 365 * 4, by = 365),
                     labels = 0:4) +
  
  scale_y_continuous(expand = expansion()) +
  
  coord_cartesian(xlim = c(0, 1100), ylim = c(0, 8)) +
  
  xlab("Time _t_ (years)") + ylab("Correlate level _c_") +
  
  plot_theme + 
  theme(axis.text.y = element_markdown(margin = margin(l = 0.8, r = 0.2, unit = "cm"))) +
  theme(legend.position = "none")

animate(
  p_heatmap + transition_manual(t, cumulative = TRUE),
  duration = 16,
  renderer = av_renderer(file = "results/heatmap.webm"),
  width = 680,
  height = 400
)

p_density <- sol %>%
  
  ggplot() +
  
  geom_col(aes(y = c, x = S), orientation = "y") + 
  ylab(NULL) +
  xlab("Proportion") +
  
  coord_cartesian(xlim = c(0.5, 0), ylim = c(0, 8)) +
  
  scale_y_continuous(expand = expansion(), position = "right") +
  
  plot_theme

p_density

cowplot::plot_grid(p_density, p_heatmap, ncol = 2)


animate(
  p_density + transition_states(t),
  duration = 16,
  renderer = av_renderer(file = "results/density.webm"),
  width = 300,
  height = 400
) 




