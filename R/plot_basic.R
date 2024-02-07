

library(tidyverse)
library(rhdf5)
library(ggtext)

source("R/plot_theme.R")



sol_S <- h5read("data/anziam2024/basic.jld2", "sol_S")
sol_I <- h5read("data/anziam2024/basic.jld2", "sol_I")


c_levels <- 0:(ncol(sol_I) - 1) / 16


sol <- data.table::data.table(t = 1:nrow(sol_I), c = rep(c_levels, each = nrow(sol_I)), I = c(sol_I), S = c(sol_S)) %>%
  as_tibble() %>%
  filter(t <= 2000)



p_summ <- sol %>%
  summarise(I = sum(I), .by = t) %>% 
  ggplot() +
  geom_line(aes(x = t, y = I),
            linewidth = 0.7) +
  
  coord_cartesian(xlim = c(0, 2000)) +
  
  scale_x_continuous(breaks = seq(0, 365 * 8, by = 365),
                     labels = 0:8) +
  
  xlab("Time (years)") + ylab("Total Infected I") +
  
  plot_theme_paper

p_summ

p_heatmap <- sol %>%
  
  mutate(S = pmin(S, 0.1)) %>% 
  
  ggplot() +
  geom_tile(aes(x = t, y = c + 1 /32, fill = S, colour = S)) +
  
  scale_fill_viridis_c(option = "B", name = "Proportion",
                       breaks = c(0.001, 0.05, 0.1), labels = c("0.00", "0.05", "≥0.10"))  +
  
  scale_colour_viridis_c(option = "B", name = "Proportion",
                       breaks = c(0.001, 0.05, 0.1), labels = c("0.00", "0.05", "≥0.10")) +
  
  scale_x_continuous(breaks = seq(0, 365 * 8, by = 365),
                     labels = 0:8) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     sec.axis = sec_axis(~ . * 16, name = expression("Strata"~italic(i)), breaks = c(3.5, 7.5, 11.5, 15.5), c(4, 8, 12, 16))) +
  
  coord_cartesian(xlim = c(0, 2000), ylim = c(0, 1)) +
  
  xlab("Time (years)") + ylab("Correlate level _c_<sub>_i_</sub>") +
  
  plot_theme_paper + 
  theme(axis.text.y = element_markdown(margin = margin(l = 0.8, r = 0.2, unit = "cm"))) +
  theme(legend.position = "right") +
  
  ggtitle("<b>B</b>")

p_heatmap

cowplot::plot_grid(p_summ, p_heatmap, ncol = 1, align = "v", axis = "lr")


ggsave(
  "results/results_basic_SkI.png",
  scale = 10 / 16,
  dpi = 300,
  width = 14, height = 9,
  bg = "white"
)

