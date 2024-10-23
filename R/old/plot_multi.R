
library(tidyverse)
library(rhdf5)

source("R/plot_theme.R")



I_t_1 <- h5read("data/paper/_multi_results_3.jld2", "sums_t")
I_t_2 <- h5read("data/paper/_multi_results_2.jld2", "sums_t")




plot_data <- bind_rows(
  tibble(
    scenario = 1,
    S1 = I_t_1[2,],
    S2 = I_t_1[3,],
    t = 1:ncol(I_t_1)
  ),
  tibble(
    scenario = 2,
    S1 = I_t_2[2,],
    S2 = I_t_2[3,],
    t = 1:ncol(I_t_2)
  )
) %>%
  pivot_longer(c("S1", "S2"), names_prefix = "S", names_to = "strain", values_to = "I")


plot_data_sum <- plot_data %>%
  group_by(t, scenario) %>%
  summarise(I = sum(I))

p1 <- plot_data %>%
  filter(scenario == 1) %>% 
  ggplot() +
  
  # geom_line(aes(x = t, y = I, colour = "Sum"),
  #           plot_data_sum %>% filter(scenario == 1),
  #           linewidth = 0.7) +
  geom_line(aes(x = t, y = I, colour = strain),
            linewidth = 0.7) +
  
  ggokabeito::scale_colour_okabe_ito(name = "Strain", order = c(9, 5, 8)) +
  
  coord_cartesian(xlim = c(0, 1100), ylim = c(0, 0.3)) +
  
  scale_x_continuous(breaks = seq(0, 365 * 8, by = 365),
                     labels = 0:8) +
  
  
  xlab("Time (years)") + ylab("Total Infected I") +
  
  plot_theme_paper


p2 <- plot_data %>%
  filter(scenario == 2) %>% 
  ggplot() +
  
  # geom_line(aes(x = t, y = I, colour = "Sum"),
  #           plot_data_sum %>% filter(scenario == 2), 
  #           linewidth = 0.7)  +
  geom_line(aes(x = t, y = I, colour = strain),
            linewidth = 0.7)+
  
  ggokabeito::scale_colour_okabe_ito(name = "Strain", order = c(9, 5, 8)) +
  
  coord_cartesian(xlim = c(0, 1100), ylim = c(0, 0.3)) +
  
  scale_x_continuous(breaks = seq(0, 365 * 8, by = 365),
                     labels = 0:8) +
  
  
  xlab("Time (years)") + ylab("Total Infected I") +
  
  plot_theme_paper


cowplot::plot_grid(p1, p2, ncol = 1, align = "v", axis = "lr")


ggsave(
  "results/_multi.png",
  scale = 10 / 16,
  dpi = 300,
  width = 16, height = 9,
  bg = "white"
)
