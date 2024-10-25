

library(tidyverse)
library(rhdf5)


source("../ode_immunity_multi/R/plot_theme.R")

sol_t <- h5read("data/paper/basic_boosting.jld2", "sol_t")
c_levels <- h5read("data/paper/basic_boosting.jld2", "c_levels")

plot_data <- sol_t %>%
  reshape2::melt(varnames = c("scenario", "t", "class", "strata"), value.name = "prevalence") %>% 
  mutate(class = c("S", "I")[class], c = c_levels[strata]) %>%
  filter(scenario == 1, t < 2000)


plot_data_summ_inf <- plot_data %>%
  group_by(t, class) %>%
  summarise(prevalence = sum(prevalence), .groups = "drop") %>%
  filter(class == "I")

plot_data_means <- plot_data %>%
  group_by(t, class) %>%
  mutate(p = prevalence / sum(prevalence)) %>% 
  summarise(c = sum(p * c), .groups = "drop") %>%
  filter(class == "S")




p_summ <- plot_data_summ_inf %>% 
  ggplot() +
  geom_line(aes(x = t, y = prevalence),
            linewidth = 0.7) +
  
  coord_cartesian(xlim = c(0, 2000)) +
  
  scale_x_continuous(breaks = scales::breaks_extended(),
                     labels = scales::label_comma()) +
  
  xlab("Time *t* (days)") + ylab("Infection<br>prevalence") +
  
  plot_theme_paper +
  theme(panel.grid.major = element_gridline) +
  
  # scale_y_log10() +
  
  ggtitle(NULL, "Infection prevalence")

p_heatmap <- plot_data %>%
  
  summarise(prevalence = sum(prevalence), .by = c("c", "t")) %>% 
  mutate(prevalence = pmin(prevalence, 0.05)) %>% 
  
  ggplot() +
  geom_tile(aes(x = t, y = c, fill = prevalence)) +
  
  geom_line(aes(x = t, y = c), plot_data_means, colour = "black", linewidth = 1.0) +
  geom_line(aes(x = t, y = c), plot_data_means, colour = "white", linewidth = 0.3) +
  
  scale_fill_viridis_c(option = "B", name = "Proportion",
                       breaks = c(0.001, 0.025, 0.05), labels = c("0.00", "0.025", "≥0.05"))  +
  
  scale_x_continuous(breaks = scales::breaks_extended(),
                     labels = scales::label_comma()) +
  
  scale_y_continuous(trans = "log2", breaks = 2^c(0, 2, 5, 8),
                     sec.axis = sec_axis(
                       transform = "identity", 
                       labels = function(x) log2(x) + 1,
                       breaks = 2^c(0, 2, 5, 8),
                       name = "Strata *i*"
                     )) +
  
  coord_cartesian(xlim = c(0, 2000), ylim = c(2^0, 2^8)) +
  
  xlab("Time *t* (days)") + ylab("Antibody<br>concentration *c~i~*") +
  
  plot_theme_paper + 
  theme(legend.position = "bottom",
        legend.title = element_text(margin = margin(b = 15, r = 20))) +
  
  ggtitle(NULL, "Population distribution of antibodies")

(p_summ / p_heatmap) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))



ggsave(
  "results/results_basic.png",
  width = 10, height = 8,
  bg = "white"
)


