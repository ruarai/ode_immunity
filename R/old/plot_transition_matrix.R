
library(tidyverse)
library(rhdf5)
library(ggtext)


mat_jump <- h5read("data/anziam2024/basic.jld2", "mat_jump")
c_levels <- h5read("data/anziam2024/basic.jld2", "c_levels")
mat_jump_no_boost <- h5read("data/anziam2024/basic.jld2", "mat_jump_no_boost")

cowplot::plot_grid(
  mat_jump_no_boost %>%
    reshape2::melt(varnames = c("i", "j")) %>%
    mutate(c_i = c_levels[i], c_j = c_levels[j],
           value = pmin(value, 0.3)) %>% 
    
    ggplot() +
    geom_tile(aes(x = c_j, y = c_i, fill = value)) +
    
    xlab("Pre-infection<br/>correlate level _j_") +
    
    ylab("Post-infection<br/>correlate level _i_") +
    
    scale_fill_viridis_c("Probability", breaks = c(0.0001, 0.15, 0.3), labels = c("0.00", "0.15", ">0.30"), limit = c(0, 0.3)) +
    coord_fixed()  +
    
    plot_theme +
    theme(legend.position = "none"),
  mat_jump %>%
    reshape2::melt(varnames = c("i", "j")) %>%
    mutate(c_i = c_levels[i], c_j = c_levels[j],
           value = pmin(value, 0.3)) %>% 
    
    ggplot() +
    geom_tile(aes(x = c_j, y = c_i, fill = value)) +
    
    xlab("Pre-infection<br/>correlate level _j_") +
    
    ylab("Post-infection<br/>correlate level _i_") +
    
    scale_fill_viridis_c("Probability", breaks = c(0.0001, 0.15, 0.3), labels = c("0.00", "0.15", ">0.30"), limit = c(0, 0.3)) +
    coord_fixed() +
    
    plot_theme +
    theme(legend.position = "none"),
  ncol = 1
)



mat_jump_no_boost %>%
  reshape2::melt(varnames = c("i", "j")) %>%
  mutate(c_i = c_levels[i], c_j = c_levels[j],
         value = pmin(value, 0.3)) %>% 
  
  ggplot() +
  geom_tile(aes(x = c_j, y = c_i, fill = value)) +
  
  xlab("Pre-infection<br/>correlate level _j_") +
  
  ylab("Post-infection<br/>correlate level _i_") +
  
  scale_fill_viridis_c("Probability", breaks = c(0.0001, 0.15, 0.3), labels = c("0.00", "0.15", ">0.30"), limit = c(0, 0.3)) +
  
  plot_theme +
  theme(legend.position = "none")
