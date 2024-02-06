
library(tidyverse)
library(rhdf5)
library(ggtext)

source("R/plot_theme.R")



sol_S_basic <- h5read("data/paper/stochastic_basic.jld2", "sol_S")
sol_I_basic <- h5read("data/paper/stochastic_basic.jld2", "sol_I")


c_levels <- 1:ncol(sol_I_basic) / 2


sol_basic <- data.table::data.table(
  t = 1:nrow(sol_I_basic), c = rep(c_levels, each = nrow(sol_I_basic)), I = c(sol_I_basic), S = c(sol_S_basic)
) %>%
  as_tibble() 


p1 <- sol_basic %>%
  summarise(I = sum(I), .by = t) %>% 
  ggplot() +
  geom_line(aes(x = t, y = I),
            linewidth = 0.5) +
  
  coord_cartesian(xlim = c(0, 365 * 15)) +
  
  scale_x_continuous(breaks = seq(0, 365 * 35, by = 365 * 2),
                     labels = seq(0, 35, by = 2)) +
  
  xlab("Time (years)") + ylab("Total Infected I") +
  
  plot_theme_paper +
  
  ggtitle("<b>A</b>")
p1



sol_S_eq <- h5read("data/paper/stochastic_around_equilbrium.jld2", "sol_S")
sol_I_eq <- h5read("data/paper/stochastic_around_equilbrium.jld2", "sol_I")



sol_eq <- data.table::data.table(
  t = 1:nrow(sol_I_eq), c = rep(c_levels, each = nrow(sol_I_eq)), I = c(sol_I_eq), S = c(sol_S_eq)
) %>%
  as_tibble() 


p2 <- sol_eq %>%
  summarise(I = sum(I), .by = t) %>% 
  ggplot() +
  geom_line(aes(x = t, y = I),
            linewidth = 0.5) +
  
  coord_cartesian(xlim = c(0, 365 * 15)) +
  
  scale_x_continuous(breaks = seq(0, 365 * 35, by = 365 * 2),
                     labels = seq(0, 35, by = 2)) +
  
  xlab("Time (years)") + ylab("Total Infected I") +
  
  plot_theme_paper +
  
  ggtitle("<b>B</b>")



sol_S_extinct <- h5read("data/paper/stochastic_extinction.jld2", "sol_S")
sol_I_extinct <- h5read("data/paper/stochastic_extinction.jld2", "sol_I")


sol_extinct <- data.table::data.table(
  t = 1:nrow(sol_I_extinct), c = rep(c_levels, each = nrow(sol_I_extinct)), I = c(sol_I_extinct), S = c(sol_S_extinct)
) %>%
  as_tibble() 


p3 <- sol_extinct %>%
  summarise(I = sum(I), .by = t) %>% 
  ggplot() +
  geom_line(aes(x = t, y = I),
            linewidth = 0.5) +
  
  coord_cartesian(xlim = c(0, 365 * 4)) +
  
  scale_x_continuous(breaks = seq(0, 365 * 4, by = 365 * 1),
                     labels = 0:4) +
  
  # scale_x_continuous(breaks = seq(0, 365 * 1, by = 365 * 2),
  #                    labels = seq(0, 34, by = 2)) +
  
  xlab("Time (years)") + ylab("Total Infected I") +
  
  plot_theme_paper +
  
  ggtitle("<b>C</b>")







inf_sims <- h5read("data/paper/stochastic_extinction_sims.jld2", "inf_sims")
ix_sim <- 1:ncol(inf_sims)


sol_survival <- data.table::data.table(t = 1:nrow(inf_sims), ix_sim = rep(ix_sim, each = nrow(inf_sims)), I = c(inf_sims)) %>%
  as_tibble() 



p4 <- sol_survival %>%
  mutate(extinct = I < 1) %>% 
  filter(t < 2000) %>% 
  group_by(t) %>%
  summarise(p_extinct = sum(extinct) / n()) %>% 
  ggplot() +
  
  annotate("rect", xmin = 0, xmax = 30, ymin = -Inf, ymax = Inf,
           fill = colour_B, alpha = 0.3) +
  
  annotate("rect", xmin = 140, xmax = 250, ymin = -Inf, ymax = Inf,
           fill = colour_A, alpha = 0.3) +
  
  annotate("rect", xmin = 380, xmax = 440, ymin = -Inf, ymax = Inf,
           fill = colour_A, alpha = 0.3) +
  
  geom_line(aes(x = t, y = p_extinct)) +
  
  scale_x_continuous(breaks = seq(0, 365 * 4, by = 365 * 1),
                     labels = 0:4) +
  
  coord_cartesian(xlim = c(0, 365 * 4)) +
  
  xlab("Time (years)") +
  ylab("Extinction probability") +
  
  plot_theme_paper +
  
  ggtitle("<b>D</b>")
p4


p_blank <- ggplot() + geom_blank() + theme_void()

cowplot::plot_grid(
  p1, p2, p_blank, p_blank, p3, p4,
  rel_widths = c(1, 0.1, 1),
  byrow = FALSE,
  ncol = 3
)



ggsave(
  "results/results_stochasticity.png",
  scale = 10 / 16,
  dpi = 300,
  width = 16, height = 9,
  bg = "white"
)







