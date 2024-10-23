

library(tidyverse)
library(rhdf5)


source("../ode_immunity_multi/R/plot_theme.R")



sol_t <- h5read("data/paper/basic_boosting.jld2", "sol_t")
c_levels <- h5read("data/paper/basic_boosting.jld2", "c_levels")

k <- 16

plot_data <- sol_t %>%
  reshape2::melt(varnames = c("scenario", "t", "class", "strata"), value.name = "prevalence") %>% 
  mutate(class = c("S", "I")[class],
         scenario = c("none", "linear", "loglinear")[scenario],
         scenario = factor(scenario, c("loglinear", "linear", "none"), labels = c("Log-linear boosting", "Linear boosting", "No boosting")),
         c = c_levels[strata])


plot_data_summ_inf <- plot_data %>%
  group_by(scenario, t, class) %>%
  summarise(prevalence = sum(prevalence), .groups = "drop") %>%
  filter(class == "I")







plot_data_means <- plot_data %>%
  group_by(t, scenario, class) %>%
  mutate(p = prevalence / sum(prevalence)) %>% 
  summarise(c = sum(p * c), .groups = "drop") %>%
  filter(class == "S")


cowplot::plot_grid(
  ggplot() +
    geom_line(aes(x = t, y = prevalence, colour = scenario),
              plot_data_summ_inf %>% filter(scenario != "None"),
              linewidth = 0.7) +
    geom_line(aes(x = t, y = prevalence, colour = scenario),
              plot_data_summ_inf %>% filter(scenario == "None"),
              linewidth = 1) +
    
    coord_cartesian(xlim = c(0, 2000)) +
    
    scale_x_continuous(breaks = scales::breaks_extended(),
                       labels = scales::label_comma()) +
    
    ggokabeito::scale_colour_okabe_ito(name = "Boosting", order = c(5, 3, 9)) +
    
    xlab("Time *t* (days)") + ylab("Infection prevalence") +
    
    plot_theme_paper +
    theme(panel.grid.major = element_gridline,
          legend.position = "none") +
    
    ggtitle("A — Infection prevalence"),
  ggplot() +
    geom_line(aes(x = t, y = c, colour = scenario),
              linewidth = 0.7,
              plot_data_means %>% filter(scenario != "None")) +
    geom_line(aes(x = t, y = c, colour = scenario),
              linewidth = 1.0,
              plot_data_means %>% filter(scenario == "None")) +
    
    coord_cartesian(xlim = c(0, 2000), ylim = c(2^-0.5, 2^6)) +
    
    xlab("Time *t* (days)") + ylab("Concentration") +
    
    ggokabeito::scale_colour_okabe_ito(name = "Scenario", order = c(5, 3, 9)) +
    
    scale_y_continuous(trans = "log2", labels = scales::label_log(base = 2), breaks = c(2^0, 2^2, 2^4, 2^6)) +
    
    scale_x_continuous(breaks = scales::breaks_extended(),
                       labels = scales::label_comma()) +
    
    plot_theme_paper +
    theme(panel.grid.major = element_gridline,
          legend.position = "bottom") +
    ggtitle("B — Population mean antibody concentration"),
  
  ncol = 1, rel_heights = c(1, 1.2)
)
  


ggsave(
  "results/results_basic_boosting.png",
  width = 10, height = 7,
  bg = "white"
)

boosting_matrices <- h5read("data/paper/basic_boosting.jld2", "boosting_matrices")

plot_data_matrices <- boosting_matrices %>%
  reshape2::melt(varnames = c("scenario", "i", "j"), value.name = "p") %>%
  mutate(
    scenario = c("none", "linear", "loglinear")[scenario],
    scenario = factor(scenario, c("loglinear", "linear", "none"), labels = c("Log-linear boosting", "Linear boosting", "No boosting"))
  ) %>%
  mutate(scenario = fct_rev(scenario)) %>%
  mutate(p = pmin(p, 0.2))


ggplot() +
  geom_tile(aes(x = j - 0.5, y = i - 0.5, fill = p),
            plot_data_matrices) +
  
  geom_abline(intercept = 0, slope = 1,
              colour = "white", linewidth = 0.7, linetype = "44") +
  
  facet_wrap(~scenario, ncol = 3) +
  
  geom_hline(yintercept = 32 * (3/8), colour = "white") +
  
  coord_fixed() +
  
  xlab("*I~i~*") + ylab("*S~j~*") +
  
  scale_fill_viridis_c(option = "A", name = "P(*I~i~* → *S~j~*)",
                       breaks = c(0.001, 0.1, 0.2), labels = c("0.0", "0.1", "≥0.2")) +
  
  plot_theme_paper +
  theme(legend.title = element_markdown()) +
  ggtitle("Transition probabilities under different boosting scenarios")


ggsave(
  "results/results_boosting_matrices.png",
  width = 10, height = 4.5,
  bg = "white"
)

