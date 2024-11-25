
library(tidyverse)
library(rhdf5)

source("R/plot_theme.R")



I_SI <- h5read("data/paper/classical.jld2", "I_SI")
I_SIR <- h5read("data/paper/classical.jld2", "I_SIR")
I_SIRS <- h5read("data/paper/classical.jld2", "I_SIRS")

I_structured <- h5read("data/paper/classical.jld2", "I_structured")
I_dde <- h5read("data/paper/classical.jld2", "I_dde")


plot_data <- tibble(
  SI = I_SI, SIR = I_SIR, SIRS = I_SIRS,
  
  structured = I_structured, dde = I_dde,
  
  t = 1:length(I_SI)
) %>%
  pivot_longer(-c(t))



plot_data %>%
  filter(name %in% c("SIR")) %>% 
  ggplot() +
  geom_line(aes(x = t, y = value, group = name),
            linewidth = 0.7) +
  
  coord_cartesian(xlim = c(0, 450)) +
  
  
  xlab("Time (days)") + ylab("Total Infected I") +
  
  plot_theme_paper +
  theme(legend.position = "bottom")



plot_data %>%
  filter(name %in% c("SIRS")) %>% 
  ggplot() +
  geom_line(aes(x = t, y = value, group = name),
            linewidth = 0.7) +
  
  coord_cartesian(xlim = c(0, 1000),
                  ylim = c(0, 0.08)) +
  
  
  xlab("Time (days)") + ylab("Total Infected I") +
  
  plot_theme_paper +
  theme(legend.position = "bottom")


plot_data %>%
  filter(name %in% c("dde")) %>% 
  ggplot() +
  geom_line(aes(x = t, y = value, group = name),
            linewidth = 0.7) +
  
  coord_cartesian(xlim = c(0, 1500),
                  ylim = c(0, 0.08)) +
  
  
  xlab("Time (days)") + ylab("Total Infected I") +
  
  plot_theme_paper +
  theme(legend.position = "bottom")


plot_data %>%
  filter(name %in% c("structured")) %>% 
  ggplot() +
  geom_line(aes(x = t, y = value, group = name),
            linewidth = 0.7) +
  
  coord_cartesian(xlim = c(0, 1500),
                  ylim = c(0, 0.08)) +
  
  
  xlab("Time (days)") + ylab("Total Infected I") +
  
  plot_theme_paper +
  theme(legend.position = "bottom")


plot_data %>%
  filter(name %in% c("structured", "dde")) %>% 
  ggplot() +
  geom_line(aes(x = t, y = value, group = name, colour = name),
            linewidth = 0.7) +
  
  ggokabeito::scale_colour_okabe_ito(order = c(3, 9)) +
  
  coord_cartesian(xlim = c(0, 1500),
                  ylim = c(0, 0.08)) +
  
  
  xlab("Time (days)") + ylab("Total Infected I") +
  
  plot_theme_paper +
  theme(legend.position = "none")




plot_data %>%
  filter(name %in% c("structured", "dde", "SIR", "SIRS")) %>% 
  mutate(name = factor(name,
                       c("SIR", "SIRS", "dde", "structured"),
                       c("<b>A</b> \u2013 SIR",
                         "<b>B</b> \u2013 SIRS",
                         "<b>C</b> \u2013 Delay differential equation",
                         "<b>D</b> \u2013 Immunity-structured"))) %>% 
  ggplot() +
  geom_line(aes(x = t, y = value),
            linewidth = 0.7) +
  
  coord_cartesian(xlim = c(0, 365 * 8)) +
  
  scale_x_continuous(breaks = seq(0, 365 * 8, by = 365),
                     labels = 0:8) +
  
  facet_wrap(~name, dir = "v", scales = "free_x") +
  

  xlab("Time (years)") + ylab("Total Infected I") +
  
  plot_theme_paper +
  theme(legend.position = "bottom",
        strip.text = element_markdown(hjust = 0))

ggsave(
  "results/results_classical.png",
  scale = 10 / 16,
  dpi = 300,
  device = png,
  width = 18, height = 9,
  bg = "white"
)


ggsave(
  "results/results_classical.pdf",
  scale = 10 / 16,
  dpi = 300,
  width = 18, height = 9,
  bg = "white"
)





