

library(tidyverse)
library(rhdf5)

source("../ode_immunity_multi/R/plot_theme.R")

y_prop <- h5read("data/paper/binding_approx.jld2", "y_prop")
x_c1 <- h5read("data/paper/binding_approx.jld2", "x_c1")
x_c2 <- h5read("data/paper/binding_approx.jld2", "x_c2")

label_c2 <- function(c2) {
  str_c("*c<sub>y'</sub>* = 10<sup>", log10(c2), "</sup>") %>%
    factor(levels = str_c("*c<sub>y'</sub>* = 10<sup>", seq(-5,5, by = 0.5), "</sup>"))
}

plot_data <- y_prop %>%
  reshape2::melt(varnames = c("c1", "c2"), value.name = "p_bound") %>%
  mutate(c1 = x_c1[c1],
         c2 = x_c2[c2],
         c2_label = label_c2(c2)) %>%
  
  filter(c2 %in% x_c2[seq(1,12, by = 2)])

plot_data_c2 <- tibble(
  c2 = x_c2
) %>%
  mutate(
    c2_label = label_c2(c2)) %>%
  
  filter(c2 %in% x_c2[seq(1,12, by = 2)])

plot_data_approx <- tibble(
  c1 = 10^seq(-3,3, by = 0.1)
) %>%
  mutate(p_bound = c1 / (1 + c1))

ggplot() +
  geom_line(aes(x = c1, y = p_bound),
            linewidth = 1.0,
            plot_data) +
  geom_line(aes(x = c1, y = p_bound),
            linewidth = 0.7, colour = colour_C,
            plot_data_approx) +
  
  geom_rug(aes(x = c2), plot_data_c2, 
           length = unit(0.15, "npc")) +
  
  xlab("*c<sub>x</sub>*") +
  ylab("*p*<sub>bound</sub>") +
  
  scale_x_log10(labels = scales::label_log()) +
  
  facet_wrap(~c2_label, ncol = 3, scales = "free_x", dir = "h") +
  
  # ggokabeito::scale_colour_okabe_ito(name = "Antigen concentration") +
  # scale_color_brewer(palette = "OrRd") +
  
  plot_theme_paper +
  theme(strip.text = element_markdown(),
        panel.grid.major = element_gridline)
