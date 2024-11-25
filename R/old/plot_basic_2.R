

library(tidyverse)
library(rhdf5)


source("R/plot_theme.R")



sol_S <- h5read("data/paper/basic.jld2", "sol_S")
sol_I <- h5read("data/paper/basic.jld2", "sol_I")


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
  
  scale_x_continuous(breaks = scales::breaks_extended(),
                     labels = scales::label_comma()) +
  
  xlab("Time (days)") + ylab("Infection prevalence") +
  
  plot_theme_paper +
  
  ggtitle("<b>A</b>")

p_heatmap <- sol %>%
  
  mutate(S = pmin(S, 0.1)) %>% 
  
  ggplot() +
  geom_tile(aes(x = t, y = c * 16, fill = S, colour = S)) +
  
  scale_fill_viridis_c(option = "B", name = "Proportion",
                       breaks = c(0.001, 0.05, 0.1), labels = c("0.00", "0.05", "≥0.10"))  +
  
  scale_colour_viridis_c(option = "B", name = "Proportion",
                         breaks = c(0.001, 0.05, 0.1), labels = c("0.00", "0.05", "≥0.10")) +
  
  scale_x_continuous(breaks = scales::breaks_extended(),
                     labels = scales::label_comma()) +
  
  coord_cartesian(xlim = c(0, 2000), ylim = c(0, 16)) +
  
  xlab("Time (days)") + ylab("Strata i") +
  
  plot_theme_paper + 
  theme(legend.position = "bottom",
        legend.title = element_text(margin = margin(b = 15, r = 20))) +
  
  ggtitle("<b>B</b>")

cowplot::plot_grid(p_summ, p_heatmap, ncol = 1, rel_heights = c(1, 1.3),
                   align = "v", axis = "lr")


ggsave(
  "results/results_singlestrain_basic.pdf",
  device = cairo_pdf,
  width = 10, height = 8,
  bg = "white"
)


