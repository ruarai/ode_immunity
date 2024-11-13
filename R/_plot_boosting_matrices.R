


library(tidyverse)
library(rhdf5)
library(patchwork)


source("../ode_immunity_multi/R/plot_theme.R")

boosting_matrices <- h5read("data/paper/basic_boosting.jld2", "boosting_matrices")

plot_data_matrices <- boosting_matrices %>%
  reshape2::melt(varnames = c("scenario", "i", "j"), value.name = "p") %>%
  mutate(
    scenario = c("none", "linear", "loglinear")[scenario],
    scenario = factor(scenario, c("loglinear", "linear", "none"), labels = c("Multiplicative boosting", "Additive boosting", "No boosting"))
  ) %>%
  mutate(scenario = fct_rev(scenario)) %>%
  mutate(p = pmin(p, 0.2))


ggplot() +
  geom_tile(aes(x = j - 0.5, y = i - 0.5, fill = p),
            plot_data_matrices) +
  
  geom_abline(intercept = 0, slope = 1,
              colour = "white", linewidth = 0.7, linetype = "44") +
  
  facet_wrap(~scenario, ncol = 3) +
  
  coord_fixed() +
  
  xlab("*I~i~*") + ylab("*S~j~*") +
  
  scale_fill_viridis_c(option = "A", name = "P(*I~i~* → *S~j~*)",
                       breaks = c(0.001, 0.1, 0.2), labels = c("0.0", "0.1", "≥0.2")) +
  
  plot_theme_paper +
  theme(legend.title = element_markdown(margin = margin(r = 0.0, b = 0.3, unit = "cm")),
        legend.position = "right")


ggsave(
  "results/results_boosting_matrices.png",
  device = png,
  width = 12, height = 4.5,
  bg = "white"
)


tibble(
  t = 1:1000
) %>%
  mutate(inf = t %in% c(50, 500, 750),
         c = 0) %>%
  
  mutate(c = c * lag(c, default ))






