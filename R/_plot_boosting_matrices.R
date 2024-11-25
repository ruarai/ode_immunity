


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


p_axes_x <- ggplot() +
  geom_blank(aes(x = 0), plot_data_matrices) +
  
  facet_wrap(~scenario, ncol = 3) +
  
  scale_x_continuous(breaks = seq(0, 32, by = 8),
                     labels = scales::label_log(base = 2)(2 ^ (8 * seq(0, 32, by = 8) / 32))) +
  
  xlab("Antibody concentration <i>c<sub>i</sub></i>") +
  
  coord_cartesian(xlim = c(0, 33)) +
  
  plot_theme_paper +
  theme(strip.text = element_markdown(size = 0))


p_axes_y <- ggplot() +
  geom_blank(aes(y = 0), plot_data_matrices) +
  
  scale_y_continuous(breaks = seq(0, 32, by = 8),
                     labels = scales::label_log(base = 2)(2 ^ (8 * seq(0, 32, by = 8) / 32))) +
  
  ylab("Antibody concentration <i>c<sub>j</sub></i>") +
  
  coord_cartesian(ylim = c(0, 33)) +
  
  plot_theme_paper +
  theme(strip.text = element_markdown(size = 0))


p_heatmap <- ggplot() +
  geom_tile(aes(x = j - 0.5, y = i - 0.5, fill = p),
            plot_data_matrices) +
  
  geom_abline(intercept = 0, slope = 1,
              colour = "white", linewidth = 0.7, linetype = "44") +
  
  facet_wrap(~scenario, ncol = 3) +
  
  scale_x_continuous(breaks = seq(0, 32, by = 8)) +
  scale_y_continuous(breaks = seq(0, 32, by = 8)) +
  
  coord_fixed(xlim = c(0, 33), ylim = c(0, 33)) +
  
  xlab("*I~i~*") + ylab("*S~j~*") +
  
  scale_fill_viridis_c(option = "A", name = "P(*I~i~* → *S~j~*)",
                       breaks = c(0.001, 0.1, 0.2), labels = c("0.0", "0.1", "≥0.2")) +
  
  plot_theme_paper +
  theme(legend.title = element_markdown(margin = margin(r = 0.0, b = 0.3, unit = "cm")),
        legend.position = "right")

# (
#   ((p_axes_y | p_heatmap) + plot_layout(widths = c(1, 40))) /
# ((plot_spacer() | p_axes_x) + plot_layout(widths = c(1, 40))) 
# ) +
#   plot_layout(heights = c(20, 1))
# 
# 
# 
# cowplot::plot_grid(
#   p_axes_y, p_heatmap, plot_spacer(), p_axes_x, ncol = 2,
#   rel_widths = c(1, 40),
#   rel_heights = c(20, 1),
#   align
# )

p_heatmap

ggsave(
  "results/results_boosting_matrices.png",
  device = png,
  width = 12, height = 4.5,
  bg = "white"
)

tibble(
  c = 2^seq(0, 8, by = 0.01)
) %>%
  # mutate(c_add = pmin(2^8, pmax(c, c + 2^4)),
  #        c_mult = pmin(2^8, pmax(c, c * 2^4)),
  #        c_null = pmin(2^8, pmax(c, 2^4))) %>% 
  
  expand_grid(scenario = 1:4) %>%
  
  mutate(
    c_post = case_when(
      scenario == 1 ~ 2^4, 
      scenario == 2 ~ c + 2^4,
      scenario == 3 ~ c * 2^4,
      scenario == 4 ~ c + 2^4 + 0.001 * c * pmax(1, 2^8 - c)
    ),
    c_post = pmin(2^8, pmax(c, c_post))
  ) %>% 
  
  # mutate(scenario = factor(scenario, 1:3, labels = c("No boosting", "Additive boosting", "Multiplicative boosting"))) %>% 
  
  ggplot() +
  
  
  annotate("linerange", y = 2^8, xmin = 2^0, xmax = 2^8,
           colour = "grey60", linetype = "44") +
  
  annotate("segment", x = 2^0, y = 2^0, xend = 2^8, yend = 2^8,
           colour = "grey60", linetype = "44") +
  
  geom_line(aes(x = c, y = c_post)) +
  
  facet_wrap(~scenario, ncol = 3) +
  
  scale_x_continuous(trans = "log2", labels = scales::label_log(base = 2)) +
  scale_y_continuous(trans = "log2", labels = scales::label_log(base = 2)) +
  
  xlab("Pre-infection concentration *c~i~*") + ylab("Post-infection<br>concentration *c~j~*") +
  
  coord_fixed(xlim = c(2^0, 2^8),
              ylim = c(2^0, 2^8)) +
  
  plot_theme_paper +
  
  theme(panel.grid = element_gridline)

