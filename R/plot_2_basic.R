


library(tidyverse)
library(rhdf5)


source("R/plot_theme.R")

sus <- h5read("data/paper/basic_boosting.jld2", "sus")
inf <- h5read("data/paper/basic_boosting.jld2", "inf")
seq_t <- h5read("data/paper/basic_boosting.jld2", "seq_t")
c_levels <- h5read("data/paper/basic_boosting.jld2", "c_levels")
p_acq <- h5read("data/paper/basic_boosting.jld2", "p_acq")

t_max <- 1000

blues <- colorspace::sequential_hcl(n = 33, h = c(265, 80), c = c(60, NA, 10), l = c(25, 95), power = c(0.7, 2))


plot_data_sus <- sus %>%
  reshape2::melt(varnames = c("t", "strata"), value.name = "prevalence") %>%
  mutate(c = c_levels[strata], p = p_acq[strata], t = seq_t[t],
         strata = strata - 1) %>% 
  filter(t < t_max)

plot_data_inf <- tibble(t = seq_t, prevalence = inf) %>%
  filter(t < t_max)


plot_data_means <- plot_data_sus %>%
  filter(t < t_max) %>% 
  group_by(t) %>%
  mutate(prevalence = prevalence / sum(prevalence)) %>% 
  summarise(c = sum(prevalence * c),
            p = sum(prevalence * p), .groups = "drop")

p_inf <- ggplot() +
  
  # geom_ribbon(aes(x = t, ymin = 0, ymax = prevalence),
  #             plot_data_inf,
  #             linewidth = 0.4, fill = colour_C, colour = "black") +
  
  geom_line(aes(x = t, y = prevalence),
            plot_data_inf, linewidth = 0.7) +
  
  coord_cartesian(xlim = c(-30, t_max), ylim = c(0, 0.075)) +
  
  scale_x_continuous(breaks = scales::breaks_extended(),
                     labels = scales::label_comma()) +
  
  xlab(NULL) + ylab("Prevalence") +
  
  plot_theme_paper +
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_gridline) +
  
  ggtitle(NULL, "<b>A</b> — Infection prevalence")

p_0 <- ggplot() +
  # geom_ribbon(aes(x = t, ymin = 0, ymax = prevalence),
  #             plot_data_sus %>% filter(strata == 0),
  #             fill = blues[33],
  #             outline.type = "upper",
  #             colour = "black", linewidth = 0.4) +
  
  geom_line(aes(x = t, y = prevalence),
            linewidth = 0.6,
            plot_data_sus %>% filter(strata == 0)) +
  
  # 
  # annotate("richtext", x = -50, y = 1.0, label = "<i>S</i><sub>0</sub>",
  #          label.r = unit(0, "cm"), label.size = 0, label.colour = "white") +
  # 
  # annotate("linerange", xmin = -10, xmax = 0, y = 1.0, linewidth = 0.4) +
  
  
  xlab(NULL) + ylab("Prevalence") +
  
  plot_theme_paper +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid.major = element_gridline) +
  
  coord_cartesian(xlim = c(-50, NA), ylim = c(0, 1.05)) +
  
  ggtitle(NULL, "<b>B</b> — Susceptible prevalence <i>S</i><sub>0</sub>")






wave_scale <- 0.01
plot_data_waves <- plot_data_sus %>%
  filter(strata > 0, t < t_max) %>%
  mutate(strata_group = factor(strata),
         offset = if_else(strata == 0, 0.0, -strata * wave_scale))

plot_data_labels <- plot_data_waves %>% filter(t == 1) %>%
  mutate(offset = if_else(strata == 0, offset + 1, offset)) %>%
  filter(strata == 1 | strata %% 8 == 0)

p_waves <- ggplot() +
  
  geom_richtext(aes(x = -50, y = offset, label = str_c("<i>S</i><sub>", strata,"</sub>")),
                label.r = unit(0, "cm"), label.size = 0, label.colour = "white",
                plot_data_labels) +
  geom_ribbon(aes(x = t, ymin = offset,
                  ymax = prevalence + offset,
                  fill = strata_group,
                  group = strata_group),
              plot_data_waves, linewidth = 0.4,
              outline.type = "upper",
              colour = "black") +
  
  geom_linerange(aes(xmin = -10, xmax = 0, y = offset), 
                 linewidth = 0.4, plot_data_labels) +
  
  xlab("Time (days)") + ylab("Prevalence") + 
  
  scale_fill_manual(values = rev(blues[2:33]) ) +
  
  scale_y_continuous(breaks = -wave_scale * 21 + c(0, 0.1), labels = c(0, 0.1)) +
  # scale_y_continuous(breaks = NULL) +
  
  coord_cartesian(ylim = c(-32 * wave_scale, 0.05),
                  xlim = c(-50, NA)) +
  
  plot_theme_paper +
  
  theme(legend.position = "none")+
  
  ggtitle(NULL, "<b>C</b> — Susceptible prevalence <i>S</i><sub>1</sub> – <i>S</i><sub><i>k</i></sub>")

p_mean_antibody <- ggplot() +
  geom_line(aes(x = t, y = c),
            plot_data_means,
            linewidth = 0.7) +
  
  scale_y_continuous(trans = "log10",
                     breaks = 10^c(0, 2, 4, 6, 8),
                     labels = scales::label_log(base = 10))  +
  
  scale_x_continuous(breaks = scales::breaks_extended(),
                     labels = scales::label_comma()) +
  
  coord_cartesian(ylim = c(10^0, 10^7)) +
  
  xlab("Time *t* (days)") + ylab("Concentration") +
  
  plot_theme_paper +
  
  ggtitle(NULL, "<b>D</b> — Mean antibody concentration") +
  
  theme(panel.grid.major = element_gridline)


(p_inf / p_waves) | (p_0 / p_mean_antibody)


ggsave(
  "results/results_basic_template.pdf",
  device = pdf,
  width = 13, height = 7,
  bg = "white"
)
 
