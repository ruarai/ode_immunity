

library(tidyverse)
library(rhdf5)


source("R/plot_theme.R")

sus <- h5read("data/paper/basic_boosting.jld2", "sus")
inf <- h5read("data/paper/basic_boosting.jld2", "inf")
seq_t <- h5read("data/paper/basic_boosting.jld2", "seq_t")
c_levels <- h5read("data/paper/basic_boosting.jld2", "c_levels")
p_acq <- h5read("data/paper/basic_boosting.jld2", "p_acq")

plot_data_sus <- sus %>%
  reshape2::melt(varnames = c("t", "strata"), value.name = "prevalence") %>%
  mutate(c = c_levels[strata], p = p_acq[strata], t = seq_t[t],
         strata = strata - 1) %>% 
  filter(t < 1500)

plot_data_inf <- tibble(t = seq_t, prevalence = inf)


plot_data_means <- plot_data_sus %>%
  group_by(t) %>%
  mutate(prevalence = prevalence / sum(prevalence)) %>% 
  summarise(c = sum(prevalence * c),
            p = sum(prevalence * p), .groups = "drop")

p_summ <- ggplot() +
  geom_line(aes(x = t, y = prevalence),
            plot_data_inf,
            linewidth = 0.7) +
  
  coord_cartesian(xlim = c(0, 1500)) +
  
  scale_x_continuous(breaks = scales::breaks_extended(),
                     labels = scales::label_comma()) +
  
  xlab(NULL) + ylab("Infection<br>prevalence") +
  
  plot_theme_paper +
  theme(panel.grid.major = element_gridline,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  
  ggtitle(NULL, "<b>A</b> — Infection prevalence")
p_summ

col_heatmap <- colorspace::sequential_hcl(n = 128, h = c(262, 136), c = c(39, 72, 0), l = c(13, 98), power = c(1.65, 1.1))

rescale_density <- function(x, range) {
  plogis(0.6 * qlogis(x) + 2)
}

p_heatmap <- plot_data_sus %>% 
  
  summarise(prevalence = sum(prevalence), .by = c("c", "t")) %>% 
  mutate(density = rescale_density(prevalence)) %>%
  
  ggplot() +
  geom_tile(aes(x = t, y = c, fill = density)) +
  
  geom_hline(aes(yintercept = 10^6),
             colour = "white", linetype = "44") +
  
  scale_x_continuous(breaks = scales::breaks_extended(),
                     labels = scales::label_comma()) +
  
  scale_y_continuous(trans = "log10", breaks = 10^c(0, 2, 4, 6, 8),
                     labels = scales::label_log(base = 10),
                     sec.axis = sec_axis(
                       transform = "identity",
                       labels = function(x) log10(x) * (32 / 8),
                       breaks = 10^c(0, 2, 4, 6, 8),
                       name = "Strata *i*"
                     )) +

  scale_fill_viridis_c(option = 8, name = "Proportion",
                       breaks = rescale_density(c(1e-8, 0.01, 0.05, 0.2, 1 - 1e-3)),
                       labels = c(0, 0.01, 0.05, 0.2, 1.0))  +
  
  # scale_fill_continuous_divergingx(palette = "PuOr", mid = rescale_density(0.1),
  #                                  breaks = rescale_density(c(1e-8, 0.01, 0.05, 0.2, 1 - 1e-3)),
  #                                  labels = c(0, 0.01, 0.05, 0.2, 1.0)) +

  coord_cartesian(xlim = c(0, 1500), ylim = c(10^0, 10^8)) +
  
  # guides(fill = guide_colourbar(barwidth = 15)) +
  
  xlab(NULL) + ylab("Antibody<br>concentration") +
  
  plot_theme_paper + 
  theme(legend.position = "none",
        legend.title = element_text(margin = margin(b = 0.6, unit = "cm")),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  
  ggtitle(NULL, "<b>B</b> — Population distribution of antibodies")

p_legend <- tibble(
  x = seq(-0.1, 1.1, by = 0.001)
) %>%
  mutate(density = rescale_density(x)) %>% 
  fill(density, .direction = "updown") %>% 
  ggplot() +
  geom_tile(aes(x = 1, y = x, fill = density)) +
  
  geom_hline(aes(yintercept = y), tibble(y = seq(0, 1, by = 0.25)),
             colour = "white") +
  
  scale_fill_viridis_c(option = 8,
                       breaks = rescale_density(c(1e-8, 0.01, 0.05, 0.2, 1 - 1e-3)),
                       labels = c(0, 0.01, 0.05, 0.2, 1.0)) +
  
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  
  plot_theme_paper + xlab(NULL) + ylab("Proportion") +
  
  theme(legend.position = "none",
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.ticks.y = element_blank())


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
  
  xlab("Time *t* (days)") + ylab("Antibody<br>concentration") +
  
  plot_theme_paper +
  
  ggtitle(NULL, "<b>C</b> — Population mean antibody concentration") +
  
  theme(panel.grid.major = element_gridline)

p_mean_antibody


p_legend_space <- (plot_spacer() / p_legend / plot_spacer()) +
  plot_layout(heights = c(1, 5, 1))

p_heatmap_with_legend <- (p_heatmap | p_legend_space) +
  plot_layout(widths = c(15, 1))

((p_summ / p_heatmap) | (plot_spacer() / p_legend_space)) +
  plot_layout(widths = c(15, 1))

(p_summ / p_heatmap_with_legend)



ggsave(
  "results/results_basic.png",
  device = png,
  width = 10, height = 7,
  bg = "white"
)


# 
# protect_colours <- colorspace::sequential_hcl(
#   n = 3, h = c(135, 100), c = c(35, 70, 5), l = c(25, 98), power = c(1, 1.5)
# )
# 
# p_protect <- plot_data %>%
#   filter(class == "S") %>%
#   group_by(t) %>%
#   summarise(p_protected = sum(prevalence * (p > 0.95)),
#             p_unprotected = sum(prevalence * (p < 0.05))) %>%
#   
#   ggplot() +
#   geom_ribbon(aes(x = t, ymin = 0, ymax = p_protected,
#                   fill = ">95%", colour = ">95%"),
#               linewidth = 0.3) +
#   geom_ribbon(aes(x = t, ymin = p_protected, ymax = 1 - p_unprotected,
#                   fill = "5-95%", colour = "5-95%"),
#               linewidth = 0.3) +
#   geom_ribbon(aes(x = t, ymin = 1 - p_unprotected, ymax = 1,
#                   fill = "<5%", colour = "<5%"),
#               linewidth = 0.3)  +
#   
#   scale_x_continuous(breaks = scales::breaks_extended(),
#                      labels = scales::label_comma()) +
#   
#   scale_fill_manual(values = protect_colours, breaks = c(">95%", "5-95%", "<5%"), name = "Protection<br>against<br>infection <i>ω</i>") +
#   scale_colour_manual(values = protect_colours, breaks = c(">95%", "5-95%", "<5%"), name = "Protection<br>against<br>infection <i>ω</i>") +
#   
#   xlab("Time *t* (days)") + ylab("Proportion<br>of population") +
#   
#   coord_cartesian(xlim = c(0, 1500), ylim = c(0, 1.0)) +
#   
#   plot_theme_paper +
#   theme(panel.grid.major = element_gridline,
#         legend.title = element_markdown()) +
#   
#   ggtitle(NULL, "Population protection")


