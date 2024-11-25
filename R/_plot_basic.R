

library(tidyverse)
library(rhdf5)


source("../ode_immunity_multi/R/plot_theme.R")

sol_t <- h5read("data/paper/basic_boosting.jld2", "sol_t")
seq_t <- h5read("data/paper/basic_boosting.jld2", "seq_t")
c_levels <- h5read("data/paper/basic_boosting.jld2", "c_levels")
p_acq <- h5read("data/paper/basic_boosting.jld2", "p_acq")

fn_p_acq <- function(c, b = 10^3, h = 3) { (c ^ h) / (b ^ h + c ^ h)}

c_05 <- uniroot(rlang::as_function(~ fn_p_acq(.x) - 0.05), c(0, max(c_levels)))$root
c_95 <- uniroot(rlang::as_function(~ fn_p_acq(.x) - 0.95), c(0, max(c_levels)))$root

plot_data <- sol_t %>%
  reshape2::melt(varnames = c("scenario", "t", "class", "strata"), value.name = "prevalence") %>% 
  mutate(class = c("S", "I")[class], c = c_levels[strata], p = p_acq[strata], t = seq_t[t]) %>%
  filter(scenario == 1, t < 1500)


plot_data_summ_inf <- plot_data %>%
  group_by(t, class) %>%
  summarise(prevalence = sum(prevalence), .groups = "drop") %>%
  filter(class == "I")

plot_data_means <- plot_data %>%
  group_by(t) %>%
  mutate(prevalence = prevalence / sum(prevalence)) %>% 
  summarise(c = sum(prevalence * c),
            p = sum(prevalence * p), .groups = "drop")




p_summ <- plot_data_summ_inf %>% 
  ggplot() +
  geom_line(aes(x = t, y = prevalence),
            linewidth = 0.7) +
  
  coord_cartesian(xlim = c(0, 1500)) +
  
  scale_x_continuous(breaks = scales::breaks_extended(),
                     labels = scales::label_comma()) +
  
  xlab(NULL) + ylab("Infection<br>prevalence") +
  
  plot_theme_paper +
  theme(panel.grid.major = element_gridline,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  
  ggtitle(NULL, "Infection prevalence")

p_heatmap <- plot_data %>% 
  
  summarise(prevalence = sum(prevalence), .by = c("c", "t")) %>% 
  mutate(prevalence = pmin(prevalence, 0.05)) %>% 
  
  ggplot() +
  geom_tile(aes(x = t, y = c, fill = prevalence)) +
  
  scale_fill_viridis_c(option = "B", name = "Proportion",
                       breaks = c(0.001, 0.025, 0.05), labels = c("0.00", "0.025", "≥0.05"))  +

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
  
  coord_cartesian(xlim = c(0, 1500), ylim = c(10^0, 10^8)) +
  
  # guides(fill = guide_colourbar(barwidth = 15)) +
  
  xlab(NULL) + ylab("Antibody<br>concentration") +
  
  plot_theme_paper + 
  theme(legend.position = "right",
        legend.title = element_text(margin = margin(b = 0.6, unit = "cm")),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  
  ggtitle(NULL, "Population distribution of antibodies")


p_heatmap

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
  
  ggtitle(NULL, "Population mean antibody concentration") +
  
  theme(panel.grid.major = element_gridline)


(p_summ / p_heatmap / p_mean_antibody) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))



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


