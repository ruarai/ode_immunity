
sou

library(tidyverse)
library(rhdf5)
library(patchwork)

source("../ode_immunity_multi/R/plot_theme.R")

x_rho <- h5read("data/paper/bifurcations.jld2", "x_rho")
y_I_sol <- h5read("data/paper/bifurcations.jld2", "y_I_sol")
y_fixed_I <- h5read("data/paper/bifurcations.jld2", "y_fixed_I")



maximums <- apply(y_I_sol[,28000:32000], 1, FUN = max)
minimums <- apply(y_I_sol[,28000:32000], 1, FUN = min)


rhos <- c(0.001, 0.003, 0.005)
y_fixed_I_rhos <- y_fixed_I[!is.na(match(x_rho, rhos))]


rho_to_halflife <- function(x) {1 / (8 * x)}

data_I_sol <- y_I_sol %>%
  reshape2::melt(varnames = c("rho", "t"), value.name = "prev") %>% 
  mutate(rho = x_rho[rho])


plot_data_bifurcation <- tibble(
  rho = x_rho,
  fixed_I = y_fixed_I,
  maximum = maximums,
  minimum = minimums
) %>%
  mutate(halflife = rho_to_halflife(rho))

bifur_point <- c(0.00393, 0.0147)

p_bifurcation <- ggplot() +
  geom_vline(aes(xintercept = rho), tibble(rho = rhos), colour = "grey80", linewidth = 1.0, alpha = 0.5) +
  
  geom_line(aes(x = rho, y = maximum),
            linewidth = 1.0,
            colour = colour_A,
            plot_data_bifurcation %>% filter(rho > 0.0002)) +
  geom_line(aes(x = rho, y = fixed_I),
            linewidth = 1.0,
            plot_data_bifurcation %>% filter(rho > bifur_point[1])) +
  geom_line(aes(x = rho, y = fixed_I),
            linewidth = 1.0, linetype = "11",
            plot_data_bifurcation %>% filter(rho <= bifur_point[1])) +
  
  annotate("point", x = bifur_point[1], y = bifur_point[2], size = 3, colour = "black") +
  annotate("point", x = bifur_point[1], y = bifur_point[2], size = 1.5, colour = "white") +
  
  xlab("Waning constant <i>ρ</i>") +
  ylab("Infection prevalence") +
  
  coord_cartesian(xlim = c(0, 0.0075),
                  ylim = c(0, 0.06),
                  expand = FALSE) +
  
  # scale_x_continuous(
  #   sec.axis = sec_axis(
  #     breaks = 8/(2^seq(10,14, by = 1)), 
  #     labels = rho_to_halflife, 
  #     transform = identity, 
  #     name = "Antibody half-life (days)"
  #   )
  # ) +
  
  
  plot_theme_paper +
  theme(panel.grid.major.y = element_gridline) +
  
  theme(axis.text.x.top = element_text(margin = margin(b = 0.3, unit = "cm"))) +
  
  ggtitle("A", "Bifurcation over <i>ρ</i>") +
  theme(plot.subtitle = element_markdown())

p_bifurcation

p_examples <- data_I_sol %>%
  filter(rho %in% rhos, t < 4000) %>% 
  mutate(rho_label = str_c("<i>ρ </i>  = ", rho)) %>% 
  ggplot() +
  
  geom_line(aes(x = t, y = prev),
            linewidth = 0.7) +
  
  geom_hline(aes(yintercept = fixed_prev, linetype = stable), 
             tibble(rho = rhos, fixed_prev = y_fixed_I_rhos) %>% 
               mutate(stable = rho > bifur_point[1],
                      rho_label = str_c("<i>ρ </i>  = ", rho)),
             linewidth = 0.7, colour = "black") +
  
  scale_x_continuous(breaks = scales::breaks_extended(),
                     labels = scales::label_comma()) +
  
  scale_y_continuous(breaks = c(0.0, 0.05, 0.1)) +
  
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "22")) +
  
  facet_wrap(~rho_label, ncol = 1, scales = "free_x") +
  
  xlab("Time (days)") + ylab("Infection prevalence") +
  
  coord_cartesian(xlim = c(0, 3400),
                  ylim = c(0, 0.08)) +
  
  plot_theme_paper +
  
  ggtitle("C") +
  theme(strip.text = element_markdown(),
        legend.position = "none",
        panel.grid.major.x = element_gridline)


p_bifurcation_min <- ggplot() +
  geom_vline(aes(xintercept = rho), tibble(rho = rhos), colour = "grey80", linewidth = 1.0, alpha = 0.5) +

  geom_line(aes(x = rho, y = minimum, colour = "min"),
            linewidth = 1.0,
            colour = colour_C,
            plot_data_bifurcation %>% filter(rho > 0.0002)) +
  geom_line(aes(x = rho, y = fixed_I),
            linewidth = 1.0,
            plot_data_bifurcation %>% filter(rho > bifur_point[1])) +
  geom_line(aes(x = rho, y = fixed_I),
            linewidth = 1.0, linetype = "11",
            plot_data_bifurcation %>% filter(rho <= bifur_point[1])) +
  annotate("point", x = bifur_point[1], y = bifur_point[2], size = 3, colour = "black") +
  annotate("point", x = bifur_point[1], y = bifur_point[2], size = 1.5, colour = "white") +
  
  xlab("Waning constant <i>ρ</i>") +
  ylab("Infection prevalence") +
  
  coord_cartesian(xlim = c(0, 0.0075),
                  ylim = c(1e-9, 0.3),
                  expand = FALSE) +
  
  scale_y_log10(labels = scales::label_log()) +
  
  plot_theme_paper +
  theme(panel.grid.major.y = element_gridline) +
  
  ggtitle("B", "Minimum infection prevalence across solution") +
  theme(axis.text.x.top = element_text(margin = margin(b = 0.3, unit = "cm")))


free(p_bifurcation / p_bifurcation_min) | (p_examples)


ggsave(
  "results/results_bifurcation.pdf",
  device = cairo_pdf,
  width = 10, height = 7,
  bg = "white"
)
