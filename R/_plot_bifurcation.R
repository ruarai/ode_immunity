
sou

library(tidyverse)
library(rhdf5)

source("../ode_immunity_multi/R/plot_theme.R")

x_rho <- h5read("data/paper/bifurcations.jld2", "x_rho")
y_I_sol <- h5read("data/paper/bifurcations.jld2", "y_I_sol")
y_fixed_I <- h5read("data/paper/bifurcations.jld2", "y_fixed_I")



maximums <- apply(y_I_sol[,28000:32000], 1, FUN = max)
minimums <- apply(y_I_sol[,28000:32000], 1, FUN = min)


rho_A <- 0.003
rho_B <- 0.005


rho_to_halflife <- function(x) {1 / (8 * x)}


plot_data_bifurcation <- tibble(
  rho = x_rho,
  fixed_I = y_fixed_I,
  maximum = maximums,
  minimum = minimums
) %>%
  mutate(halflife = rho_to_halflife(rho))

bifur_point <- c(0.00393, 0.0147)

p_bifurcation <- ggplot() +
  geom_vline(xintercept = rho_B,
             colour = "grey80", linewidth = 1.0, alpha = 0.5) +

  geom_vline(xintercept = rho_A,
             colour = "grey80", linewidth = 1.0, alpha = 0.5) +
  
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
  
  xlab("Waning constant <i> ρ</i>") +
  ylab("Infection prevalence") +
  
  coord_cartesian(xlim = c(0, 0.0075),
                  ylim = c(0, 0.06),
                  expand = FALSE) +
  
  scale_x_continuous(sec.axis = sec_axis(breaks = 8/(2^seq(10,14, by = 1)), labels = rho_to_halflife, transform = identity, name = "Antibody half-life")) +
  
  
  plot_theme_paper +
  
  theme(axis.text.x.top = element_text(margin = margin(b = 0.3, unit = "cm"))) +
  
  ggtitle("A")

p_bifurcation



p_examples <- cowplot::plot_grid(
  tibble(I = y_I_sol[which(x_rho == rho_A),]) %>% 
    mutate(t = row_number()) %>%
    
    ggplot() +
    
    geom_line(aes(x = t, y = I),
              linewidth = 0.7) +
    
    geom_hline(yintercept = y_fixed_I[which(x_rho == rho_A)],
               linewidth = 0.7, colour = colour_C) +
    
    scale_x_continuous(breaks = scales::breaks_extended(),
                       labels = scales::label_comma())  +
    
    scale_y_continuous(breaks = c(0.0, 0.05, 0.1)) +
    
    xlab("Time (days)") + ylab("Infection prevalence") +
    
    coord_cartesian(xlim = c(0, 3400),
                    ylim = c(0, 0.08)) +
    
    plot_theme_paper +
    
    ggtitle("B", str_c("<i>ρ </i>  = ", rho_A)) +
    theme(plot.subtitle = element_markdown()),
  
  tibble(I = y_I_sol[which(x_rho == rho_B),]) %>% 
    mutate(t = row_number()) %>%
    
    ggplot() +
    
    geom_line(aes(x = t, y = I),
              linewidth = 0.7) +
    
    geom_hline(yintercept = maximums[which(x_rho == rho_B)],
               linewidth = 0.7, colour = colour_C) +
    
    scale_x_continuous(breaks = scales::breaks_extended(),
                       labels = scales::label_comma()) +
    
    scale_y_continuous(breaks = c(0.0, 0.05, 0.1)) +
    
    xlab("Time (days)") + ylab("Infection prevalence") +
    
    coord_cartesian(xlim = c(0, 3400),
                    ylim = c(0, 0.08)) +
    
    plot_theme_paper +
    
    ggtitle(" ", str_c("<i>ρ </i> = ", rho_B)) +
    theme(plot.subtitle = element_markdown()),
  
  ncol = 2, align = "h", axis = "tb"
)


cowplot::plot_grid(
  p_bifurcation,
  p_examples,
  ncol = 1
)

ggsave(
  "results/results_bifurcation.pdf",
  device = cairo_pdf,
  width = 10, height = 7,
  bg = "white"
)

ggplot() +
  geom_vline(xintercept = 0.008,
             colour = "grey80", linewidth = 1.0, alpha = 0.5) +
  
  geom_vline(xintercept = 0.005,
             colour = "grey80", linewidth = 1.0, alpha = 0.5) +
  
  geom_line(aes(x = rho, y = minimum, colour = "min"),
            linewidth = 1.0,
            colour = colour_C,
            plot_data_bifurcation %>% filter(rho > 0.0002)) +
  geom_line(aes(x = rho, y = maximum, colour = "max"),
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
  
  coord_cartesian(xlim = c(0, 0.005),
                  ylim = c(1e-8, 0.5),
                  expand = FALSE) +
  
  scale_y_log10(labels = scales::label_log()) +
  
  scale_x_continuous(sec.axis = sec_axis(breaks = 8/(2^seq(10,14, by = 1)), labels = rho_to_halflife, transform = identity, name = "Antibody half-life")) +
  
  plot_theme_paper +
  
  theme(axis.text.x.top = element_text(margin = margin(b = 0.3, unit = "cm")))



