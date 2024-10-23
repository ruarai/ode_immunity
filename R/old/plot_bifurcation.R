


library(tidyverse)
library(rhdf5)


x_lambda <- h5read("data/paper/bifurcations.jld2", "x_lambda")
y_I_sol <- h5read("data/paper/bifurcations.jld2", "y_I_sol")
y_fixed_I <- h5read("data/paper/bifurcations.jld2", "y_fixed_I")



maximums <- apply(y_I_sol[,28000:32000], 1, FUN = max)




plot_data_bifurcation <- tibble(
  lambda = x_lambda,
  fixed_I = y_fixed_I,
  maximum = maximums
) 

p_bifurcation <- ggplot() +
  
  # annotate("rug", x = 0.008,
  #          colour = "grey80", linewidth = 1.0) +
  # 
  # annotate("rug", x = 0.006,
  #          colour = "grey80", linewidth = 1.0) +
  # 
  geom_vline(xintercept = 0.008,
             colour = "grey80", linewidth = 1.0, alpha = 0.5) +

  geom_vline(xintercept = 0.005,
             colour = "grey80", linewidth = 1.0, alpha = 0.5) +
  
  geom_line(aes(x = lambda, y = maximum),
            linewidth = 1.0,
            colour = colour_A,
            plot_data_bifurcation %>% filter(lambda > 0.0002)) +
  geom_line(aes(x = lambda, y = fixed_I),
            linewidth = 1.0,
            plot_data_bifurcation %>% filter(lambda > 0.0065)) +
  geom_line(aes(x = lambda, y = fixed_I),
            linewidth = 1.0, linetype = "11",
            plot_data_bifurcation %>% filter(lambda <= 0.033)) +
  
  # annotate("point", x = 0.008, y = 0.0167, size = 3, colour = "white") +
  # annotate("point", x = 0.008, y = 0.0167, size = 2, colour = colour_C) +
  # 
  # annotate("point", x = 0.006, y = 0.0127, size = 3, colour = "white") +
  # annotate("point", x = 0.006, y = 0.0127, size = 2, colour = colour_C) +
  
  annotate("point", x = 0.0065, y = 0.015, size = 4) +
  annotate("point", x = 0.0065, y = 0.015, size = 2, colour = "white") +
  
  xlab("Waning constant <i> ρ</i>") +
  ylab("Infection prevalence") +
  
  coord_cartesian(xlim = c(0, 0.015),
                  ylim = c(0, 0.07),
                  expand = FALSE) +
  
  plot_theme_paper +
  
  ggtitle("<b>A</b>")



p_examples <- cowplot::plot_grid(
  tibble(I = y_I_sol[which(x_lambda == 0.005),]) %>% 
    mutate(t = row_number()) %>%
    
    ggplot() +
    
    geom_line(aes(x = t, y = I),
              linewidth = 0.7) +
    
    geom_hline(yintercept = y_fixed_I[which(x_lambda == 0.005)],
               linewidth = 0.7, colour = colour_C) +
    
    scale_x_continuous(breaks = scales::breaks_extended(),
                       labels = scales::label_comma())  +
    
    scale_y_continuous(breaks = c(0.0, 0.05, 0.1)) +
    
    xlab("Time (days)") + ylab("Infection prevalence") +
    
    coord_cartesian(xlim = c(0, 3400),
                    ylim = c(0, 0.1)) +
    
    plot_theme_paper +
    
    ggtitle("<b>B</b>", "<i>ρ </i>  = 0.005") +
    theme(plot.subtitle = element_markdown()),
  
  tibble(I = y_I_sol[which(x_lambda == 0.008),]) %>% 
    mutate(t = row_number()) %>%
    
    ggplot() +
    
    geom_line(aes(x = t, y = I),
              linewidth = 0.7) +
    
    geom_hline(yintercept = maximums[which(x_lambda == 0.008)],
               linewidth = 0.7, colour = colour_C) +
    
    scale_x_continuous(breaks = scales::breaks_extended(),
                       labels = scales::label_comma()) +
    
    scale_y_continuous(breaks = c(0.0, 0.05, 0.1)) +
    
    xlab("Time (days)") + ylab("Infection prevalence") +
    
    coord_cartesian(xlim = c(0, 3400),
                    ylim = c(0, 0.1)) +
    
    plot_theme_paper +
    
    ggtitle("<b></b>", "<i>ρ </i> = 0.008") +
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


