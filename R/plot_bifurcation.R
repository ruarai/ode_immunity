


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

  geom_vline(xintercept = 0.006,
             colour = "grey80", linewidth = 1.0, alpha = 0.5) +
  
  geom_line(aes(x = lambda, y = maximum),
            linewidth = 1.0,
            colour = colour_A,
            plot_data_bifurcation %>% filter(lambda > 0.0002)) +
  geom_line(aes(x = lambda, y = fixed_I),
            linewidth = 1.0,
            plot_data_bifurcation %>% filter(lambda > 0.0071)) +
  geom_line(aes(x = lambda, y = fixed_I),
            linewidth = 1.0, linetype = "11",
            plot_data_bifurcation %>% filter(lambda <= 0.033)) +
  
  annotate("point", x = 0.008, y = 0.0167, size = 3, colour = "white") +
  annotate("point", x = 0.008, y = 0.0167, size = 2, colour = colour_C) +
  
  annotate("point", x = 0.006, y = 0.0127, size = 3, colour = "white") +
  annotate("point", x = 0.006, y = 0.0127, size = 2, colour = colour_C) +
  
  annotate("point", x = 0.0071, y = 0.015, size = 4) +
  annotate("point", x = 0.0071, y = 0.015, size = 2, colour = "white") +
  
  xlab("Decay rate λ") +
  ylab("Infected I") +
  
  coord_cartesian(xlim = c(0, 0.015),
                  ylim = c(0, 0.07),
                  expand = FALSE) +
  
  plot_theme_paper +
  
  ggtitle("<b>A</b>")



p_bifurcation




p_examples <- cowplot::plot_grid(
  tibble(I = y_I_sol[which(x_lambda == 0.006),]) %>% 
    mutate(t = row_number()) %>%
    
    ggplot() +
    
    geom_line(aes(x = t, y = I),
              linewidth = 0.7) +
    
    geom_hline(yintercept = y_fixed_I[which(x_lambda == 0.006)],
               linewidth = 0.7, colour = colour_C) +
    
    scale_x_continuous(breaks = seq(0, 365 * 12, by = 365),
                       labels = 0:12) +
    
    xlab("Time (years)") + ylab("Infected I") +
    
    coord_cartesian(xlim = c(0, 3400),
                    ylim = c(0, 0.1)) +
    
    plot_theme_paper +
    
    ggtitle("<b>B</b>", "λ = 0.006"),
  
  tibble(I = y_I_sol[which(x_lambda == 0.008),]) %>% 
    mutate(t = row_number()) %>%
    
    ggplot() +
    
    geom_line(aes(x = t, y = I),
              linewidth = 0.7) +
    
    geom_hline(yintercept = maximums[which(x_lambda == 0.008)],
               linewidth = 0.7, colour = colour_C) +
    
    scale_x_continuous(breaks = seq(0, 365 * 12, by = 365),
                       labels = 0:12) +
    
    xlab("Time (years)") + ylab("Infected I") +
    
    coord_cartesian(xlim = c(0, 3400),
                    ylim = c(0, 0.1)) +
    
    plot_theme_paper +
    
    ggtitle("<b>C</b>", "λ = 0.008"),
  
  ncol = 2, align = "h", axis = "tb"
)


cowplot::plot_grid(
  p_bifurcation + theme(plot.margin = margin(t = 0.3, b = 0.3, l = 5, r = 5, unit = "cm")),
  p_examples,
  ncol = 1
)

ggsave(
  "results/results_bifurcation.png",
  device = png,
  scale = 10 / 16,
  dpi = 300,
  width = 16, height = 10,
  bg = "white"
)


