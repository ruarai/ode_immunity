


library(tidyverse)
library(rhdf5)


x_lambda <- h5read("data/anziam2024/bifurcations.jld2", "x_lambda")
y_I_sol <- h5read("data/anziam2024/bifurcations.jld2", "y_I_sol")
y_fixed_I <- h5read("data/anziam2024/bifurcations.jld2", "y_fixed_I")



plot(y_I_sol[4,])

maximums <- apply(y_I_sol[,28000:32000], 1, FUN = max)



plot_data_bifurcation <- tibble(
  lambda = x_lambda,
  fixed_I = y_fixed_I,
  maximum = maximums
) 



ggplot() +
  
  geom_vline(xintercept = 0.02,
             colour = colour_B, linetype = "24", linewidth = 1.0) +
  
  geom_vline(xintercept = 0.04,
             colour = colour_B, linetype = "24", linewidth = 1.0) +

  geom_line(aes(x = lambda, y = maximum),
            linewidth = 1.0,
            colour = colour_A,
            plot_data_bifurcation %>% filter(lambda > 0.002)) +
  geom_line(aes(x = lambda, y = fixed_I),
            linewidth = 1.0,
            plot_data_bifurcation %>% filter(lambda > 0.032)) +
  geom_line(aes(x = lambda, y = fixed_I),
            linewidth = 1.0, linetype = "11",
            plot_data_bifurcation %>% filter(lambda <= 0.033)) +
  
  annotate("point", x = 0.0325, y = 0.0115, size = 4) +
  annotate("point", x = 0.0325, y = 0.0115, size = 2, colour = "white") +
  
  annotate("label", hjust = 1, x = 0.01, y = 0.02, label = "Unstable\nfixed\npoint", size = 7, lineheight = 0.8, label.size = 0) +
  annotate("label", hjust = 0, x = 0.031, y = 0.025, label = "Stable fixed\npoint", size = 7, lineheight = 0.8, label.size = 0) +
  annotate("label", hjust = 0, x = 0.021, y = 0.04, label = "Limit cycle", size = 7, lineheight = 0.8, label.size = 0) +

  annotate("segment", x = 0.01, y = 0.017, xend = 0.017, yend = 0.008, linewidth = 0.7) +
  annotate("segment", x = 0.025, y = 0.037, xend = 0.024, yend = 0.033, linewidth = 0.7) +
  annotate("segment", x = 0.036, y = 0.02, xend = 0.038, yend = 0.015, linewidth = 0.7) +
  
  xlab("Decay rate λ") +
  ylab("Infected I") +
  
  coord_cartesian(xlim = c(0, 0.05),
                  ylim = c(0, 0.05),
                  expand = FALSE) +
  
  plot_theme



library(ggtext)

cowplot::plot_grid(
  
  
  tibble(I = y_I_sol[which(x_lambda == 0.04),]) %>% 
    mutate(t = row_number()) %>%
    
    ggplot() +
    
    geom_line(aes(x = t, y = I),
              linewidth = 0.7,
              colour = colour_B) +
    
    geom_hline(yintercept = y_fixed_I[which(x_lambda == 0.04)],
               linetype = "24", linewidth = 0.7) +
    
    scale_x_continuous(breaks = seq(0, 365 * 12, by = 365),
                       labels = 0:12) +
    
    xlab("Time (years)") + ylab("Infected I") +
    
    coord_cartesian(xlim = c(0, 3400),
                    ylim = c(0, 0.1)) +
    
    plot_theme,
  
  tibble(I = y_I_sol[which(x_lambda == 0.02),]) %>% 
    mutate(t = row_number()) %>%
    
    ggplot() +
    
    geom_hline(yintercept = maximums[which(x_lambda == 0.02)],
               linewidth = 0.7, colour = colour_A) +
    
    geom_line(aes(x = t, y = I),
              linewidth = 0.7,
              colour = colour_B) +
    
    geom_hline(yintercept = y_fixed_I[which(x_lambda == 0.02)],
               linetype = "24", linewidth = 0.7) +
    
    scale_x_continuous(breaks = seq(0, 365 * 12, by = 365),
                       labels = 0:12) +
    
    xlab("Time (years)") + ylab("Infected I") +
    
    coord_cartesian(xlim = c(0, 3400),
                    ylim = c(0, 0.1)) +
    
    plot_theme,
  
  ncol = 1, align = "v"
)

tibble(I = y_I_sol[which(x_lambda == 0.02),]) %>% 
  mutate(t = row_number()) %>%
  
  ggplot() +
  
  geom_hline(yintercept = maximums[which(x_lambda == 0.02)],
             linewidth = 0.7, colour = colour_A) +
  
  geom_line(aes(x = t, y = I),
            linewidth = 0.7,
            colour = colour_B) +
  
  geom_hline(yintercept = y_fixed_I[which(x_lambda == 0.02)],
             linetype = "24", linewidth = 0.7) +
  
  scale_x_continuous(breaks = seq(0, 365 * 12, by = 365),
                     labels = 0:12) +
  
  xlab("Time (years)") + ylab("Infected I") +
  
  coord_cartesian(xlim = c(0, 3400),
                  ylim = c(0, 0.1)) +
  
  plot_theme



tibble(I = y_I_sol[which(x_lambda == 0.006),]) %>% mutate(t = row_number()) %>%
  
  ggplot() +
  
  geom_line(aes(x = t, y = I),
            linewidth = 0.7,
            colour = colour_B) +
  
  geom_hline(yintercept = y_fixed_I[which(x_lambda == 0.01)],
             linetype = "24", linewidth = 0.7) +
  
  xlab("Time _t_ (days)") + ylab("Infected _I_") +
  
  plot_theme




