


library(tidyverse)
library(rhdf5)


x_lambda <- h5read("data/bifurcations.jld2", "x_lambda")
y_I_sol <- h5read("data/bifurcations.jld2", "y_I_sol")
y_fixed_I <- h5read("data/bifurcations.jld2", "y_fixed_I")



plot(y_I_sol[4,])

maximums <- apply(y_I_sol[,20000:30000], 1, FUN = max)



plot_data_bifurcation <- tibble(
  lambda = x_lambda,
  fixed_I = y_fixed_I,
  maximum = maximums
) %>%
  filter(lambda > 0.006)


colour_A <- "#7BB5CA"
colour_B <- "#CC9D57"

ggplot() +
  
  annotate("segment", x = 0.03, y = -0.05, xend = 0.03, yend = 0.25,
           linewidth = 1.0, colour = colour_B, linetype = "24") +
  
  annotate("segment", x = 0.35, y = -0.05, xend = 0.35, yend = 0.25,
           linewidth = 1.0, colour = colour_B, linetype = "24") +
  
  geom_vline(xintercept = 0.35,
             linewidth = 1.0, 
             colour = colour_B, linetype = "24") +
  
  geom_line(aes(x = lambda, y = maximum),
            linewidth = 1.0,
            colour = colour_A,
            plot_data_bifurcation) +
  geom_line(aes(x = lambda, y = fixed_I),
            linewidth = 1.0,
            plot_data_bifurcation %>% filter(lambda < 0.02)) +
  geom_line(aes(x = lambda, y = fixed_I),
            linewidth = 1.0,
            plot_data_bifurcation %>% filter(lambda > 0.3)) +
  geom_line(aes(x = lambda, y = fixed_I),
            linewidth = 1.0, linetype = "11",
            plot_data_bifurcation %>% filter(lambda <= 0.31)) +
  
  annotate("point", x = 0.3, y = 0.075, size = 4) +
  annotate("point", x = 0.3, y = 0.075, size = 2, colour = "white") +
  
  annotate("point", x = 0.02, y = 0.006, size = 4) +
  annotate("point", x = 0.02, y = 0.006, size = 2, colour = "white") +
  
  annotate("label", hjust = 1, x = 0.15, y = 0.12, label = "Unstable\nfixed\npoint", size = 7, lineheight = 0.8, label.size = 0) +
  annotate("label", hjust = 0, x = 0.3, y = 0.15, label = "Stable fixed\npoint", size = 7, lineheight = 0.8, label.size = 0) +
  annotate("label", hjust = 0, x = 0.24, y = 0.22, label = "Limit cycle", size = 7, lineheight = 0.8, label.size = 0) +
  
  annotate("label", hjust = 0.5, x = 0.03, y = -0.05, label = "λ = 0.03", size = 7, lineheight = 0.8, label.size = 0.8, label.r = unit(0, "cm")) +
  annotate("label", hjust = 0.5, x = 0.35, y = -0.05, label = "λ = 0.35", size = 7, lineheight = 0.8, label.size = 0.8, label.r = unit(0, "cm")) +
  
  annotate("segment", x = 0.37, y = 0.14, xend = 0.39, yend = 0.105, linewidth = 0.7) +
  annotate("segment", x = 0.3, y = 0.2, xend = 0.27, yend = 0.17, linewidth = 0.7) +
  annotate("segment", x = 0.15, y = 0.09, xend = 0.2, yend = 0.06, linewidth = 0.7) +
  
  xlab("Decay rate λ") +
  ylab("Infected _I_") +
  
  coord_cartesian(xlim = c(0, 0.5),
                  ylim = c(0, 0.25),
                  expand = FALSE,
                  clip = "off") +
  
  plot_theme +
  
  theme(axis.title.x = element_markdown(margin = margin(b = 0, t = 1.3, unit = "cm")))



library(ggtext)

cowplot::plot_grid(
  
  tibble(I = y_I_sol[which(x_lambda == 0.03),]) %>% 
    mutate(t = row_number()) %>%
    filter(t < 1200) %>%
    
    ggplot() +
    
    geom_line(aes(x = t, y = I),
              linewidth = 0.7,
              colour = colour_B) +
    
    geom_hline(yintercept = y_fixed_I[which(x_lambda == 0.03)],
               linetype = "24", linewidth = 0.7) +
    
    xlab("Time _t_ (days)") + ylab("Infected _I_") +
    
    coord_cartesian(xlim = c(0, 1000),
                    ylim = c(0, 0.3)) +
    
    plot_theme,
    
  tibble(I = y_I_sol[which(x_lambda == 0.35),]) %>% 
    mutate(t = row_number()) %>%
    filter(t < 1200) %>%
    
    ggplot() +
    
    geom_line(aes(x = t, y = I),
              linewidth = 0.7,
              colour = colour_B) +
    
    geom_hline(yintercept = y_fixed_I[which(x_lambda == 0.35)],
               linetype = "24", linewidth = 0.7) +
    
    xlab("Time _t_ (days)") + ylab("Infected _I_") +
    
    coord_cartesian(xlim = c(0, 1000),
                    ylim = c(0, 0.3)) +
    
    plot_theme,
  
  ncol = 1, align = "v"


)



tibble(I = y_I_sol[which(x_lambda == 0.006),]) %>% mutate(t = row_number()) %>%
  
  ggplot() +
  
  geom_line(aes(x = t, y = I),
            linewidth = 0.7,
            colour = colour_B) +
  
  geom_hline(yintercept = y_fixed_I[which(x_lambda == 0.01)],
             linetype = "24", linewidth = 0.7) +
  
  xlab("Time _t_ (days)") + ylab("Infected _I_") +
  
  plot_theme




