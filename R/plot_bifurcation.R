


library(tidyverse)
library(rhdf5)


x_lambda <- h5read("data/bifurcations.jld2", "x_lambda")
y_I_sol <- h5read("data/bifurcations.jld2", "y_I_sol")
y_fixed_I <- h5read("data/bifurcations.jld2", "y_fixed_I")



plot(y_I_sol[4,])

maximums <- apply(y_I_sol[,10000:12000], 1, FUN = max)



plot_data_bifurcation <- tibble(
  lambda = x_lambda,
  fixed_I = y_fixed_I,
  maximum = maximums
)


colour_A <- "#7BB5CA"
colour_B <- "#CC9D57"

ggplot() +
  
  annotate("segment", x = 0.23, y = -0.05, xend = 0.23, yend = 0.25,
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
            plot_data_bifurcation %>% filter(lambda > 0.3)) +
  geom_line(aes(x = lambda, y = fixed_I),
            linewidth = 1.0, linetype = "11",
            plot_data_bifurcation %>% filter(lambda <= 0.31)) +
  
  annotate("point", x = 0.3, y = 0.075, size = 4) +
  annotate("point", x = 0.3, y = 0.075, size = 2, colour = "white") +
  
  annotate("label", hjust = 1, x = 0.15, y = 0.12, label = "Unstable\nfixed\npoint", size = 7, lineheight = 0.8, label.size = 0) +
  annotate("label", hjust = 0, x = 0.3, y = 0.15, label = "Stable fixed\npoint", size = 7, lineheight = 0.8, label.size = 0) +
  annotate("label", hjust = 0, x = 0.24, y = 0.22, label = "Limit cycle", size = 7, lineheight = 0.8, label.size = 0) +
  
  annotate("label", hjust = 0.5, x = 0.23, y = -0.05, label = "λ = 0.23", size = 7, lineheight = 0.8, label.size = 0.8, label.r = unit(0, "cm")) +
  annotate("label", hjust = 0.5, x = 0.35, y = -0.05, label = "λ = 0.35", size = 7, lineheight = 0.8, label.size = 0.8, label.r = unit(0, "cm")) +
  
  annotate("segment", x = 0.37, y = 0.14, xend = 0.39, yend = 0.105, linewidth = 0.7) +
  annotate("segment", x = 0.3, y = 0.2, xend = 0.27, yend = 0.17, linewidth = 0.7) +
  annotate("segment", x = 0.15, y = 0.09, xend = 0.2, yend = 0.06, linewidth = 0.7) +
  
  xlab("Decay rate λ") +
  ylab("Infected I") +
  
  coord_cartesian(xlim = c(0, 0.5),
                  ylim = c(0, 0.25),
                  expand = FALSE,
                  clip = "off") +
  
  theme_minimal() +
  
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(),
        text = element_text(family = "Helvetica", colour = "black", size = 20),
        line = element_line(linewidth = 0.7),
        axis.title.x = element_text(margin = margin(b = 0, t = 1.3, unit = "cm")),
        axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm")),
        axis.text.y = element_text(margin = margin(r = 0.3, unit = "cm")),
        plot.margin = margin(t = 0.5, l = 0.5, b = 0.5, r = 0.5, unit = "cm"),
        axis.ticks.length=unit(-0.1, "cm"))



