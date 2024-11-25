

library(tidyverse)
library(rhdf5)

library(patchwork)

source("R/plot_theme.R")

trans <- function(p, log_a, b) {
  return(plogis(b * qlogis(p) - b * log_a))
}

KD <- 1
b1 <- 2.5
log_a1 <- -2

neut_threshold <- 0.95
c_neut_threshold <- KD * exp(log_a1) * exp(1) ^ (qlogis(neut_threshold) / b1)


plot_data <- tibble(
  c = 10^seq(-2,3,by = 0.05)
) %>%
  mutate(p_bound = c / (KD + c),
         p_neut = trans(p_bound, log_a = log_a1, b = b1),
         p_protect = trans(p_neut, log_a = 4, b = 2))

p_common <- list(
  plot_theme_paper,
  theme(panel.grid.major.x = element_line(colour = "grey80", linetype = "48", linewidth = 0.5),
        plot.margin = margin(r = 1, l = 1, unit = "cm"))
)



p_A <- ggplot() +
  geom_line(aes(x = c, y = p_neut),
            plot_data) +
  
  geom_vline(xintercept = c_neut_threshold, colour = colour_C) +
  annotate("point", x = c_neut_threshold, y = 0.95, colour = colour_C) +
  
  scale_x_log10(labels = scales::label_log(), 
                breaks = scales::breaks_log(n = 8)) +
  xlab("Concentration *c<sub>x</sub>*") +
  ylab("*p*<sub>neut</sub>") +
  
  plot_theme_paper +
  theme(panel.grid.major.x = element_line(colour = "grey70", linetype = "48", linewidth = 0.5),
        axis.text.x.top = element_text(margin = margin(b = 0.3, unit = "cm")),
        plot.margin = margin(t = 1))

p_B <- ggplot() +
  
  geom_vline(aes(xintercept = c / c_neut_threshold), 
             tibble(c = 10^seq(-2, 3, by = 1)), 
             colour = "grey70", linetype = "48") +
  geom_line(aes(x = c / c_neut_threshold, y = p_neut),
            plot_data) +
  geom_vline(xintercept = 1, colour = colour_C) +
  annotate("point", x = 1, y = 0.95, colour = colour_C) +
  
  scale_x_log10(labels = scales::label_comma(), 
                breaks = 10^(0:3)) +
  xlab("Viral neutralisation titre mu") +
  ylab("*p*<sub>neut</sub>") +
  
  plot_theme_paper +
  theme(panel.grid.major.x = element_line(colour = shades::opacity(colour_C, 0.95), linetype = "84", linewidth = 0.5),
        axis.text.x.top = element_text(margin = margin(b = 0.3, unit = "cm")),
        plot.margin = margin(t = -2, b = 1))


p_A / p_B


p_C <- ggplot() +
  geom_line(aes(x = c, y = p_neut),
            plot_data) +
  
  geom_vline(aes(xintercept = c), 
             tibble(c = c_neut_threshold * 10^seq(-2, 3, by = 1)), 
             colour = colour_C, alpha = 0.5, linetype = "84") +
  geom_vline(xintercept = c_neut_threshold, colour = colour_C) +
  annotate("point", x = c_neut_threshold, y = 0.95, colour = colour_C) +
  
  scale_x_log10(labels = scales::label_log(), 
                breaks = scales::breaks_log(n = 8),
                sec.axis = sec_axis(transform = function(x) {x / c_neut_threshold},
                                    breaks = scales::breaks_log(8),
                                    labels = scales::label_log(),
                                    name = "Viral neutralisation titre")) +
  xlab("Concentration *c<sub>x</sub>*") +
  
  plot_theme_paper +
  theme(panel.grid.major.x = element_line(colour = "grey70", linetype = "48", linewidth = 0.5),
        axis.text.x.top = element_text(margin = margin(b = 0.3, unit = "cm")),
        plot.margin = margin(t = 10, b = 1))

p_A / p_C
