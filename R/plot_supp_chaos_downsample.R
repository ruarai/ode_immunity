
library(tidyverse)
library(rhdf5)
library(patchwork)


source("R/plot_theme.R")

hdf5_file <- "data/chaos_downsample.jld2"

t_post_burn_in <- h5read(hdf5_file, "t_post_burn_in")
x_eta <- h5read(hdf5_file, "x_eta")
y_inc <- h5read(hdf5_file, "y_inc")

downsample_rate <- h5read(hdf5_file, "downsample_rate")
chaos_result <- h5read(hdf5_file, "chaos_result")




y_tbl <- y_inc %>%
  reshape2::melt(varnames = c("ix", "t"), value.name = "incidence") %>%
  mutate(t = t_post_burn_in[t],
         eta = x_eta[ix])


p_bifur <- ggplot() +
  
  geom_vline(aes(xintercept = eta),
             chaos_tbl %>% filter(chaos > 0.5, downsample == 80),
             alpha = 0.5, colour = ggokabeito::palette_okabe_ito(5),
             linewidth = 1.0) +
  
  geom_point(aes(x = eta, y = incidence),
             size = 0.2, stroke = 0.1,
             y_tbl %>% filter(t %% 365 == 180))  +
  
  coord_cartesian(ylim = c(0, 0.001), xlim = c(0.25, 0.45)) +
  
  xlab("Seasonal forcing strength <i>η</i>") + ylab("Infection incidence") +
  plot_theme_paper +
  
  ggtitle("<b>A</b> — Bifurcation diagram over seasonal forcing strength")

p_bifur

chaos_tbl <- chaos_result %>%
  reshape2::melt(varnames = c("ix", "ix_d"), value.name = "chaos") %>%
  mutate(downsample = downsample_rate[ix_d],
         eta = x_eta[ix])


p_chaos_test <- ggplot() +
  
  geom_hline(yintercept = 80, linewidth = 2.0, colour = "grey80") +
  
  geom_tile(aes(x = eta, y = downsample),
            fill = ggokabeito::palette_okabe_ito(5), alpha = 0.5,
            chaos_tbl %>% filter(chaos > 0.5)) +
  
  coord_cartesian(xlim = c(0.25, 0.45), ylim = c(0, NA)) +
  
  xlab("Seasonal forcing strength <i>η</i>") + ylab("Downsampling rate") +
  
  plot_theme_paper +
  
  ggtitle("<b>B</b> — Chaos classification by downsampling rate")

p_bifur / p_chaos_test


ggsave(
  "results/results_supp_chaos_downsample.png",
  width = 8, height = 6.75,
  bg = "white"
)


