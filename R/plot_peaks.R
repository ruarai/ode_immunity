library(tidyverse)
library(rhdf5)
library(patchwork)

source("R/plot_theme.R")
source("R/read_seasonality_data.R")

season_duration <- 365 / 4
half_season <- season_duration / 2

y_peaks <- h5read("data/paper/period_over_grid.jld2", "y_peaks")


seasonality_data <- read_seasonality_data("data/paper/period_over_grid.jld2")



peaks_data <- tibble(
  eta = y_peaks[1,],
  r = y_peaks[2,],
  t = y_peaks[3,],
  prevalence = y_peaks[4,]
) %>%
  mutate(t_mod = (t + 180) %% 365 - 180,
         t_mod_10 = t %% (365 * 8),
         t_start = t - 400000) %>%
  mutate(r_label = str_c("Decay rate <i>r</i> = ", r)) %>%
  
  left_join(seasonality_data, by = c("eta", "r"))


peaks_summ <- peaks_data %>%
  mutate(t_mod = t %% 365,
         theta = t_mod / 365 * 2 * pi,
         x = cos(theta)) %>% 
  
  group_by(eta, r, r_label) %>%
  summarise(mean_x = weighted.mean(x, prevalence)) 


ggplot() +
  geom_tile(aes(x = eta, y = r, fill = mean_x),
            peaks_summ) +
  
  scale_fill_distiller(type = "div",
                       palette = "RdBu",
                       direction = 1,
                       limits = c(-1, 1))
  


peaks_subset <- peaks_data %>%
  # filter(r %in% c(0.01, 0.02, 0.03, 0.04))
  filter(r %in% c(0.01, 0.03, 0.05, 0.07, 0.09)) %>% 
  mutate(period = period / 365,
         period = pmin(period, 8),
         period = factor(round(period)))


years_span <- 10
peaks_data %>%
  filter(eta %in% c(0.1, 0.2, 0.3, 0.4), t_start < 365 * years_span) %>%
  ggplot() +
  geom_rect(aes(xmin = t - half_season, xmax = t + half_season, ymin = -Inf, ymax = Inf),
            tibble(t = seq(1, years_span - 1, by = 1) * 365),
            alpha = 0.1, fill = "black") +
  
  geom_point(aes(x = t_start, y = r, colour = prevalence),
             size = 0.5) +
  
  scale_x_continuous(breaks = c(0:years_span) * 365,
                     labels = function(x) {x / 365}) +
  
  scale_colour_gradientn(name = "Peak infection\nprevalence",
                         colours = RColorBrewer::brewer.pal(9, "Blues")[3:9],
                         limits = c(0, NA)) +
  
  facet_wrap(~(factor(eta)), ncol = 1) +
  
  # facet_wrap(~fct_rev(factor(r_label)), ncol = 1, dir = "v", scales = "free_y") +
  
  plot_theme_paper +
  
  theme(strip.text = element_markdown()) +
  
  ggtitle(NULL, "Infection peak prevalence occurrence by time of year")

peaks_data %>%
  filter(r == 0.09, eta == 0.2)


# \sin(\frac{2 \pi}{365} t + \frac{\pi}{2})
data_seasonality <- tibble(
  t = 0:365,
  rel_seasonality = sin(2 * pi / 365 * t + pi / 2)
)

ggplot() +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = -30, ymax = 30 * 2, fill = "black",
           alpha = 0.2) +
  
  geom_point(aes(x = eta, y = t_mod, colour = prevalence),
             peaks_subset,
             size = 0.3, pch = 15) +
  
  plot_theme_paper +
  
  scale_colour_gradientn(name = "Peak infection\nprevalence",
                         colours = RColorBrewer::brewer.pal(9, "Blues")[3:9],
                         limits = c(0, NA)) +
  
  xlab("Seasonality constant <i>η</i>") +
  ylab("Time of year") +
  
  scale_y_continuous(breaks = c(-180, -90, 0, 90, 180),
                     labels = c("Jul", "Oct", "Jan", "Apr", "Jul")) +
  
  coord_cartesian(ylim = c(-180, 180)) +
  
  scale_x_continuous(breaks = c(0, 0.25, 0.5)) +
  
  facet_wrap(~fct_rev(factor(r_label)), ncol = 1, dir = "v", scales = "free_y") +
  
  theme(strip.text = element_markdown()) +
  
  ggtitle(NULL, "Infection peak prevalence occurrence by time of year")



ggplot() +
  # annotate("rect", xmin = -half_season, xmax = half_season, ymin = 0, ymax = 0.1, fill = "black",
  #          alpha = 0.2) +
  
  geom_point(aes(x = t_mod, y = r, colour = prevalence),
             peaks_data %>% filter(eta %in% c(0.1, 0.25, 0.5)),
             size = 0.3, pch = 15) +
  
  plot_theme_paper +
  
  scale_colour_gradientn(name = "Peak infection\nprevalence",
                         colours = RColorBrewer::brewer.pal(9, "Blues")[3:9],
                         limits = c(0, NA)) +
  
  # xlab("Seasonality constant <i>η</i>") +
  xlab("Time of year") +
  
  scale_x_continuous(breaks = c(-180, -90, 0, 90, 180),
                     labels = c("Jul", "Oct", "Jan", "Apr", "Jul")) +
  
  coord_cartesian(xlim = c(-180, 180)) +
  
  # scale_x_continuous(breaks = c(0, 0.25, 0.5)) +
  
  facet_wrap(~(factor(eta)), ncol = 3, dir = "v", scales = "free_y") +
  
  theme(strip.text = element_markdown(),
        panel.grid.major.x = element_gridline) +
  
  ggtitle(NULL, "Infection peak prevalence occurrence by time of year")
