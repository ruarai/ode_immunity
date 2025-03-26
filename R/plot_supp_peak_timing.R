
library(tidyverse)
library(rhdf5)
library(patchwork)

source("R/plot_theme.R")


x_eta <- h5read("data/period_over_grid_examples.jld2", "x_eta")
x_r <- h5read("data/period_over_grid_examples.jld2", "x_r")
y_inf <- h5read("data/period_over_grid_examples.jld2", "y_inf")
y_sus <- h5read("data/period_over_grid_examples.jld2", "y_sus")


t_ex_start <- 365 * 100
t_ex_end <- 365 * (100 + 70)


plot_data_ex_inf <- y_inf %>%
  reshape2::melt(c("i", "t"), value.name = "prevalence") %>%
  mutate(eta = x_eta[i], r = x_r[i], label = x_labels[i]) %>%
  filter(t >= t_ex_start, t < t_ex_end) %>%
  group_by(label, eta, r, t) %>%
  summarise(prevalence = sum(prevalence))


peaks <- plot_data_ex_inf %>%
  group_by(label, eta, r) %>%
  filter(prevalence > lag(prevalence) & prevalence > lead(prevalence)) %>%
  mutate(year = (t) / 365,
         group = floor(year))

peaks_end_of_year <- peaks %>% 
  group_by(label) %>%
  mutate(leading_year = lead(year), leading_t = lead(t)) %>%
  group_by(label, group) %>% 
  slice(n())




ggplot() +
  geom_path(aes(x = year, y = t %% 365, group = group),
            linewidth = 0.3, alpha = 0.3,
            peaks) +

  geom_segment(aes(x = year, y = t %% 365, xend = leading_year, yend = leading_t %% 365 + 365),
               linewidth = 0.3, alpha = 0.3,
               peaks_end_of_year) +

  geom_segment(aes(x = year, y = t %% 365 - 365, xend = leading_year, yend = leading_t %% 365),
               linewidth = 0.3, alpha = 0.3,
               peaks_end_of_year) +
  
  geom_point(aes(x = year, y = t %% 365),
             size = 0.5,
             peaks) +
  
  facet_wrap(~label, ncol = 1) +
  
  scale_y_continuous(breaks = c(0, 180, 365)) +
  
  coord_cartesian(ylim = c(0, 365)) +
  
  plot_theme_paper +
  
  xlab("Time (years)") + ylab("Time of year (days)") +
  
  theme(strip.text = element_markdown())




ggplot() +
  
  geom_segment(aes(x = year, y = t %% 365, xend = year + 3 / 365, yend = t %% 365 + 3, colour = log10(prevalence) ),
               linewidth = 2.5,
               plot_data_ex_inf %>% mutate(year = t / 365)) +
  
  facet_wrap(~label, ncol = 1) +
  
  scale_y_continuous(breaks = c(0, 180, 365)) +
  
  coord_cartesian(ylim = c(0, 365)) +
  
  scale_colour_viridis_c() +
  
  plot_theme_paper +
  
  xlab("Time (years)") + ylab("Time of year (days)") +
  
  theme(strip.text = element_markdown())







