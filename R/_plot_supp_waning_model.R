

library(tidyverse)
library(patchwork)

source("R/plot_theme.R")




k <- 128
rho <- 0.005
C <- 8
t_max <- 200

num_data <- expand_grid(
  t = seq(0, t_max, by = 0.1),
  i = 0:k
) %>%
  mutate(
    p_t = dpois(i, rho * k * sqrt(t)),
    f = 10^(-C * i / k)
  )

closed_data <- expand_grid(
  t = 0:t_max
) %>%
  mutate(mean = exp(-C * rho * log(10) * t),
         var = (C^2 * rho * log(10)^2)/(k * t),
         upper = exp(-(C * rho * log(10) - sqrt(var)) * t),
         lower = exp(-(C * rho * log(10) + sqrt(var)) * t)
         )

ggplot() +
  geom_tile(aes(x = t, y = f, fill = pmin(0.2, p_t)),
            num_data) +
  
  geom_line(aes(x = t, y = mean), closed_data) +
  geom_line(aes(x = t, y = upper), closed_data) +
  geom_line(aes(x = t, y = lower), closed_data) +
  
  scale_fill_viridis_c() +
  
  scale_y_log10()


ggplot() +
  geom_line(aes(x = f, y = p_t),
            num_data %>% filter(t %in% c(0, 5, 10, 15))) +
  
  facet_wrap(~t, ncol = 1, scales = "free_y") +
  
  coord_cartesian(xlim = c(10^-4, NA)) +
  
  scale_x_log10()


r_data <- tibble(
  t = 1:t_max
) %>%
  rowwise() %>%
  # mutate(r = list(rpois(1000, rho * k * t) * C * log(10) / (k * t))) %>%
  mutate(r = list(rpois(1000, rho * k * t) * C * log(10) / (k * sqrt(t)))) %>%
  unnest(r) %>%
  
  group_by(t) %>%
  summarise(mu = mean(r),
            var = var(r))

ggplot() +
  
  geom_line(aes(x = t, y = var),
            r_data) 




 expand_grid(
  t = 1:t_max,
  k = 10:100
) %>%
  mutate(mean = exp(-C * rho * log(10) * t),
         var = (C^2 * rho * log(10)^2)/(k * t),
         upper = exp(-(C * rho * log(10) - sqrt(var)) * t),
         lower = exp(-(C * rho * log(10) + sqrt(var)) * t)
  ) %>%
   
   ggplot() +
   geom_tile(aes(x = t, y = k, fill = var))
  
  






