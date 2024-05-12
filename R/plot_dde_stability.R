



library(tidyverse)

plot_data <- expand_grid(
  R0 = seq(0, 4, by = 0.05),
  lambda = seq(0, 0.2, by = 0.001),
  mu = 1,
  sigma = 0.05,
  gamma = c(0.5, 1.0, 1.5)
) %>%
  mutate(stable_lower = 2 * lambda * (lambda + mu * gamma) /  (gamma ^ 2 * (sigma ^ 2 + mu ^ 2)) + 1 <= R0,
         stable_upper = R0 <= 2 * mu * gamma / lambda + 3,
         stable_2 = mu / lambda > 1 / gamma,
         stable = stable_lower & stable_upper & stable_2)




ggplot(plot_data) +
  #geom_tile(aes(x = R0, y = lambda, fill = stable_lower)) +
  geom_tile(aes(x = R0, y = lambda, fill = stable)) +
  
  annotate("point", x = 1.5, y = 0.05) +
  
  facet_wrap(~gamma)












