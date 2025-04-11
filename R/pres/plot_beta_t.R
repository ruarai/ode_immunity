

eta <- 0.5
beta <- 0.375

tibble(
  t = seq(0, 365 * 5)
) %>%
  mutate(beta_t = beta * (1 + eta * cos(2 * pi / 365 * t))) %>%
  
  ggplot() +
  geom_line(aes(x = t, y = beta_t)) +
  
  geom_vline(aes(xintercept = d), tibble(d = seq(0, 365 * 5, by = 365)),
             linewidth = 0.6, linetype = "42") +
  
  
  scale_x_continuous(breaks = seq(0, 365 * 5, by = 365),
                     labels = function(x) round(x/365)) +
  
  coord_cartesian(ylim = c(0, 0.8)) +
  
  ylab("<i>β</i>(<i>t</i>)") +
  
  xlab("Time (years)") +
  
  plot_theme_paper
