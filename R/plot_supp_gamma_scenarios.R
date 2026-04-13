
# \frac{c_i^b}{c_\text{mid}^b+c_i^b}

p_acq <- tibble(
  c = 2^seq(0, 8, length.out = 33),
  b = 8,
  c_mid = 2^ 3
) %>% 
  mutate(
    p_acq = (c ^ b) / (c_mid ^ b + c ^ b)
  ) %>% 
  
  
  ggplot() +
  
  geom_vline(xintercept = 2^3, linetype = "44") +
  
  geom_line(aes(x = c, y = p_acq), linewidth = 0.7) +
  
  scale_x_continuous(trans = "log2", breaks = scales::log_breaks(base = 2, n = 8)) +
  
  ylab("Protection against infection") +
  
  xlab("Effective antibody level") +
  
  plot_theme_paper


p_rec <- expand_grid(
  c = 2^seq(0, 8, length.out = 33),
  b = 8,
  scenario = c("Scenario 1", "Scenario 2", "Baseline")
) %>% 
  mutate(
    c_mid = case_when(
      scenario == "Scenario 1" ~ 2^2,
      scenario == "Scenario 2" ~ 2^4,
      TRUE ~ NA
    )
  ) %>% 
  mutate(
    gamma = 0.8 + (c ^ b) / (c_mid ^ b + c ^ b),
    gamma = replace_na(gamma, 0.8)
  ) %>% 
  
  
  ggplot() +
  
  geom_line(aes(x = c, y = gamma, colour = scenario), linewidth = 0.7) +
  
  geom_vline(
    aes(xintercept = c_mid, colour = scenario), linetype = "44",
    tibble(scenario = c("Scenario 1", "Scenario 2"), c_mid = c(2^2, 2^4))
  ) +
  
  scale_colour_brewer(name = "Scenario", type = "qual", palette = 2) +
  
  scale_x_continuous(trans = "log2", breaks = scales::log_breaks(base = 2, n = 8)) +
  
  ylab("Recovery rate") +
  
  xlab("Effective antibody level") +
  
  coord_cartesian(ylim = c(0, 2.3)) +
  
  plot_theme_paper +
  theme(legend.position = "inside",
        legend.position.inside = c(0.1, 0.7))

library(patchwork)
(p_acq + ggtitle("<b>A</b> — Protection against infection <i>ω</i><sub><i>i</i></sub>")) / 
  (p_rec + ggtitle("<b>B</b> — Recovery rate <i>γ</i><sub><i>i</i></sub>"))


ggsave(
  "results/results_supp_gamma_scenarios.png",
  device = png, bg = "white",
  width = 10, height = 7
)
