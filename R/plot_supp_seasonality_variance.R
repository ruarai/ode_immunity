

plot_data %>%
  mutate(
    class = case_when(
      quasiperiodic ~ "Quasiperiodic",
      periodic ~ "Periodic",
      chaotic ~ "Chaotic",
      TRUE ~ "Unclassified"
    )
  ) %>%
  
  mutate(
    class = factor(
      class, 
      levels = c("Periodic", "Quasiperiodic", "Chaotic", "Unclassified")
    )
  ) %>% 
  
  filter(class != "Unclassified") %>% 
  
  ggplot() +
  geom_histogram(aes(x = season_var), binwidth = 0.01,
                 fill = "grey5") +
  
  facet_wrap(~class, ncol = 1, scales = "free_y") +
  
  xlab("Circular variance") +
  ylab("Frequency") +
  
  scale_y_continuous(labels = scales::label_comma()) +
  
  plot_theme_paper


ggsave(
  "results/results_supp_variance_class.pdf",
  device = cairo_pdf,
  width = 8, height = 6,
  bg = "white"
)



