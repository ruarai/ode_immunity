

library(tidyverse)
library(rhdf5)
library(patchwork)


source("R/plot_theme.R")

scenario_labels <- c(
  "Baseline model", 
  "No boosting", 
  "Multiplicative boosting"
)


seq_t <- h5read("data/paper/supp_boosting_basic.jld2", "seq_t")
c_levels <- h5read("data/paper/supp_boosting_basic.jld2", "c_levels")

plot_data_inf <- h5read("data/paper/supp_boosting_basic.jld2", "results_inf") %>%
  reshape2::melt(varnames = c("t", "scenario"), value.name = "prevalence") %>%
  as_tibble() %>% 
  mutate(
    t = seq_t[t],
    scenario = c("independent", "none", "multiplicative")[scenario],
    scenario = factor(scenario, c("independent", "none", "multiplicative"), labels = scenario_labels)
  ) %>%
  
  filter(t < 1500)


plot_data_sus <- h5read("data/paper/supp_boosting_basic.jld2", "results_sus") %>%
  reshape2::melt(varnames = c("t", "strata", "scenario"), value.name = "prevalence") %>%
  as_tibble() %>% 
  mutate(
    t = seq_t[t],
    c = c_levels[strata],
    scenario = c("independent", "none", "multiplicative")[scenario],
    scenario = factor(scenario, c("independent", "none", "multiplicative"), labels = c("Baseline model", "No boosting", "Multiplicative boosting"))
  )


plot_data_means <- plot_data_sus %>%
  group_by(t, scenario) %>%
  mutate(p = prevalence / sum(prevalence)) %>% 
  summarise(c = sum(p * c), .groups = "drop") %>%
  
  filter(t < 1500)


p_example_prevalence <- ggplot() +
  geom_line(aes(x = t, y = prevalence, colour = scenario),
            plot_data_inf %>% filter(scenario == "Multiplicative boosting"),
            linewidth = 0.7) +
  geom_line(aes(x = t, y = prevalence, colour = scenario),
            plot_data_inf %>% filter(scenario == "No boosting"),
            linewidth = 0.9) +
  
  scale_x_continuous(breaks = scales::breaks_extended(),
                     labels = scales::label_comma()) +
  
  ggokabeito::scale_colour_okabe_ito(name = "Scenario", order = c(5, 9)) +
  guides(colour = guide_legend(reverse = TRUE)) +
  
  xlab("Time (days)") + ylab("Prevalence") +
  
  plot_theme_paper +
  theme(panel.grid.major = element_gridline,
        legend.position = "bottom",
        legend.direction = "horizontal") +
  
  ggtitle(NULL, "<b>D</b> — Infection prevalence")

p_example_mean_antibodies <- ggplot() +
  geom_line(aes(x = t, y = c, colour = scenario),
            linewidth = 0.7,
            plot_data_means %>% filter(scenario == "Multiplicative boosting")) +
  geom_line(aes(x = t, y = c, colour = scenario),
            linewidth = 0.9,
            plot_data_means %>% filter(scenario == "No boosting")) +
  
  coord_cartesian(ylim = c(10^-0.5, 10^7)) +
  
  xlab("Time (days)") + ylab("Concentration") +
                                   
  ggokabeito::scale_colour_okabe_ito(name = "Scenario", order = c(5, 9)) +
  guides(colour = guide_legend(reverse = TRUE)) +
  
  scale_y_continuous(trans = "log10", labels = scales::label_log(base = 10), breaks = 10^c(0, 2, 4, 6, 8)) +
  
  scale_x_continuous(breaks = scales::breaks_extended(),
                     labels = scales::label_comma()) +
  
  plot_theme_paper +
  theme(panel.grid.major = element_gridline,
        legend.position = "bottom") +
  ggtitle(NULL, "<b>E</b> — Population mean antibody concentration")

(p_example_prevalence / p_example_mean_antibodies) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")



x_r <- h5read("data/paper/bifurcations.jld2", "x_r")


data_I_sol <- bind_rows(
  reshape2::melt(h5read("data/paper/bifurcations.jld2", "y_I_sol"), 
                 varnames = c("r", "t"), value.name = "prev") %>% 
    mutate(r = x_r[r],
           scenario = "independent") %>%
    filter(r > 0),
  
  reshape2::melt(h5read("data/paper/bifurcations_boosting.jld2", "y_I_sol"),
                 varnames = c("r", "t"), value.name = "prev") %>% 
    mutate(r = x_r[r],
           scenario = "multiplicative") %>%
    filter(r > 0)
) %>%
  mutate(scenario = factor(scenario, c("independent", "multiplicative"), labels = c("Baseline model", "Multiplicative boosting")))
  

data_inc_sol <- bind_rows(
  reshape2::melt(h5read("data/paper/bifurcations.jld2", "y_inc_sol"),
                 varnames = c("r", "t"), value.name = "inc") %>% 
    mutate(r = x_r[r],
           scenario = "independent"),
  
  reshape2::melt(h5read("data/paper/bifurcations_boosting.jld2", "y_inc_sol"),
                 varnames = c("r", "t"), value.name = "inc") %>% 
    mutate(r = x_r[r],
           scenario = "multiplicative")
) %>%
  mutate(scenario = factor(scenario, c("independent", "multiplicative"), labels = c("Baseline model", "Multiplicative boosting")))

data_mean_incidence <- data_inc_sol %>%
  filter(t > 10000) %>% 
  group_by(r, scenario) %>%
  summarise(mean_inc = mean(inc))

maxmins <- data_I_sol %>%
  filter(t > days_burn_in, r > 0.002) %>% 
  group_by(r, scenario) %>%
  summarise(max = max(prev), min = min(prev))


data_fixed <- bind_rows(
  reshape2::melt(h5read("data/paper/bifurcations.jld2", "y_fixed_I"),
                 varnames = c("r"), value.name = "prev") %>% 
    mutate(r = x_r[r],
           scenario = "independent"),
  reshape2::melt(h5read("data/paper/bifurcations_boosting.jld2", "y_fixed_I"),
                 varnames = c("r"), value.name = "prev") %>% 
    mutate(r = x_r[r],
           scenario = "multiplicative")
) %>%
  mutate(scenario = factor(scenario, c("independent", "multiplicative"), labels = c("Baseline model", "Multiplicative boosting"))) %>%
  filter(r > 0)

bifur_points <- maxmins %>% 
  left_join(data_fixed) %>%
  mutate(diff = max - min) %>%
  ungroup() %>% 
  group_by(scenario) %>% 
  filter(diff > 1e-4) %>%
  slice(n()) %>%
  mutate(mean = max / 2 + min / 2) %>% 
  select(r, scenario, prev = mean)


p_bifurcation <- ggplot() +
  geom_vline(xintercept = 0.05,
             colour = "grey80", linewidth = 1.0, alpha = 0.5) +
  geom_line(aes(x = r, y = max, colour = scenario, linewidth = scenario),
            maxmins) +
  geom_line(aes(x = r, y = min, colour = scenario, linewidth = scenario),
            maxmins) +
  geom_line(aes(x = r, y = prev, colour = scenario, linewidth = scenario),
            linetype = "44",
            data_fixed) +
  
  
  geom_point(aes(x = r, y = prev, colour = scenario), 
             bifur_points,
             size = 3) +
  
  geom_point(aes(x = r, y = prev), 
             bifur_points,
             size = 1.5, colour = "white") +
  
  ggokabeito::scale_colour_okabe_ito(name = "Scenario", order = c(5, 9)) +
  
  scale_linewidth_manual(values = c("Multiplicative boosting" = 0.7, "Baseline model" = 0.9)) +
  
  xlab("Mean antibody decay rate <i>r</i>") +
  ylab("Prevalence") +
  
  # coord_cartesian(xlim = c(-0.0002, 0.0075),
  #                 ylim = c(-0.005, 0.06),
  #                 expand = FALSE) +
  
  plot_theme_paper +
  theme(legend.position = "none",
        axis.text.x.top = element_text(margin = margin(b = 0.3, unit = "cm")),
        panel.grid.major = element_gridline,
        plot.subtitle = element_markdown()) +
  
  ggtitle(NULL, "<b>A</b> — Bifurcation over antibody decay rate <i>r</i>")

p_bifurcation


data_period <- bind_rows(
  reshape2::melt(h5read("data/paper/bifurcations.jld2", "period"),
                 varnames = c("r", "name"), value.name = "value") %>%
    mutate(r = x_r[r],
           scenario = "independent"),
  reshape2::melt(h5read("data/paper/bifurcations_boosting.jld2", "period"),
                 varnames = c("r", "name"), value.name = "value") %>%
    mutate(r = x_r[r],
           scenario = "multiplicative"),
) %>%
  mutate(
    scenario = factor(scenario, c("independent", "multiplicative"), labels = c("Baseline model", "Multiplicative boosting")),
    name = c("period", "period_sd", "period_n")[name]
  ) %>%
  
  pivot_wider() %>% 
  
  left_join(bifur_points %>% select(r_bifur = r, scenario)) %>%
  filter(r < r_bifur, period_sd < 1, period_n > 1)



min_periods <- data_period %>% group_by(scenario) %>% slice(1) %>%
  select(r_min = r, scenario)

p_period <- ggplot() +
  geom_vline(xintercept = 0.05,
             colour = "grey80", linewidth = 1.0, alpha = 0.5) +
  geom_line(aes(x = r, y = 365 / period, colour = scenario, linewidth = scenario),
            data_period) +
  
  geom_point(aes(x = r, y = 365 / period, colour = scenario), 
             data_period %>% group_by(scenario) %>% filter(r == max(r)),
             size = 3) +
  
  geom_point(aes(x = r, y = 365 / period), 
             data_period %>% group_by(scenario) %>% filter(r == max(r)),
             size = 1.5, colour = "white") +
  
  ggokabeito::scale_colour_okabe_ito(name = "Scenario", order = c(5, 9)) +
  
  scale_linewidth_manual(values = c("Multiplicative boosting" = 0.7, "Baseline model" = 0.9)) +
  
  xlab("Mean antibody decay rate <i>r</i>") +
  ylab("Frequency (years<sup>-1</sup>)") +
  
  coord_cartesian(xlim = c(0, 0.15),
                  ylim = c(-0.1, 4)) +
  
  plot_theme_paper +
  theme(legend.position = "none",
        panel.grid.major = element_gridline) +
  
  ggtitle(NULL,"<b>B</b> — Periodic solution frequency")

p_period


p_attack_rate <- ggplot() +
  geom_vline(xintercept = 0.05,
             colour = "grey80", linewidth = 1.0, alpha = 0.5) +
  geom_line(aes(x = r, y = mean_inc * 365, colour = scenario, linewidth = scenario),
            data_mean_incidence) +
  
  xlab("Mean antibody decay rate <i>r</i>") +
  ylab("Infection incidence") +
  
  ggokabeito::scale_colour_okabe_ito(name = "Scenario", order = c(5, 9)) +
  
  scale_linewidth_manual(values = c("Multiplicative boosting" = 0.7, "Baseline model" = 0.9)) +
  
  
  plot_theme_paper +
  theme(legend.position = "none",
        panel.grid.major = element_gridline) +
  
  ggtitle(NULL, "<b>C</b> — Average annual infection incidence at solution")
p_attack_rate

p_left <- (p_bifurcation / p_period / p_attack_rate) +
  plot_layout(tag_level = "new")

p_right <- (p_example_prevalence / p_example_mean_antibodies) +
  plot_layout(tag_level = "keep")

(p_left | p_right) / guide_area() + 
  plot_layout(guides = "collect", heights = c(1, 0.1)) +
  theme(plot.tag = element_text(face = "bold", size = 15))

ggsave(
  "results/results_supp_boosting.pdf",
  device = cairo_pdf,
  width = 14, height = 10,
  bg = "white"
)

