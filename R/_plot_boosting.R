

library(tidyverse)
library(rhdf5)
library(patchwork)


source("R/plot_theme.R")



sol_t <- h5read("data/paper/basic_boosting.jld2", "sol_t")
c_levels <- h5read("data/paper/basic_boosting.jld2", "c_levels")


plot_data <- sol_t %>%
  reshape2::melt(varnames = c("scenario", "t", "class", "strata"), value.name = "prevalence") %>% 
  mutate(class = c("S", "I")[class],
         scenario = c("none", "loglinear")[scenario],
         scenario = factor(scenario, c("loglinear", "linear", "none"), labels = c("Multiplicative boosting", "Additive boosting", "No boosting")),
         c = c_levels[strata])


plot_data_summ_inf <- plot_data %>%
  group_by(scenario, t, class) %>%
  summarise(prevalence = sum(prevalence), .groups = "drop") %>%
  filter(class == "I")

plot_data_means <- plot_data %>%
  group_by(t, scenario, class) %>%
  mutate(p = prevalence / sum(prevalence)) %>% 
  summarise(c = sum(p * c), .groups = "drop") %>%
  filter(class == "S")


p_example_prevalence <- ggplot() +
  geom_line(aes(x = t, y = prevalence, colour = scenario),
            plot_data_summ_inf %>% filter(scenario != "No boosting"),
            linewidth = 0.7) +
  geom_line(aes(x = t, y = prevalence, colour = scenario),
            plot_data_summ_inf %>% filter(scenario == "No boosting"),
            linewidth = 0.9) +
  
  coord_cartesian(xlim = c(0, 2000)) +
  
  scale_x_continuous(breaks = scales::breaks_extended(),
                     labels = scales::label_comma()) +
  
  ggokabeito::scale_colour_okabe_ito(name = "Scenario", order = c(5, 9)) +
  guides(colour = guide_legend(reverse = TRUE)) +
  
  xlab("Time *t* (days)") + ylab("Prevalence") +
  
  plot_theme_paper +
  theme(panel.grid.major = element_gridline,
        legend.position = "right",
        legend.direction = "horizontal") +
  
  ggtitle(NULL, "Infection prevalence")

p_example_mean_antibodies <- ggplot() +
  geom_line(aes(x = t, y = c, colour = scenario),
            linewidth = 0.7,
            plot_data_means %>% filter(scenario != "No boosting")) +
  geom_line(aes(x = t, y = c, colour = scenario),
            linewidth = 0.9,
            plot_data_means %>% filter(scenario == "No boosting")) +
  
  coord_cartesian(xlim = c(0, 2000), ylim = c(10^-0.5, 10^7)) +
  
  xlab("Time *t* (days)") + ylab("Concentration") +
                                   
  ggokabeito::scale_colour_okabe_ito(name = "Scenario", order = c(5, 9)) +
  guides(colour = guide_legend(reverse = TRUE)) +
  
  scale_y_continuous(trans = "log10", labels = scales::label_log(base = 10), breaks = 10^c(0, 2, 4, 6, 8)) +
  
  scale_x_continuous(breaks = scales::breaks_extended(),
                     labels = scales::label_comma()) +
  
  plot_theme_paper +
  theme(panel.grid.major = element_gridline,
        legend.position = "none") +
  ggtitle(NULL, "Population mean antibody concentration")



x_r <- h5read("data/paper/bifurcations_w_boost.jld2", "x_r")
y_I_sol <- h5read("data/paper/bifurcations_w_boost.jld2", "y_I_sol")
y_fixed_I <- h5read("data/paper/bifurcations_w_boost.jld2", "y_fixed_I")
y_inc_sol <- h5read("data/paper/bifurcations_w_boost.jld2", "y_inc_sol")


data_I_sol <- y_I_sol %>%
  reshape2::melt(varnames = c("r", "scenario", "t"), value.name = "prev") %>% 
  mutate(scenario = c("none", "linear", "loglinear")[scenario],
         scenario = factor(scenario, c("loglinear", "linear", "none"), labels = c("Multiplicative boosting", "Additive boosting", "No boosting")),
         r = x_r[r])

maxmins <- data_I_sol %>%
  filter(t > 20000, r > 0.003) %>% 
  group_by(r, scenario) %>%
  summarise(max = max(prev), min = min(prev))

data_fixed <- y_fixed_I %>%
  reshape2::melt(varnames = c("r", "scenario"), value.name = "prev") %>% 
  mutate(scenario = c("none", "linear", "loglinear")[scenario],
         scenario = factor(scenario, c("loglinear", "linear", "none"), labels = c("Multiplicative boosting", "Additive boosting", "No boosting")),
         r = x_r[r])


data_inc_sol <- y_inc_sol %>%
  reshape2::melt(varnames = c("r", "scenario", "t"), value.name = "inc") %>% 
  mutate(scenario = c("none", "linear", "loglinear")[scenario],
         scenario = factor(scenario, c("loglinear", "linear", "none"), labels = c("Multiplicative boosting", "Additive boosting", "No boosting")),
         r = x_r[r])

data_mean_incidence <- data_inc_sol %>%
  filter(t > 10000) %>% 
  group_by(scenario, r) %>%
  summarise(mean_inc = mean(inc))

r_to_halflife <- function(x) {1 / (8 * x)}

bifur_points <- maxmins %>% 
  left_join(data_fixed) %>%
  mutate(diff = max - prev) %>%
  ungroup() %>% 
  group_by(scenario) %>% 
  filter(diff > 1e-4) %>%
  slice(n()) %>%
  select(r, scenario, prev = prev)


p_bifurcation <- ggplot() +
  # geom_vline(xintercept = 0.0025,
  #            colour = "grey80", linewidth = 1.0, alpha = 0.5) +
  geom_line(aes(x = r, y = max, colour = scenario, linewidth = scenario),
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
  
  scale_linewidth_manual(values = c("Multiplicative boosting" = 0.7, "Additive boosting" = 0.7, "No boosting" = 0.9)) +
  
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
  
  ggtitle(NULL, "Bifurcation over <i>ρ</i>")

p_bifurcation


period <- h5read("data/paper/bifurcations_w_boost.jld2", "period")
attack_rate <- h5read("data/paper/bifurcations_w_boost.jld2", "attack_rate")

data_period <- period %>%
  reshape2::melt(varnames = c("r", "scenario", "name"), value.name = "value") %>% 
  mutate(scenario = c("none", "linear", "loglinear")[scenario],
         name = c("period", "period_sd", "period_n")[name],
         scenario = factor(scenario, c("loglinear", "linear", "none"), labels = c("Multiplicative boosting", "Additive boosting", "No boosting")),
         r = x_r[r]) %>%
  
  pivot_wider() %>% 
  
  left_join(bifur_points %>% select(r_bifur = r, scenario)) %>%
  filter(r < r_bifur, period_sd < 1, period_n > 1)

min_periods <- data_period %>% group_by(scenario) %>% slice(1) %>%
  select(r_min = r, scenario)

data_attack_rate <- attack_rate %>%
  reshape2::melt(varnames = c("r", "scenario"), value.name = "attack_rate") %>% 
  mutate(scenario = c("none", "linear", "loglinear")[scenario],
         scenario = factor(scenario, c("loglinear", "linear", "none"), labels = c("Multiplicative boosting", "Additive boosting", "No boosting")),
         r = x_r[r]) %>%
  
  left_join(bifur_points %>% select(r_bifur = r, scenario)) %>%
  left_join(min_periods) %>% 
  filter(r < r_bifur, r >= r_min) %>%
  
  left_join(data_period %>% select(r, scenario, period))


p_period <- ggplot() +
  # geom_vline(xintercept = 0.0025,
  #            colour = "grey80", linewidth = 1.0, alpha = 0.5) +
  geom_line(aes(x = r, y = 365 / period, colour = scenario, linewidth = scenario),
            data_period) +
  
  geom_point(aes(x = r, y = 365 / period, colour = scenario), 
             data_period %>% group_by(scenario) %>% filter(r == max(r)),
             size = 3) +
  
  geom_point(aes(x = r, y = 365 / period), 
             data_period %>% group_by(scenario) %>% filter(r == max(r)),
             size = 1.5, colour = "white") +
  
  ggokabeito::scale_colour_okabe_ito(name = "Scenario", order = c(5, 9)) +
  
  scale_linewidth_manual(values = c("Multiplicative boosting" = 0.7, "Additive boosting" = 0.7, "No boosting" = 0.9)) +
  
  xlab("Mean antibody decay rate <i>r</i>") +
  ylab("Frequency (years<sup>-1</sup>)") +
  
  coord_cartesian(xlim = c(0, 0.1),
                  ylim = c(-0.1, 4)) +
  
  plot_theme_paper +
  theme(legend.position = "none",
        panel.grid.major = element_gridline) +
  
  ggtitle(NULL,"Periodic solution frequency")


p_attack_rate <- ggplot() +
  # geom_vline(xintercept = 0.0025,
  #            colour = "grey80", linewidth = 1.0, alpha = 0.5) +
  geom_line(aes(x = r, y = mean_inc * 365, colour = scenario, linewidth = scenario),
            data_mean_incidence) +
  
  xlab("Mean antibody decay rate <i>r</i>") +
  ylab("Attack rate") +
  
  # coord_cartesian(xlim = c(-0.0002, 0.0075),
  #                 ylim = c(-0.3, 0.007 * 365),
  #                 expand = FALSE) +
  
  ggokabeito::scale_colour_okabe_ito(name = "Scenario", order = c(5, 9)) +
  
  scale_linewidth_manual(values = c("Multiplicative boosting" = 0.7, "Additive boosting" = 0.7, "No boosting" = 0.9)) +
  
  
  plot_theme_paper +
  theme(legend.position = "none",
        panel.grid.major = element_gridline) +
  
  ggtitle(NULL, "Yearly infection attack rate at solution")
p_attack_rate

p_left <- (p_bifurcation / p_period / p_attack_rate) +
  plot_layout(tag_level = "new")

p_right <- (p_example_prevalence / p_example_mean_antibodies) +
  plot_layout(tag_level = "keep")

(p_left | p_right) / guide_area() + 
  plot_layout(guides = "collect", heights = c(1, 0.1)) +
  plot_annotation(tag_levels = list(c("A", "B", " "), c("i.", "ii.", "iii.")), tag_sep = " ") &
  theme(plot.tag = element_text(face = "bold", size = 15))

ggsave(
  "results/results_boosting.pdf",
  device = cairo_pdf,
  width = 14, height = 10,
  bg = "white"
)

