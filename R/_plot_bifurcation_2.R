

library(tidyverse)
library(rhdf5)
library(patchwork)


source("../ode_immunity_multi/R/plot_theme.R")

rhos <- c(0.001, 0.003, 0.005)



rho_to_halflife <- function(x) {1 / (8 * x)}

sec_x_axis <- list(
  scale_x_continuous(
    labels = NULL,
    sec.axis = sec_axis(
      identity,
      name = "Antibody half-life (days)",
      labels = rho_to_halflife,
      breaks = rho_to_halflife(c(5, 10, 15, 20, 30, 50, 100, 400))
  ))
)


x_rho <- h5read("data/paper/bifurcations_w_boost.jld2", "x_rho")
y_I_sol <- h5read("data/paper/bifurcations_w_boost.jld2", "y_I_sol")
y_inc_sol <- h5read("data/paper/bifurcations_w_boost.jld2", "y_inc_sol")
y_fixed_I <- h5read("data/paper/bifurcations_w_boost.jld2", "y_fixed_I")

hide_x_axis <- list(
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
)

days_burn_in <- 20000


data_I_sol <- y_I_sol %>%
  reshape2::melt(varnames = c("rho", "scenario", "t"), value.name = "prev") %>% 
  mutate(rho = x_rho[rho]) %>%
  filter(scenario == 1)


data_inc_sol <- y_inc_sol %>%
  reshape2::melt(varnames = c("rho", "scenario", "t"), value.name = "inc") %>% 
  mutate(rho = x_rho[rho]) %>%
  filter(scenario == 1)

data_mean_incidence <- data_inc_sol %>%
  filter(t > 10000) %>% 
  group_by(rho) %>%
  summarise(mean_inc = mean(inc))

maxmins <- data_I_sol %>%
  filter(t > days_burn_in, rho > 0.0003) %>% 
  group_by(rho, scenario) %>%
  summarise(max = max(prev), min = min(prev))

data_fixed <- y_fixed_I %>%
  reshape2::melt(varnames = c("rho", "scenario"), value.name = "prev") %>% 
  mutate(rho = x_rho[rho]) %>%
  filter(scenario == 1)

bifur_points <- maxmins %>% 
  left_join(data_fixed) %>%
  mutate(diff = max - prev) %>%
  ungroup() %>% 
  group_by(scenario) %>% 
  filter(diff > 1e-4) %>%
  slice(n()) %>%
  select(rho, scenario, prev = prev)


p_bifurcation <- ggplot() +
  geom_vline(aes(xintercept = rho),
             tibble(rho = rhos),
             colour = "grey80", linewidth = 1.0, alpha = 0.3) +
  
  geom_line(aes(x = rho, y = max),
            linewidth = 1.0,
            colour = colour_A,
            maxmins %>% filter(rho <= bifur_points$rho[1])) +
  geom_line(aes(x = rho, y = prev),
            linewidth = 1.0,
            data_fixed %>% filter(rho >= bifur_points$rho[1])) +
  geom_line(aes(x = rho, y = prev),
            linewidth = 1.0, linetype = "44",
            data_fixed) +
  
  
  geom_point(aes(x = rho, y = prev), 
             colour = "black",
             bifur_points,
             size = 3) +
  
  geom_point(aes(x = rho, y = prev), 
             bifur_points,
             size = 1.5, colour = "white") +
  
  xlab(NULL) +
  ylab("Prevalence") +
  
  coord_cartesian(xlim = c(-0.0002, 0.0075),
                  ylim = c(-0.005, 0.06),
                  expand = FALSE) +
  
  # sec_x_axis +
  
  plot_theme_paper +
  
  theme(legend.position = "none",
        axis.text.x.top = element_text(margin = margin(b = 0.3, unit = "cm")),
        panel.grid.major = element_gridline,
        plot.subtitle = element_markdown()) +
  
  ggtitle(NULL, "Bifurcation over <i>ρ</i>")

p_bifurcation

p_bifurcation_min <- ggplot() +
  geom_vline(aes(xintercept = rho),
             tibble(rho = rhos),
             colour = "grey80", linewidth = 1.0, alpha = 0.3) +
  geom_line(aes(x = rho, y = min),
            linewidth = 1.0,
            colour = colour_C,
            maxmins %>% filter(rho <= bifur_points$rho[1])) +
  geom_line(aes(x = rho, y = prev),
            linewidth = 1.0,
            data_fixed %>% filter(rho >= bifur_points$rho[1])) +
  geom_line(aes(x = rho, y = prev),
            linewidth = 1.0, linetype = "44",
            data_fixed) +
  
  
  geom_point(aes(x = rho, y = prev), 
             colour = "black",
             bifur_points,
             size = 3) +
  
  geom_point(aes(x = rho, y = prev), 
             bifur_points,
             size = 1.5, colour = "white") +
  
  xlab("Waning constant <i>ρ</i>") +
  ylab("Prevalence") +
  
  coord_cartesian(xlim = c(-0.0002, 0.0075),
                  ylim = c(1e-10, 1.0),
                  expand = FALSE) +
  
  scale_y_log10(labels = scales::label_log(),
                breaks = scales::breaks_log(n = 5)) +
  
  plot_theme_paper +
  theme(legend.position = "none",
        panel.grid.major = element_gridline,
        plot.subtitle = element_markdown()) +
  
  ggtitle(NULL, "Minimum infection prevalence across solution")

p_bifurcation_min


period <- h5read("data/paper/bifurcations_w_boost.jld2", "period")
attack_rate <- h5read("data/paper/bifurcations_w_boost.jld2", "attack_rate")

data_period <- period %>%
  reshape2::melt(varnames = c("rho", "scenario", "name"), value.name = "value") %>% 
  mutate(name = c("period", "period_sd", "period_n")[name],
         rho = x_rho[rho]) %>%
  filter(scenario == 1) %>%
  
  pivot_wider() %>% 
  
  left_join(bifur_points %>% select(rho_bifur = rho, scenario)) %>%
  filter(rho < rho_bifur, period_sd < 1, period_n > 1)

min_periods <- data_period %>% group_by(scenario) %>% slice(1) %>%
  select(rho_min = rho, scenario)

data_attack_rate <- attack_rate %>%
  reshape2::melt(varnames = c("rho", "scenario"), value.name = "attack_rate") %>% 
  mutate(rho = x_rho[rho]) %>%
  filter(scenario == 1) %>%
  
  left_join(bifur_points %>% select(rho_bifur = rho, scenario)) %>%
  left_join(min_periods) %>% 
  filter(rho < rho_bifur, rho >= rho_min) %>%
  left_join(data_period)


p_period <- ggplot() +
  geom_vline(aes(xintercept = rho),
             tibble(rho = rhos),
             colour = "grey80", linewidth = 1.0, alpha = 0.3) +
  geom_line(aes(x = rho, y = 365 / period),
            linewidth = 1.0,
            data_period) +
  
  xlab(NULL) +
  ylab("Frequency (years<sup>-1</sup>)") +
  
  coord_cartesian(xlim = c(-0.0002, 0.0075),
                  ylim = c(-0.1, 2.5),
                  expand = FALSE) +
  
  scale_y_continuous(labels = scales::label_comma()) +
  # sec_x_axis +
  
  plot_theme_paper +
  
  theme(legend.position = "none",
        axis.text.x.top = element_text(margin = margin(b = 0.3, unit = "cm")),
        panel.grid.major = element_gridline) +
  
  ggtitle(NULL,"Periodic solution frequency")

p_period

p_attack_rate <- ggplot() +
  geom_vline(aes(xintercept = rho),
             tibble(rho = rhos),
             colour = "grey80", linewidth = 1.0, alpha = 0.3) +
  geom_line(aes(x = rho, y = mean_inc * 365),
            linewidth = 1.0,
            data_mean_incidence) +
  
  xlab("Waning constant <i>ρ</i>") +
  ylab("Attack rate") +
  
  coord_cartesian(xlim = c(-0.0002, 0.0075),
                  ylim = c(-0.3, 0.007 * 365),
                  expand = FALSE) +
  
  
  plot_theme_paper +
  theme(legend.position = "none",
        panel.grid.major = element_gridline) +
  
  ggtitle(NULL, "Yearly infection attack rate")

plot_data_ex_fixed_points <- data_fixed %>%
  filter(rho %in% rhos) %>%
  mutate(stable = rho >= bifur_points$rho[[1]]) %>% 
  mutate(rho_label = str_c("<i>ρ </i>  = ", rho))


p_examples <- data_I_sol %>%
  filter(rho %in% rhos, t < 4000) %>% 
  mutate(rho_label = str_c("<i>ρ </i>  = ", rho)) %>% 
  ggplot() +
  
  geom_line(aes(x = t, y = prev),
            linewidth = 0.7) +

  geom_hline(aes(yintercept = prev, linetype = stable),
             plot_data_ex_fixed_points,
             linewidth = 0.7, colour = "black") +
  
  scale_x_continuous(breaks = scales::breaks_extended(),
                     labels = scales::label_comma()) +
  
  scale_y_continuous(breaks = c(0.0, 0.05, 0.1)) +
  
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "22")) +
  
  facet_wrap(~rho_label, ncol = 3, scales = "free") +
  
  xlab("Time (days)") + ylab("Prevalence") +
  
  coord_cartesian(xlim = c(0, 3400),
                  ylim = c(0, 0.08)) +
  
  plot_theme_paper +
  theme(strip.text = element_markdown(),
        legend.position = "none",
        panel.grid.major.x = element_gridline)


p_top <- (
  (p_bifurcation / p_bifurcation_min) |
  (p_period / p_attack_rate)
) +
  plot_layout(tag_level = "new")

p_top

(p_top / p_examples) +
  plot_layout(heights = c(2, 1)) + 
  plot_annotation(tag_levels = list("A", c("i.", "ii.", "iii.", "iv.")), tag_sep = " ") &
  theme(plot.tag = element_text(face = "bold", size = 15))




ggsave(
  "results/results_bifurcation.pdf",
  device = cairo_pdf,
  width = 14, height = 9,
  bg = "white"
)




data_I_sol %>%
  filter(t > 20000, rho > 0.0003) %>% 
  group_by(rho, scenario) %>%
  mutate(inc = prev)
