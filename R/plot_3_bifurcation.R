

library(tidyverse)
library(rhdf5)
library(patchwork)


source("R/plot_theme.R")

rs <- c(0.015, 0.03)

x_r <- h5read("data/paper/bifurcations.jld2", "x_r")
y_I_sol <- h5read("data/paper/bifurcations.jld2", "y_I_sol")
y_inc_sol <- h5read("data/paper/bifurcations.jld2", "y_inc_sol")
y_fixed_I <- h5read("data/paper/bifurcations.jld2", "y_fixed_I")
y_means <- h5read("data/paper/bifurcations.jld2", "y_means")
y_real_eigs <- h5read("data/paper/bifurcations.jld2", "y_real_eigs")
y_imag_eigs <- h5read("data/paper/bifurcations.jld2", "y_imag_eigs")

hide_x_axis <- list(
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
)

days_burn_in <- 30000


data_I_sol <- y_I_sol %>%
  reshape2::melt(varnames = c("r", "t"), value.name = "prev") %>% 
  mutate(r = x_r[r])

data_y_means <- y_means %>%
  reshape2::melt(varnames = c("r", "t"), value.name = "mean") %>% 
  mutate(r = x_r[r])

data_inc_sol <- y_inc_sol %>%
  reshape2::melt(varnames = c("r", "t"), value.name = "inc") %>% 
  mutate(r = x_r[r])

data_mean_incidence <- data_inc_sol %>%
  filter(t > 10000) %>% 
  group_by(r) %>%
  summarise(mean_inc = mean(inc))

data_eigs <- left_join(
  reshape2::melt(y_real_eigs, varnames = c("i", "r"), value.name = "e"),
  reshape2::melt(y_imag_eigs, varnames = c("i", "r"), value.name = "e_i")
) %>% 
  mutate(r = x_r[r]) %>%
  as_tibble()

eigs_stable <- data_eigs %>%
  group_by(r) %>%
  summarise(stable = all(e < 1e-10),
            maximum = max(e))

data_fixed <- y_fixed_I %>%
  reshape2::melt(varnames = c("r"), value.name = "prev") %>% 
  mutate(r = x_r[r]) %>%
  left_join(eigs_stable) %>%
  mutate(group = data.table::rleid(stable))

bifur_points <- data_fixed %>% 
  filter((stable & lead(!stable)) | (lag(!stable) & stable))

maxmins <- data_I_sol %>%
  filter(t > days_burn_in, r > 0.0003) %>% 
  group_by(r) %>%
  summarise(max = max(prev), min = min(prev))


p_bifurcation_min <- ggplot() +
  geom_vline(aes(xintercept = r),
             tibble(r = rs),
             colour = "grey80", linewidth = 1.0, alpha = 0.3) +
  geom_line(aes(x = r, y = max),
            colour = colour_A,
            linewidth = 1.0,
            maxmins %>% filter(r <= bifur_points$r[2])) +
  geom_line(aes(x = r, y = min),
            colour = colour_A,
            linewidth = 1.0,
            maxmins %>% filter(r <= bifur_points$r[2], min > 0)) +
  geom_line(aes(x = r, y = prev, group = group),
            linewidth = 1.0,
            data_fixed %>% filter(stable)) +
  geom_line(aes(x = r, y = prev),
            linewidth = 1.0, linetype = "44",
            data_fixed %>% filter(!stable)) +
  
  
  geom_point(aes(x = r, y = prev), 
             colour = "black",
             bifur_points,
             size = 3) +
  
  geom_point(aes(x = r, y = prev), 
             bifur_points,
             size = 1.5, colour = "white") +
  
  xlab("Antibody decay rate <i>r</i>") +
  ylab("Infection prevalence") +
  
  coord_cartesian(ylim = c(1e-10, 2.0)) +
  
  scale_y_log10(labels = scales::label_log(),
                breaks = 10^c(-10, -8, -6, -4, -2, 0)) +
  
  plot_theme_paper +
  theme(legend.position = "none",
        panel.grid.major = element_gridline,
        plot.subtitle = element_markdown()) +
  
  ggtitle(NULL, "<b>A</b> — Bifurcation over antibody<br> decay rate <i>r</i>")

p_bifurcation_min


period <- h5read("data/paper/bifurcations.jld2", "period")
# attack_rate <- h5read("data/paper/bifurcations.jld2", "attack_rate")

data_period <- period %>%
  reshape2::melt(varnames = c("r", "name"), value.name = "value") %>% 
  mutate(name = c("period", "period_sd", "period_n")[name],
         r = x_r[r]) %>%
  
  pivot_wider() %>% 
  
  filter(r < bifur_points$r[2], period_sd < 1, period_n > 1)

min_periods <- data_period %>% slice(1) %>%
  select(r_min = r)


p_period <- ggplot() +
  geom_vline(aes(xintercept = r),
             tibble(r = rs),
             colour = "grey80", linewidth = 1.0, alpha = 0.3) +
  geom_line(aes(x = r, y = 365 / period),
            linewidth = 1.0,
            data_period) +
  
  
  geom_point(aes(x = r, y = 365 / period), 
             colour = "black",
             data_period %>% filter(r == max(r)),
             size = 3) +
  
  geom_point(aes(x = r, y = 365 / period), 
             data_period %>% filter(r == max(r)),
             size = 1.5, colour = "white") +
  
  xlab("Antibody decay rate <i>r</i>") +
  ylab("Frequency (years<sup>-1</sup>)") +
  
  coord_cartesian(xlim = c(0, 0.05), ylim = c(0, 4)) +
  
  scale_y_continuous(labels = scales::label_comma()) +
  # sec_x_axis +
  
  plot_theme_paper +
  
  theme(legend.position = "none",
        axis.text.x.top = element_text(margin = margin(b = 0.3, unit = "cm")),
        panel.grid.major = element_gridline,
        plot.subtitle = element_markdown()) +
  
  ggtitle(NULL,"<b>B</b> — Periodic solution frequency")

p_period

p_attack_rate <- ggplot() +
  geom_vline(aes(xintercept = r),
             tibble(r = rs),
             colour = "grey80", linewidth = 1.0, alpha = 0.3) +
  geom_line(aes(x = r, y = mean_inc * 365),
            linewidth = 1.0,
            data_mean_incidence) +
  
  xlab("Antibody decay rate <i>r</i>") +
  ylab("Infection incidence") +
  
  coord_cartesian(ylim = c(NA, 3.0)) +
  
  
  plot_theme_paper +
  theme(legend.position = "none",
        panel.grid.major = element_gridline,
        plot.subtitle = element_markdown()) +
  
  ggtitle(NULL, "<b>C</b> — Average annual infection incidence<br> at solution")


plot_data_ex_fixed_points <- data_fixed %>%
  filter(r %in% rs) %>%
  mutate(stable = r >= bifur_points$r[[2]]) %>% 
  mutate(r_label = str_c("Antibody decay rate <i>r </i>  = ", r))


p_examples <- data_I_sol %>%
  filter(r %in% rs, t < 4000) %>% 
  mutate(r_label = str_c("Antibody decay rate <i>r </i>  = ", r)) %>% 
  ggplot() +
  
  geom_line(aes(x = t, y = prev),
            linewidth = 0.7) +

  geom_hline(aes(yintercept = prev, linetype = stable),
             plot_data_ex_fixed_points,
             linewidth = 0.7, colour = "black") +
  
  scale_x_continuous(breaks = scales::breaks_extended(),
                     labels = scales::label_comma()) +
  
  scale_y_continuous(breaks = c(0.0, 0.05, 0.1)) +
  
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "44")) +
  
  facet_wrap(~r_label, ncol = 3, scales = "free") +
  
  xlab("Time (days)") + ylab("Infection prevalence") +
  
  coord_cartesian(xlim = c(0, 3400),
                  ylim = c(0, 0.08)) +
  
  plot_theme_paper +
  theme(strip.text = element_markdown(),
        plot.subtitle = element_markdown(),
        legend.position = "none",
        panel.grid.major.x = element_gridline) +
  
  ggtitle(NULL, "<b>D</b> — Exemplar infection prevalence")

p_examples


p_examples_antibodies <- data_y_means %>%
  filter(r %in% rs, t < 4000) %>% 
  mutate(r_label = str_c("Antibody decay rate <i>r </i>  = ", r)) %>% 
  ggplot() +
  
  geom_line(aes(x = t, y = mean),
            linewidth = 0.7) +
  
  scale_x_continuous(breaks = scales::breaks_extended(),
                     labels = scales::label_comma()) +
  
  scale_y_log10(breaks = scales::breaks_log(),
                labels = scales::label_log()) +
  
  facet_wrap(~r_label, ncol = 3, scales = "free") +
  
  xlab("Time (days)") + ylab("Concentration") +
  
  coord_cartesian(xlim = c(0, 3400),
                  ylim = c(10^2, 10^6.5)) +
  
  plot_theme_paper +
  theme(strip.text = element_markdown(),
        plot.subtitle = element_markdown(),
        legend.position = "none",
        panel.grid.major.x = element_gridline) +
  
  ggtitle(NULL, "<b>E</b> — Exemplar mean antibody concentration")


# p_top <- (
#   (p_bifurcation / p_bifurcation_min) |
#   (p_period / p_attack_rate)
# ) +
#   plot_layout(tag_level = "new")


p_top <- (p_bifurcation_min | p_period | p_attack_rate) +
  plot_layout(tag_level = "new")




# p_top

(p_top / p_examples) +
  plot_layout(heights = c(1.5, 1)) + 
  # plot_annotation(tag_levels = list("A", c("1", "2", "3")), tag_sep = " ") &
  theme(plot.tag = element_text(face = "bold", size = 15))




ggsave(
  "results/results_bifurcation.pdf",
  device = cairo_pdf,
  width = 14, height = 8,
  bg = "white"
)




p_supp_A <- ggplot() +
  geom_path(aes(x = e, y = e_i, group = i),
            linewidth = 0.5,
            data_eigs) +
  
  geom_hline(yintercept = 0, linetype = "44") +
  geom_vline(xintercept = 0, linetype = "44") +
  
  xlab("Re(λ<sub><i>i</i></sub>)") +
  ylab("Im(λ<sub><i>i</i></sub>)") +
  
  coord_cartesian(xlim = c(-0.01,0.01),
                  ylim = c(-0.1, 0.1)) +
  
  ggtitle(NULL, "<b>A</b> — Eigenvalues for varying<br>antibody decay rate <i>r</i>") +
  
  plot_theme_paper


p_supp_B <- ggplot() +
  geom_line(aes(x = r, y = maximum),
            linewidth = 0.8,
            eigs_stable) +
  
  geom_hline(yintercept = 0, linetype = "44") +
  
  ylab("max<sub><i>i</i></sub> Re(λ<sub><i>i</i></sub>)") +
  xlab("Antibody decay rate <i>r</i>") +
  
  coord_cartesian(ylim = c(-0.005, 0.005)) +
  
  ggtitle(NULL, "<b>B</b> — Maximum real part<br>of eigenvalue") +
  
  plot_theme_paper


p_supp_A | p_supp_B



ggsave(
  "results/results_supp_eigenvalues.pdf",
  device = cairo_pdf,
  width = 8, height = 4,
  bg = "white"
)
 