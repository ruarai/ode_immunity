

library(tidyverse)
library(rhdf5)
library(patchwork)


source("R/plot_theme.R")

rs <- c(0.003, 0.015, 0.03)

x_r <- h5read("data/bifurcations.jld2", "x_r")
y_I_sol <- h5read("data/bifurcations.jld2", "y_I_sol")
y_inc_sol <- h5read("data/bifurcations.jld2", "y_inc_sol")
y_fixed_I <- h5read("data/bifurcations.jld2", "y_fixed_I")
y_means <- h5read("data/bifurcations.jld2", "y_means")
y_real_eigs <- h5read("data/bifurcations.jld2", "y_real_eigs")
y_imag_eigs <- h5read("data/bifurcations.jld2", "y_imag_eigs")

period <- h5read("data/bifurcations.jld2", "period")
# attack_rate <- h5read("data/bifurcations.jld2", "attack_rate")

hide_x_axis <- list(
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
)

data_I_sol <- y_I_sol %>%
  reshape2::melt(varnames = c("r", "t"), value.name = "prev") %>% 
  mutate(r = x_r[r])

data_y_means <- y_means %>%
  reshape2::melt(varnames = c("r", "t"), value.name = "mean") %>% 
  mutate(r = x_r[r])

data_inc_sol <- y_inc_sol %>%
  reshape2::melt(varnames = c("r", "t"), value.name = "inc") %>% 
  mutate(r = x_r[r])

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


data_period <- period %>%
  reshape2::melt(varnames = c("r", "name"), value.name = "value") %>% 
  mutate(name = c("period", "period_sd", "period_n")[name],
         r = x_r[r]) %>%
  
  pivot_wider() %>% 
  
  filter(r < bifur_points$r[2], period_sd < 1, period_n > 1)

min_periods <- data_period %>% slice(1) %>%
  select(r_min = r)


ggplot() +
  geom_line(aes(x = r, y = period / 365),
            linewidth = 1.0, colour = colour_A,
            data_period) +
  
  
  geom_point(aes(x = r, y = period / 365), 
             colour = "black",
             data_period %>% filter(r == max(r)),
             size = 3) +
  
  geom_point(aes(x = r, y = period / 365), 
             data_period %>% filter(r == max(r)),
             size = 1.5, colour = "white") +
  
  xlab("Effective antibody decay rate <i>r</i> (days<sup>-1</sup>)") +
  ylab("Periods (years)") +
  
  coord_cartesian(xlim = c(0, 0.03), ylim = c(0, 4)) +
  
  scale_y_continuous(labels = scales::label_comma()) +
  # sec_x_axis +
  
  plot_theme_paper +
  
  theme(legend.position = "none",
        axis.text.x.top = element_text(margin = margin(b = 0.3, unit = "cm")),
        panel.grid.major.y = element_gridline,
        plot.subtitle = element_markdown()) +
  
  ggtitle(NULL, "Period of the periodic solution")

