
library(tidyverse)

library(rhdf5)


source("../ode_immunity_multi/R/plot_theme.R")


x_vals <- h5read("data/paper/period_over_grid.jld2", "x_vals")
y_inf_summary <- h5read("data/paper/period_over_grid.jld2", "y_inf_summary")
y_period <- h5read("data/paper/period_over_grid.jld2", "y_period")
y_attack_rate <- h5read("data/paper/period_over_grid.jld2", "y_attack_rate")

plot_data <- tibble(
  eta = x_vals[1, ], rho = x_vals[2, ],
  inf_min = y_inf_summary[, 1], inf_max = y_inf_summary[, 2],  inf_mean = y_inf_summary[, 3],
  inc_min = y_inf_summary[, 4], inc_max = y_inf_summary[, 5],  inc_mean = y_inf_summary[, 6],
  period = y_period[,1], period_sd = y_period[,2], period_n = y_period[,3],
  attack_rate = y_attack_rate
) %>%
  mutate(inf_diff = inf_max - inf_min,
         periodic = (period_sd < 1) & (period_n >= 5),
         quasiperiodic = (period_sd >= 1) & (period_n >= 5))




plot_data_periodic <- plot_data %>%
  filter(eta > 0) %>% 
  filter(periodic) %>% 
  mutate(period = period / 365,
         period = pmin(period, 8),
         period = factor(round(period)))

plot_data_quasiperiodic <- plot_data %>% filter(quasiperiodic)

plot_data_eta_zero <- plot_data %>% filter(eta == 0)

bifur_zero <- plot_data_eta_zero %>% filter(inf_diff < 1e-3) %>% pull(rho) %>% head(1)

year_stops <- c(1/2, 2/3, 1, 3/2, 2, 3, 4)
year_marks <- approxfun(plot_data_eta_zero$period, plot_data_eta_zero$rho)(365 * year_stops)
plot_data_year_marks <- tibble(rho_0 = year_marks, year = year_stops) %>%
  mutate(year_label = str_c(scales::label_comma()(year), " yr"))

plot_annotations <- list(
  geom_point(aes(x = -0.01, y = rho_0), plot_data_year_marks, pch = "-", size = 6),
  
  annotate("linerange", x = -0.01, ymin = bifur_zero, ymax = 0.005),
  annotate("point", x = -0.01, y = bifur_zero, pch = "-", size = 6),
  geom_text(aes(x = -0.07, y = rho_0, label = year_label), hjust = 0, plot_data_year_marks),
  annotate("text", x = -0.07, y = 0.00465, label = "Fixed\npoint", hjust = 0),
  geom_linerange(aes(ymin = 0, ymax = 0.005, x = eta), tibble(eta = seq(0, 0.5, by = 0.1)), colour = "black", alpha = 0.25, linetype = "88"),
  geom_linerange(aes(xmin = 0, xmax = 0.5, y = rho), tibble(rho = seq(0, 0.005, by = 0.001)), colour = "black", alpha = 0.25, linetype = "88")
)


period_cols <- viridis::inferno(n = 8, direction = -1, begin = 0.1)
p_period <- ggplot() +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.005, fill = "white") +
  geom_tile(aes(x = eta, y = rho, fill = period),
            plot_data_periodic) +
  # geom_tile(aes(x = eta, y = rho, fill = factor(9)),
  #           alpha = 0,
  #           plot_data_quasiperiodic) +
  geom_tile(aes(x = eta, y = rho, fill = factor(4.5),
            plot_data_quasiperiodic) +
  
  plot_annotations +
  
  
  scale_fill_manual(name = "Period",
                    values = c(period_cols[1:4], "#2260BE", period_cols[5:8]) %>% `names<-`(c(1:4, "4.5", 5:8)),
                    
                    labels = c("1 yr", str_c(2:4, "yrs"), "Quasiperiodic", str_c(5:7, "yrs"), "≥8 yrs") %>% `names<-`(c(1:4, "4.5", 5:8)),
                    breaks = c(1:4, 4.5, 5:8)
                    ) +

  
  coord_fixed(ratio = 100) +
  xlab("Seasonality constant <i>η</i>") + ylab("Waning constant <i>ρ</i>") +
  guides(fill = guide_legend(nrow = 2, ncol = 5)) +
  
  plot_theme_paper +
  theme(legend.position = "bottom", legend.byrow = TRUE)

# p_period

p_attack_rate <- ggplot()  +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.005, fill = "white")+
  geom_tile(aes(x = eta, y = rho, fill = inc_mean * 365),
            plot_data) +
  
  plot_annotations + 
  
  scale_fill_stepsn(
    colours = colorspace::sequential_hcl(n = 20, h = c(300, 75), c = c(40, NA, 95), l = c(15, 90), power = c(1, 1.1)),
    name = "Yearly infection\nattack rate",
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, 0.05),
    labels = c("0", "", "0.1", "", "0.2", "", "0.3", "", "0.4", "", "0.5")
  ) +
  
  coord_fixed(ratio = 100) +
  xlab("Seasonality constant <i>η</i>") + ylab("Waning constant <i>ρ</i>")  +
  
  plot_theme_paper +
  guides(fill = guide_colourbar(barwidth = 15)) +
  theme(legend.position = "bottom")

p_max <- ggplot()  +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.005, fill = "white")+
  geom_tile(aes(x = eta, y = rho, fill = inf_max),
            plot_data) +
  
  plot_annotations + 
  
  scale_fill_viridis_b(
    name = "Peak\ninfection\nprevalence", option = "mako",
    limits = c(0, 0.25),
    breaks = seq(0, 0.25, 0.025),
    labels = c("0", "", "0.05", "", "0.1", "", "0.15", "", "0.2", "", "0.25")
  ) +
  coord_fixed(ratio = 100) +
  xlab("Seasonality constant <i>η</i>") + ylab("Waning constant <i>ρ</i>")  +
  
  plot_theme_paper +
  guides(fill = guide_colourbar(barwidth = 15)) +
  theme(legend.position = "bottom")

# p_max


p_min <- ggplot()  +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.005, fill = "white")+
  geom_tile(aes(x = eta, y = rho, fill = pmax(-14, log10(inf_min))),
            plot_data) +
  
  plot_annotations + 
  
  scale_fill_stepsn(
    colours = rev(colorspace::diverging_hcl(n = 20, h = c(240, 15), c = c(60, 80), l = c(75, 5), power = c(1.2, 1.5))),
    name = "Minimum\ninfection\nprevalence (log10)",
    limits = c(-15, -1),
    breaks = seq(-15, -1, 0.5),
    labels = c("-15", "", "", "", "-13", "", "", "", "-11", "", "", "", "-9", 
               "", "", "", "-7", "", "", "", "-5", "", "", "", "-3", "", "", 
               "", "-1")
  ) +
  
  coord_fixed(ratio = 100) +
  xlab("Seasonality constant <i>η</i>") + ylab("Waning constant <i>ρ</i>")  +
  
  plot_theme_paper +
  guides(fill = guide_colourbar(barwidth = 15)) +
  theme(legend.position = "bottom")

# p_min



p_full <- (p_period | p_attack_rate) / (p_max | p_min) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))







ggsave(
  "results/results_supp_seasonality.png",
  p_full,
  device = png,
  width = 14, height = 14,
  bg = "white"
)