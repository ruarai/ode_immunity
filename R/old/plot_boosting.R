
library(tidyverse)
library(rhdf5)


source("R/plot_theme.R")

track_state_no_boost <- h5read("data/paper/agent_based_no_boost.jld2", "track_state")
I_t_no_boost <- h5read("data/paper/agent_based_no_boost.jld2", "I_t")

track_state_with_boost <- h5read("data/paper/agent_based_with_boost.jld2", "track_state")
I_t_with_boost <- h5read("data/paper/agent_based_with_boost.jld2", "I_t")



make_abm_plot <- function(track_state, I_t, i_focus, title) {
  
  tbl_state <- track_state %>%
    `colnames<-`(1:ncol(track_state)) %>%
    `rownames<-`(1:nrow(track_state)) %>%
    as_tibble(rownames = "t") %>%
    pivot_longer(-t, names_to = "i", values_to = "state") %>%
    mutate(i = as.numeric(i), t = as.numeric(t)) %>% 
    
    mutate(class = if_else(state > 16, "I", "S"),
           strata = (state - 1) %% 16 + 1)  %>%
    mutate(corr = strata / 16 - 1 / 16)
  
  plot_track_data <- tbl_state %>%
    filter(i == i_focus, t < 500) %>%
    mutate(plot_group = data.table::rleid(class))
  
  plot_infection_events <- plot_track_data %>%
    filter(class == "S", lead(class) == "I")
  
  plot_recovery_transitions <- bind_rows(
    plot_track_data %>% filter(class == "I", lead(class) == "S") %>% mutate(t2 = t),
    plot_track_data %>% filter(class == "S", lag(class) == "I") %>% mutate(t2 = t - 1)
  ) %>%
    arrange(t2) %>%
    mutate(plot_group = data.table::rleid(t2))
  
  plot_I_t <- tibble(t = 1:length(I_t), I = I_t) %>%
    filter(t < 500)
  
  cowplot::plot_grid(
    ggplot(plot_I_t) +
      geom_line(aes(x = t, y = I)) +
      
      xlab(NULL) + ylab("Total Infected I") +
      
      scale_y_continuous(breaks = scales::breaks_extended(3),
                         labels = scales::label_comma()) +
      
      plot_theme_paper +
      ggtitle(title),
    ggplot() +
      geom_line(aes(x = t, y = corr, group = plot_group),
                plot_track_data %>% filter(class == "S")) +
      geom_line(aes(x = t, y = corr, group = plot_group),
                colour = colour_B,
                plot_track_data %>% filter(class == "I")) +
      geom_line(aes(x = t, y = corr, group = plot_group),
                plot_recovery_transitions,
                linetype = "dotted") +
      
      geom_point(aes(x = t, y = corr), colour = colour_B,
                 plot_infection_events) +
      
      coord_cartesian(ylim = c(0, 1)) +
      
      scale_y_continuous(breaks = c(0, 0.5, 1.0),
                         sec.axis = sec_axis(~ . * 16, name = expression("Strata"~italic(i)), breaks = c(3.5, 7.5, 11.5, 15.5), c(4, 8, 12, 16))) +
      
      xlab("Time (days)") + ylab("Correlate level _c_<sub>_i_</sub>") +
      
      plot_theme_paper,
    
    ncol = 1, align = "v", axis = "lr",
    rel_heights = c(1, 1)
  )
}

p_abm <- cowplot::plot_grid(
  
  make_plot(track_state_no_boost, I_t_no_boost, 43, "<b>C</b>"),
  make_plot(track_state_with_boost, I_t_with_boost, 62, "<b>D</b>"),
  
  ncol = 2
)
p_abm

x_rho <- h5read("data/paper/bifurcations.jld2", "x_rho")
y_I_sol <- h5read("data/paper/bifurcations.jld2", "y_I_sol")
y_fixed_I <- h5read("data/paper/bifurcations.jld2", "y_fixed_I")



x_rho_boost <- h5read("data/paper/bifurcations.jld2", "x_rho")
y_I_sol_boost <- h5read("data/paper/bifurcations.jld2", "y_I_sol")
y_fixed_I_boost <- h5read("data/paper/bifurcations.jld2", "y_fixed_I")

plot(y_I_sol_boost[x_rho_boost == 0.006,28000:30000])
plot(y_I_sol[x_rho_boost == 0.006,28000:30000])



maximums <- apply(y_I_sol[,28000:32000], 1, FUN = max)
maximums_boost <- apply(y_I_sol_boost[,28000:32000], 1, FUN = max)



plot_data_bifurcation <- tibble(
  rho = x_rho,
  fixed_I = y_fixed_I,
  maximum = maximums
) 


plot_data_bifurcation_boost <- tibble(
  rho = x_rho_boost,
  fixed_I = y_fixed_I_boost,
  maximum = maximums_boost
) 


p_bifurcation <- ggplot() +
  
  geom_line(aes(x = rho, y = maximum),
            linewidth = 1.0,
            colour = ggokabeito::palette_okabe_ito(3), alpha = 0.5,
            plot_data_bifurcation_boost %>% filter(rho > 0.0002)) +
  geom_line(aes(x = rho, y = fixed_I),
            linewidth = 1.0, colour = ggokabeito::palette_okabe_ito(3),
            plot_data_bifurcation_boost %>% filter(rho > 0.0071)) +
  geom_line(aes(x = rho, y = fixed_I),
            linewidth = 1.0, linetype = "11",
            plot_data_bifurcation_boost %>% filter(rho <= 0.0071)) +
  
  geom_line(aes(x = rho, y = maximum),
            linewidth = 0.5, 
            plot_data_bifurcation %>% filter(rho > 0.0002, rho <= 0.0071)) +
  geom_line(aes(x = rho, y = fixed_I),
            linewidth = 0.5,
            plot_data_bifurcation %>% filter(rho > 0.0071)) +
  geom_line(aes(x = rho, y = fixed_I),
            linewidth = 0.5, linetype = "11", alpha = 0.5,
            plot_data_bifurcation %>% filter(rho <= 0.0071)) +
  
  xlab("Decay rate λ") +
  ylab("Infected I") +
  
  coord_cartesian(xlim = c(0, 0.015),
                  ylim = c(0, 0.07)) +
  
  plot_theme_paper +
  
  ggtitle("<b>B</b>")



p_bifurcation


sol_I_no_boost <- h5read("data/paper/basic.jld2", "sol_I")
sol_I_boosting <- h5read("data/paper/basic_boosting.jld2", "sol_I")


c_levels <- 0:(ncol(sol_I) - 1) / 16


sol_no_boost <- data.table::data.table(t = 1:nrow(sol_I_no_boost), c = rep(c_levels, each = nrow(sol_I_no_boost)), I = c(sol_I_no_boost)) %>%
  as_tibble() %>%
  mutate(group = "No boosting")

sol_boosting <- data.table::data.table(t = 1:nrow(sol_I_boosting), c = rep(c_levels, each = nrow(sol_I_boosting)), I = c(sol_I_boosting)) %>%
  as_tibble() %>%
  mutate(group = "Boosting")



plot_data <- bind_rows(
  sol_no_boost, sol_boosting
) %>%
  summarise(I = sum(I), .by = c(t, group)) %>%
  mutate(group = factor(group, c("Boosting", "No boosting")))

 p_basic <- ggplot() +
   geom_line(aes(x = t, y = I),
             plot_data %>% filter(group == "Boosting"),
             colour = ggokabeito::palette_okabe_ito(3),
             linewidth = 0.7) +
   geom_line(aes(x = t, y = I),
             plot_data %>% filter(group == "No boosting"),
             linewidth = 0.7) +
  
   coord_cartesian(xlim = c(0, 4000),
                   ylim = c(0, 0.07)) +
  
  scale_x_continuous(breaks = seq(0, 365 * 12, by = 365),
                     labels = 0:12) +
   
  
  xlab("Time (years)") + ylab("Infected I") +
  
  plot_theme_paper +
  
  ggtitle("<b>A</b>") +
   
  theme(legend.position = "bottom")
p_basic


p_legend <- ggplot() +
  geom_line(aes(x = 1, y = 1, colour = group),
            tibble(group = c("With boosting",
                             "Without boosting"))) +
  
  ggokabeito::scale_colour_okabe_ito(name = "Model", order = c(3, 9)) +
  
  plot_theme_paper +
  guides(colour = guide_legend(keywidth = 3, override.aes = list(linewidth = 0.7))) +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        legend.text = element_markdown())
p_legend


cowplot::plot_grid(
  
  cowplot::plot_grid(
    p_basic, p_bifurcation,
    ncol = 2,
    align = "h", axis = "tb"
  ),
  cowplot::get_legend(p_legend),
  ncol = 1, rel_heights = c(10, 2)
)


ggsave(
  "results/results_SkIk.png",
  scale = 10 / 16,
  device = png,
  dpi = 300,
  width = 16, height = 8,
  bg = "white"
)

    


