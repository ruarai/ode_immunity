library(tidyverse)
library(rhdf5)
library(patchwork)

source("R/plot_theme.R")
source("R/read_seasonality_data.R")

plot_data_raw <- read_seasonality_data("data/period_over_grid_supp.jld2")

plot_data <- plot_data_raw %>% 
  mutate(model_label = case_when(
    model_version == 1 ~ "<b>A</b>. Baseline",
    model_version == 2 ~ "<b>B</b>. Lower <i>h</i> (<i>h</i> = 4)",
    model_version == 3 ~ "<b>C</b>. Higher <i>h</i> (<i>h</i> = 12)",
    model_version == 4 ~ "<b>D</b>. Narrower jump dist. (<i>&sigma;</i> = 0.1)",
    model_version == 5 ~ "<b>E</b>. Wider jump dist. (<i>&sigma;</i> = 1.0)",
    model_version == 6 ~ "<b>F</b>. Lower <i>c</i><sub>mid</sub> (<i>c</i><sub>mid</sub> = 4)",
    model_version == 7 ~ "<b>G</b>. Lower <i>c</i><sub>mid</sub> (<i>c</i><sub>mid</sub> = 16)",
    model_version == 8 ~ "<b>H</b>. Higher <i>k</i> (<i>k</i> = 64)",
    model_version == 9 ~ "<b>I</b>. Higher <i>k</i> (<i>k</i> = 128)",
    model_version == 10 ~ "<b>J</b>. Including importations"
  ))


plot_data_periodic <- plot_data %>%
  filter(eta > 0) %>% 
  filter(periodic) %>% 
  mutate(period = period / 365,
         period = pmin(period, 8),
         period = factor(round(period)))


plot_data_chaotic <- plot_data %>% filter(chaotic)

plot_data_empty <- plot_data %>% filter(!periodic, !quasiperiodic, !chaotic)

plot_data_quasiperiodic <- plot_data %>% filter(quasiperiodic)


period_cols <- viridis::inferno(n = 8, direction = -1, begin = 0.1)

p_period <- ggplot() +
  
  geom_tile(aes(x = eta, y = r, fill = period),
            plot_data_periodic) +
  geom_tile(aes(x = eta, y = r, fill = factor(9)),
            plot_data_quasiperiodic) +
  geom_tile(aes(x = eta, y = r, fill = factor(10)),
            plot_data_chaotic) +
  geom_tile(aes(x = eta, y = r, fill = factor(11)),
            plot_data_empty) +
  
  scale_fill_manual(
    name = "Period",
    values = c(period_cols, "#0076BC", "#9ED7F3", "white") %>% 
      `names<-`(1:11),
    
    labels = c("1 year", str_c(2:7, " years"), "≥8 years", "Quasiperiodic", "Chaotic", "Unclassified") %>%
      `names<-`(1:11),
    
    breaks = c(9, 1:4, 10, 5:8, 11)
  ) +
  
  facet_wrap(~model_label) +
  
  
  coord_fixed(ratio = 16.66, ylim = c(0, 0.03)) +
  xlab("Seasonal forcing strength <i>η</i>") + ylab("Effective antibody decay rate <i>r</i>") +
  guides(fill = guide_legend(nrow = 3, ncol = 5),
         colour = guide_none()) +
  
  plot_theme_paper +
  theme(legend.position = "bottom", legend.byrow = TRUE,
        legend.key = element_rect(colour = "grey80", linewidth = 0.5),
        strip.text = element_markdown())

p_period

color_tbl <- read_delim("data/lipari.txt", col_names = c("R", "G", "B")) %>%
  mutate(c = rgb(R, G, B))

min_cols <- color_tbl$c[seq(1, nrow(color_tbl), length.out = 128)]

p_min <- ggplot(plot_data) +
  geom_tile(aes(x = eta, y = r, fill = log10(inf_min))) +
  
  scale_fill_stepsn(
    # colours = rev(colorspace::diverging_hcl(n = 20, h = c(240, 15), c = c(60, 80), l = c(75, 5), power = c(1.2, 1.5))),
    colours = min_cols,
    name = "Minimum\ninfection\nprevalence (log10)",
    limits = c(-15, -1),
    breaks = seq(-15, -1, 0.5),
    labels = c("<-15", "", "", "", "-13", "", "", "", "-11", "", "", "", "-9", 
               "", "", "", "-7", "", "", "", "-5", "", "", "", "-3", "", "", 
               "", "-1")
  ) +
  
  facet_wrap(~model_label) +
  
  plot_theme_paper +
  
  coord_fixed(ratio = 16.66, ylim = c(0, 0.03)) +
  xlab("Seasonal forcing strength <i>η</i>") + ylab("Effective antibody decay rate <i>r</i>")  +
  
  plot_theme_paper +
  guides(fill = guide_colourbar(barwidth = 15),
         colour = guide_none()) +
  theme(legend.position = "bottom",
        legend.title = element_text(margin = margin(r = 0.7, unit = "cm")),
        strip.text = element_markdown())

p_period / p_min


ggsave(
  "results/results_supp_seasonality_sensitivity.png",
  device = png, bg = "white",
  width = 14, height = 18
)






