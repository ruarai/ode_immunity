

library(tidyverse)

source("R/plot_theme.R")

ili_data <- read_csv("data_old/flu-ili-byregion-fluseason.csv")


plot_data <- ili_data %>%
  mutate(date = mdy(weekending)) %>%
  group_by(date) %>%
  summarise(count = sum(Total_ILI)) %>%
  
  filter(date >= ymd("2010-06-01"),
         date < ymd("2015-06-01"))


ggplot() +
  geom_linerange(aes(x = date, ymin = 0, ymax = count),
                 linewidth = 1, colour = colour_A,
                 plot_data) +
  
  geom_vline(aes(xintercept = d), tibble(d = seq(ymd("2011-01-01"), ymd("2015-01-01"), by = "year")),
             linewidth = 0.8, linetype = "42") +
  
  scale_x_date(
    date_breaks = "year",
    labels = scales::label_date_short()
  ) +
  
  xlab("Date") + ylab("Count") +
  
  coord_cartesian(ylim = c(0, NA)) +
  
  plot_theme_paper +
  
  ggtitle("Counts of influenza-like illness")
