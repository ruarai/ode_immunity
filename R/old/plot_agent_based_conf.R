
library(tidyverse)
library(rhdf5)


source("R/plot_theme.R")

track_state <- h5read("data/paper/agent_based_transition.jld2", "track_state")
track_acq <- h5read("data/paper/agent_based_transition.jld2", "track_acq")
I_t <- h5read("data/paper/agent_based_transition.jld2", "I_t")


tbl_state <- track_state %>%
  reshape2::melt(varnames = c("t", "i"), value.name = "state") %>%
  as_tibble() %>%
  mutate(i = as.numeric(i), t = as.numeric(t)) %>% 
  
  mutate(class = if_else(state > k + 1, "I", "S"),
         strata = (state - 1) %% (k + 1) + 1)  %>%
  mutate(corr = strata / k)


tbl_acq <- track_acq %>%
  reshape2::melt(varnames = c("t", "i"), value.name = "acq") %>%
  as_tibble() %>%
  mutate(i = as.numeric(i), t = as.numeric(t))

k <- 128


tbl_acq %>%
  count(t, acq) %>%
  complete(t, acq, fill = list(n = 0)) %>% 
  mutate(n = pmin(n, 200)) %>% 
  ggplot() +
  geom_point(aes(x = t, y = acq, colour = n))


tbl_state %>%
  count(t, corr) %>%
  complete(t, corr, fill = list(n = 0)) %>% 
  mutate(n = pmin(n, 200)) %>% 
  ggplot() +
  geom_point(aes(x = t, y = corr, colour = n))


tbl_state_corr <- tbl_state %>%
  mutate(corr = 1 / (1 + exp(-40 * (corr - 0.5)))) %>% 
  filter(i > 20, t < 800)

plot_track_data <- tbl_state_corr %>%
  group_by(t) %>%
  summarise(corr = mean(corr))



plot_I_t <- tibble(t = 1:length(I_t), I = I_t) %>%
  filter(t < 800) %>%
  mutate(I = zoo::rollmean(I, k = 1, fill = NA))

cowplot::plot_grid(
  ggplot(plot_I_t) +
    geom_line(aes(x = t, y = I)) +
    
    scale_y_continuous(breaks = scales::breaks_extended(3),
                       labels = scales::label_comma()) +
    
    xlab(NULL) + ylab(NULL) +
    
    plot_theme_paper +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank()),
  
  ggplot() +
    geom_line(aes(x = t, y = corr),
              plot_track_data) +
    
    xlab(NULL) + ylab(NULL) +
    
    plot_theme_paper +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank()),
  
  ncol = 1, align = "v", axis = "lr",
  rel_heights = c(1, 1)
)



tbl_acq_mean <- tbl_acq %>%
  group_by(t) %>%
  summarise(mean = mean(acq))


tbl_acq %>%
  filter(i > 10, i < 100) %>%
  group_by(i) %>% 
  mutate(newinf = lag(acq) < acq,
         newinf_id = data.table::rleid(newinf),
         acq = acq + rnorm(1, mean = 0.01, sd = 0.01)) %>% 
  ggplot() +
  geom_line(aes(x = t, y = acq, group = interaction(i, newinf_id)),
            alpha = 0.25) +
  
  geom_line(aes(x = t, y = mean),
            tbl_acq_mean)




library(distributional)



hdr_to_tbl <- function(h) {
  lowers <- as.vector(vctrs::field(h, "lower")[[1]])
  uppers <- as.vector(vctrs::field(h, "upper")[[1]])
  
  min_len <- min(length(lowers), length(uppers))
  
  tibble(
    lower = lowers[1:min_len],
    upper = uppers[1:min_len],
    level = vctrs::field(h, "level")
  )
}




x_dist <- distributional::dist_sample(list(c(1,1,2)))

x_hdr <- hdr(x_dist, size = 90, n = 512)

hdr_to_tbl(x_hdr)

library(vctrs)

new_hdr <- function(lower = list_of(.ptype = double()),
                    upper = list_of(.ptype = double()),
                    size = double()) {
  lower <- as_list_of(lower)
  upper <- as_list_of(upper)
  vec_assert(lower, list_of(.ptype = double()))
  vec_assert(upper, list_of(.ptype = double()))
  vec_assert(size, double())
  if (any(size < 0 | size > 100, na.rm = TRUE))
    abort("'size' must be between [0, 100].")
  
  
  out <- vec_recycle_common(lower = lower, upper = upper)
  mapply(
    function(l,u) #if (any(u<l, na.rm = TRUE)) abort("`upper` can't be lower than `lower`."),
    l = out[["lower"]], u = out[["upper"]]
  )
  out[["level"]] <- vctrs::vec_recycle(size, vec_size(out[[1]]))
  
  vctrs::new_rcrd(out, class = "hdr")
}

assignInNamespace("new_hdr", new_hdr,ns="distributional")



tbl_interp <- tbl_state %>%
  filter(i < 5000) %>% 
  expand_grid(t_add = seq(0.0, 0.9, by = 0.1)) %>%
  mutate(t = t + t_add) %>%
  select(-t_add) %>% 
  arrange(t) %>%
  mutate(corr = if_else(floor(t) != t, NA_real_, corr)) %>%
  group_by(i) %>% 
  mutate(corr = zoo::na.approx(corr, na.rm = FALSE))


tbl_interp %>%
  count(t = floor(t),
        corr = floor(corr * 50) / 50) %>% 
  group_by(t) %>%
  mutate(n = n / sum(n)) %>% 
  ggplot() +
  geom_tile(aes(x = t, y = corr, fill = n))


hdr_list <- tbl_state %>%
  filter(i > 10, i <= 10000) %>%
  
  filter(t < 100) %>%
  group_by(t) %>%
  summarise(dist = dist_sample(list(corr))) %>%
  expand_grid(level = seq(20, 90, by = 10)) %>%
  group_by(t, level) %>% 
  mutate(hdr = hdr(dist, n = 100, size = level))


x <- pmap(list(hdr_list$t, hdr_list$hdr), function(t, hdr) hdr_to_tbl(hdr) %>% mutate(t = t) ) %>%
  bind_rows() %>%
  arrange(desc(level))  %>%
  mutate(level = factor(level))


ggplot(x) +
  geom_linerange(aes(x = t, ymin = lower, ymax = upper, colour = level),
                 size = 2,
                 alpha = 0.5) +
  
  scale_colour_brewer(palette = "Blues", direction = -1)


ggplot(tbl_state %>% filter(i < 100)) +
  geom_point(aes(x = t, y = corr),
             alpha = 0.1)










