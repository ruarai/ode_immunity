


data <- read_csv("antibody_binding_neutralisation.csv")


data_A <- data %>%
  filter(mAbs %in% c("P2C-1F11", "P2B-2F6", "P2C-1A3"))

# 600kDa = 600,000 Da
approx_molar_mass <- 600 * 1e3 #g/mol = Da


ugmL_to_nM <- function(x_ugmL) {
  x_gL <- x_ugmL * 1e-6 * 1e3
  
  (x_gL / approx_molar_mass) * 1e9
}

data_binding <- data_A %>%
  expand_grid(c_nM = 10^seq(-5,5,by = 0.01)) %>%
  mutate(binding = (c_nM / RBD_KD) / (1 + c_nM / RBD_KD),
         c_M = c_nM * 1e-9)

hill_fn <- function(x, h, ic50) {
  return((x ^ h) / (ic50^h + x^h) )
}

get_h <- function(ic50, ic80) {
  hill_err <- function(h) {
    abs(hill_fn(ic50, h, ic50) - 0.5) + abs(hill_fn(ic80, h, ic50) - 0.8)
  }
  
  optimize(hill_err, c(0.1, 10))$minimum
}

data_neut <- data_A %>%
  rowwise() %>%
  mutate(h = get_h(ic50_livevirus, ic80_livevirus)) %>%
  
  expand_grid(c_ugmL = 10^seq(-6,10,by = 0.01)) %>%
  mutate(neut = hill_fn(c_ugmL, h, ic50_livevirus),
         c_nM = ugmL_to_nM(c_ugmL),
         c_M = c_nM * 1e-9,
         neut = if_else(mAbs == "P2C-1C10", neut * 0.85, neut),
         neut = if_else(mAbs == "P2A-1A10", neut * 0.8, neut),
         
         ic50 = ugmL_to_nM(ic50_livevirus) * 1e-9,
         ic80 = ugmL_to_nM(ic80_livevirus) * 1e-9) %>%
  
  select(c_M, mAbs, neut, ic50, ic80) %>% 
  
  left_join(data_binding %>% select(c_M, mAbs, binding), by = join_by(mAbs, closest(x$c_M < y$c_M)))

ggplot() +
  geom_line(aes(x = c_M.x, y = neut, colour = mAbs, linetype = "neut"),
            data_neut) +
  geom_line(aes(x = c_M.y, y = binding, colour = mAbs, linetype = "binding"),
            data_neut) +
  
  geom_point(aes(x = ic50, y = 0.5, colour = mAbs),
             data_neut %>% group_by(mAbs) %>% slice(1)) +
  
  geom_point(aes(x = ic80, y = 0.8, colour = mAbs),
             data_neut %>% group_by(mAbs) %>% slice(1)) +
  
  scale_x_log10()

data_neut %>%
  ggplot() +
  geom_line(aes(x = binding, y = neut, colour = mAbs)) +
  coord_fixed()

data_neut %>%
  ggplot() +
  geom_line(aes(x = binding, y = neut, colour = mAbs)) +
  
  geom_abline(slope = 1, intercept = 0) +
  
  scale_y_continuous(trans = "logit", breaks = seq(0.1, 0.9, by = 0.2)) +
  scale_x_continuous(trans = "logit", breaks = seq(0.1, 0.9, by = 0.2)) +
  
  coord_fixed(xlim = c(0.05, 0.95), ylim = c(0.05, 0.95))







