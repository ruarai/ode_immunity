

read_seasonality_data <- function(file) {
  
  x_vals <- h5read(file, "x_vals")
  y_inf_summary <- h5read(file, "y_inf_summary")
  y_period <- h5read(file, "y_period")
  y_seasonality <- h5read(file, "y_seasonality")
  
  
  plot_data <- tibble(
    eta = x_vals[1, ], r = x_vals[2, ],
    inf_min = y_inf_summary[, 1], inf_max = y_inf_summary[, 2],
    inf_mean = y_inf_summary[, 3], inf_chaos = y_inf_summary[, 4],
    
    inc_min = y_inf_summary[, 5], inc_max = y_inf_summary[, 6],  
    inc_mean = y_inf_summary[, 7], inc_chaos = y_inf_summary[ , 8],
    lyapunov = y_inf_summary[, 9],
    
    period = y_period[,1], period_sd = y_period[,2], period_n = y_period[,3],
    
    season_x = y_seasonality[,1], season_y = y_seasonality[,2]
  ) %>%
    mutate(
      inf_diff = inf_max - inf_min,
      period_error = pmin(abs(period %% 365 - 365), period %% 365),
      
      periodic = (period_error < 1) & (period_n > 1),
      chaotic = (inf_chaos > 0) & (!periodic),
      
      quasiperiodic = (period_n > 1) & (!periodic) & (eta > 0)
    )
  
  return(plot_data)
}

