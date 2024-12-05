

global const n_compartments = 3

global const c_sus = 1
global const c_inf = 2
global const c_count = 3

global const n_inf_0 = 0.001


global const n_days_burn_in = 30000
global const n_days = 100000


global const baseline_k = 32
global const baseline_C = 8.0

global const baseline_R = 1.5
global const baseline_gamma = 0.25
global const baseline_beta = baseline_R * baseline_gamma

global const baseline_r = 0.05
global const baseline_b = 10^3
global const baseline_h = 3

global const baseline_c_jump_dist = Normal(6, 0.5)


global const periodic_ϵ = 10^-6