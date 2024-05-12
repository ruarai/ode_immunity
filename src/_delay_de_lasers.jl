

n_days = 1000
tau = 50
τ = Normal(tau, 0.001)


beta = 0.75
gamma = 0.5
R0 = 0.75 / 0.5

t_days, S_t, I_t, R_t, dI_t = sim_dde(;
    R0 = R0,
    γ = gamma,
    τ = τ,
    dt = 0.01,
    n_days, 
    n_inf_0 = 0.001
)


plot(t_days, I_t)

I_c = ((1 / beta) * (R0 - 1)) / (1 - exp(-))

eps_a = beta * I_c

η = a - b


eps = sqrt(beta / (gamma ^ 2))

x_f = (2 * tau -)