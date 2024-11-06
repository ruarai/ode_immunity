
include("dependencies.jl")

using JLD2

k = 32
C = 8.0
R = 1.5
gamma = 0.25
beta = R * gamma
rho_0 = 0.1
b = 2^3
h = 8
c_jump_dist = Normal(2^6, 2^5)

model_params_0 = make_model_parameters(
    k = k, beta = beta, gamma = gamma, C = C, rho = rho_0,
    b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none"
)

ode_sparsity = ode_get_sparsity(model_params_0)

ϵ = 10^-6

Δt = 0.25
n_inf_0 = 0.01
burn_in_days = 20000
n_days = 40000

x_rho = collect(0.00120:0.00005:0.0018)

y_fixed_I = zeros(length(x_rho))
y_I_sol = zeros(length(x_rho), n_days)


period = zeros(length(x_rho), 3)
attack_rate = zeros(length(x_rho))

@showprogress Threads.@threads for i in eachindex(x_rho)
    model_params = make_model_parameters(
        k = k, beta = beta, gamma = gamma, C = C, rho = x_rho[i],
        b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none", eta = 0.2
    )

    # Calculate the fixed point/steady state solution
    sus_vec, inf_vec = get_steady_state(model_params)
    y_fixed_I[i] = sum(inf_vec)


    # Calculate a (not-necessarily-stable) solution of the ODE over time
    ode_solution = ode_solve(model_params, n_days, n_inf_0, ode_sparsity, saveat = Δt)

    for d in 1:n_days
        y_I_sol[i, d] = sum(ode_solution(d)[ode_ix(c_inf, 1:model_params.S, model_params.S)])
    end


    period_mean, period_sd, period_n = get_period(ode_solution, model_params, burn_in_days, n_days, Δt, ϵ)

    period[i, :] = [period_mean period_sd period_n]
    attack_rate[i] = get_periodic_attack_rate(ode_solution, model_params, burn_in_days, n_days, period_mean)


end


plot(x_rho, maximum(y_I_sol[:, burn_in_days:end], dims = 2))

plot(y_I_sol[:, (end - 365 * 10):end]')
plot(y_I_sol[7,(end - 365 * 10):end ])


plot(x_rho, period[:,1])
plot(x_rho, period[:,2])
plot(x_rho, attack_rate)