#include("dependencies.jl")

k = 32
C = 8.0
R = 1.5
gamma = 0.25
beta = R * gamma
rho = 0.003
b = 2^3
h = 8
c_jump_dist = Normal(2^6, 2^5)

x_eta = 0:0.01:1.0


n_days = 40000
burn_in_days = 10000

p_max = 10000
x_period = 50:5:p_max
err_periods = zeros(length(x_period), length(x_eta))

function get_period_err(ode_sol, p, Δt, p_max, burn_in_days, model_params)
    y1 = ode_sol(burn_in_days:Δt:(burn_in_days + p_max))[1:(model_params.S * 2), :]'
    y2 = ode_sol((burn_in_days + p):Δt:(burn_in_days + p_max + p + (Δt * 0.4)))[1:(model_params.S * 2), :]'

    return sum(abs.(y1 .- y2)) 
end

for i_eta in eachindex(x_eta)
    println(i_eta)
    model_params = make_model_parameters(
        k = k, beta = beta, gamma = gamma, C = C, rho = rho,
        b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none", eta = x_eta[i_eta]
    )

    # ode_sparsity = ode_get_sparsity(model_params)

    n_inf_0 = 0.0001
    ode_solution = @time ode_solve(model_params, n_days, n_inf_0, ode_sparsity, saveat = 0.1)

    @showprogress Threads.@threads for i in eachindex(x_period)
        p = x_period[i]
        err_periods[i, i_eta] = get_period_err(
            ode_solution, p, 1.0, p_max, burn_in_days, model_params
        )
    end
end


heatmap(x_eta, x_period, log.(err_periods), ylim = (0, 4000), xlim = (0, 0.5))
plot(x_period, log.(err_periods), xlim = (0, 1000))
plot!(x_period, x_period .* 0.01)

