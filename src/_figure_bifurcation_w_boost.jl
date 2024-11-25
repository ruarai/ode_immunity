
println("Bifurcation (with boost) plots")
include("dependencies.jl")

k = 32
C = 8.0
R = 1.5
gamma = 0.25
beta = R * gamma
rho_0 = 0.003
b = 10^3
h = 3
c_jump_dist = Normal(10^6, 10^7)

boost_scenarios = ["none" "linear" "loglinear"]
jump_dist_scenarios = [Normal(10^6, 10^7) Normal(10^6, 10^7) Normal(10, 7/8)]

model_params_0 = make_model_parameters(
    k = k, beta = beta, gamma = gamma, C = C, rho = rho_0,
    b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none"
)

ode_sparsity = ode_get_sparsity(model_params_0)

ϵ = 10^-6

Δt = 0.25
n_inf_0 = 0.01

x_rho = collect(0.0:0.00005:0.007)

y_fixed_I = zeros(length(x_rho), length(boost_scenarios))
y_I_sol = zeros(length(x_rho), length(boost_scenarios), n_days)
y_inc_sol = zeros(length(x_rho), length(boost_scenarios), n_days)


period = zeros(length(x_rho), 3, 3)
attack_rate = zeros(length(x_rho), 3)

@showprogress Threads.@threads for i in eachindex(x_rho)
    for j in eachindex(boost_scenarios)
        model_params = make_model_parameters(
            k = k, beta = beta, gamma = gamma, C = C, rho = x_rho[i],
            b = b, h = h, c_jump_dist = jump_dist_scenarios[j]; boosting = boost_scenarios[j]
        )
    
        # Calculate the fixed point/steady state solution
        sus_vec, inf_vec = get_steady_state(model_params)
        y_fixed_I[i, j] = sum(inf_vec)
    
    
        # Calculate a (not-necessarily-stable) solution of the ODE over time
        ode_solution = ode_solve(model_params, n_days, n_inf_0, ode_sparsity, saveat = Δt)
    
        for d in 1:n_days
            y_I_sol[i, j, d] = sum(ode_solution(d)[ode_ix(c_inf, 1:model_params.S, model_params.S)])
            y_inc_sol[i, j, d] = sum(ode_solution(d)[ode_ix(c_count, 1:model_params.S, model_params.S)])
        end


        period_mean, period_sd, period_n = get_period(ode_solution, model_params, burn_in_days, n_days, Δt, ϵ)

        period[i, j, :] = [period_mean period_sd period_n]
        attack_rate[i, j] = get_periodic_attack_rate(ode_solution, model_params, burn_in_days, n_days, period_mean)


    end
end

y_inc_sol = diff(y_inc_sol, dims = 3)

ix_good = period[:,:,2] .< 1

plot(x_rho, period[:,:,1])
plot(x_rho, period[:,:,2])
plot(x_rho, attack_rate)


jldsave("data/paper/bifurcations_w_boost.jld2"; x_rho, y_fixed_I, y_I_sol, y_inc_sol, period, attack_rate)