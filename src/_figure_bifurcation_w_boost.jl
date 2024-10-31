
println("Bifurcation (with boost) plots")
include("dependencies.jl")

using JLD2

k = 32
C = 8.0
R = 1.5
gamma = 0.25
beta = R * gamma
rho_0 = 0.003
b = 2^3
h = 8
c_jump_dist = Normal(2^6, 2^5)

boost_scenarios = ["none" "linear" "loglinear"]
jump_dist_scenarios = [Normal(2^6, 2^5) Normal(2^6, 2^5) Normal(6, 5/8)]

model_params_0 = make_model_parameters(
    k = k, beta = beta, gamma = gamma, C = C, rho = rho_0,
    b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none"
)

ode_sparsity = ode_get_sparsity(model_params_0)

n_inf_0 = 0.01
n_days = 64000

x_rho = collect(0.0:0.00005:0.007)

y_fixed_I = zeros(length(x_rho), length(boost_scenarios))
y_I_sol = zeros(length(x_rho), length(boost_scenarios), n_days)
y_count_sol = zeros(length(x_rho), length(boost_scenarios), n_days)

Threads.@threads for i in eachindex(x_rho)
    for j in eachindex(boost_scenarios)
        model_params = make_model_parameters(
            k = k, beta = beta, gamma = gamma, C = C, rho = x_rho[i],
            b = b, h = h, c_jump_dist = jump_dist_scenarios[j]; boosting = boost_scenarios[j]
        )
    
        # Calculate the fixed point/steady state solution
        sus_vec, inf_vec = get_steady_state(model_params)
        y_fixed_I[i, j] = sum(inf_vec)
    
    
        # Calculate a (not-necessarily-stable) solution of the ODE over time
        ode_solution = ode_solve(model_params, n_days, n_inf_0, ode_sparsity)
    
        for d in 1:n_days
            y_I_sol[i, j, d] = sum(ode_solution(d)[ode_ix(c_inf, 1:model_params.S, model_params.S)])
            y_count_sol[i, j, d] = sum(ode_solution(d)[ode_ix(c_count, 1:model_params.S, model_params.S)])
        end

    end
    println(i)
end




plot(y_I_sol[20,1,1:1000])





y_inc_sol = diff(y_count_sol, dims = 3)

burn_in_days = 40000


period = zeros(Int, length(x_rho), 3)
attack_rate = zeros(length(x_rho), 3)

@showprogress Threads.@threads for i in eachindex(x_rho)
    for j in 1:3
        incidence = y_inc_sol[i, j, burn_in_days:end]
        period[i, j] = get_period(incidence)
        attack_rate[i, j] = get_periodic_attack_rate(incidence, period[i, j])
    end
end


plot(x_rho, period)
plot(x_rho, attack_rate)


jldsave("data/paper/bifurcations_w_boost.jld2"; x_rho, y_fixed_I, y_I_sol, y_inc_sol, period, attack_rate)