
println("Bifurcation (with boost) plots")
include("dependencies.jl")

model_params_0 = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    C = baseline_C, r = 0.05,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist
)

boost_scenarios = ["none" "loglinear"]

ode_sparsity = ode_get_sparsity(model_params_0)

ϵ = 10^-6

Δt = 0.25
n_inf_0 = 0.01
x_r = collect(0:0.001:0.1)

y_fixed_I = zeros(length(x_r), length(boost_scenarios))
y_I_sol = zeros(length(x_r), length(boost_scenarios), n_days)
y_inc_sol = zeros(length(x_r), length(boost_scenarios), n_days)


period = zeros(length(x_r), 3, 3)
attack_rate = zeros(length(x_r), 3)

@showprogress Threads.@threads for i in eachindex(x_r)
    for j in eachindex(boost_scenarios)
        model_params = make_model_parameters(
            k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
            C = baseline_C, r = x_r[i],
            b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist;
            boosting = boost_scenarios[j]
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


        period_mean, period_sd, period_n = get_period(ode_solution, model_params, n_days_burn_in, n_days, Δt, ϵ)

        period[i, j, :] = [period_mean period_sd period_n]
        attack_rate[i, j] = get_periodic_attack_rate(ode_solution, model_params, n_days_burn_in, n_days, period_mean)
    end
end



jldsave(
    "data/paper/bifurcations_w_boost.jld2";
    x_r, y_fixed_I, y_I_sol, y_inc_sol, period, attack_rate
)