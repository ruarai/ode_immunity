include("dependencies.jl")


model_params_0 = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    C = baseline_C, r = baseline_r,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist
)


n_days_burn_in = 50000
n_days = 100000
t_seq = 0:n_days

periodic_Δt = 0.25
x_r = collect(0:0.002:0.15)

y_fixed_I = zeros(length(x_r))
y_I_sol = zeros(length(x_r), length(t_seq))
y_inc_sol = zeros(length(x_r), length(t_seq))


period = zeros(length(x_r), 3)

@showprogress Threads.@threads for i in eachindex(x_r)
    model_params = make_model_parameters(
        k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
        C = baseline_C, r = x_r[i],
        b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist
    )

    # Calculate the fixed point/steady state solution
    u_steady = get_steady_state(model_params)
    y_fixed_I[i] = u_steady[model_params.S + 1]


    # Calculate a (not-necessarily-stable) solution of the ODE over time
    ode_solution = ode_solve(model_params, n_days, n_inf_0, saveat_step = periodic_Δt)

    y_I_sol[i, :] = get_inf(ode_solution, t_seq, model_params)
    y_inc_sol[i, :] = get_inc(ode_solution, t_seq, model_params)

    period_mean, period_sd, period_n = get_period(ode_solution, model_params, n_days_burn_in, n_days, periodic_Δt, periodic_ϵ)

    period[i, :] = [period_mean period_sd period_n]
end

jldsave(
    "data/paper/bifurcations_w_boost.jld2";
    x_r, y_fixed_I, y_I_sol, y_inc_sol, period
)