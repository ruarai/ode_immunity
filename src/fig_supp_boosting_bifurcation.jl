include("dependencies.jl")

n_days_burn_in = 100 * 365
n_days = 200 * 365
t_seq = 0:n_days

periodic_Δt = 0.25
x_r = collect(0:0.002:0.15)

y_fixed_I = zeros(length(x_r))
y_I_sol = zeros(length(x_r), length(t_seq))
y_inc_sol = zeros(length(x_r), length(t_seq))


period = zeros(length(x_r), 3)

@showprogress Threads.@threads for i in eachindex(x_r)
    model_params_boosting = make_model_parameters(
        k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
        a = baseline_a, r = x_r[i],
        b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist;

        boosting = "multiplicative"
    )

    # Calculate the fixed point/steady state solution
    u_steady = get_steady_state_boosting(model_params_boosting)
    y_fixed_I[i] = sum(u_steady[ode_ix_boosting(c_inf, model_params_boosting.S, 1:model_params_boosting.S)])


    # Calculate a (not-necessarily-stable) solution of the ODE over time
    ode_solution = ode_solve_boosting(model_params_boosting, n_days, n_inf_0, saveat_step = periodic_Δt)

    sus, inf, inc = get_summ_boosting(ode_solution, t_seq, model_params_boosting)

    y_I_sol[i, :] = inf
    y_inc_sol[i, :] = inc

    period_mean, period_sd, period_n = get_period(ode_solution, model_params_boosting, n_days_burn_in, n_days, periodic_Δt, periodic_ϵ)

    period[i, :] = [period_mean period_sd period_n]
end


plot(y_fixed_I, ylim = (0, 0.1))
plot!(mean(y_I_sol[:, 50000:end], dims = 2))

jldsave(
    "data/paper/bifurcations_boosting.jld2";
    x_r, y_fixed_I, y_I_sol, y_inc_sol, period
)