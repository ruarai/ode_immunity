include("dependencies.jl")
using ForwardDiff

n_days_burn_in = 100 * 365
n_days = 200 * 365
t_seq = 0:n_days

periodic_Δt = 0.25
x_r = collect(0.0005:0.0005:0.05)

y_fixed_I = zeros(length(x_r))
y_I_sol = zeros(length(x_r), length(t_seq))
y_inc_sol = zeros(length(x_r), length(t_seq))

y_eigs = Vector{Vector{ComplexF64}}(undef, length(x_r))

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


    du0 = zeros(Float64, length(u_steady))
    u0 = copy(u_steady)
    J = ForwardDiff.jacobian((du, u) -> ode_step_no_count_boosting!(du, u, model_params_boosting, 0.0), du0, u0)

    y_eigs[i] = eigvals(J)
end


y_eigs_stack = stack(y_eigs)
y_real_eigs = real.(y_eigs_stack)
y_imag_eigs = imag.(y_eigs_stack)

jldsave(
    "data/paper/bifurcations_boosting.jld2";
    x_r, y_fixed_I, y_I_sol, y_inc_sol, period, y_real_eigs
)