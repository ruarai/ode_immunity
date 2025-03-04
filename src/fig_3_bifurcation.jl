include("dependencies.jl")
using ForwardDiff

n_days_burn_in = 50000
n_days = 100000
t_seq = 0:n_days

periodic_Δt = 0.25
x_r = collect(0.002:0.002:0.15)

y_fixed_I = zeros(length(x_r))
y_I_sol = zeros(length(x_r), length(t_seq))
y_inc_sol = zeros(length(x_r), length(t_seq))
y_means = zeros(length(x_r), length(t_seq))


y_eigs = Vector{Vector{ComplexF64}}(undef, length(x_r))

period = zeros(length(x_r), 3)

@showprogress Threads.@threads for i in eachindex(x_r)
    model_params = make_model_parameters(
        k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
        a = baseline_a, r = x_r[i],
        b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist
    )

    # Calculate the fixed point/steady state solution
    u_steady = get_steady_state(model_params)
    y_fixed_I[i] = u_steady[model_params.S + 1]


    # Calculate a (not-necessarily-stable) solution of the ODE over time
    ode_solution = ode_solve(model_params, n_days, n_inf_0, saveat_step = periodic_Δt)

    y_I_sol[i, :] = get_inf(ode_solution, t_seq, model_params)
    y_inc_sol[i, :] = get_inc(ode_solution, t_seq, model_params)



    S_sol = get_sus(ode_solution, t_seq, model_params)
    y_means[i, :] = [sum(S_sol[t, :] / (1 - y_I_sol[i, t]) .* model_params.c_levels) for t in eachindex(t_seq)]

    period_mean, period_sd, period_n = get_period(ode_solution, model_params, n_days_burn_in, n_days, periodic_Δt, periodic_ϵ)

    period[i, :] = [period_mean period_sd period_n]

    du0 = zeros(Float64, model_params.S)
    u0 = u_steady[1:33]
    J = ForwardDiff.jacobian((du, u) -> ode_step_minimal!(du, u, model_params, 0.0), du0, u0)

    y_eigs[i] = eigvals(J)
end

y_eigs_stack = stack(y_eigs)
y_real_eigs = real.(y_eigs_stack)
y_imag_eigs = imag.(y_eigs_stack)

plot(y_real_eigs', y_imag_eigs', seriestype = :path, legend = false, xlim = (-0.01, 0.01))

scatter(y_eigs_stack[:, 3], xlim = (-0.1, 0.1), ylim = (-0.1, 0.1))

jldsave(
    "data/paper/bifurcations.jld2";
    x_r, y_fixed_I, y_I_sol, y_inc_sol, y_means, period, y_real_eigs, y_imag_eigs
)