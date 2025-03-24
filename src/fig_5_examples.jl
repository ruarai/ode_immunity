include("dependencies.jl")

n_days = 150000

t_seq = collect(1:n_days)

# 0.0, 0.018, "i", factor(0),
# 0.07, 0.018, "ii", factor(4.5),
# 0.2, 0.018, "iii", factor(1),
# 0.27, 0.018, "iv", factor(2),
# 0.37, 0.018, "v", factor(9)
x_eta = [0 0.07 0.2 0.27 0.37]
x_r = [0.018 0.018 0.018 0.018 0.018]

y_inf = zeros(length(x_eta), length(t_seq))
y_sus = zeros(length(x_eta), length(t_seq), baseline_k + 1)

@showprogress Threads.@threads for i in eachindex(x_eta)
    model_params = make_model_parameters(
        k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
        a = baseline_a, r = x_r[i],
        b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist;
        eta = x_eta[i]
    )

    ode_solution = ode_solve(model_params, n_days, n_inf_0)

    y_inf[i, :] = get_inf(ode_solution, t_seq, model_params)
    y_sus[i, :, :] = get_sus(ode_solution, t_seq, model_params)
end

jldsave("data/paper/period_over_grid_examples.jld2"; x_eta, x_r, y_inf, y_sus)