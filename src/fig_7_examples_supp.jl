include("dependencies.jl")

n_days = 150000

t_seq = collect(1:n_days)

# 0.015, 0.015
# 0.22, 0.0292
x_eta = [0.015 0.022]
x_r = [0.015 0.0292]

y_sus = zeros(length(x_eta), length(t_seq), baseline_k + 1)
y_inf = zeros(length(x_eta), length(t_seq))
y_inc = zeros(length(x_eta), length(t_seq))


@showprogress Threads.@threads for i in eachindex(x_eta)
    model_params = make_model_parameters(
        k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
        a = baseline_a, r = x_r[i],
        b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist;
        eta = x_eta[i]
    )

    ode_solution = ode_solve(model_params, n_days, n_inf_0)

    sus, inf, inc = get_results(ode_solution, t_seq, model_params)

    y_sus[i, :, :] = sus
    y_inf[i, :] = inf
    y_inc[i, :] = inc
end



jldsave("data/paper/seasonality_bias_examples_supp.jld2"; x_eta, x_r, y_sus, y_inf, y_inc)