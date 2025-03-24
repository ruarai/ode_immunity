include("dependencies.jl")

n_days_burn_in = 1000 * 365
n_days = n_days_burn_in + 250 * 365

t_seq = 1:n_days
t_post_burn_in = n_days_burn_in:n_days

# 0.25, 0.01, "iv",
# 0.25, 0.013, "iii",
# 0.25, 0.016, "ii",
# 0.05, 0.016, "i"
x_eta = [0.25 0.25 0.25 0.05]
x_r = [0.01 0.013 0.0165 0.0165]

y_sus = zeros(length(x_eta), length(t_seq), baseline_k + 1)
y_inf = zeros(length(x_eta), length(t_seq))
y_inc = zeros(length(x_eta), length(t_seq))
y_seasonality = zeros(length(x_eta), 2)

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


    seasonality_x, seasonality_y = get_seasonality_coordinates(ode_solution, t_post_burn_in, model_params)

    y_seasonality[i, :] = [seasonality_x, seasonality_y]
end

jldsave("data/seasonality_bias_examples.jld2"; x_eta, x_r, y_sus, y_inf, y_inc, y_seasonality)
