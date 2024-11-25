include("dependencies.jl")


model_params_0 = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    C = baseline_C, r = baseline_r,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist
)

ode_sparsity = ode_get_sparsity(model_params_0)

n_inf_0 = 0.01
t_seq = collect(1:n_days)

# 0.3, 0.0028, "D",
# 0.02, 0.0028, "C",
# 0.3, 0.0023, "A",
# 0.3, 0.0033, "B"
x_eta = [0.3 0.02 0.3 0.3]
x_rho = [0.0028 0.0028 0.0023 0.0033] # TODO

y_inf = zeros(length(x_eta), n_days)
sol_t = zeros(length(x_eta), n_days, 2, model_params_0.S)

@showprogress Threads.@threads for i in eachindex(x_eta)

    model_params = make_model_parameters(
        k = k, beta = beta, gamma = gamma, C = C, rho = x_rho[i],
        b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none", eta = x_eta[i]
    )

    # Calculate a (not-necessarily-stable) solution of the ODE over time
    ode_solution = ode_solve(model_params, n_days, n_inf_0, ode_sparsity)

    sol_t[i, :, 1, :] = ode_solution(t_seq)[1:model_params_0.S, :]'
    sol_t[i, :, 2, :] = ode_solution(t_seq)[(model_params_0.S + 1):(model_params_0.S * 2), :]'

    println(i)
end

jldsave("data/paper/period_over_grid_examples.jld2"; x_eta, x_rho, sol_t, t_seq)