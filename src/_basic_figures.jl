
include("dependencies.jl")

model_params = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    C = baseline_C, r = baseline_r,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist
)

ode_sparsity = ode_get_sparsity(model_params)

n_days_short = 10000
n_inf_0 = 0.0001
t = 0:n_days_short

sol_t = zeros(2, n_days_short + 1, 2, model_params.S)

ode_solution = @time ode_solve(model_params, n_days_short, n_inf_0, ode_sparsity, saveat = 0.25)

sol_t[1, :, 1, :] .= ode_solution(t)[ode_ix(c_sus, 1:model_params.S, model_params.S), :]'
sol_t[1, :, 2, :] .= ode_solution(t)[ode_ix(c_inf, 1:model_params.S, model_params.S), :]'

model_params_boosting = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    C = baseline_C, r = baseline_r,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist;
    boosting = "loglinear"
)

ode_solution_boosting = @time ode_solve(model_params_boosting, n_days_short, n_inf_0, ode_sparsity)

sol_t[2, :, 1, :] .= ode_solution_boosting(t)[ode_ix(c_sus, 1:model_params.S, model_params.S), :]'
sol_t[2, :, 2, :] .= ode_solution_boosting(t)[ode_ix(c_inf, 1:model_params.S, model_params.S), :]'

c_levels = model_params.c_levels

boosting_matrices = zeros(2, model_params.S, model_params.S)
boosting_matrices[1,:,:] = model_params.M
boosting_matrices[2,:,:] = model_params_boosting.M


p_acq = model_params.p_acq
seq_t = collect(t)
jldsave("data/paper/basic_boosting.jld2"; seq_t, sol_t, c_levels, boosting_matrices, p_acq)
