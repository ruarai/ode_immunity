
include("dependencies.jl")

model_params = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    C = baseline_C, r = baseline_r,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist
)

ode_sparsity = ode_get_sparsity(model_params)

n_days_short = 10000
seq_t = collect(0:n_days_short)

ode_solution = @time ode_solve(model_params, n_days_short, n_inf_0, ode_sparsity, saveat = 0.25)

sus, inf, inc = get_results(ode_solution, seq_t, model_params)

c_levels = model_params.c_levels
p_acq = model_params.p_acq
jldsave("data/paper/basic_boosting.jld2"; seq_t, sus, inf, inc, c_levels, p_acq)