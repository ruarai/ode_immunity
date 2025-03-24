include("dependencies.jl")

model_params_no_boosting = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    a = baseline_a, r = baseline_r,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist;
    boosting = "independent"
)

model_params_boosting_A = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    a = baseline_a, r = baseline_r,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist;
    boosting = "none"
)


model_params_boosting_B = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    a = baseline_a, r = baseline_r,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist;
    boosting = "multiplicative"
)


n_days = 8000
seq_t = collect(0:n_days)
c_levels = model_params_no_boosting.c_levels


ode_solution_no_boosting = @time ode_solve(model_params_no_boosting, n_days, n_inf_0)
sus, inf, inc = get_results(ode_solution_no_boosting, seq_t, model_params_no_boosting)

ode_solution_boosting_A = @time ode_solve_boosting(model_params_boosting_A, n_days, n_inf_0)
sus_A, inf_A, inc_A = get_summ_boosting(ode_solution_boosting_A, seq_t, model_params_boosting_A)

ode_solution_boosting_B = @time ode_solve_boosting(model_params_boosting_B, n_days, n_inf_0)
sus_B, inf_B, inc_B = get_summ_boosting(ode_solution_boosting_B, seq_t, model_params_boosting_B)

results_sus = reduce((x,y) -> cat(x, y, dims = 3), [sus, sus_A, sus_B])
results_inf = reduce(hcat, [inf, inf_A, inf_B])
results_inc = reduce(hcat, [inc, inc_A, inc_B])


jldsave("data/supp_boosting_basic.jld2"; seq_t, c_levels, results_sus, results_inf, results_inc)
