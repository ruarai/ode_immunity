
include("dependencies.jl")

model_params = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    C = baseline_C, r = baseline_r,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist,
    eta = 0.0
)


ode_sparsity = ode_get_sparsity(model_params)

Δt = 0.25
ode_solution = @time ode_solve(model_params, n_days_short, n_inf_0, ode_sparsity, saveat = Δt)

t_seq = collect(0:Δt:n_days_short)

y = ode_solution(t_seq)[1:(model_params.S + 1), :]


jldsave("data/paper/ex_period.jld2"; t_seq, y)
