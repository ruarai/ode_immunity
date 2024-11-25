
include("dependencies.jl")

model_params = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    C = baseline_C, r = baseline_r,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist,
    eta = 0.2
)

ode_sparsity = ode_get_sparsity(model_params)
n_inf_0 = 0.0001


Δt = 0.25
ode_solution = @time ode_solve(model_params, n_days, n_inf_0, ode_sparsity, saveat = Δt)

t = collect(0:Δt:n_days)
y = ode_solution(t)[1:(model_params.S * 2), :]


jldsave("data/paper/ex_period.jld2"; t, y)
