
# include("dependencies.jl")

model_params = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    C = baseline_C, r = 0.05,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist;
    eta = 0.5
)


u0 = zeros(Float64, n_compartments * model_params.S)
u0[ode_ix(c_sus, 1, model_params.S)] = 1.0 - 0.01
u0[ode_ix(c_inf, 1, model_params.S)] = 0.01

n_days_short = 40000
n_inf_0 = 0.0001
t = 0:n_days_short

ode_sparsity = ode_get_sparsity(model_params)
ode_solution = @time ode_solve(model_params, n_days_short, n_inf_0, ode_sparsity, saveat = 1.0)

inf = vec(sum(ode_solution(t)[ode_ix(c_inf, 1:model_params.S, model_params.S), :], dims = 1))[20000:1:end]

testchaos01(inf[1:20:end])