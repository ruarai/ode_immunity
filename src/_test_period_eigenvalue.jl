
# include("dependencies.jl")

model_params = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    C = baseline_C, r = 0.04,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist
)

S_steady, I_steady = get_steady_state(model_params, true)

using ForwardDiff


u0 = zeros(Double64,  2 * model_params.S)
u0[ode_ix(c_sus, 1:model_params.S, model_params.S)] = S_steady
u0[ode_ix(c_inf, 1:model_params.S, model_params.S)] = I_steady

du0 = copy(u0)

n_days_short = 40000
n_inf_0 = 0.0001
t = 0:n_days_short

ode_sparsity = ode_get_sparsity(model_params)
ode_solution = @time ode_solve(model_params, n_days_short, n_inf_0, ode_sparsity, saveat = 0.25)

inf = vec(sum(ode_solution(t)[ode_ix(c_inf, 1:model_params.S, model_params.S), :], dims = 1))

plot(inf[(end-5000):end])

period = get_period(ode_solution, model_params, 30000, n_days_short, 0.25, 1e-6)


J = ForwardDiff.jacobian((du, u) -> ode_step_no_count!(du, u, model_params, 0.0), du0, u0)

# Why are there zeros?
heatmap(J)

eigs_J = eigvals(J, maxiters = 10000)

eigs_J = eigvals(J)

period_eig = 2 * π ./ imag.(eigs_J)[end]


period[1] ./ period_eig