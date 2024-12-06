
include("dependencies.jl")

model_params = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    C = baseline_C, r = 0.04,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist
)

u_sol = get_steady_state(model_params, true)

using ForwardDiff


u0 = convert.(Double64, copy(u_sol))
du0 = copy(u0)

n_days_short = 40000
n_inf_0 = 0.0001
t = 0:n_days_short

ode_sparsity = ode_get_sparsity(model_params)
ode_solution = @time ode_solve(model_params, n_days_short, n_inf_0, ode_sparsity, saveat = 0.25)

inf = get_inf(ode_solution, t, model_params)

plot(inf[(end-5000):end])

period = get_period(ode_solution, model_params, n_days_burn_in, n_days_short, 0.25, 1e-6)

using ForwardDiff
J = ForwardDiff.jacobian((du, u) -> ode_step_no_count!(du, u, model_params, 0.0), du0, u0)

# Why are there zeros?
heatmap(J)

eigs_J = eigvals(J)
scatter(eigs_J)

scatter(2 * π ./ abs.(eigs_J), ylim = (0, 500))

period_eig = 2 * π ./ imag.(eigs_J)
scatter(period_eig)
hline!([period[1]])

period_eig[end]
period[1]


period[1] ./ period_eig

using LinearAlgebra

eigfact