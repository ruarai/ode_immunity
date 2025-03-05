include("dependencies.jl")

model_params = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    a = baseline_a, r = 0.003,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist
)

n_days = 10000
seq_t = collect(0:n_days)


u0_fixed = get_steady_state(model_params)

using Random
rng = MersenneTwister(3)

u0_noisy = u0_fixed .+ rand(rng, LogNormal(-4.8, 1), model_params.S + 1)
u0_noisy = u0_noisy / sum(u0_noisy)

ode_solution = @time ode_solve(model_params, n_days, 0, u0 = u0_noisy);

sus, inf, inc = get_results(ode_solution, seq_t, model_params)

jldsave("data/paper/bifurcation_stable_fixed.jld2"; seq_t, inf)