include("dependencies.jl")


model_params = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    a = baseline_a, r = baseline_r,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist,
    eta = 0.5
)

using Random
r = Random.Xoshiro(0)

n_samples = 10000

u = zeros(model_params.S + 2, n_samples)
t = rand(r, n_samples) .* 365

res_du = zeros(model_params.S + 2, n_samples)

for i in 1:n_samples
    u[:, i] = rand(r, model_params.S + 2)
    u[:, i] = u[:, i] / sum(u[:, i])

    du = zeros(model_params.S + 2)

    ode_step!(du, u[:, i], model_params, t[i])

    res_du[:, i] = du
end

heatmap(res_du)