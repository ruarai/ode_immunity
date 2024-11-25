
include("dependencies.jl")

model_params = make_model_parameters(
    k = baseline_k, beta = baseline_beta, gamma = baseline_gamma,
    C = baseline_C, r = 0.02,
    b = baseline_b, h = baseline_h, c_jump_dist = baseline_c_jump_dist;
    eta = 0.5
)

ode_sparsity = ode_get_sparsity(model_params)

n_days_short = 80000
n_inf_0 = 0.0001
t = 0:n_days_short

sol_t = zeros(2, n_days_short + 1, 2, model_params.S)

ode_solution = @time ode_solve(model_params, n_days_short, n_inf_0, ode_sparsity, saveat = 0.25)



inf_t = ode_solution(t)[ode_ix(c_inf, 1:model_params.S, model_params.S), :]
inf_t = vec(sum(inf_t, dims = 1))

x_t = ode_solution(t)[ode_ix(c_sus, 1, model_params.S), :]

plot(inf_t)

# x_t = log10.(x_t)


pq = zeros(2, n_days_short)

c = pi * 0.15
for i in 2:n_days_short
    pq[1, i] = pq[1, i - 1] + x_t[i] * cos(i * c)
    pq[2, i] = pq[1, i - 1] + x_t[i] * sin(i * c)
end
plot(pq[1, 5000:end], pq[2, 5000:end], linetype = :path)


plot(pq')


