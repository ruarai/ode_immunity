
include("src/dependencies.jl")

c_max = 8
k = 2
R = 1.5
gamma = 0.25
beta = R * gamma
lambda = 0.02
b = 2
m = 5
c_jump_dist = Normal(4.5, 0.5)

model_params = make_model_parameters(
    c_max = c_max, k = k, beta = beta, gamma = gamma, lambda = lambda,
    b = b, m = m, c_jump_dist = c_jump_dist
)
ode_sparsity = ode_get_sparsity(model_params)


plot(model_params.c_levels, model_params.p_acq)
heatmap(model_params.M)


n_inf_0 = 0.0001
n_days = 365*10


ode_solution = @time ode_solve(model_params, n_days, n_inf_0, ode_sparsity)

sol_I = zeros(n_days, model_params.N)
sol_S = zeros(n_days, model_params.N)

for d in 1:n_days, i in 1:model_params.N
    sol_S[d, :] = ode_solution(d)[ode_ix(c_sus, 1:model_params.N, model_params.N)]
    sol_I[d, :] = ode_solution(d)[ode_ix(c_inf, 1:model_params.N, model_params.N)]
end

plot(sum(sol_I, dims = 2))
heatmap(min.(sol_S[1:500,:], 0.1)')

jldsave("data/anziam2024/basic.jld2"; model_params, sol_I, sol_S)

