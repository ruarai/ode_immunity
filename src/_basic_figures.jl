
include("dependencies.jl")

k = 16
R = 1.5
gamma = 0.25
beta = R * gamma
lambda = 0.002
b = 0.25
m = 40
c_jump_dist = Normal(0.8, 0.05)

model_params = make_model_parameters(
    k = k, beta = beta, gamma = gamma, lambda = lambda,
    b = b, m = m, c_jump_dist = c_jump_dist; boosting = false
)

ode_sparsity = ode_get_sparsity(model_params)


plot(model_params.c_levels, model_params.p_acq)
heatmap(model_params.c_levels, model_params.c_levels, model_params.M)

plot(model_params.c_levels, model_params.M[:,1])

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
heatmap(min.(sol_S[:,:], 0.1)')
heatmap(min.(sol_S[1:500,:], 0.1)')


mat_jump = model_params.M
c_levels = model_params.c_levels

mat_jump_no_boost = build_immunity_matrix_no_boost(model_params.N, model_params.c_levels, model_params.c_jump_dist)

jldsave("data/paper/basic.jld2"; mat_jump, c_levels, mat_jump_no_boost, sol_I, sol_S)

