
include("dependencies.jl")

using JLD2

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

n_inf_0 = 0.0001
n_days = 365*20

ode_solution = @time ode_solve(model_params, n_days, n_inf_0, ode_sparsity)

sol_I = zeros(n_days, model_params.S)
sol_S = zeros(n_days, model_params.S)

for d in 1:n_days, i in 1:model_params.S
    sol_S[d, :] = ode_solution(d)[ode_ix(c_sus, 1:model_params.S, model_params.S)]
    sol_I[d, :] = ode_solution(d)[ode_ix(c_inf, 1:model_params.S, model_params.S)]
end


c_levels = model_params.c_levels

jldsave("data/paper/basic.jld2"; c_levels, sol_I, sol_S)



model_params_boosting = make_model_parameters(
    k = k, beta = beta, gamma = gamma, lambda = lambda,
    b = b, m = m, c_jump_dist = c_jump_dist; boosting = true
)

n_inf_0 = 0.0001
n_days = 365*20

ode_solution = @time ode_solve(model_params_boosting, n_days, n_inf_0, ode_sparsity)

sol_I = zeros(n_days, model_params.S)
sol_S = zeros(n_days, model_params.S)

for d in 1:n_days, i in 1:model_params.S
    sol_S[d, :] = ode_solution(d)[ode_ix(c_sus, 1:model_params.S, model_params.S)]
    sol_I[d, :] = ode_solution(d)[ode_ix(c_inf, 1:model_params.S, model_params.S)]
end


jldsave("data/paper/basic_boosting.jld2"; sol_I, sol_S)