
include("dependencies.jl")

using JLD2

k = 32
C = 8.0
R = 1.5
gamma = 0.25
beta = R * gamma
rho = 0.003
b = 2^6
h = 8
c_jump_dist = Normal(2^9, 0.0)

model_params = make_model_parameters(
    k = k, beta = beta, gamma = gamma, C = C, rho = rho,
    b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none", eta = 0.0
)

heatmap(model_params.M)
ode_sparsity = ode_get_sparsity(model_params)

n_inf_0 = 0.0001
n_days = 365*20


Δt = 1
ode_solution = @time ode_solve(model_params, n_days, n_inf_0, ode_sparsity, saveat = Δt)


t = collect(0:Δt:n_days)

y = ode_solution(t)[ode_ix(c_sus, 1:model_params.S, model_params.S), (365*15):end ]

heatmap(min.(y,0.05))

plot(y', legend = false)


y = ode_solution(t)[ode_ix(c_inf, 1:model_params.S, model_params.S), :]
I_y = vec(sum(y, dims = 1))

plot(I_y)