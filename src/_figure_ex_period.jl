
include("dependencies.jl")

using JLD2

k = 32
C = 8.0
R = 1.5
gamma = 0.25
beta = R * gamma
rho = 0.003
b = 2^3
h = 8
c_jump_dist = Normal(2^6, 2^5)

model_params = make_model_parameters(
    k = k, beta = beta, gamma = gamma, C = C, rho = rho,
    b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none", eta = 0.2
)

ode_sparsity = ode_get_sparsity(model_params)

n_inf_0 = 0.0001
n_days = 365*40


Δt = 0.25
ode_solution = @time ode_solve(model_params, n_days, n_inf_0, ode_sparsity, saveat = Δt)


t = collect(0:Δt:n_days)
y = ode_solution(t)[1:(model_params.S * 2), :]


jldsave("data/paper/ex_period.jld2"; t, y)
