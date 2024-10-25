
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
    b = b, h = h, c_jump_dist = c_jump_dist; boosting = "none"
)

plot(model_params.c_levels)

p_acq = (model_params.c_levels .^ h) ./ (b ^ h .+ model_params.c_levels .^ h)
plot(log2.(model_params.c_levels), p_acq)


heatmap(min.(model_params.M, 0.05))

ode_sparsity = ode_get_sparsity(model_params)

n_inf_0 = 0.0001
n_days = 365*20


sol_t = zeros(3, n_days, 2, model_params.S)

ode_solution = @time ode_solve(model_params, n_days, n_inf_0, ode_sparsity)

for d in 1:n_days, i in 1:model_params.S
    sol_t[1, d, 1, :] = ode_solution(d)[ode_ix(c_sus, 1:model_params.S, model_params.S)]
    sol_t[1, d, 2, :] = ode_solution(d)[ode_ix(c_inf, 1:model_params.S, model_params.S)]
end

model_params_boosting = make_model_parameters(
    k = k, beta = beta, gamma = gamma, C = C, rho = rho,
    b = b, h = h, c_jump_dist = c_jump_dist; boosting = "linear"
)


heatmap(model_params_boosting.M)

ode_solution = @time ode_solve(model_params_boosting, n_days, n_inf_0, ode_sparsity)

for d in 1:n_days, i in 1:model_params.S
    sol_t[2, d, 1, :] = ode_solution(d)[ode_ix(c_sus, 1:model_params.S, model_params.S)]
    sol_t[2, d, 2, :] = ode_solution(d)[ode_ix(c_inf, 1:model_params.S, model_params.S)]
end




model_params_boosting_loglinear = make_model_parameters(
    k = k, beta = beta, gamma = gamma, C = C, rho = rho,
    b = b, h = h, c_jump_dist = Normal(6, 5/8); boosting = "loglinear"
)
heatmap(model_params_boosting_loglinear.M)
heatmap(model_params.M)


ode_solution = @time ode_solve(model_params_boosting_loglinear, n_days, n_inf_0, ode_sparsity)

for d in 1:n_days, i in 1:model_params.S
    sol_t[3, d, 1, :] = ode_solution(d)[ode_ix(c_sus, 1:model_params.S, model_params.S)]
    sol_t[3, d, 2, :] = ode_solution(d)[ode_ix(c_inf, 1:model_params.S, model_params.S)]
end

plot(sol_t[1,:,2,:])

c_levels = model_params.c_levels

boosting_matrices = zeros(3, model_params.S, model_params.S)
boosting_matrices[1,:,:] = model_params.M
boosting_matrices[2,:,:] = model_params_boosting.M
boosting_matrices[3,:,:] = model_params_boosting_loglinear.M


jldsave("data/paper/basic_boosting.jld2"; sol_t, c_levels, boosting_matrices)
