
include("dependencies.jl")

using JLD2

k = 16
R = 1.5
gamma = 0.25
beta = R * gamma
lambda = 0.003
b = 0.25
m = 40
c_jump_dist = Normal(0.8, 0.05)

model_params = make_model_parameters(
    k = k, beta = beta, gamma = gamma, lambda = lambda,
    b = b, m = m, c_jump_dist = c_jump_dist; boosting = false
)

ode_sparsity = ode_get_sparsity(model_params)

n_inf_0 = 0.0001
n_days = 365*10

ode_solution = @time ode_solve(model_params, n_days, n_inf_0, ode_sparsity)


# Reducing mean to match mean protection threshold
# But we also need to increase variance here to capture the increased
# uncertainty from waning process + possibly also the p_acq function

c_jump_dist_dde = Normal(0.8 - 0.25, 0.22)
lambda = 0.003

τ = c_jump_dist_dde / lambda

t_days, S_t, I_t, R_t, dI_t = sim_dde(;
    R0 = model_params.beta / model_params.gamma ,
    γ = model_params.gamma,
    τ = τ,
    dt = 0.01,
    n_days, 
    n_inf_0 = n_inf_0,
    return_each_step = false
)


I_structured = zeros(n_days)
I_dde = zeros(n_days)
for d in 1:n_days
    I_structured[d] = sum(ode_solution(d)[ode_ix(c_inf, 1:model_params.S, model_params.S)])
    I_dde[d] = I_t[d]
end

# Plot demonstrating reduction in subsequent peak sizes
plot(I_structured)
plot!(I_dde)


# Now, ? plot demonstrating why this occurs - change in net susceptibility

S_t_ode = zeros(n_days, model_params.S)
pop_protection_ode = zeros(n_days)

S_t_dde = zeros(n_days)

for d in 1:n_days
    S_t_ode[d, :] = ode_solution(d)[ode_ix(c_sus, 1:model_params.S, model_params.S)]

    pop_protection_ode[d] = sum(S_t_ode[d, :] .* model_params.p_acq)

    S_t_dde[d] = S_t[d]
end

plot(pop_protection)
plot!(1 .- S_t_dde)

plot(S_t_ode, legend = false)