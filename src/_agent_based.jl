


include("dependencies.jl")

using JLD2

k = 16
R = 1.5
gamma = 0.25
beta = R * gamma
b = 0.25
lambda = 0.008
m = 40
c_jump_dist = Normal(0.8, 0.05)

# Temporary model_parameters so sparsity can be generated
model_params = make_model_parameters(
    k = k, beta = beta, gamma = gamma, lambda = lambda,
    b = b, m = m, c_jump_dist = c_jump_dist; boosting = true
)


n_days = 1000

n_pop_size = 100000
n_track = 100


# track_state, I_t = simulate_agent_based(n_days, n_track, n_pop_size, model_params, 2)

# jldsave("data/paper/agent_based_with_boost.jld2"; track_state, I_t)


# Temporary model_parameters so sparsity can be generated
model_params_no_boost = make_model_parameters(
    k = k, beta = beta, gamma = gamma, lambda = lambda,
    b = b, m = m, c_jump_dist = c_jump_dist; boosting = false
)

track_state, I_t = simulate_agent_based(n_days, n_track, n_pop_size, model_params_no_boost, 2)

jldsave("data/paper/agent_based_no_boost.jld2"; track_state, I_t)