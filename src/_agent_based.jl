


include("dependencies.jl")

using JLD2


model_params = make_model_parameters(
    k = 128, beta = 1.5 * 0.25, gamma = 0.25, lambda = 0.002,
    b = 0.6, m = 20, c_jump_dist = Normal(0.8, 0.1); boosting = false
)

plot(model_params.p_acq)

heatmap(model_params.M)

n_days = 2000

n_pop_size = 100000
n_track = 10000


track_state, I_t = simulate_agent_based(n_days, n_track, n_pop_size, model_params, 1)

plot(I_t)
track_corr = (track_state .- 1) .% model_params.S

track_acq = 1 ./ (1 .+ exp.(-model_params.m .* (track_corr / model_params.k .- model_params.b)))

plot(track_corr[:,10:300], legend = false, lc = "black", linealpha = 0.1)
heatmap(track_corr[:,10:1000]')

plot(track_acq[:,10:20], legend = false)
heatmap(track_acq[:,10:1000]')


jldsave("data/paper/agent_based_transition.jld2"; track_state, track_acq, I_t)